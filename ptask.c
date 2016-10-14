
/**
 * @file ptask.c
 * @brief parallel task manager
 */

/* import unittest */
#define UNITTEST_UNIQUE_ID			100
#define UNITTEST 					1

#include  "unittest.h"

#include <stdint.h>
#include <pthread.h>
#include "ptask.h"
#include "queue.h"
#include "log.h"


/**
 * @macro _force_inline
 * @brief inline directive for gcc-compatible compilers
 */
#define _force_inline	inline

/* constants */
#define PTASK_DISPATCHER_EXIT		( (void *)((int64_t)-1) )

/**
 * context container
 */
struct ptask_container_s {
	pthread_t th;
	queue_t *inq, *outq;
	void *worker_arg;
	void *(*worker)(void *worker_arg, void *item);
};
struct ptask_context_s {
	int64_t num_threads;
	int64_t queue_size;

	/* worker for the single thread mode */
	void *worker_arg;
	void *(*worker)(void *worker_arg, void *item);

	/* threads */
	struct ptask_container_s c[];
};

/**
 * @fn ptask_dispatcher
 */
static
void *ptask_dispatcher(void *s)
{
	struct ptask_container_s *c = (struct ptask_container_s *)s;
	while(1) {
		void *item;		/* pointer to item */
		if(queue_get_wait(c->inq, (void **)&item) == 0) {
			if(item == PTASK_DISPATCHER_EXIT) { break; }
			debug("worker_arg(%p), item(%p)", c->worker_arg, item);
			void *result = c->worker(c->worker_arg, item);
			queue_put_wait(c->outq, result);
		}
	}
	return(NULL);
}

/**
 * @fn ptask_get_num_cores
 */
int64_t ptask_get_num_cores(
	void)
{
	/** !! fixme !! */
	return(4);
}

/**
 * @fn ptask_init
 */
ptask_t *ptask_init(
	void *(*worker)(void *arg, void *item),
	void *worker_arg[],
	int64_t num_threads,
	int64_t queue_size)
{
	/* negative thread count is not allowed */
	if(num_threads < 0) {
		return(NULL);
	}

	/* worker must not be NULL */
	if(worker == NULL) {
		return(NULL);
	}

	/* malloc context */
	uint64_t malloc_size = sizeof(struct ptask_context_s)
		+ num_threads * sizeof(struct ptask_container_s);
	struct ptask_context_s *ctx = (struct ptask_context_s *)malloc(malloc_size);
	memset(ctx, 0, malloc_size);

	/* store params */
	ctx->num_threads = num_threads;
	ctx->queue_size = (queue_size <= 0) ? 1024 : queue_size;

	/* for single thread mode */
	ctx->worker_arg = (worker_arg != NULL) ? worker_arg[0] : NULL;
	ctx->worker = worker;

	/* create threads */
	for(int64_t i = 0; i < num_threads; i++) {
		debug("%lld\n", i);
		ctx->c[i].inq = queue_create_limited(ctx->queue_size);
		ctx->c[i].outq = queue_create_limited(ctx->queue_size);
		ctx->c[i].worker_arg = (worker_arg != NULL) ? worker_arg[i] : NULL;
		ctx->c[i].worker = worker;
		pthread_create(&ctx->c[i].th, NULL, ptask_dispatcher, (void *)&ctx->c[i]);
	}
	debug("ctx(%p)", ctx);
	return((ptask_t *)ctx);
}

/**
 * @fn ptask_clean
 */
void ptask_clean(
	ptask_t *_ctx)
{
	struct ptask_context_s *ctx = (struct ptask_context_s *)_ctx;
	if(ctx == NULL) {
		return;
	}

	debug("ctx(%p), num_threads(%lld)", ctx, ctx->num_threads);

	for(int64_t i = 0; i < ctx->num_threads; i++) {
		debug("%lld\n", i);

		/* send destroy signal */
		queue_put(ctx->c[i].inq, PTASK_DISPATCHER_EXIT);

		/* wait for the thread terminates */
		void *status;
		pthread_join(ctx->c[i].th, &status);

		/* cleanup queues */
		queue_flush(ctx->c[i].inq);
		queue_flush(ctx->c[i].outq);
		queue_destroy(ctx->c[i].inq);
		queue_destroy(ctx->c[i].outq);
	}
	free(ctx);
	return;
}

/**
 * @fn ptask_parallel
 */
int ptask_parallel(
	ptask_t *_ctx,
	void *items[],
	void *results[])
{
	/* check context */
	struct ptask_context_s *ctx = (struct ptask_context_s *)_ctx;
	if(ctx == NULL) {
		return(PTASK_ERROR);
	}

	if(ctx->num_threads == 0) {
		/* single thread mode */
		void *result = ctx->worker(ctx->worker_arg,
			(items != NULL) ? items[0] : NULL);
		debug("single thread mode item(%p), result(%p)", items[0], result);
		if(results != NULL) {
			results[0] = result;
		}
		return(PTASK_SUCCESS);
	}

	/* push item to queues */
	for(int64_t i = 0; i < ctx->num_threads; i++) {
		void *item = (items != NULL) ? items[i] : NULL;
		debug("put queue i(%lld), item(%p)", i, item);
		queue_put_wait(ctx->c[i].inq, item);
	}

	/* wait... */
	for(int64_t i = 0; i < ctx->num_threads; i++) {
		void *result;
		queue_get_wait(ctx->c[i].outq, (void **)&result);
		debug("get queue i(%lld), result(%p)", i, result);
		
		/* store the result */
		if(results != NULL) {
			results[i] = result;
		}
	}

	return(PTASK_SUCCESS);
}

/**
 * @struct ptask_stream_status_s
 */
#define CONT 			( 0 )
#define TERM 			( 1 )
struct ptask_stream_status_s {
	int64_t status;
	int64_t cnt;
};

/**
 * @fn ptask_stream
 * @brief get an item from source, throw it to worker, and gather the results into drain.
 */
static _force_inline
int ptask_stream_single(
	struct ptask_context_s *ctx,
	void *(*source)(void *arg), void *source_arg,
	void (*drain)(void *arg, void *result), void *drain_arg)
{
	void *item = NULL;
	while((item = source(source_arg)) != NULL) {
		// log("item(%p)", item);
		void *result = ctx->worker(ctx->worker_arg, item);
		drain(drain_arg, result);
	}
	return(PTASK_SUCCESS);
}
static _force_inline
struct ptask_stream_status_s ptask_stream_bulk_source(
	struct ptask_context_s *ctx,
	void *(*source)(void *arg), void *source_arg,
	int64_t bulk_elems_per_thread)
{
	int64_t rem_size[ctx->num_threads];

	// log("enter ptask_stream_bulk_source");

	/* get number of elems to be pushed to the queue */
	int64_t min_rem_size = bulk_elems_per_thread;
	for(int64_t j = 0; j < ctx->num_threads; j++) {
		rem_size[j] = ctx->queue_size - queue_elements(ctx->c[j].inq);
		min_rem_size = (rem_size[j] < min_rem_size) ? rem_size[j] : min_rem_size;
	}

	int64_t num_items = 0;
	/* read */ {
		/* bulk read */
		int64_t i;
		for(i = 0; i < min_rem_size; i++) {
			for(int64_t j = 0; j < ctx->num_threads; j++) {
				void *item = source(source_arg);
				if(item == NULL) {
					/* source reached the end, calc number of items read */
					num_items += i * ctx->num_threads + j;
					return((struct ptask_stream_status_s){
						.status = TERM,
						.cnt = num_items
					});
				}
				queue_put(ctx->c[j].inq, item);
			}
		}

		/* adjust num_items and rem_size */
		num_items += ctx->num_threads * min_rem_size;
		for(int64_t j = 0; j < ctx->num_threads; j++) {
			rem_size[j] -= min_rem_size;
		}

		/* push remainings */
		for(; i < bulk_elems_per_thread; i++) {
			for(int64_t j = 0; j < ctx->num_threads; j++) {
				/* check queue has the room */
				if(rem_size[j] <= 0) { continue; }
				rem_size[j]--;

				void *item = source(source_arg);
				if(item == NULL) {
					return((struct ptask_stream_status_s){
						.status = TERM,
						.cnt = num_items
					});
				}
				num_items++;
				queue_put(ctx->c[j].inq, item);
			}
		}
	}
	return((struct ptask_stream_status_s){
		.status = CONT,
		.cnt = num_items
	});
}
static _force_inline
struct ptask_stream_status_s ptask_stream_bulk_drain(
	struct ptask_context_s *ctx,
	void (*drain)(void *arg, void *result), void *drain_arg,
	int64_t bulk_elems_per_thread)
{
	int64_t res_cnt[ctx->num_threads];

	// log("enter ptask_stream_bulk_drain");

	/* get number of elems in the result queue */
	int64_t num_items = 0;
	int64_t min_res_cnt = bulk_elems_per_thread;
	for(int64_t j = 0; j < ctx->num_threads; j++) {
		res_cnt[j] = queue_elements(ctx->c[j].outq);
		num_items += res_cnt[j];
		min_res_cnt = (res_cnt[j] < min_res_cnt) ? res_cnt[j] : min_res_cnt;
		debug("j(%lld), res_cnt(%lld)", j, res_cnt[j]);
	}

	/** write */ {
		/* bulk write */
		int64_t i;
		for(i = 0; i < min_res_cnt; i++) {
			for(int64_t j = 0; j < ctx->num_threads; j++) {
				void *result;
				queue_get_wait(ctx->c[j].outq, (void **)&result);
				drain(drain_arg, result);
			}
		}

		/* adjust cnt */
		for(int64_t j = 0; j < ctx->num_threads; j++) {
			res_cnt[j] -= min_res_cnt;
			debug("j(%lld), res_cnt(%lld)", j, res_cnt[j]);
		}

		/* read remainings */
		for(int64_t j = 0; j < ctx->num_threads; j++) {
			for(int64_t k = 0; k < res_cnt[j]; k++) {
				debug("j(%lld), res_cnt(%lld)", j, res_cnt[j]);

				void *result;
				queue_get_wait(ctx->c[j].outq, (void **)&result);
				drain(drain_arg, result);
			}
		}
	}
	return((struct ptask_stream_status_s){
		.status = TERM,
		.cnt = num_items
	});
}
int ptask_stream(
	ptask_t *_ctx,
	void *(*source)(void *arg), void *source_arg,
	void (*drain)(void *arg, void *result), void *drain_arg,
	int64_t bulk_elems)
{
	/* check context sanity */
	struct ptask_context_s *ctx = (struct ptask_context_s *)_ctx;
	if(ctx == NULL || source == NULL || drain == NULL) {
		return(PTASK_ERROR);
	}

	/* branch to single thread mode if num_threads == 0 */
	if(ctx->num_threads == 0) {
		return(ptask_stream_single(ctx,
			source, source_arg,
			drain, drain_arg));
	}

	/* multi thread mode */
	int64_t icnt = 0, ocnt = 0;

	/* restore defaults */
	bulk_elems = (bulk_elems <= 0) ? 512 : bulk_elems;
	int64_t bulk_elems_per_thread = bulk_elems / ctx->num_threads;

	/* scatter and gather */
	// int64_t blk_icnt = ctx->num_threads * bulk_elems_per_thread;
	while(1) {
		/* bulk source */
		struct ptask_stream_status_s s = ptask_stream_bulk_source(ctx,
			source, source_arg,
			bulk_elems_per_thread);
		icnt += s.cnt;
		debug("icnt(%lld), ocnt(%lld)", icnt, ocnt);

		/* check source termination */
		if(s.status != CONT) { break; }

		/* bulk drain */
		ocnt += ptask_stream_bulk_drain(ctx,
			drain, drain_arg,
			bulk_elems_per_thread).cnt;
		debug("icnt(%lld), ocnt(%lld)", icnt, ocnt);
	}
	debug("icnt(%lld), ocnt(%lld)", icnt, ocnt);
	while(ocnt < icnt) {
		// printf("out %lld\n", ocnt);
		for(int64_t j = 0; j < ctx->num_threads; j++) {
			if(queue_elements(ctx->c[j].outq) <= 0) { continue; }

			ocnt++;
			void *result;
			queue_get_wait(ctx->c[j].outq, (void **)&result);
			drain(drain_arg, result);
			debug("icnt(%lld), ocnt(%lld)", icnt, ocnt);
		}
	}

	// printf("%d, %d\n", queue_elements(c.inq), queue_elements(c.outq));

	debug("end");
	return(PTASK_SUCCESS);
}


/* unittests */

/**
 * @fn unittest_init_working_arr
 */
#define UNITTEST_WORKING_ARR_LEN		( 65536 )
#define _ptr(_p)		( (int64_t **)(_p) )
#define _src(_p)		( (int64_t *)((int64_t **)(_p) + UNITTEST_WORKING_ARR_LEN) )
#define _dst(_p)		( (int64_t *)((int64_t **)(_p) + 2 * UNITTEST_WORKING_ARR_LEN) )
static
void *unittest_init_working_arr(
	void *params)
{
	int64_t const len = UNITTEST_WORKING_ARR_LEN;

	/* malloc mem */
	void *p = malloc(len * (2 * sizeof(int64_t) + sizeof(int64_t *)));
	int64_t **ptr = _ptr(p);
	int64_t *src = _src(p);
	int64_t *dst = _dst(p);

	/* init mem */
	for(int64_t i = 0; i < len; i++) {
		ptr[i] = &src[i];
		src[i] = (int64_t)params;
		dst[i] = 0xcafe;
	}
	return((void *)p);
}

/**
 * @fn unittest_check_working_arr
 */
static
int unittest_check_working_arr(
	void *p,
	int64_t num,
	int64_t check_len)
{
	int64_t **ptr = _ptr(p);
	int64_t *src = _src(p);
	int64_t *dst = _dst(p);
	for(int64_t i = 0; i < check_len; i++) {
		if(ptr[i] != &src[i]) {
			debug("ptr check failed i(%lld), ptr[i](%p), &src[i](%p)", i, ptr[i], &src[i]);
			return(0);
		}
		if(src[i] != 0) {
			debug("src check failed i(%lld), src[i](%lld), num(%lld)", i, src[i], num);
			return(0);
		}
		if(dst[i] != num) {
			debug("dst check failed i(%lld), dst[i](%lld), num(%lld)", i, dst[i], num);
			return(0);
		}
	}
	return(1);
}

/**
 * @macro format_array
 */
#define format_array(_p) ({ \
	int64_t *dst = _dst(_p); \
	char *base = alloca(8 * UNITTEST_WORKING_ARR_LEN); \
	char *p = base; \
	for(int64_t i = 0; i < UNITTEST_WORKING_ARR_LEN; i++) { \
		p += sprintf(p, "%" PRId64 ", ", dst[i]); \
	} \
	p[-2] = '\0'; \
	base; \
})

unittest_config(
	.name = "ptask",
	.params = (void *)0,
	.init = unittest_init_working_arr,
	.clean = free
);

/**
 * @fn unittest_init_const_arr
 */
#define UNITTEST_CONST_ARR_LEN		( 32 )
static
void *unittest_init_const_arr(
	void *params)
{
	int64_t const len = UNITTEST_CONST_ARR_LEN;

	/* malloc mem */
	int64_t *arr = (int64_t *)malloc(len * sizeof(int64_t));

	/* init mem */
	for(int64_t i = 0; i < len; i++) {
		arr[i] = (int64_t)params;
	}
	return((void *)arr);
}

/**
 * @macro with_arr
 */
#define with_arr(_num) \
	.params = (void *)(_num), \
	.init = unittest_init_const_arr, \
	.clean = free

/**
 * @fn unittest_source
 * @brief take a pointer to an array generated by unittest_init_arr,
 * returns a pointer to an element in the array
 */
static int64_t unittest_source_cnt = 0;
static
void *unittest_source(
	void *arg)
{
	/* arg must be a pointer to an array */
	debug("source_cnt(%lld)", unittest_source_cnt);
	if(unittest_source_cnt >= UNITTEST_WORKING_ARR_LEN) {
		return(NULL);
	}
	int64_t *arr = _src(arg);
	return(&arr[unittest_source_cnt++]);
}

/**
 * @fn unittest_worker
 * @brief take a constant (arg) and pointer to int64_t (item), return the sum of arg and item.
 */
static
void *unittest_worker(
	void *arg,
	void *item)
{
	int64_t *val = (int64_t *)item;
	int64_t add = (int64_t)arg;

	debug("worker val(%p), add(%lld), *val(%lld), res(%lld)",
		val, add, *val, *val + add);
	return((void *)(*val + add));
}

/**
 * @fn unittest_drain
 * @brief take a pointer to an array generated by unittest_init_arr,
 * push a result to the array.
 */
static int64_t unittest_drain_cnt = 0;
static
void unittest_drain(
	void *arg,
	void *result)
{
	debug("drain_cnt(%lld)", unittest_drain_cnt);
	if(unittest_drain_cnt >= UNITTEST_WORKING_ARR_LEN) {
		return;
	}
	int64_t *dst = _dst(arg);
	dst[unittest_drain_cnt++] = (uint64_t)result;
	return;
}

/**
 * @fn unittest_init_stream
 */
static
void unittest_init_stream(
	void)
{
	unittest_source_cnt = 0;
	unittest_drain_cnt = 0;
	return;
}

/* check if ptask_init returns valid context */
unittest()
{
	/* single-thread, default queue size */
	struct ptask_context_s *p = ptask_init(
		unittest_worker, gctx, 0, 0);
	assert(p != NULL);
	ptask_clean(p);

	/* 4-threads, default queue size */
	p = ptask_init(unittest_worker, gctx, 4, 0);
	assert(p != NULL);
	ptask_clean(p);

	/* single-thread, queue size == 2048 */
	p = ptask_init(unittest_worker, gctx, 0, 2048);
	assert(p != NULL);
	ptask_clean(p);

	/* 4-threads, queue size == 2048 */
	p = ptask_init(unittest_worker, gctx, 4, 2048);
	assert(p != NULL);
	ptask_clean(p);
}

unittest(with_arr(100))
{
	struct ptask_context_s *p = ptask_init(
		unittest_worker, ctx, 0, 0);

	ptask_parallel(p, gctx, (void *)_dst(gctx));
	assert(unittest_check_working_arr(gctx, 100, 1), "%s", format_array(gctx));

	ptask_clean(p);
}

unittest(with_arr(100))
{
	struct ptask_context_s *p = ptask_init(
		unittest_worker, ctx, 4, 0);

	ptask_parallel(p, gctx, (void *)_dst(gctx));
	assert(unittest_check_working_arr(gctx, 100, 4), "%s", format_array(gctx));

	ptask_clean(p);
}

unittest(with_arr(100))
{
	struct ptask_context_s *p = ptask_init(
		unittest_worker, ctx, 0, 2048);

	ptask_parallel(p, gctx, (void *)_dst(gctx));
	assert(unittest_check_working_arr(gctx, 100, 1), "%s", format_array(gctx));

	ptask_clean(p);
}

unittest(with_arr(100))
{
	struct ptask_context_s *p = ptask_init(
		unittest_worker, ctx, 4, 2048);

	ptask_parallel(p, gctx, (void *)_dst(gctx));
	assert(unittest_check_working_arr(gctx, 100, 4), "%s", format_array(gctx));

	ptask_clean(p);
}

unittest(with_arr(100))
{
	struct ptask_context_s *p = ptask_init(
		unittest_worker, ctx, 0, 0);
	unittest_init_stream();

	ptask_stream(p,
		unittest_source, gctx,
		unittest_drain, gctx,
		0);
	assert(unittest_check_working_arr(gctx, 100, UNITTEST_WORKING_ARR_LEN),
		"%s", format_array(gctx));

	ptask_clean(p);
}

unittest(with_arr(100))
{
	struct ptask_context_s *p = ptask_init(
		unittest_worker, ctx, 4, 0);
	unittest_init_stream();

	ptask_stream(p,
		unittest_source, gctx,
		unittest_drain, gctx,
		0);
	assert(unittest_check_working_arr(gctx, 100, UNITTEST_WORKING_ARR_LEN),
		"%s", format_array(gctx));

	ptask_clean(p);
}

unittest(with_arr(100))
{
	struct ptask_context_s *p = ptask_init(
		unittest_worker, ctx, 0, 16);
	unittest_init_stream();

	ptask_stream(p,
		unittest_source, gctx,
		unittest_drain, gctx,
		0);
	assert(unittest_check_working_arr(gctx, 100, UNITTEST_WORKING_ARR_LEN),
		"%s", format_array(gctx));

	ptask_clean(p);
}

unittest(with_arr(100))
{
	struct ptask_context_s *p = ptask_init(
		unittest_worker, ctx, 4, 16);
	unittest_init_stream();

	ptask_stream(p,
		unittest_source, gctx,
		unittest_drain, gctx,
		0);
	assert(unittest_check_working_arr(gctx, 100, UNITTEST_WORKING_ARR_LEN),
		"%s", format_array(gctx));

	ptask_clean(p);
}


unittest(with_arr(100))
{
	struct ptask_context_s *p = ptask_init(
		unittest_worker, ctx, 0, 2048);
	unittest_init_stream();

	ptask_stream(p,
		unittest_source, gctx,
		unittest_drain, gctx,
		0);
	assert(unittest_check_working_arr(gctx, 100, UNITTEST_WORKING_ARR_LEN),
		"%s", format_array(gctx));

	ptask_clean(p);
}

unittest(with_arr(100))
{
	struct ptask_context_s *p = ptask_init(
		unittest_worker, ctx, 4, 2048);
	unittest_init_stream();

	ptask_stream(p,
		unittest_source, gctx,
		unittest_drain, gctx,
		0);
	assert(unittest_check_working_arr(gctx, 100, UNITTEST_WORKING_ARR_LEN),
		"%s", format_array(gctx));

	ptask_clean(p);
}

/**
 * end of ptask.c
 */
