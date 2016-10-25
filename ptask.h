
/**
 * @file ptask.h
 *
 * @brief parallel task manager
 *
 * @author Hajime Suzuki
 * @date 2016/3/22
 * @license MIT
 */
#ifndef _PTASK_H_INCLUDED
#define _PTASK_H_INCLUDED

#include <stdint.h>

/**
 * @enum ptask_status
 */
enum ptask_status {
	PTASK_SUCCESS	= 0,
	PTASK_ERROR		= 1
};

/**
 * @type ptask_t
 * @brief context container
 */
typedef struct ptask_context_s ptask_t;

/**
 * @fn ptask_init
 * @brief initialize threads
 */
ptask_t *ptask_init(
	void *(*worker)(void *arg, void *item),
	void *worker_arg[],
	int64_t num_threads,
	int64_t queue_size);

/**
 * @fn ptask_clean
 */
void ptask_clean(
	ptask_t *ctx);

/**
 * @fn ptask_parallel
 */
int ptask_parallel(
	ptask_t *ctx,
	void *items[],
	void *results[]);

/**
 * @fn ptask_stream
 * @brief get an item from source, throw it to worker, and gather the results into drain.
 */
int ptask_stream(
	ptask_t *ctx,
	void *(*source)(void *arg),
	void *source_arg,
	void (*drain)(void *arg, void *result),
	void *drain_arg,
	int64_t bulk_elems);

#endif /* _PTASK_H_INCLUDED */
/**
 * end of ptask.h
 */
