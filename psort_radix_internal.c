
/**
 * @file psort_radix_intl.c
 *
 * @brief write-combining radix sort implementation. included from psort.c
 *
 * @author Hajime Suzuki
 * @date 2016/3/20
 * @license MIT
 */

/* check if element type and read / write macros are defined */
#if !defined(elem_t) || !defined(rd) || !defined(wr) || !defined(ex)
#  error "elem_t, rd, wr, and ex macros must be defined."
#endif

/* constant */
#define WCR_BUF_ELEM_COUNT		( WCR_BUF_SIZE / sizeof(elem_t) )

/**
 * @fn psort_count_occ
 */
static
void join(psort_count_occ_, SUFFIX)(
	struct psort_thread_context_s *ctx)
{
	/* initialize occ and counter array */
	uint64_t *occ = (uint64_t *)ctx->occ;
	memset(occ, 0, (sizeof(uint64_t) + sizeof(uint8_t)) * WCR_OCC_SIZE);	

	/* count occurrences */
	elem_t *sp = (elem_t *)ctx->src + ctx->from;
	elem_t *ep = (elem_t *)ctx->src + ctx->to;
	uint64_t digit = ctx->digit;
	for(elem_t *p = sp; p < ep; p++) {
		occ[ex(rd(p), digit)]++;
	}
	return;
}

/**
 * @fn psort_gather_occ
 * @brief gather occurrence table. item is ignored.
 */
static
void join(psort_gather_occ_, SUFFIX)(
	struct psort_thread_context_s *ctx)
{
	uint64_t threads = ctx->num_threads;
	elem_t *ptr = ctx->dst;
	for(uint64_t i = 0; i < WCR_OCC_SIZE; i++) {
		for(uint64_t j = 0; j < threads; j++) {
			/* store base pointer to occ[j].occ[i] */
			uint64_t curr_occ = ctx[j].occ[i];
			ctx[j].occ[i] = (uint64_t)ptr;
			ptr += curr_occ;
		}
	}
	return;
}

/**
 * @fn psort_scatter
 */
static
void join(psort_scatter_, SUFFIX)(
	struct psort_thread_context_s *ctx)
{
	/* extract pointers */
	elem_t **base = (elem_t **)ctx->occ;
	uint8_t *cnt = (uint8_t *)ctx->cnt;
	elem_t (*buf)[WCR_BUF_ELEM_COUNT] = (elem_t (*)[WCR_BUF_ELEM_COUNT])ctx->buf;

	/* extract elem pointer */
	elem_t *sp = (elem_t *)ctx->src + ctx->from;
	elem_t *ep = (elem_t *)ctx->src + ctx->to;

	/* scatter */
	uint64_t digit = ctx->digit;
	for(elem_t *p = sp; p < ep; p++) {
		/* load an element, store to the buffer */
		register elem_t e = rd(p);
		register uint64_t n = ex(e, digit);
		register elem_t *q = base[n];

		if(((uint64_t)q & (WCR_BUF_SIZE - 1)) != 0) {
			/* p is not aligned to cache boundary, write directly */
			wr(q, e);
			base[n] = q + 1;
			continue;
		}

		/* read counter */
		register uint8_t c = cnt[n];
		wr(&buf[n][c], e);

		/** check if flush is needed */
		if(++c == WCR_BUF_ELEM_COUNT) {
			/* bulk copy and update pointer */
			memcpy_buf(q, buf[n]);
			base[n] = q + WCR_BUF_ELEM_COUNT;
			c = 0;
		}

		/* write back cnt */
		cnt[n] = c;
	}
	return;
}

/**
 * @fn psort_flush
 */
static
void join(psort_flush_, SUFFIX)(
	struct psort_thread_context_s *ctx)
{
	/** flush the remaining content */
	uint64_t threads = ctx->num_threads;
	for(uint64_t i = 0; i < WCR_OCC_SIZE; i++) {
		for(uint64_t j = 0; j < threads; j++) {
			if(ctx[j].cnt[i] == 0) { continue; }
			elem_t *dp = (elem_t *)ctx[j].occ[i];
			elem_t *sp = (elem_t *)ctx[j].buf[i];
			for(uint64_t k = 0; k < ctx[j].cnt[i]; k++) {
				wr(&dp[k], rd(&sp[k]));
			}
		}
	}
	return;
}

/**
 * @fn psort_copyback
 */
static
void join(psort_copyback_, SUFFIX)(
	struct psort_thread_context_s *ctx)
{
	/* extract elem pointer */
	elem_t *src = (elem_t *)ctx->src;
	elem_t *dst = (elem_t *)ctx->dst;
	uint64_t from = ctx->from;
	uint64_t to = ctx->to;

	/* copy back */
	memcpy(dst + from, src + from, sizeof(elem_t) * (to - from));
	return;
}

/**
 * @fn psort_partialsort_parallel
 */
static
void join(psort_partialsort_parallel_, SUFFIX)(
	void *src,
	uint64_t len,
	uint64_t num_threads,
	uint64_t lower_digit,
	uint64_t higher_digit)
{
	uint64_t nt = (num_threads == 0) ? 1 : num_threads;

	/* malloc buffer */
	void *ptr = aligned_malloc(
		  nt * (
		  	  sizeof(struct psort_thread_context_s *)
		  	+ 3 * sizeof(void *)		/* function pointers */
			+ sizeof(struct psort_thread_context_s))
		+ sizeof(elem_t) * len,
		32);

	/* array of pointers to thread contexts */
	struct psort_thread_context_s **pth = (struct psort_thread_context_s **)ptr;
	
	/* array of pointers to functions */
	void **count_occ = (void **)&pth[nt];
	void **scatter = (void **)&count_occ[nt];
	void **copyback = (void **)&scatter[nt];

	/* thread contexts and working buffers */
	struct psort_thread_context_s *th = (struct psort_thread_context_s *)&copyback[nt];
	void *dst = (void *)&th[nt];

	/* initialize thread contexts */
	for(uint64_t i = 0; i < nt; i++) {
		/* pointer to thread context */
		pth[i] = &th[i];

		/* pointer to functions */
		count_occ[i] = (void *)join(psort_count_occ_, SUFFIX);
		scatter[i] = (void *)join(psort_scatter_, SUFFIX);
		copyback[i] = (void *)join(psort_copyback_, SUFFIX);

		/* num_threads */
		th[i].num_threads = nt;

		/* array */
		th[i].from = i * len / nt;
		th[i].to = (i + 1) * len / nt;
	}

	/* initialize ptask object */
	ptask_t *pt = ptask_init(psort_dispatcher, (void **)pth, num_threads, 1024);

	/* LSB first radixsort */
	for(uint64_t i = lower_digit; i < higher_digit; i++) {
		/* set digit and pointers */
		for(uint64_t j = 0; j < nt; j++) {
			th[j].digit = i;
			th[j].src = src;
			th[j].dst = dst;
		}

		/* count occ */
		ptask_parallel(pt, count_occ, NULL);

		/* gather occ */
		join(psort_gather_occ_, SUFFIX)(th);

		/* scatter */
		ptask_parallel(pt, scatter, NULL);

		/* flush */
		join(psort_flush_, SUFFIX)(th);

		/* swap */
		void *tmp = src; src = dst; dst = tmp;
	}

	/* copyback */
	if((higher_digit - lower_digit) & 0x01) {
		for(uint64_t j = 0; j < nt; j++) {
			th[j].src = src;
			th[j].dst = dst;
		}
		ptask_parallel(pt, copyback, NULL);
	}

	/* cleanup ptask object */
	ptask_clean(pt);

	/* cleanup working memory */
	free(ptr);
	return;
}

/* unittests */
unittest_config(
	.name = "psort_radix_intl",
	.depends_on = { "ptask" }
);

/* small integer, single thread */
unittest()
{
	uint64_t raw[] =    { 1, 0, 2, 1, 0, 2, 0, 0, 1, 1 };
	uint64_t sorted[] = { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2 };

	elem_t *arr = (elem_t *)aligned_malloc(sizeof(elem_t) * 10, 16);
	for(uint64_t i = 0; i < 10; i++) {
		wr(arr + i, p(raw[i]));
	}

	/* sort */
	join(psort_partialsort_parallel_, SUFFIX)(
		arr, 10, 0, 0, sizeof(elem_t));

	/* check */
	for(uint64_t i = 0; i < 10; i++) {
		assert(e(rd(arr + i)) == sorted[i],
			"%llu, %llu", e(rd(arr + i)), (uint64_t)sorted[i]);
	}

	/* cleanup */
	free(arr);
}

/* small integer, 4-thread */
unittest()
{
	uint64_t raw[] =    { 1, 0, 2, 1, 0, 2, 0, 0, 1, 1 };
	uint64_t sorted[] = { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2 };

	elem_t *arr = (elem_t *)aligned_malloc(sizeof(elem_t) * 10, 16);
	for(uint64_t i = 0; i < 10; i++) {
		wr(arr + i, p(raw[i]));
	}

	/* sort */
	join(psort_partialsort_parallel_, SUFFIX)(
		arr, 10, 4, 0, sizeof(elem_t));

	/* check */
	for(uint64_t i = 0; i < 10; i++) {
		assert(e(rd(arr + i)) == sorted[i],
			"%llu, %llu", e(rd(arr + i)), (uint64_t)sorted[i]);
	}

	/* cleanup */
	free(arr);
}

/* middle integer, single thread */
unittest()
{
	uint64_t raw[] =    { 1000, 0, 2000, 1000, 0, 2000, 0, 0, 1000, 1000 };
	uint64_t sorted[] = { 0, 0, 0, 0, 1000, 1000, 1000, 1000, 2000, 2000 };

	elem_t *arr = (elem_t *)aligned_malloc(sizeof(elem_t) * 10, 16);
	for(uint64_t i = 0; i < 10; i++) {
		wr(arr + i, p(raw[i]));
	}

	/* sort */
	join(psort_partialsort_parallel_, SUFFIX)(
		arr, 10, 0, 0, sizeof(elem_t));

	/* check */
	for(uint64_t i = 0; i < 10; i++) {
		assert(e(rd(arr + i)) == sorted[i],
			"%llu, %llu", e(rd(arr + i)), (uint64_t)sorted[i]);
	}

	/* cleanup */
	free(arr);
}

/* middle integer, 4-thread */
unittest()
{
	uint64_t raw[] =    { 1000, 0, 2000, 1000, 0, 2000, 0, 0, 1000, 1000 };
	uint64_t sorted[] = { 0, 0, 0, 0, 1000, 1000, 1000, 1000, 2000, 2000 };

	elem_t *arr = (elem_t *)aligned_malloc(sizeof(elem_t) * 10, 16);
	for(uint64_t i = 0; i < 10; i++) {
		wr(arr + i, p(raw[i]));
	}

	/* sort */
	join(psort_partialsort_parallel_, SUFFIX)(
		arr, 10, 4, 0, sizeof(elem_t));

	/* check */
	for(uint64_t i = 0; i < 10; i++) {
		assert(e(rd(arr + i)) == sorted[i],
			"%llu, %llu", e(rd(arr + i)), (uint64_t)sorted[i]);
	}

	/* cleanup */
	free(arr);
}

/* inverse, long array, single thread */
unittest()
{
	uint64_t const len = 10000;
	elem_t *arr = (elem_t *)aligned_malloc(
		sizeof(elem_t) * len,
		16);

	for(uint64_t i = 0; i < len; i++) {
		wr(arr + i, p((len - 1) - i));
	}

	/* sort */
	join(psort_partialsort_parallel_, SUFFIX)(
		arr, len, 0, 0, sizeof(elem_t));

	/* check */
	for(uint64_t i = 0; i < len; i++) {
		assert(e(rd(arr + i)) == i, "%llu", e(rd(arr + i)));
	}

	/* cleanup */
	free(arr);
	return;
}

/* inverse, long array, 4-thread */
unittest()
{
	uint64_t const len = 10000;
	elem_t *arr = (elem_t *)aligned_malloc(
		sizeof(elem_t) * len,
		16);

	for(uint64_t i = 0; i < len; i++) {
		wr(arr + i, p((len - 1) - i));
	}

	/* sort */
	join(psort_partialsort_parallel_, SUFFIX)(
		arr, len, 4, 0, sizeof(elem_t));

	/* check */
	for(uint64_t i = 0; i < len; i++) {
		assert(e(rd(arr + i)) == i, "%llu", e(rd(arr + i)));
	}

	/* cleanup */
	free(arr);
	return;
}
#if 0
/* benchmark */
#include <sys/time.h>
unittest()
{
	uint64_t const len = 200000000;
	elem_t *arr = (elem_t *)aligned_malloc(
		sizeof(elem_t) * len,
		16);

	for(uint64_t i = 0; i < 5; i++) {
		/* init array */
		for(uint64_t j = 0; j < len; j++) {
			wr(arr + j, p((len - 1) - j));
		}

		struct timeval ts, te;
		gettimeofday(&ts, NULL);

		/* sort */
		join(psort_partialsort_parallel_, SUFFIX)(
			arr, len, 4, 0, sizeof(elem_t));

		gettimeofday(&te, NULL);

		/* check */
		for(uint64_t j = 1; j < len; j++) {
			assert(e(rd(arr + j - 1)) <= e(rd(arr + j)),
				"%llu, %llu",
				e(rd(arr + j - 1)), e(rd(arr + j)));
		}

		fprintf(stderr, "%lu us\n",
			(te.tv_sec - ts.tv_sec) * 1000000 + (te.tv_usec - ts.tv_usec));
	}

	/* cleanup */
	free(arr);
	return;
}
#endif

/* cleanup macros */
#undef WCR_BUF_ELEM_COUNT
#undef elem_t
#undef SUFFIX
#undef rd
#undef wr
#undef ex
#undef p
#undef e
#undef UNITTEST_UNIQUE_ID

/**
 * end of psort_radix_intl.c
 */
