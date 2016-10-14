
/**
 * @file psort.c
 *
 * @brief parallel sort library frontend implementation
 *
 * @author Hajime Suzuki
 * @date 2016/3/20
 * @license MIT
 */

/* import unittest */
#define UNITTEST_UNIQUE_ID			200
#define UNITTEST 					1

#include  "unittest.h"

#include <stdint.h>
#include <string.h>
#include "arch/arch.h"
#include "psort.h"
#include "ptask.h"		/* pthread parallel task execution library */
#include "log.h"
#include "sassert.h"


/* constants */
#define WCR_OCC_SIZE				( 1<<8 )	/** 8bit */

/**
 * @macro _likely, _unlikely
 * @brief branch prediction hint for gcc-compatible compilers
 */
#define _likely(x)		__builtin_expect(!!(x), 1)
#define _unlikely(x)	__builtin_expect(!!(x), 0)

/**
 * @macro _force_inline
 * @brief inline directive for gcc-compatible compilers
 */
#define _force_inline	inline
// #define _force_inline

/* max / min */
#define MAX2(x, y)					( (x) < (y) ? (y) : (x) )
#define MIN2(x, y)					( (x) > (y) ? (y) : (x) )

/**
 * @struct psort_occ_s
 */
struct psort_occ_s {
	uint64_t occ[WCR_OCC_SIZE];
};

/**
 * @struct psort_buffer_counter_s
 */
struct psort_buffer_counter_s {
	uint8_t cnt[WCR_OCC_SIZE];
};

/**
 * @struct psort_buffer_s
 */
struct psort_buffer_s {
	uint8_t buf[WCR_OCC_SIZE][WCR_BUF_SIZE];
};

/**
 * @struct psort_thread_context_s
 */
struct psort_thread_context_s {
	// struct psort_occ_s *occ;
	// struct psort_buffer_counter_s *cnt;

	/* aligned on 64byte boundary */
	uint64_t occ[WCR_OCC_SIZE];			/* (2048) */
	uint8_t cnt[WCR_OCC_SIZE];			/* (256) */

	/* aligned on 64byte boundary (2304) */
	uint8_t buf[WCR_OCC_SIZE][WCR_BUF_SIZE];/* (8192) */

	/* aligned on 64byte boundary (10496) */
	uint32_t digit;						/* (4) */
	uint32_t num_threads;				/* (4) */
	void *src;							/* (8) */
	void *dst;							/* (8) */
	uint64_t from, to;					/* (16) */
	uint64_t _pad[3];					/* (24) */

	/* aligned on 64byte boundary (10560) */
};
_static_assert(sizeof(struct psort_thread_context_s) == 2368 + WCR_OCC_SIZE * WCR_BUF_SIZE);

/**
 * @fn aligned_malloc
 *
 * @brief an wrapper of posix_memalign function
 */
static _force_inline
void *aligned_malloc(size_t size, size_t align)
{
	void *ptr;
	if(posix_memalign(&ptr, align, size) != 0) {
		return(NULL);
	}
	return(ptr);
}

/**
 * @fn psort_dispatcher
 */
static
void *psort_dispatcher(
	void *arg,
	void *item)
{
	debug("arg(%p), item(%p)", arg, item);
	((void (*)(struct psort_thread_context_s *))item)(
		(struct psort_thread_context_s *)arg);
	return(NULL);
}

/* instanciate radixsort and quicksort */
#define _join(a, b)			a##b
#define join(a, b)			_join(a, b)

#undef UNITTEST_UNIQUE_ID

/* 16bit */
#define elem_t				uint16_t
#define SUFFIX				16
#define rd 					_rd
#define wr 					_wr
#define ex 					_ex
#define p 					_p
#define e 					_e
#define UNITTEST_UNIQUE_ID	SUFFIX
#include "psort_radix_internal.c"
// #include "psort_quick_intl.c"

/* 32bit */
#define elem_t				uint32_t
#define SUFFIX				32
#define rd 					_rd
#define wr 					_wr
#define ex 					_ex
#define p 					_p
#define e 					_e
#define UNITTEST_UNIQUE_ID	SUFFIX
#include "psort_radix_internal.c"
// #include "psort_quick_intl.c"

/* 64bit */
#define elem_t				uint64_t
#define SUFFIX				64
#define rd 					_rd
#define wr 					_wr
#define ex 					_ex
#define p 					_p
#define e 					_e
#define UNITTEST_UNIQUE_ID	SUFFIX
#include "psort_radix_internal.c"
// #include "psort_quick_intl.c"

/* 128bit */
#define elem_t				elem_128_t
#define SUFFIX				128
#define rd 					rd_128
#define wr 					wr_128
#define ex 					ex_128
#define p 					p_128
#define e 					e_128
#define UNITTEST_UNIQUE_ID	SUFFIX
#include "psort_radix_internal.c"
// #include "psort_quick_intl.c"

/**
 * @fn psort_full
 */
int psort_full(
	void *arr,
	uint64_t len,
	uint64_t elem_size,
	uint64_t num_threads)
{
	switch(elem_size) {
		case 2: psort_partialsort_parallel_16(arr, len, num_threads, 0, 2); return(0);
		case 4: psort_partialsort_parallel_32(arr, len, num_threads, 0, 4); return(0);
		case 8: psort_partialsort_parallel_64(arr, len, num_threads, 0, 8); return(0);
		case 16: psort_partialsort_parallel_128(arr, len, num_threads, 0, 16); return(0);
		default: return(-1);
	}
	return(-1);
}

/**
 * @fn psort_half
 */
int psort_half(
	void *arr,
	uint64_t len,
	uint64_t elem_size,
	uint64_t num_threads)
{
	switch(elem_size) {
		case 2: psort_partialsort_parallel_16(arr, len, num_threads, 0, 1); return(0);
		case 4: psort_partialsort_parallel_32(arr, len, num_threads, 0, 2); return(0);
		case 8: psort_partialsort_parallel_64(arr, len, num_threads, 0, 4); return(0);
		case 16: psort_partialsort_parallel_128(arr, len, num_threads, 0, 8); return(0);
		default: return(-1);
	}
	return(-1);
}

/**
 * @fn psort_partial
 */
int psort_partial(
	void *arr,
	uint64_t len,
	uint64_t elem_size,
	uint64_t num_threads,
	uint64_t from,
	uint64_t to)
{
	switch(elem_size) {
		case 2: psort_partialsort_parallel_16(arr, len, num_threads, from, to); return(0);
		case 4: psort_partialsort_parallel_32(arr, len, num_threads, from, to); return(0);
		case 8: psort_partialsort_parallel_64(arr, len, num_threads, from, to); return(0);
		case 16: psort_partialsort_parallel_128(arr, len, num_threads, from, to); return(0);
		default: return(-1);
	}
	return(-1);
}

/* unittest */
#include <time.h>

#define UNITTEST_UNIQUE_ID 		61
unittest_config(
	.name = "psort",
	.depends_on = { "psort_radix_internal" }
);

#define UNITTEST_ARR_LEN		10000

/* srand */
unittest()
{
	srand(time(NULL));
}

/* full sort 16bit */
unittest()
{
	/* init */
	uint16_t *arr = (uint16_t *)malloc(sizeof(uint16_t) * UNITTEST_ARR_LEN);
	for(int64_t i = 0; i < UNITTEST_ARR_LEN; i++) {
		arr[i] = rand() % UINT16_MAX;
	}

	/* sort */
	psort_full(arr, UNITTEST_ARR_LEN, sizeof(uint16_t), 4);

	/* check */
	for(int64_t i = 1; i < UNITTEST_ARR_LEN; i++) {
		assert(arr[i - 1] <= arr[i], "%lld, %d, %d", i, arr[i - 1], arr[i]);
	}
	free(arr);
}

/* full sort 32bit */
unittest()
{
	/* init */
	uint32_t *arr = (uint32_t *)malloc(sizeof(uint32_t) * UNITTEST_ARR_LEN);
	for(int64_t i = 0; i < UNITTEST_ARR_LEN; i++) {
		arr[i] = rand() % UINT32_MAX;
	}

	/* sort */
	psort_full(arr, UNITTEST_ARR_LEN, sizeof(uint32_t), 4);

	/* check */
	for(int64_t i = 1; i < UNITTEST_ARR_LEN; i++) {
		assert(arr[i - 1] <= arr[i], "%lld, %d, %d", i, arr[i - 1], arr[i]);
	}
	free(arr);
}

/* full sort 64bit */
unittest()
{
	/* init */
	uint64_t *arr = (uint64_t *)malloc(sizeof(uint64_t) * UNITTEST_ARR_LEN);
	for(int64_t i = 0; i < UNITTEST_ARR_LEN; i++) {
		arr[i] = (uint64_t)rand() + ((uint64_t)rand()<<32);
	}

	/* sort */
	psort_full(arr, UNITTEST_ARR_LEN, sizeof(uint64_t), 4);

	/* check */
	for(int64_t i = 1; i < UNITTEST_ARR_LEN; i++) {
		assert(arr[i - 1] <= arr[i], "%lld, %lld, %lld", i, arr[i - 1], arr[i]);
	}
	free(arr);
}

/* half sort 16bit */
unittest()
{
	/* init */
	uint16_t *arr = (uint16_t *)malloc(sizeof(uint16_t) * UNITTEST_ARR_LEN);
	for(int64_t i = 0; i < UNITTEST_ARR_LEN; i++) {
		arr[i] = rand() % UINT16_MAX;
	}

	/* sort */
	psort_half(arr, UNITTEST_ARR_LEN, sizeof(uint16_t), 4);

	/* check */
	for(int64_t i = 1; i < UNITTEST_ARR_LEN; i++) {
		assert((0xff & arr[i - 1]) <= (0xff & arr[i]),
			"%lld, %d, %d", i, 0xff & arr[i - 1], 0xff & arr[i]);
	}
	free(arr);
}

/* half sort 32bit */
unittest()
{
	/* init */
	uint32_t *arr = (uint32_t *)malloc(sizeof(uint32_t) * UNITTEST_ARR_LEN);
	for(int64_t i = 0; i < UNITTEST_ARR_LEN; i++) {
		arr[i] = rand() % UINT32_MAX;
	}

	/* sort */
	psort_half(arr, UNITTEST_ARR_LEN, sizeof(uint32_t), 4);

	/* check */
	for(int64_t i = 1; i < UNITTEST_ARR_LEN; i++) {
		assert((0xffff & arr[i - 1]) <= (0xffff & arr[i]),
			"%lld, %d, %d", i, 0xffff & arr[i - 1], 0xffff & arr[i]);
	}
	free(arr);
}

/* half sort 64bit */
unittest()
{
	/* init */
	uint64_t *arr = (uint64_t *)malloc(sizeof(uint64_t) * UNITTEST_ARR_LEN);
	for(int64_t i = 0; i < UNITTEST_ARR_LEN; i++) {
		arr[i] = (uint64_t)rand() + ((uint64_t)rand()<<32);
	}

	/* sort */
	psort_half(arr, UNITTEST_ARR_LEN, sizeof(uint64_t), 4);

	/* check */
	for(int64_t i = 1; i < UNITTEST_ARR_LEN; i++) {
		assert((0xffffffff & arr[i - 1]) <= (0xffffffff & arr[i]),
			"%lld, %lld, %lld", i, 0xffffffff & arr[i - 1], 0xffffffff & arr[i]);
	}
	free(arr);
}

/* partial sort 16bit */
unittest()
{
	/* init */
	uint16_t *arr = (uint16_t *)malloc(sizeof(uint16_t) * UNITTEST_ARR_LEN);
	for(int64_t i = 0; i < UNITTEST_ARR_LEN; i++) {
		arr[i] = rand() % UINT16_MAX;
	}

	/* sort */
	psort_partial(arr, UNITTEST_ARR_LEN, sizeof(uint16_t), 4, 1, 2);

	/* check */
	for(int64_t i = 1; i < UNITTEST_ARR_LEN; i++) {
		assert((0xff00 & arr[i - 1]) <= (0xff00 & arr[i]),
			"%lld, %d, %d", i, 0xff00 & arr[i - 1], 0xff00 & arr[i]);
	}
	free(arr);
}

/* partial sort 32bit */
unittest()
{
	/* init */
	uint32_t *arr = (uint32_t *)malloc(sizeof(uint32_t) * UNITTEST_ARR_LEN);
	for(int64_t i = 0; i < UNITTEST_ARR_LEN; i++) {
		arr[i] = rand() % UINT32_MAX;
	}

	/* sort */
	psort_partial(arr, UNITTEST_ARR_LEN, sizeof(uint32_t), 4, 2, 4);

	/* check */
	for(int64_t i = 1; i < UNITTEST_ARR_LEN; i++) {
		assert((0xffff0000 & arr[i - 1]) <= (0xffff0000 & arr[i]),
			"%lld, %d, %d", i, 0xffff0000 & arr[i - 1], 0xffff0000 & arr[i]);
	}
	free(arr);
}

/* partial sort 64bit */
unittest()
{
	/* init */
	uint64_t *arr = (uint64_t *)malloc(sizeof(uint64_t) * UNITTEST_ARR_LEN);
	for(int64_t i = 0; i < UNITTEST_ARR_LEN; i++) {
		arr[i] = (uint64_t)rand() + ((uint64_t)rand()<<32);
	}

	/* sort */
	psort_partial(arr, UNITTEST_ARR_LEN, sizeof(uint64_t), 4, 4, 8);

	/* check */
	uint64_t const mask = 0xffffffff00000000;
	for(int64_t i = 1; i < UNITTEST_ARR_LEN; i++) {
		assert((mask & arr[i - 1]) <= (mask & arr[i]),
			"%lld, %lld, %lld", i, mask & arr[i - 1], mask & arr[i]);
	}
	free(arr);
}

/**
 * end of psort.c
 */
