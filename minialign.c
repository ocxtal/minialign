
/**
 * @file minialign.c
 *
 * @brief minialign core implementation
 *
 * @author Hajime Suzuki (original files by Heng Li)
 * @license MIT
 */
// #define DEBUG
/* configurations */
/**
 * @macro MM_VERSION
 */
#ifndef MM_VERSION
#  define MM_VERSION		"minialign-0.6.0-devel"
#endif

/**
 * @macro STRICT_STREAM_ORDERING
 * @brief set non-zero value to ensure the order of output records
 */
#define STRICT_STREAM_ORDERING		( 1 )

/**
 * @macro COLLECT_SUPPLEMENTARY
 * @brief collect suppementary alignments, set zero for compatibility with 0.4.x
 */
#define COLLECT_SUPPLEMENTARY		( 1 )

/**
 * @macro USE_CRC32_HASH
 * @brief use crc32 for hash64, set zero for compatibility with 0.4.x
 */
#define USE_CRC32_HASH				( 1 )

/**
 * @macro MAX_THREADS
 * @brief max #threads, set larger value if needed
 */
#define MAX_THREADS					( 128 )

/**
 * @macro MAX_FRQ_CNT
 * @brief max #frq in seed filtering
 */
#define MAX_FRQ_CNT					( 7 )

/**
 * @macro UNITTEST
 * @brief set zero to disable unittests
 */
#ifndef UNITTEST
#  define UNITTEST 					( 0 )
#endif

/* make sure POSIX APIs are properly activated */
#if defined(__linux__) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE		200112L
#endif
#if defined(__darwin__) && !defined(_DARWIN_C_FULL)
#  define _DARWIN_C_SOURCE	_DARWIN_C_FULL
#endif

#include <getopt.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/types.h>
#include "sassert.h"


/* max, min */
#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )

/* _likely, _unlikely */
#define _likely(x)		__builtin_expect(!!(x), 1)
#define _unlikely(x)	__builtin_expect(!!(x), 0)

/* _force_inline */
#define _force_inline	inline

/* add namespace */
#ifndef _export
#  ifdef NAMESPACE
#    define _export_cat(x, y)		x##_##y
#    define _export_cat2(x, y)		_export_cat(x, y)
#    define _export(_base)			_export_cat2(NAMESPACE, _base)
#  else
#    define _export(_base)			_base
#  endif
#endif


/**
 * @fn cputime
 * @brief returns consumed CPU time in seconds (double)
 */
static _force_inline
double cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);

	return(
		  r.ru_utime.tv_sec
		+ r.ru_stime.tv_sec
		+ 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec)
	);
}

/**
 * @fn realtime
 * @brief returns elapsed real time in seconds (double)
 */
static _force_inline
double realtime(void)
{
	struct timeval tp;
	gettimeofday(&tp, NULL);

	return(tp.tv_sec + tp.tv_usec * 1e-6);
}

/**
 * @fn version
 * @brief returns pointer to version string
 */
static _force_inline
char const *version(void)
{
	char const *prefix = "minialign-";
	uint64_t spos = strncmp(
		MM_VERSION,
		prefix,
		MIN2(strlen(MM_VERSION), strlen(prefix)))
			? 0
			: strlen(prefix);
	return(&MM_VERSION[spos]);
}

/* malloc wrappers */
/**
 * @struct mm_info_t
 * @brief per-thread information container: dumped when SEGV or other faults are caught.
 */
typedef struct {
	uint64_t enabled;
	char const *msg;
	void *_pad[14];
} mm_info_t;
_static_assert(sizeof(mm_info_t) == 128);	/* fits in a cache line */
mm_info_t info[MAX_THREADS+1] __attribute__(( aligned(64) ));

#define enable_info(t)		{ info[t].enabled = 1; }
#define disable_info(t)		{ info[t].enabled = 0; }
#define set_info(t, x)		{ info[t].msg = (char const *)(x); }

/**
 * @fn oom_abort
 * @brief called when malloc / realloc failed. dump thread information and exit (noreturn).
 */
static
void oom_abort(
	char const *name,
	size_t req)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);

	fprintf(stderr,
		"[E::%s] Out of memory. (required: %zu MB, maxrss: %ld MB)\n",
		name, req / 1024, r.ru_maxrss);

	for(uint64_t i = 0; i < MAX_THREADS + 1; i++) {
		if(info[i].enabled) {
			fprintf(stderr,
				"[E::%s]  thread %" PRIu64 ": %s\n",
				name, i,
				info[i].msg ? info[i].msg : "No information available.");
		}
	}
	exit(128);	/* 128 reserved for out of memory */
}

/**
 * @macro mm_malloc
 */
#define mm_malloc(x) ({ \
	void *_ptr = malloc((size_t)(x)); \
	if(_unlikely(_ptr == NULL)) { \
		oom_abort(__func__, (x)); \
	} \
	_ptr; \
})
#define malloc(x)		mm_malloc(x)

/**
 * @macro mm_realloc
 */
#define mm_realloc(x, y) ({ \
	void *_ptr = realloc((void *)(x), (size_t)(y)); \
	if(_unlikely(_ptr == NULL)) { \
		oom_abort(__func__, (y)); \
	} \
	_ptr; \
})
#define realloc(x, y)	mm_realloc(x, y)

/**
 * @macro mm_calloc
 */
#define mm_calloc(x, y) ({ \
	void *_ptr = calloc((size_t)(x), (size_t)(y)); \
	if(_unlikely(_ptr == NULL)) { \
		oom_abort(__func__, (y)); \
	} \
	_ptr; \
})
#define calloc(x, y)	mm_calloc(x, y)

/**
 * @macro _pnum, _pstr
 * @brief dump number and string to buffer
 */
#define _pnum(type, _buf, _n) ({ \
	uint8_t _b[16] = {0}; \
	type _m = (type)(_n); int64_t _i = 0; \
	while(_m) { \
		_b[_i++] = _m % 10; \
		_m /= 10; \
	} \
	uint64_t _len = _i + (_i == 0); \
	for(int64_t _j = 0; _j < _len; _j++) { \
		_buf[_j] = _b[_len - _j - 1] + '0'; \
	} \
	_len; \
})
#define _pstr(_buf, _str) ({ \
	char const *_p = (_str); \
	char *_q = (_buf); \
	while((*_q++ = *_p++)) {} \
	_q - _buf - 1; \
})

/**
 * @fn mm_rand64
 * @brief 64bit random number generator (for unittest)
 */
static _force_inline
uint64_t mm_rand64(void)
{
	uint64_t bits = 31, n = 0;
	for(uint64_t acc = 0; acc < 64; acc += bits) {
		n <<= bits; n ^= (uint64_t)rand();
	}
	return n;
}

/* end of misc.c */

#define UNITTEST_UNIQUE_ID		1
#include "unittest.h"

#include "lmm.h"
#include "log.h"
#include "gaba_wrap.h"
#include "arch/arch.h"
#include "kvec.h"

unittest_config( .name = "minialign" );

/* test for timers and mallocs */
unittest( .name = "misc.timer" ) {
	double cpu = cputime(), real = realtime();	/* make sure they do not raise segv */
	assert(isnan(cpu) == 0);
	assert(isnan(real) == 0);
}

unittest( .name = "misc.malloc" ) {
	uint64_t const size = 1024 * 1024 * 1024;
	uint8_t *p = malloc(size);
	assert(p != NULL);

	memset(p, 0, size);							/* make sure we can touch this area */

	p = realloc(p, 2*size);
	assert(p != NULL);

	p = realloc(p, size/2);
	assert(p != NULL);
	free(p);
}
/* end of unittest */

/* minimap.h */
typedef union {
	uint64_t u64[1];
	uint32_t u32[2];
} v2u32_t;
_static_assert(sizeof(v2u32_t) == 8);

typedef union {
	uint64_t u64[2];
	uint32_t u32[4];
} v4u32_t;
_static_assert(sizeof(v4u32_t) == 16);

typedef struct { size_t n, m; v4u32_t *a; } v4u32_v;
typedef struct { size_t n, m; v2u32_t *a; } v2u32_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; uint32_t *a; } uint32_v;
typedef struct { size_t n, m; uint16_t *a; } uint16_v;
typedef struct { size_t n, m; uint8_t *a; } uint8_v;
typedef struct { size_t n, m; void **a; } ptr_v;

#include "ksort.h"
#define sort_key_128x(a)		( (a).u64[0] )
KRADIX_SORT_INIT(128x, v4u32_t, sort_key_128x, 8) 
#define sort_key_64x(a)			( (a).u32[0] )
KRADIX_SORT_INIT(64x, v2u32_t, sort_key_64x, 4)
KSORT_INIT_GENERIC(uint32_t)


/* hash.c */

/**
 * @struct kh_s
 * @brief 64-bit key 64-bit value Robinhood hash table
 */
typedef struct kh_s {
	uint32_t mask, max, cnt, ub;		/* size of the array, element count, upper bound */
	v4u32_t *a;
} kh_t;
_static_assert(sizeof(kh_t) == 24);

#define KH_SIZE			( 256 )			/* initial table size */
#define KH_THRESH		( 0.4 )			/* max occupancy */
#define KH_DST_MAX		( 16 )			/* max distance from base pos */
#define KH_INIT_VAL		( UINT64_MAX )	/* initial vals */

#define kh_size(h)		( (h)->mask + 1 )
#define kh_cnt(h)		( (h)->cnt )
#define kh_exist(h, i)	( (h)->a[i].u64[0] + 2 >= 2 )
#define kh_key(h, i)	( (h)->a[i].u64[0] )
#define kh_val(h, i)	( (h)->a[i].u64[1] )

/**
 * @fn kh_init
 * @brief initialize hash with *size* table
 */
static _force_inline
void kh_init_static(kh_t *h, uint64_t size)
{
	/* roundup to power of two */
	size = 0x8000000000000000>>(lzcnt(size - 1) - 1);
	size = MAX2(size, KH_SIZE);

	/* initialize kh object */
	*h = (kh_t){
		.mask = size - 1,		/* in-use table size */
		.max = size,			/* malloc'd table size */
		.cnt = 0,
		.ub = size * KH_THRESH,
		.a = malloc(sizeof(v4u32_t) * size)
	};

	/* init table with invalid key-val pairs */
	for(uint64_t i = 0; i < size; i++) {
		h->a[i].u64[0] = UINT64_MAX;
		h->a[i].u64[1] = KH_INIT_VAL;
	}
	return;
}
#define kh_init(_size)				({ kh_t *h = calloc(1, sizeof(kh_t)); kh_init_static(h, _size); h; })

/**
 * @fn kh_destroy
 */
static _force_inline
void kh_destroy(kh_t *h)
{
	if(h == NULL) { return; }

	free(h->a);
	free(h);
	return;
}
#define kh_destroy_static(_h)	{ free(((kh_t *)_h)->a); }

/**
 * @type khread_t, khwrite_t
 * @brief hash table serializer / deserializer I/O function types
 */
typedef uint64_t (*khread_t)(void *fp, void *buf, uint64_t size);
typedef uint64_t (*khwrite_t)(void *fp, void *buf, uint64_t size);

/**
 * @struct kh_hdr_t
 */
typedef struct kh_hdr_s {
	uint32_t size, cnt;
} kh_hdr_t;
_static_assert(sizeof(kh_hdr_t) == 8);

/**
 * @fn kh_dump
 * @brief binary serializer, fp is an opaque pointer passed to the first argument of wfp
 */
static _force_inline
void kh_dump(kh_t *h, void *fp, khwrite_t wfp)
{
	kh_hdr_t hdr = { 0 };

	/* dump a mark of zero-sized table */
	if(h == NULL || h->a == NULL) {
		wfp(fp, &hdr, sizeof(kh_hdr_t));
		return;
	}

	/* dump size */
	hdr = (kh_hdr_t){
		.size = h->mask + 1,	/* table size */
		.cnt = h->cnt			/* occupancy */
	};
	wfp(fp, &hdr, sizeof(kh_hdr_t));

	/* dump table */
	wfp(fp, h->a, sizeof(v4u32_t) * hdr.size);
	return;
}

/**
 * @fn kh_load_static
 * @brief binary deserializer, fp is an opaque pointer
 */
static _force_inline
void kh_load_static(kh_t *h, void *fp, khread_t rfp)
{
	kh_hdr_t hdr = { 0 };

	/* read sizes, return NULL if table size is zero */
	if((rfp(fp, &hdr, sizeof(kh_hdr_t))) != sizeof(kh_hdr_t) || hdr.size == 0) {
		*h = (kh_t){ 0 };
		return;
	}

	/* create hash object */
	*h = (kh_t){
		.mask = hdr.size - 1,		/* in-use table size */
		.max = hdr.size,			/* malloc'd table size */
		.cnt = hdr.cnt,
		.ub = hdr.size * KH_THRESH,
		.a = malloc(sizeof(v4u32_t) * hdr.size)
	};

	/* read table */
	if((rfp(fp, h->a, sizeof(v4u32_t) * hdr.size)) != sizeof(v4u32_t) * hdr.size) {
		free(h->a);
		free(h);
		return;
	}
	return;
}
#define kh_load(_fp, _rfp)			({ kh_t *h = calloc(1, sizeof(kh_t)); kh_load_static(h, _fp, _rfp); h; })

/**
 * @fn kh_clear
 * @brief flush the content of the hash table
 */
static _force_inline
void kh_clear(kh_t *h)
{
	if(h == 0) { return; }

	/* don't clear max since it holds the malloc'd table size */
	h->mask = KH_SIZE - 1;
	h->cnt = 0;
	h->ub = KH_SIZE * KH_THRESH;

	for(uint64_t i = 0; i < KH_SIZE; i++) {
		h->a[i].u64[0] = UINT64_MAX;
		h->a[i].u64[1] = KH_INIT_VAL;
	}
	return;
}

#if 0
/**
 * @fn kh_put_intl
 * @brief insert key-val pair or replace val if key already exists, returns #bin newly used
 */
static _force_inline
uint64_t kh_put_intl(v4u32_t *a, uint64_t k, uint64_t v, uint64_t mask)
{
	uint64_t kc = k, pos = k & mask;		/* current key and candidate pos */
	while(1) {
		uint64_t kt = a[pos].u64[0];		/* load key onto register */

		/* test if the current key-val pair can be pushed to the pos */
		if(kc <= kt) {						/* k < kt || is_empty(kt) || is_moved(kt) */
			uint64_t vt = a[pos].u64[1];	/* save existing val onto register */
			a[pos] = (v4u32_t){
				.u64 = { k, v }
			};
			if(kt + 2 < 2) { return(1); }	/* is_empty(kt) || is_moved(kt) */
			if(k == kt) { return(0); }		/* k == kt, val was replaced by the new one */
			k = kc = kt; v = vt;			/* robinhood swap */
		}

		/* adjust key when rounded, update pos for the next itr */
		kc -= (mask + 1) & (pos + 1);
		pos = mask & (pos + 1);
	}
	return(0);
}
#endif

/**
 * @fn kh_allocate
 * @brief allocate a bucket for key, returns the index of allocated bucket.
 */
typedef struct { uint64_t idx, n; } kh_bidx_t;
static _force_inline
kh_bidx_t kh_allocate(v4u32_t *a, uint64_t k, uint64_t v, uint64_t mask)
{
	#define _poll_bucket(_i, _b0) ({ \
		int64_t _b = (_b0); \
		uint64_t _k1; \
		while(1) { \
			/* test: is_empty(kt) || is_moved(kt) */ \
			_k1 = a[_i].u64[0]; \
			if(_b <= (int64_t)(_k1 & mask) + (_k1 + 2 < 2)) { break; } \
			_b -= (_i + 1) & (mask + 1); \
			_i = (_i + 1) & mask; \
		} \
		_k1;	/* return key */ \
	})


	uint64_t i = k & mask, k0 = k, v0 = v;
	uint64_t k1 = _poll_bucket(i, i & mask);		/* search first bucket */

	if(k0 == k1) {									/* duplicated key */
		// a[i].u64[0] = k0;
		return((kh_bidx_t){ i, 0 });
	}

	uint64_t j = i;
	a[i].u64[0] = k0;								/* update key */
	while(k1 + 2 >= 2) {
		uint64_t v1 = a[i].u64[1];					/* load for swap */
		a[i].u64[1] = v0;							/* save previous value */
		k0 = k1; v0 = v1;							/* robinhood swap */
		i = (i + 1) & mask;							/* advance index */
		k1 = _poll_bucket(i, k0 & mask);			/* search next bucket */
		a[i].u64[0] = k0;							/* update key */
	}
	a[i].u64[1] = v0;								/* save previous value */
	return((kh_bidx_t){ j, 1 });

	#undef _poll_bucket
}

/**
 * @fn kh_extend
 * @brief extend hash table
 */
static _force_inline
void kh_extend(kh_t *h)
{
	uint64_t prev_size = h->mask + 1, size = 2 * prev_size, mask = size - 1;

	debug("kh_extend called, new_size(%llx), cnt(%u), ub(%u)", size, h->cnt, h->ub);

	/* update size */
	h->mask = mask;
	h->ub = size * KH_THRESH;

	/* double the table if needed */
	if(size > h->max) {
		h->a = realloc(h->a, sizeof(v4u32_t) * size);
		h->max = size;
	}

	/* clear the extended area */
	for(uint64_t i = 0; i < prev_size; i++) {
		h->a[i + prev_size].u64[0] = UINT64_MAX;
		h->a[i + prev_size].u64[1] = KH_INIT_VAL;
	}

	/* reallocate bins */
	for(uint64_t i = 0; i < size; i++) {
		uint64_t k = h->a[i].u64[0];	/* load key */

		/* test if reallocate is required */
		if(k + 2 < 2 || (k & mask) == i) { continue; }

		uint64_t v = h->a[i].u64[1];	/* load val */
		h->a[i] = (v4u32_t){			/* mark the current bin moved */
			.u64 = { UINT64_MAX-1, KH_INIT_VAL }
		};
		kh_allocate(h->a, k, v, mask);
	}
	return;
}

/**
 * @fn kh_put
 */
static _force_inline
void kh_put(kh_t *h, uint64_t key, uint64_t val)
{
	if(h->cnt >= h->ub) {
		debug("trigger extend, cnt(%u), ub(%u)", h->cnt, h->ub);
		kh_extend(h);
	}
	kh_bidx_t b = kh_allocate(h->a, key, val, h->mask);
	h->cnt += b.n;
	// h->ub = ((b.idx - key) & h->mask) > KH_DST_MAX ? 0 : h->ub;
	h->a[b.idx] = (v4u32_t){
		.u64 = { key, val }
	};
	return;
}

/**
 * @fn kh_put_ptr
 */
static _force_inline
uint64_t *kh_put_ptr(kh_t *h, uint64_t key, uint64_t extend)
{
	if(extend != 0 && h->cnt >= h->ub) {
		debug("trigger extend, cnt(%u), ub(%u)", h->cnt, h->ub);
		kh_extend(h);
	}
	kh_bidx_t b = kh_allocate(h->a, key, KH_INIT_VAL, h->mask);
	debug("allocated hash bin (%llu, %llu), (%llx, %llx), dst(%lld)", b.idx, b.n, h->a[b.idx].u64[0], h->a[b.idx].u64[1], (b.idx - key) & h->mask);

	h->cnt += b.n;
	// h->ub = ((b.idx - key) & h->mask) > KH_DST_MAX ? 0 : h->ub;
	return(&h->a[b.idx].u64[1]);
}

/**
 * @fn kh_get
 * @brief returns val or UINT64_MAX if not found
 */
static _force_inline
uint64_t kh_get(kh_t *h, uint64_t key)
{
	uint64_t mask = h->mask, pos = key & mask, k;
	do {
		if((k = h->a[pos].u64[0]) == key) {
			return(h->a[pos].u64[1]);
		}
		pos = mask & (pos + 1);
	} while(k + 1 != 0);			/* !is_empty(k) || is_moved(k) */
	return(UINT64_MAX);
}

/**
 * @fn kh_get_ptr
 * @brief returns pointer to val or NULL if not found
 */
static uint64_t *kh_get_ptr(kh_t *h, uint64_t key)
{
	uint64_t mask = h->mask, pos = key & mask, k;
	do {
		if((k = h->a[pos].u64[0]) == key) {
			return(&h->a[pos].u64[1]);
		}
		pos = mask & (pos + 1);
	} while(k + 1 != 0);			/* !is_empty(k) || is_moved(k) */
	return(NULL);
}

unittest( .name = "kh.base" ) {
	kh_t *h = kh_init(0);
	assert(h != NULL);

	uint64_t const kmask = mm_rand64(), vmask = mm_rand64(), cnt = 1024 * 1024;
	uint64_t const *p;

	assert(kh_cnt(h) == 0, "cnt(%lu)", kh_cnt(h));
	for(uint64_t i = 0; i < 2*cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p == NULL);
		assert(kh_get(h, i^kmask) == UINT64_MAX);
	}

	/* put key-val pairs */
	for(uint64_t i = 0; i < cnt; i++) kh_put(h, i^kmask, i^vmask);
	assert(kh_cnt(h) == cnt, "cnt(%lu)", kh_cnt(h));
	for(uint64_t i = 0; i < cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p != NULL, "i(%lu)", i);
		assert(*p == (i^vmask), "i(%lu), val(%lu, %lu)", i, *p, (i^vmask));
		assert(kh_get(h, i^kmask) == (i^vmask), "i(%lu), val(%lu, %lu)", i, *p, (i^vmask));
	}
	for(uint64_t i = cnt; i < 2*cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p == NULL);
		assert(kh_get(h, i^kmask) == UINT64_MAX);
	}

	/* clear */
	kh_clear(h);
	assert(kh_cnt(h) == 0, "cnt(%lu)", kh_cnt(h));
	for(uint64_t i = 0; i < cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p == NULL);
		assert(kh_get(h, i^kmask) == UINT64_MAX);
	}

	/* add another set of key-val pairs */
	for(uint64_t i = cnt; i < 2*cnt; i++) kh_put(h, i^kmask, i^vmask);
	assert(kh_cnt(h) == cnt, "cnt(%lu)", kh_cnt(h));
	for(uint64_t i = 0; i < cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p == NULL);
		assert(kh_get(h, i^kmask) == UINT64_MAX);
	}
	for(uint64_t i = cnt; i < 2*cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p != NULL, "i(%lu)", i);
		assert(*p == (i^vmask), "i(%lu), val(%lu, %lu)", i, *p, (i^vmask));
		assert(kh_get(h, i^kmask) == (i^vmask), "i(%lu), val(%lu, %lu)", i, *p, (i^vmask));
	}

	kh_destroy(h);
}

unittest( .name = "kh.io" ) {
	kh_t *h = kh_init(0);
	uint64_t const kmask = mm_rand64(), vmask = mm_rand64(), cnt = 1024 * 1024;
	for(uint64_t i = 0; i < cnt; i++) kh_put(h, i^kmask, i^vmask);

	/* dump */
	char const *filename = "./minialign.unittest.kh.tmp";
	gzFile fp = gzopen(filename, "w");
	assert((void *)fp != NULL);
	kh_dump(h, (void *)fp, (khwrite_t)gzwrite);
	gzclose(fp);
	kh_destroy(h);

	/* restore */
	fp = gzopen(filename, "r");
	assert((void *)fp != NULL);
	h = kh_load((void *)fp, (khread_t)gzread);
	assert(h != NULL);
	gzclose(fp);

	uint64_t const *p;
	assert(kh_cnt(h) == cnt, "cnt(%lu)", kh_cnt(h));
	for(uint64_t i = 0; i < cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p != NULL, "i(%lu)", i);
		assert(*p == (i^vmask), "i(%lu), val(%lu, %lu)", i, *p, (i^vmask));
		assert(kh_get(h, i^kmask) == (i^vmask), "i(%lu), val(%lu, %lu)", i, *p, (i^vmask));
	}
	for(uint64_t i = cnt; i < 2*cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p == NULL);
		assert(kh_get(h, i^kmask) == UINT64_MAX);
	}
	kh_destroy(h);
	remove(filename);
}
/* end of hash.c */

/* queue.c */
/**
 * @type pt_source_t, pt_worker_t, pt_drain_t
 * @brief callback functions types for multithreaded stream
 */
typedef void *(*pt_source_t)(uint32_t tid, void *arg);
typedef void *(*pt_worker_t)(uint32_t tid, void *arg, void *item);
typedef void (*pt_drain_t)(uint32_t tid, void *arg, void *item);

/**
 * @struct pt_q_s
 * @brief lock-based queue context
 */
typedef struct pt_q_s {
	uint64_t lock, head, tail, size;
	void **elems;
	uint64_t _pad1[3];
} pt_q_t;
_static_assert(sizeof(pt_q_t) == 64);		/* to occupy a single cache line */

/**
 * @struct pt_thread_s
 * @brief thread-local worker context
 */
typedef struct pt_thread_s {
	pthread_t th;
	uint64_t tid;
	pt_q_t *in, *out;
	volatile pt_worker_t wfp;
	volatile void *warg;
} pt_thread_t;

/**
 * @struct pt_s
 * @brief parallel task processor context
 */
typedef struct pt_s {
	pt_q_t in, out;
	uint32_t nth;
	pt_thread_t c[];	/* [0] is reserved for master */
} pt_t;

#define PT_EMPTY	( (void *)(UINT64_MAX) )
#define PT_EXIT		( (void *)(UINT64_MAX-1) )

/**
 * @fn pt_enq
 * @brief enqueue; concurrent queue is better while currently lock-based for simplicity
 */
static _force_inline
uint64_t pt_enq(pt_q_t *q, uint64_t tid, void *elem)
{
	uint64_t z, ret = UINT64_MAX;

	/* lock by thread id */
	do { z = UINT32_MAX; } while(!cas(&q->lock, &z, tid));

	/* push element to queue */
	uint64_t head = q->head, tail = q->tail, size = q->size;
	if(((head + 1) % size) != tail) {
		q->elems[head] = elem;
		q->head = (head + 1) % size;
		ret = 0;
	}

	/* release */
	do { z = tid; } while(!cas(&q->lock, &z, UINT32_MAX));
	return(ret);
}

/**
 * @fn pt_deq
 */
static _force_inline
void *pt_deq(pt_q_t *q, uint64_t tid)
{
	void *elem = PT_EMPTY;
	uint64_t z;

	/* lock by thread id */
	do { z = UINT32_MAX; } while(!cas(&q->lock, &z, tid));

	/* get element from queue */
	uint64_t head = q->head, tail = q->tail, size = q->size;
	if(head != tail) {
		elem = q->elems[tail]; q->tail = (tail + 1) % size;
	}

	/* release */
	do { z = tid; } while(!cas(&q->lock, &z, UINT32_MAX));
	return(elem);
}

/**
 * @fn pt_dispatch
 * @brief per-thread function dispatcher, with ping-pong prefetching
 */
static _force_inline
void *pt_dispatch(void *s)
{
	pt_thread_t *c = (pt_thread_t *)s;
	struct timespec tv = { .tv_nsec = 512 * 1024 };

	void *ping = PT_EMPTY, *pong = PT_EMPTY;
	while(1) {
		ping = pt_deq(c->in, c->tid);		/* prefetch */
		if(ping == PT_EMPTY && pong == PT_EMPTY) {
			nanosleep(&tv, NULL);			/* no task is available, sleep for a while */
		}
		if(pong != PT_EMPTY) {
			pt_enq(c->out, c->tid, c->wfp(c->tid, (void *)c->warg, pong));
		}
		if(ping == PT_EXIT) { break; }		/* terminate thread */

		pong = pt_deq(c->in, c->tid);		/* prefetch */
		if(ping == PT_EMPTY && pong == PT_EMPTY) {
			nanosleep(&tv, NULL);			/* no task is available, sleep for a while */
		}
		if(ping != PT_EMPTY) {
			pt_enq(c->out, c->tid, c->wfp(c->tid, (void *)c->warg, ping));
		}
		if(pong == PT_EXIT) { break; }		/* terminate thread */
	}
	return(NULL);
}

/**
 * @fn pt_destroy
 */
static _force_inline
void pt_destroy(pt_t *pt)
{
	void *status;

	/* send termination signal */
	for(uint64_t i = 1; i < pt->nth; i++) {
		pt_enq(pt->c->in, pt->c->tid, PT_EXIT);
	}

	/* wait for the threads terminate */
	for(uint64_t i = 1; i < pt->nth; i++) {
		disable_info(i);		/* disable debug information */
		pthread_join(pt->c[i].th, &status);
	}

	/* clear queues */
	while(pt_deq(pt->c->in, pt->c->tid) != PT_EMPTY) {}
	while(pt_deq(pt->c->out, pt->c->tid) != PT_EMPTY) {}

	/* clear queues and objects */
	free(pt->in.elems);
	free(pt->out.elems);
	free(pt);
	return;
}

/**
 * @fn pt_init
 */
static _force_inline
pt_t *pt_init(uint32_t nth)
{
	/* init object */
	nth = (nth == 0)? 1 : nth;
	pt_t *pt = calloc(1, sizeof(pt_t) + sizeof(pt_thread_t) * nth);

	/* init queues (note: #elems can be larger for better performance?) */
	pt->in = (pt_q_t){
		.lock = UINT32_MAX,
		.elems = calloc(8 * nth, sizeof(void *)),
		.size = 8 * nth
	};
	pt->out = (pt_q_t){
		.lock = UINT32_MAX,
		.elems = calloc(8 * nth, sizeof(void *)),
		.size = 8 * nth
	};

	/* init parent thread info */
	pt->nth = nth; pt->c[0].tid = 0;
	pt->c[0].in = &pt->in;
	pt->c[0].out = &pt->out;

	/* init children info, create children */
	for(uint64_t i = 1; i < nth; i++) {
		pt->c[i].tid = i;
		pt->c[i].in = &pt->in;
		pt->c[i].out = &pt->out;

		pthread_create(&pt->c[i].th, NULL, pt_dispatch, (void *)&pt->c[i]);
		enable_info(i);			/* enable debug information */
	}
	return pt;
}

/**
 * @fn pt_set_worker
 * @brief update worker function and argument pointers
 */
static _force_inline
int pt_set_worker(
	pt_t *pt,
	pt_worker_t wfp, void **warg)
{
	void *item;

	/* fails when unprocessed object exists in the queue */
	if((item = pt_deq(pt->c->in, pt->c->tid)) != PT_EMPTY) {
		pt_enq(pt->c->in, pt->c->tid, item);
		return(-1);
	}

	/* update pointers */
	for(uint64_t i = 0; i < pt->nth; i++) {
		pt->c[i].wfp = wfp;
		pt->c[i].warg = warg ? warg[i] : NULL;
	}

	/* syncronize */
	fence();
	return(0);
}

/**
 * @fn pt_stream
 * @brief multithreaded stream, source and drain are always called in the parent thread
 */
static _force_inline
int pt_stream(
	pt_t *pt,
	pt_source_t sfp, void *sarg,
	pt_worker_t wfp, void **warg,
	pt_drain_t dfp, void *darg)
{
	if(pt_set_worker(pt, wfp, warg)) { return(-1); }

	/* keep balancer between [lb, ub) */
	uint64_t const lb = 2 * pt->nth, ub = 4 * pt->nth;
	uint64_t bal = 0;
	void *item;

	while((item = sfp(pt->c->tid, sarg)) != NULL) {
		pt_enq(pt->c->in, pt->c->tid, item);
		if(++bal < ub) { continue; }

		while(bal > lb) {
			/* flush to drain (note: while loop is better?) */
			if((item = pt_deq(pt->c->out, pt->c->tid)) != PT_EMPTY) {
				bal--;
				dfp(pt->c->tid, darg, item);
			}

			/* process one in the master (parent) thread */
			if((item = pt_deq(pt->c->in, pt->c->tid)) != PT_EMPTY) {
				bal--;
				dfp(pt->c->tid, darg, wfp(pt->c->tid, warg ? *warg : NULL, item));
			}
		}
	}

	/* source depleted, process remainings */
	while((item = pt_deq(pt->c->in, pt->c->tid)) != PT_EMPTY) {
		pt_enq(pt->c->out, pt->c->tid, wfp(pt->c->tid, warg ? *warg : NULL, item));
	}

	/* flush results */
	while(bal > 0) {
		if((item = pt_deq(pt->c->out, pt->c->tid)) != PT_EMPTY) {
			bal--;
			dfp(pt->c->tid, darg, item);
		}
	}
	return(0);
}

/**
 * @fn pt_parallel
 */
static _force_inline
int pt_parallel(
	pt_t *pt,
	pt_worker_t wfp, void **warg,
	void **item)
{
	/* fails when unprocessed element exists */
	if(pt_set_worker(pt, wfp, warg)) { return(-1); }

	/* push items */
	for(uint64_t i = 1; i < pt->nth; i++) {
		pt_enq(pt->c->in, pt->c->tid, item? item[i] : NULL);
	}
	debug("pushed items");

	/* process the first one in the master (parent) thread */
	wfp(pt->c->tid, warg? warg[0] : NULL, item? item[0] : NULL);
	debug("finished master");

	/* wait for the children done */
	for(uint64_t i = 1; i < pt->nth; i++) {
		while(pt_deq(pt->c->out, pt->c->tid) == PT_EMPTY) {
			struct timespec tv = { .tv_nsec = 512 * 1024 };
			nanosleep(&tv, NULL);	/* wait for a while */
			/* sched_yield(); */
		}
		debug("joined i(%llu)", i);
	}
	return 0;
}

/* unittests */
static void *pt_unittest_source(uint32_t tid, void *arg)
{
	uint64_t *s = (uint64_t*)arg;
	if(*s >= 1024) return NULL;
	uint64_t *p = malloc(sizeof(uint64_t));
	*p = *s; *s = *s + 1;
	return p;
}
static void *pt_unittest_worker(uint32_t tid, void *arg, void *item)
{
	uint64_t *p = (uint64_t *)item, *a = (uint64_t *)arg;
	*p += *a;
	return p;
}
static void pt_unittest_drain(uint32_t tid, void *arg, void *item)
{
	uint64_t *d = (uint64_t*)arg, *p = (uint64_t*)item;
	*d += *p;
	free(item);
}

unittest( .name = "pt.single" ) {
	pt_t *pt = pt_init(1);
	assert(pt != NULL);

	uint64_t icnt = 0, ocnt = 0, inc = 1, *arr[1] = { &inc };
	pt_stream(pt, pt_unittest_source, (void *)&icnt, pt_unittest_worker, (void **)arr, pt_unittest_drain, (void *)&ocnt);
	assert(icnt == 1024, "icnt(%lu)", icnt);
	assert(ocnt == 512*1025, "ocnt(%lu)", ocnt);

	uint64_t s[4] = { 0 }, *sp[1] = { &s[0] };
	pt_parallel(pt, pt_unittest_worker, (void **)arr, (void **)sp);
	assert(s[0] == 1, "d[0](%lu)", s[0]);
	pt_destroy(pt);
}

unittest( .name = "pt.multi" ) {
	pt_t *pt = pt_init(4);
	assert(pt != NULL);

	uint64_t icnt = 0, ocnt = 0, inc = 1, *arr[4] = { &inc, &inc, &inc, &inc };
	pt_stream(pt, pt_unittest_source, (void *)&icnt, pt_unittest_worker, (void **)arr, pt_unittest_drain, (void *)&ocnt);
	assert(icnt == 1024, "icnt(%lu)", icnt);
	assert(ocnt == 512*1025, "ocnt(%lu)", ocnt);

	uint64_t s[4] = { 0,1,2,3 }, *sp[4] = { &s[0],&s[1],&s[2],&s[3] };
	pt_parallel(pt, pt_unittest_worker, (void **)arr, (void **)sp);
	for(uint64_t i = 0; i < 4; i++) {
		assert(s[i] == (i+1), "i(%lu), d[i](%lu)", i, s[i]);
	}
	pt_destroy(pt);
}

/* stdio stream with multithreaded compression / decompression */
/**
 * @struct pg_block_t
 * @brief compression / decompression unit block
 */
typedef struct {
	uint64_t head, len;				/* head pointer (index) and block length */
	uint32_t id;					/* block id */
	uint8_t raw, flush, _pad[2];	/* raw: 1 if compressed, flush: 1 if needed to dump */
	uint8_t buf[];
} pg_block_t;

/**
 * @struct pg_t
 * @brief context
 */
typedef struct {
	FILE *fp;
	pt_t *pt;
	pg_block_t *s;
	uint32_t ub, lb, bal, icnt, ocnt, eof, nth;
	uint64_t block_size;
	kvec_t(v4u32_t) hq;
	void *c[];
} pg_t;
#define incq_comp(a, b)		( (int64_t)(a).u64[0] - (int64_t)(b).u64[0] )

/**
 * @fn pg_deflate
 * @brief compress block
 */
static
pg_block_t *pg_deflate(pg_block_t *in, uint64_t block_size)
{
	/* create dest block */
	uint64_t buf_size = block_size * 1.2;
	pg_block_t *out = malloc(sizeof(pg_block_t) + buf_size);

	/* compress (deflate) */
	z_stream zs = {
		.next_in = in->buf,   .avail_in = in->len,
		.next_out = out->buf, .avail_out = buf_size
	};
	deflateInit2(&zs, 1, Z_DEFLATED, 15, 8, Z_DEFAULT_STRATEGY);
	deflate(&zs, Z_FINISH);
	deflateEnd(&zs);

	/* set metadata */
	out->head = 0;
	out->len = buf_size - zs.avail_out;
	out->id = in->id;
	out->raw = 0;
	out->flush = 1;

	/* cleanup input block */
	free(in);
	return(out);
}

/**
 * @fn pg_inflate
 * @brief decompress block
 */
static
pg_block_t *pg_inflate(pg_block_t *in, uint64_t block_size)
{
	/* create dest block */
	uint64_t buf_size = block_size * 1.2;
	pg_block_t *out = malloc(sizeof(pg_block_t) + buf_size);

	/* inflate */
	z_stream zs = {
		.next_in = in->buf, .avail_in = in->len,
		.next_out = out->buf, .avail_out = buf_size
	};
	inflateInit2(&zs, 15);
	inflate(&zs, Z_FINISH);
	inflateEnd(&zs);

	/* set metadata */
	out->head = 0;
	out->len = buf_size - zs.avail_out;
	out->id = in->id;
	out->raw = 1;
	out->flush = 0;

	/* cleanup input block */
	free(in);
	return(out);
}

/**
 * @fn pg_worker
 * @brief thread-local worker
 */
static
void *pg_worker(uint32_t tid, void *arg, void *item)
{
	pg_t *pg = (pg_t*)arg;
	pg_block_t *s = (pg_block_t*)item;
	if(s == NULL || s->len == 0) { return(s); }

	char buf[128], *p = buf;
	p += _pstr(p, "[pg_worker] bucket id: "); p += _pnum(uint32_t, p, s->id); *p = '\0';
	set_info(tid, buf);
	return((s->raw? pg_deflate : pg_inflate)(s, pg->block_size));
}

/**
 * @fn pg_read_block
 * @brief read compressed block from input stream
 */
static _force_inline
pg_block_t *pg_read_block(pg_t *pg)
{
	pg_block_t *s = malloc(sizeof(pg_block_t) + pg->block_size);

	/* read block */
	if(fread(&s->len, sizeof(uint64_t), 1, pg->fp) != 1 || s->len == 0) { goto _fail; }
	if(fread(s->buf, sizeof(uint8_t), s->len, pg->fp) != s->len) { goto _fail; }

	/* set metadata */
	s->head = 0;
	s->id = pg->icnt++;
	s->raw = 0;
	s->flush = 0;
	return(s);
_fail:
	free(s);
	return(NULL);
}

/**
 * @fn pg_write_block
 * @brief write compressed block to output stream
 */
static _force_inline
void pg_write_block(pg_t *pg, pg_block_t *s)
{
	if(s->len == 0) return;

	/* dump header then block */
	fwrite(&s->len, sizeof(uint64_t), 1, pg->fp);
	fwrite(s->buf, sizeof(uint8_t), s->len, pg->fp);
	free(s);
	return;
}

/**
 * @fn pg_init
 * @brief initialize stream with fp
 */
static _force_inline
pg_t *pg_init(FILE *fp, uint32_t nth)
{
	/* create context */
	pg_t *pg = calloc(1, sizeof(pg_t) + nth * sizeof(pg_t*));
	*pg = (pg_t){
		.fp = fp,
		.pt = pt_init(nth),
		.lb = nth, .ub = 3 * nth, .bal = 0, .nth = nth,
		.block_size = 1024 * 1024
	};
	kv_hq_init(pg->hq);

	/* init worker args */
	for(uint64_t i = 0; i < nth; i++) {
		pg->c[i] = (void *)pg;
	}
	pt_set_worker(pg->pt, pg_worker, pg->c);

	return(pg);
}

/**
 * @fn pg_destroy
 */
static _force_inline
void pg_destroy(pg_t *pg)
{
	/* process current working block */
	pg_block_t *s = pg->s, *t;
	if(s && s->flush == 1 && s->head != 0) {
		s->len = s->head;
		if(pg->nth == 1) {
			pg_write_block(pg, pg_deflate(s, pg->block_size));
		} else {
			pg->bal++;
			pt_enq(pg->pt->c->in, pg->pt->c->tid, s);
		}
	} else {
		free(s);
	}

	/* process remainings */
	while(pg->bal > 0) {
		if((t = pt_deq(pg->pt->c->out, pg->pt->c->tid)) == PT_EMPTY) {
			/* wait for a while(note: nanosleep is better?) */
			sched_yield(); continue;
		}

		/* item fetched, flush if needed */
		pg->bal--;
		if(t && !t->flush) { free(t); continue; }

		/* push heapqueue to sort */
		kv_hq_push(v4u32_t, incq_comp, pg->hq, ((v4u32_t){.u64 = {t->id, (uintptr_t)t}}));
	}

	/* flush heapqueue */
	while(pg->hq.n > 1) {
		pg->ocnt++;
		t = (pg_block_t*)kv_hq_pop(v4u32_t, incq_comp, pg->hq).u64[1];
		pg_write_block(pg, t);
	}

	/* write terminator */
	uint64_t z = 0;
	fwrite(&z, sizeof(uint64_t), 1, pg->fp);
	
	/* cleanup contexts */
	pt_destroy(pg->pt);
	kv_hq_destroy(pg->hq);
	free(pg);
	return;
}

/**
 * @fn pgread
 */
static _force_inline
uint64_t pgread(pg_t *pg, void *dst, uint64_t len)
{
	uint64_t rem = len;
	pg_block_t *s = pg->s, *t;
	if(pg->eof == 2) return(0);

	while(rem > 0) {
		/* check and prepare a valid inflated block */
		while(s == NULL || s->head == s->len) {
			free(s); s = NULL;
			if(pg->nth == 1) {
				/* single-threaded */
				if((t = pg_read_block(pg)) == NULL) {
					pg->eof = 2;
					return(len-rem);
				}
				pg->s = s = pg_inflate(t, pg->block_size);
			} else {
				/* multithreaded; read compressed blocks and push them to queue */
				while(!pg->eof && pg->bal < pg->ub) {
					if((t = pg_read_block(pg)) == NULL) {
						pg->eof = 1; break;
					}
					pg->bal++;
					pt_enq(pg->pt->c->in, pg->pt->c->tid, t);
				}

				/* fetch inflated blocks and push heapqueue to sort */
				while((t = pt_deq(pg->pt->c->out, pg->pt->c->tid)) != PT_EMPTY) {
					pg->bal--;
					kv_hq_push(v4u32_t, incq_comp, pg->hq, ((v4u32_t){.u64 = {t->id, (uintptr_t)t}}));
				}

				/* check if input depleted */
				if(pg->ocnt >= pg->icnt) {
					pg->eof = 2;
					return(len-rem);
				}

				/* check if inflated block is available */
				if(pg->hq.n < 2 || pg->hq.a[1].u64[0] > pg->ocnt) {
					sched_yield(); continue;	/* wait for a while */
				}

				/* block is available! */
				pg->ocnt++;
				pg->s = s = (pg_block_t*)kv_hq_pop(v4u32_t, incq_comp, pg->hq).u64[1];
			}
		}

		/* copy to dst buffer */
		uint64_t adv = MIN2(rem, s->len-s->head);
		memcpy(dst + len - rem, s->buf + s->head, adv);
		rem -= adv; s->head += adv;
	}
	return(len);
}

/**
 * @fn pgwrite
 */
static _force_inline
uint64_t pgwrite(pg_t *pg, const void *src, uint64_t len)
{
	uint64_t rem = len;
	pg_block_t *s = pg->s, *t;
	while(rem > 0) {
		/* push the current block to queue and prepare an empty one if the current one is full */
		if(s == NULL || s->head == s->len) {
			if(pg->nth == 1) {
				/* single-threaded */
				if(s != NULL) {
					pg_write_block(pg, pg_deflate(s, pg->block_size));
				}
			} else {
				/* multithreaded; fetch copressed block and push to heapqueue to sort */
				while(pg->bal > pg->lb) {
					if((t = pt_deq(pg->pt->c->out, pg->pt->c->tid)) == PT_EMPTY) {
						if(pg->bal >= pg->ub) {
							sched_yield(); continue;		/* queue full, wait for a while */
						} else {
							break;							/* unexpected error */
						}
					}
					pg->bal--;
					kv_hq_push(v4u32_t, incq_comp, pg->hq, ((v4u32_t){.u64 = {t->id, (uintptr_t)t}}));
				}

				/* flush heapqueue */
				while(pg->hq.n > 1 && pg->hq.a[1].u64[0] <= pg->ocnt) {
					pg->ocnt++;
					t = (pg_block_t*)kv_hq_pop(v4u32_t, incq_comp, pg->hq).u64[1];
					pg_write_block(pg, t);
				}

				/* push the current block to deflate queue */
				if(s != NULL) {
					pg->bal++;
					pt_enq(pg->pt->c->in, pg->pt->c->tid, s);
				}
			}

			/* create new block */
			s = malloc(sizeof(pg_block_t) + pg->block_size);
			s->head = 0;
			s->len = pg->block_size;
			s->id = pg->icnt++;
			s->raw = 1;
			s->flush = 1;
		}

		/* copy the content */
		uint64_t adv = MIN2(rem, s->len-s->head);
		memcpy(s->buf + s->head, src + len - rem, adv);
		rem -= adv; s->head += adv;
	}
	pg->s = s;
	return len;
}

unittest( .name = "pg.single" ) {
	uint64_t const size = 1024 * 1024 * 1024;
	uint8_t *p = malloc(size), *q = malloc(size);
	for(uint64_t i = 0; i < size; i++) p[i] = i % 253;

	char const *filename = "./minialign.unittest.pg.tmp";
	FILE *fp = fopen(filename, "w");
	assert(fp != NULL);
	pg_t *pg = pg_init(fp, 1);
	assert(pg != NULL);

	for(uint64_t i = 0; i < 3; i++) {
		uint64_t l = pgwrite(pg, p, size);
		assert(l == size, "l(%lu), size(%lu)", l, size);
	}
	pg_destroy(pg);
	fclose(fp);


	fp = fopen(filename, "r");
	assert(fp != NULL);
	pg = pg_init(fp, 1);
	assert(pg != NULL);

	for(uint64_t i = 0; i < 3; i++) {
		uint64_t l = pgread(pg, q, size);
		assert(l == size, "l(%lu), size(%lu)", l, size);
		assert(memcmp(p, q, size) == 0);
	}
	pg_destroy(pg);
	fclose(fp);
	free(p); free(q); remove(filename);
}

unittest( .name = "pg.multi" ) {
	uint64_t const size = 1024 * 1024 * 1024;
	uint8_t *p = malloc(size), *q = malloc(size);
	for(uint64_t i = 0; i < size; i++) p[i] = i % 253;

	char const *filename = "./minialign.unittest.pg.tmp";
	FILE *fp = fopen(filename, "w");
	assert(fp != NULL);
	pg_t *pg = pg_init(fp, 4);
	assert(pg != NULL);

	for(uint64_t i = 0; i < 3; i++) {
		uint64_t l = pgwrite(pg, p, size);
		assert(l == size, "l(%lu), size(%lu)", l, size);
	}
	pg_destroy(pg);
	fclose(fp);


	fp = fopen(filename, "r");
	assert(fp != NULL);
	pg = pg_init(fp, 4);
	assert(pg != NULL);

	for(uint64_t i = 0; i < 3; i++) {
		uint64_t l = pgread(pg, q, size);
		assert(l == size, "l(%lu), size(%lu)", l, size);
		assert(memcmp(p, q, size) == 0);
	}
	pg_destroy(pg);
	fclose(fp);
	free(p); free(q); remove(filename);
}
/* end of queue.c */

/* bamlite.c in bwa */

/**
 * @struct bam_header_t
 * @brief bam file header container
 */
typedef struct {
	int32_t n_targets;
	char **target_name;
	uint32_t *target_len;
	size_t l_text, n_text;
	char *text;
} bam_header_t;

/**
 * @struct bam_core_t
 * @brief bam record header
 */
typedef struct {
	int32_t tid;
	int32_t pos;
	uint8_t l_qname, qual;
	uint16_t bin, n_cigar, flag;
	int32_t l_qseq;
	int32_t mtid;
	int32_t mpos;
	int32_t isize;
} bam_core_t;
_static_assert(sizeof(bam_core_t) == 32);

/**
 * @fn bam_header_destroy
 */
static _force_inline
void bam_header_destroy(bam_header_t *header)
{
	if(header == 0) { return; }

	/* cleanup name array */
	if(header->target_name) {
		for(uint64_t i = 0; i < (uint64_t)header->n_targets; i++) {
			if(header->target_name[i] != NULL) {
				free(header->target_name[i]);
			}
		}
		free(header->target_name);
	}

	/* cleanup len array */
	if(header->target_len) { free(header->target_len); }

	/* raw text */
	if(header->text) { free(header->text); }
	free(header);

	return;
}

/**
 * @fn bam_read_header
 */
static _force_inline
bam_header_t *bam_read_header(gzFile fp)
{
	/* read "BAM1" (magic) */
	char magic[4];
	if(gzread(fp, magic, 4) != 4 || strncmp(magic, "BAM\001", 4) != 0) {
		return NULL;
	}

	/* create object */
	bam_header_t *header = calloc(1, sizeof(bam_header_t));
	
	/* read plain text and the number of reference sequences */
	if(gzread(fp, &header->l_text, 4) != 4) {
		goto _bam_read_header_fail;
	}
	header->text = (char*)calloc(header->l_text + 1, 1);
	if((size_t)gzread(fp, header->text, header->l_text) != header->l_text) {
		goto _bam_read_header_fail;
	}
	if(gzread(fp, &header->n_targets, 4) != 4) {
		goto _bam_read_header_fail;
	}

	/* read reference sequence names and lengths */
	header->target_name = (char**)calloc(header->n_targets, sizeof(char*));
	header->target_len = (uint32_t*)calloc(header->n_targets, 4);
	for(uint64_t i = 0; i != (uint64_t)header->n_targets; i++) {
		int32_t name_len;
		if(gzread(fp, &name_len, 4) != 4) {
			goto _bam_read_header_fail;
		}

		header->target_name[i] = (char*)calloc(name_len, 1);
		if(gzread(fp, header->target_name[i], name_len) != name_len) {
			goto _bam_read_header_fail;
		}
		if(gzread(fp, &header->target_len[i], 4) != 4) {
			goto _bam_read_header_fail;
		}
	}
	return(header);

_bam_read_header_fail:
	bam_header_destroy(header);
	return(NULL);
}

/* end of bamlite.c */

/* bseq.c */
#define BSEQ_MGN			( 64 )			/* buffer margin length */

/**
 * @struct bseq_file_t
 * @brief context container
 */
typedef struct {
	gzFile fp;
	bam_header_t *bh;
	uint8_t *p, *base, *tail;				/* input seq buffer */
	uint64_t size, acc;						/* size equal to batch_size */
	uint16_t *tags;
	uint32_t l_tags, n_seq, min_len;
	uint8_t is_eof, delim, keep_qual, keep_comment, state;
} bseq_file_t;

/**
 * @struct bseq_tag_t
 * @brief sam optional tag container
 */
typedef struct {
	uint32_t size;							/* size of data array */
	char tag[2], type;
	uint8_t flag;							/* 1 if only if contained in primary */
	uint8_t data[];
} bseq_tag_t;

/**
 * @struct bseq_t
 * @brief record container
 */
typedef struct {
	uint32_t l_seq;							/* sequence length */
	uint16_t l_name, n_tag;					/* name length, #tags */
	char *name;								/* pointers */
	uint8_t *seq, *qual, *tag;
	void *reg;								/* reserved for alignment result */
} bseq_t;
_static_assert(sizeof(bseq_t) == 48);
typedef struct { size_t n, m; bseq_t *a; } bseq_v;

/**
 * @fn bseq_search_tag
 * @brief search 2-byte sam tag (t2:t1) from tag array by linear probing
 */
static _force_inline
uint64_t bseq_search_tag(uint32_t l_tags, uint16_t const *tags, uint16_t t1, uint16_t t2)
{
	v32i16_t tv = _set_v32i16((t2<<8) | t1);
	for(uint64_t i = 0; i < ((l_tags + 0x1f) & ~0x1f); i+=0x20) {
		if(((v32_masku_t){ .mask = _mask_v32i16(_eq_v32i16(_loadu_v32i16(&tags[i]), tv)) }).all != 0) {
			return(1);
		}
	}
	return(0);
}

/**
 * @fn bseq_open
 */
static _force_inline
bseq_file_t *bseq_open(
	char const *fn,
	uint64_t batch_size,					/* buffer (block) size */
	uint32_t keep_qual,						/* 1 to keep quality string */
	uint32_t min_len,						/* minimum length cutoff (to filter out short seqs) */
	uint32_t l_tags,
	uint16_t const *tags)					/* tags to be preserved */
{
	int c;
	set_info(0, "[bseq_open] initialize bseq object");

	/* open stream */
	gzFile f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if(f == NULL) { return(NULL); }

	/* create instance */
	bseq_file_t *fp = (bseq_file_t*)calloc(1, sizeof(bseq_file_t));
	*fp = (bseq_file_t){
		.fp = f,
		.keep_qual = keep_qual,
		.min_len = min_len
	};

	/* determine file type; allow some invalid spaces at the head */
	for(uint64_t i = 0; i < 4; i++) {
		if((c = gzgetc(fp->fp)) == 'B') {	/* test bam signature */
			gzungetc(c, fp->fp); fp->bh = bam_read_header(fp->fp); break;
		} else if(c == '>' || c == '@') {	/* test fasta/q delimiter */
			gzungetc(c, fp->fp); fp->delim = c; break;
		}
	}
	if(!fp->bh && !fp->delim) {
		free(fp);
		return(NULL);
	}

	/* init buffer */
	fp->tail = fp->p = (fp->base = malloc((fp->size = batch_size) + 2*BSEQ_MGN)) + batch_size;
	memset(fp->base, 0, BSEQ_MGN);
	memset(fp->base + batch_size + BSEQ_MGN, 0, BSEQ_MGN);

	/* init tags */
	fp->tags = calloc(((fp->l_tags = l_tags) + 0x1f) & ~0x1f, sizeof(uint16_t));
	if(l_tags && tags) {
		memcpy(fp->tags, tags, l_tags * sizeof(uint16_t));
	}
	fp->keep_comment = bseq_search_tag(fp->l_tags, fp->tags, 'C', 'O');

	return(fp);
}

/**
 * @fn bseq_close
 * @brief returns #seqs read
 */
static _force_inline
uint32_t bseq_close(bseq_file_t *fp)
{
	if(fp == NULL) { return(0); }
	uint32_t n_seq = fp->n_seq;
	gzclose(fp->fp);
	bam_header_destroy(fp->bh);
	free(fp->base);
	free(fp->tags);
	free(fp);
	return(n_seq);
}

/**
 * @fn bseq_save_tags
 * @brief save bam tag to buffer (must have enough space), returns #tags saved
 */
static _force_inline
uint64_t bseq_save_tags(
	uint32_t l_tags, uint16_t const *tags,	/* tags */
	uint32_t l_arr, uint8_t const *arr,		/* src buffer pointer and len */
	uint8_v *restrict mem)					/* dst buffer pointer */
{
	static uint8_t const tag_size[32] = {
		0, 1, 0xfe, 1,  0, 0, 4, 0,  0xfe, 4, 0, 0,  0, 0, 0, 0,
		0, 0, 0, 2,     0, 0, 0, 0,  0, 0, 0xff, 0,  0, 0, 0, 0
	};
	uint64_t n_tag = 0;
	uint8_t *q = &mem->a[mem->n];				/* memory region must be already reserved */
	uint8_t const *p = arr, *tail = arr + l_arr;
	while(p < tail) {
		/* first check if the tag must be saved */
		uint64_t keep = bseq_search_tag(l_tags, tags, p[0], p[1]), size;
		n_tag += keep != 0;

		/* test the size */
		if((size = tag_size[p[2] & 0x1f]) == 0xff) {
			/* string, read until '\0' */
			if(keep != 0) {
				while((*q++ = *p++)) {}		/* save */
			} else {
				while(*p++) {}				/* skip */
			}
		} else {
			/* others */
			uint64_t len = (size == 0xfe)
				? 8 + tag_size[p[3] & 0x1f] * *((uint32_t *)&p[4])	/* array */
				: 3 + size;											/* single value */
			if(keep != 0) {
				memcpy(q, p, len); q += len;
			}
			p += len;
		}
	}
	// return(p - arr);

	mem->n = q - mem->a;					/* write back #elems */
	return(n_tag);
}

/**
 * @fn bseq_read_bam
 */
static _force_inline
uint64_t bseq_read_bam(
	bseq_file_t const *fp,
	uint64_t blk_size,
	bseq_v *restrict seq,					/* sequence metadata array */
	uint8_v *restrict mem)					/* block buffer */
{
	uint8_t *sname, *sseq, *squal, *stag;
	uint32_t l_tag;

	/* test header and length */
	bam_core_t const *c = (bam_core_t const *)fp->p;
	if(c->flag & 0x900) { return(0); }							/* skip supp / secondary */
	if((uint32_t)c->l_qseq < fp->min_len) { return(0); }		/* skip short reads */

	/* extract pointers */
	sname = fp->p + sizeof(bam_core_t);
	sseq = sname + c->l_qname + sizeof(uint32_t) * c->n_cigar;
	squal = sseq + (c->l_qseq + 1) / 2;
	stag = squal + c->l_qseq;
	l_tag = ((uint8_t *)fp->p + blk_size) - stag;				/* FIXME: tag section length */

	/* reserve memory */
	uint64_t req_size = (
		  c->l_qname + 1 + c->l_qseq + 1 						/* name and seq with '\0' */
		+ ((fp->keep_qual && *squal != 0xff)					/* quality string if available */
			? c->l_qseq
			: 0) + 1
		+ (fp->l_tags? l_tag : 0) + 1							/* tag section */
	);
	bseq_t *s = kv_pushp(bseq_t, *seq);
	kv_reserve(uint8_t, *mem, mem->n + req_size);

	/* set lengths */
	s->l_seq = c->l_qseq;
	s->l_name = c->l_qname - 1;									/* remove tail '\0' */
	// s->l_tag = fp->l_tags? l_tag : 0;

	/* copy name */
	s->name = (char *)mem->n;
	memcpy(mem->a + mem->n, sname, c->l_qname); mem->n += c->l_qname;

	/* copy seq */
	s->seq = (uint8_t *)mem->n;
	if(c->flag & 0x10) {
		for(int64_t i = c->l_qseq-1; i >= 0; i--) {
			mem->a[mem->n++] = "\x0\x8\x4\x0\x2\x0\x0\x0\x1"[0x0f & (sseq[i>>1]>>((~i&0x01)<<2))];
		}
	} else {
		for(int64_t i = 0; i < c->l_qseq; i++) {
			mem->a[mem->n++] = 0x0f & (sseq[i>>1]>>((~i&0x01)<<2));
		}
	}
	mem->a[mem->n++] = '\0';

	/* copy qual if available */
	s->qual = (uint8_t *)mem->n;
	if(fp->keep_qual && *squal != 0xff) {
		if(c->flag & 0x10) {
			for(int64_t i = c->l_qseq-1; i >= 0; i--) {
				mem->a[mem->n++] = squal[i] + 33;
			}
		} else {
			for(int64_t i = 0; i < c->l_qseq; i++) {
				mem->a[mem->n++] = squal[i] + 33;
			}
		}
	}
	mem->a[mem->n++] = '\0';

	/* save tags */
	s->tag = (uint8_t *)mem->n;
	if(fp->l_tags && l_tag) {
		// mem->n += bseq_save_tags(fp->l_tags, fp->tags, l_tag, stag, mem->a + mem->n);
		s->n_tag = bseq_save_tags(fp->l_tags, fp->tags, l_tag, stag, mem);
	}
	mem->a[mem->n++] = '\0';
	return(0);
}

/**
 * SIMD string parsing macros
 */
#define _match(_v1, _v2)	( ((v32_masku_t){ .mask = _mask_v32i8(_eq_v32i8(_v1, _v2)) }).all )
#define _strip(_p, _t, _v) ({ \
	v32i8_t _r = _loadu_v32i8(_p); \
	ZCNT_RESULT uint64_t _l = MIN2(tzcnt(~((uint64_t)_match(_r, _v))), _t - _p); \
	uint64_t _len = _l; \
	_p += _len; _len; \
})
#define _readline(_p, _t, _q, _dv, _op) ({ \
	uint64_t _m1, _m2; \
	uint64_t _len; \
	v32i8_t const _lv = _set_v32i8('\n'); \
	do { \
		v32i8_t _r = _loadu_v32i8(_p), _s = _op(_r); _storeu_v32i8(_q, _s); \
		_m1 = _match(_r, _dv); _m2 = _match(_r, _lv); \
		ZCNT_RESULT uint64_t _l = MIN2(tzcnt(_m1 | _m2), _t - _p); \
		_len = _l; _p += 32; _q += 32; \
	} while(_len >= 32); \
	_p += _len - 32; _q += _len - 32; _q -= _q[-1] == 0x0f; _m1>>_len; \
})
#define _skipline(_p, _t) ({ \
	uint64_t _m; \
	uint64_t _len; \
	uint8_t const *_b = _p; \
	v32i8_t const _lv = _set_v32i8('\n'); \
	do { \
		v32i8_t _r = _loadu_v32i8(_p); \
		_m = _match(_r, _lv); \
		ZCNT_RESULT uint64_t _l = MIN2(tzcnt(_m), _t - _p); _len = _l; _p += 32; \
	} while(_len >= 32); \
	_p += _len - 32; _p - _b; \
})
#define _beg(_q, _b)		( (uint8_t *)(_q - _b) )
#define _term(_q, _b, _ofs) ({ \
	uint64_t _len = (uint64_t)(_q - &(_b)[(uint64_t)(_ofs)]); \
	*_q++ = '\0'; _len; \
})

/**
 * @fn bseq_read_sam
 */
static _force_inline
uint64_t bseq_read_sam(
	bseq_file_t *fp,
	bseq_v *seq,
	uint8_v *mem)
{
	return(0);
}

/**
 * @fn bseq_read_fasta
 * @brief returns 0 when correctly finished, 1 when buffer starved, >=2 when broken
 */
static _force_inline
uint64_t bseq_read_fasta(
	bseq_file_t *restrict fp,
	bseq_v *restrict seq,			/* src */
	uint8_v *restrict mem)			/* dst, must have enough space (e.g. 2 * buffer) */
{
	#define _id(x)			(x)
	#define _trans(x)		( _shuf_v32i8(cv, _and_v32i8(fv, x)) )
	/* keep them on registers */
	v32i8_t const dv = _set_v32i8(fp->delim == '@'? '+' : fp->delim);
	v32i8_t const sv = _set_v32i8(' '), lv = _set_v32i8('\n'), fv = _set_v32i8(0xf);
	v32i8_t const cv = _from_v16i8_v32i8(_seta_v16i8(0,0,15,0,0,0,0,0,4,0,8,8,2,0,1,0));

	bseq_t *s = &seq->a[seq->n-1];
	uint8_t *p = fp->p, *q = &mem->a[mem->n];
	uint8_t const *t = fp->tail;
	uint64_t m, acc, lim, ret = 1;

	debug("enter, state(%u), p(%p), t(%p), eof(%u)", fp->state, p, t, fp->is_eof);
	if(p >= t) { return(1); }		/* refill needed */
	switch (fp->state) {
		case 0:						/* idle */
			debug("state(%u)", fp->state);
			if(*p++ != fp->delim) { return(2); }/* broken */
			s = kv_pushp(bseq_t, *seq);			/* create new sequence */
			fp->state = 1;						/* transition to spaces between delim and name */
		case 1:
			debug("state(%u)", fp->state);
			_strip(p, t, sv);
			if(_unlikely(p >= t)) { goto _refill; }
			s->name = (char*)_beg(q, mem->a);
			fp->state = 2;
		case 2:
			debug("state(%u)", fp->state);
			m = _readline(p, t, q, sv, _id);
			if(_unlikely(p >= t)) { goto _refill; }
			p++;								/* skip '\n' or ' ' */
			s->l_name = _term(q, mem->a, s->name);
			s->tag = _beg(q, mem->a);
			s->n_tag = m & 0x01;				/* set n_tag if comment line found */
			if((m & 0x01) == 0) { goto _seq_head; }
			fp->state = 3;
		case 3:
			debug("state(%u)", fp->state);
			_strip(p, t, sv);
			if(_unlikely(p >= t)) { goto _refill; }
			*q++ = 'C'; *q++ = 'O'; *q++ = 'Z';
			fp->state = 4;
		case 4:									/* parsing comment */
			debug("state(%u)", fp->state);
			_readline(p, t, q, lv, _id);
			if(_unlikely(p >= t)) { goto _refill; }	/* refill needed, comment continues */
			p++;								/* skip '\n' */
			while(q[-1] == ' ') { q--; }		/* strip spaces */
			if(!fp->keep_comment) q = mem->a + (uint64_t)s->tag;
		_seq_head:
			// s->l_tag = _term(q, mem->a, s->tag);
			_term(q, mem->a, s->tag);
			s->seq = _beg(q, mem->a);
			fp->state = 5;
		case 5:									/* parsing seq */
			debug("state(%u)", fp->state);
			while(1) {
				m = _readline(p, t, q, dv, _trans);
				if(_unlikely(p >= t)) { m |= fp->is_eof; break; }
				if(m & 0x01) { break; }
				p++;							/* skip '\n' */
			}
			if((m & 0x01) == 0) { goto _refill; }
			s->l_seq = _term(q, mem->a, s->seq);
			s->qual = _beg(q, mem->a);
			if(fp->delim == '>' || _unlikely(p >= t)) { goto _qual_tail; }
			fp->state = 6;
		case 6:
			debug("state(%u)", fp->state);
			_skipline(p, t);
			if(_unlikely(p >= t)) { goto _refill; }
			p++;
			fp->state = 7; fp->acc = 0;
		case 7:									/* parsing qual */
			debug("state(%u)", fp->state);
			acc = fp->acc, lim = s->l_seq;
			if(fp->keep_qual) {
				while(1) {
					uint8_t const *b = q; _readline(p, t, q, lv, _id); acc += q - b;
					if(_unlikely(p >= t)) { fp->acc = acc; goto _refill; }
					if(acc >= lim) { break; }
					p++;						/* skip '\n' */
				}
			} else {
				while(1) {
					acc += _skipline(p, t);
					if(_unlikely(p >= t)) { fp->acc = acc; goto _refill; }
					if(acc >= lim) { break; }
					p++;						/* skip '\n' */
				}
			}
			fp->state = 8;
		case 8:
			debug("state(%u)", fp->state);
			if(_unlikely(p >= t)) { goto _refill; }
			_strip(p, t, lv);
		_qual_tail:
			_term(q, mem->a, s->qual);
			break;
		default:								/* invalid state */
			return(2);							/* broken */
	}
	debug("break loop, state(%u), eof(%u)", fp->state, fp->is_eof)
	ret = p >= t; fp->state = 0;				/* back to idle */
	if((uint32_t)s->l_seq < fp->min_len) {
		seq->n--;
		q = mem->a + (uint64_t)s->name;
	}
_refill:
	debug("return, p(%p), t(%p)", p, t);
	fp->p = p; mem->n = q - mem->a;
	return(ret);

	#undef _id
	#undef _trans
}
#undef _match
#undef _strip
#undef _readline
#undef _skipline
#undef _beg
#undef _term

/**
 * @fn bseq_read
 */
static _force_inline
bseq_t *bseq_read(
	bseq_file_t *restrict fp,
	uint32_t *restrict n,
	void *restrict *base,
	uint64_t *restrict size)
{
	uint8_v mem = { 0 };
	bseq_v seq = { 0 };
	static uint8_t const margin[64] = { 0 };

	set_info(0, "[bseq_read] read sequence block from file");

	/* check the context is valid */
	if(fp->is_eof && fp->p >= fp->tail) { return(NULL); }

	/* reserve mem */
	kv_reserve(uint8_t, mem, 2*fp->size);
	kv_pushm(uint8_t, mem, margin, 64);
	if(fp->bh) {
		/* bam */
		while(mem.n < fp->size) {
			uint32_t size;
			if(gzread(fp->fp, &size, 4) != 4) {
				fp->is_eof = 1;
				break;
			}
			if(fp->size < size) {
				fp->tail = (fp->base = fp->p = realloc(fp->base, fp->size = size)) + size;
			}
			if(gzread(fp->fp, fp->p, size) != size) {
				seq.n = 0;
				fp->is_eof = 2;
				break;
			}
			bseq_read_bam(fp, size, &seq, &mem);
		}
	} else {
		/* fasta/q */
		while(mem.n < fp->size + 64) {
			uint64_t ret;
			while((ret = bseq_read_fasta(fp, &seq, &mem)) == 1) {
				if(fp->is_eof) { goto _tail; }				/* correctly reached the end */

				/* reload buffer */
				fp->p = fp->base + BSEQ_MGN;
				fp->tail = fp->p + gzread(fp->fp, fp->p, fp->size);
				memset(fp->tail, 0, BSEQ_MGN);				/* fill margin */

				if(fp->tail < fp->base + BSEQ_MGN + fp->size) {
					fp->is_eof = 1;							/* reached the tail */
				}
				if(mem.n + 2*fp->size > mem.m) {
					mem.a = realloc(mem.a, mem.m *= 2);		/* expand mem */
				}
			}
			if(ret > 1) {									/* error occurred */
				seq.n = 0;
				fp->is_eof = 2;
				goto _tail;
			}
		}
	_tail:;
	}

	kv_pushm(uint8_t, mem, margin, 64);

	/* adjust pointers */
	for(uint64_t i = 0; i < seq.n; i++) {
		seq.a[i].name += (ptrdiff_t)mem.a;
		seq.a[i].seq += (ptrdiff_t)mem.a;
		seq.a[i].qual += (ptrdiff_t)mem.a;
		seq.a[i].tag += (ptrdiff_t)mem.a;
	}

	/* seq.n == 0 indicates error */
	if(seq.n == 0) {
		free(mem.a);
		mem.a = 0;
		fp->is_eof = 1;
	}

	fp->n_seq += seq.n;
	*n = seq.n;
	*base = (void *)mem.a;
	*size = mem.n;
	return(seq.a);
}

unittest( .name = "bseq.fasta" ) {
	char const *filename = "./minialign.unittest.bseq.tmp";
	char const *content =
		">test0\nAAAA\n"
		"> test1\nATAT\nCGCG\n\n"
		">  test2\n\nAAAA\n"
		">test3 comment comment  \nACGT\n\n";

	FILE *fp = fopen(filename, "w");
	assert(fp != NULL);
	fwrite(content, 1, strlen(content), fp);
	fclose(fp);

	uint16_t const tags[1] = { ((uint16_t)'C') | (((uint16_t)'O')<<8) };
	bseq_file_t *b = bseq_open(filename, 64, 1, 0, 1, tags);
	assert(b != NULL);

	uint32_t n_seq = 0;
	void *ptr = NULL;
	uint64_t size = 0;
	bseq_t *s = bseq_read(b, &n_seq, &ptr, &size);

	assert(s != NULL);
	assert(n_seq == 4, "n_seq(%u)", n_seq);
	assert(ptr != NULL);
	assert(size > 0, "size(%lu)", size);

	assert(s[0].l_name == 5, "l_name(%u)", s[0].l_name);
	assert(s[0].l_seq == 4, "l_seq(%u)", s[0].l_seq);
	assert(s[0].n_tag == 0, "n_tag(%u)", s[0].n_tag);
	assert(strcmp((const char*)s[0].name, "test0") == 0, "name(%s)", s[0].name);
	assert(strcmp((const char*)s[0].seq, "\x1\x1\x1\x1") == 0, "seq(%s)", s[0].seq);
	assert(strcmp((const char*)s[0].qual, "") == 0, "qual(%s)", s[0].qual);
	assert(strcmp((const char*)s[0].tag, "") == 0, "tag(%s)", s[0].tag);

	assert(s[1].l_name == 5, "l_name(%u)", s[1].l_name);
	assert(s[1].l_seq == 8, "l_seq(%u)", s[1].l_seq);
	assert(s[1].n_tag == 0, "n_tag(%u)", s[1].n_tag);
	assert(strcmp((const char*)s[1].name, "test1") == 0, "name(%s)", s[1].name);
	assert(strcmp((const char*)s[1].seq, "\x1\x8\x1\x8\x2\x4\x2\x4") == 0, "seq(%s)", s[1].seq);
	assert(strcmp((const char*)s[1].qual, "") == 0, "qual(%s)", s[1].qual);
	assert(strcmp((const char*)s[1].tag, "") == 0, "tag(%s)", s[1].tag);

	assert(s[2].l_name == 5, "l_name(%u)", s[2].l_name);
	assert(s[2].l_seq == 4, "l_seq(%u)", s[2].l_seq);
	assert(s[2].n_tag == 0, "n_tag(%u)", s[2].n_tag);
	assert(strcmp((const char*)s[2].name, "test2") == 0, "name(%s)", s[2].name);
	assert(strcmp((const char*)s[2].seq, "\x1\x1\x1\x1") == 0, "seq(%s)", s[2].seq);
	assert(strcmp((const char*)s[2].qual, "") == 0, "qual(%s)", s[2].qual);
	assert(strcmp((const char*)s[2].tag, "") == 0, "tag(%s)", s[2].tag);

	assert(s[3].l_name == 5, "l_name(%u)", s[3].l_name);
	assert(s[3].l_seq == 4, "l_seq(%u)", s[3].l_seq);
	assert(s[3].n_tag == 1, "n_tag(%u)", s[3].n_tag);
	assert(strcmp((const char*)s[3].name, "test3") == 0, "name(%s)", s[3].name);
	assert(strcmp((const char*)s[3].seq, "\x1\x2\x4\x8") == 0, "seq(%s)", s[3].seq);
	assert(strcmp((const char*)s[3].qual, "") == 0, "qual(%s)", s[3].qual);
	assert(strcmp((const char*)s[3].tag, "COZcomment comment") == 0, "tag(%s)", s[3].tag);

	n_seq = bseq_close(b);
	assert(n_seq == 4);
	remove(filename);
}

unittest( .name = "bseq.fastq" ) {
	char const *filename = "./minialign.unittest.bseq.tmp";
	char const *content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\n12+3\n+123\n"
		"@  test2\n\nAAAA\n+  test2\n\n\n12@3\n\n"
		"@test3  comment comment   \nACGT\n\n+ test3\n@123";

	FILE *fp = fopen(filename, "w");
	assert(fp != NULL);
	fwrite(content, 1, strlen(content), fp);
	fclose(fp);

	uint16_t const tags[1] = { ((uint16_t)'C') | (((uint16_t)'O')<<8) };
	bseq_file_t *b = bseq_open(filename, 256, 1, 0, 1, tags);
	assert(b != NULL);

	uint32_t n_seq = 0;
	void *ptr = NULL;
	uint64_t size = 0;
	bseq_t *s = bseq_read(b, &n_seq, &ptr, &size);

	assert(s != NULL);
	assert(n_seq == 4, "n_seq(%u)", n_seq);
	assert(ptr != NULL);
	assert(size > 0, "size(%lu)", size);

	assert(s[0].l_name == 5, "l_name(%u)", s[0].l_name);
	assert(s[0].l_seq == 4, "l_seq(%u)", s[0].l_seq);
	assert(s[0].n_tag == 0, "n_tag(%u)", s[0].n_tag);
	assert(strcmp((const char*)s[0].name, "test0") == 0, "name(%s)", s[0].name);
	assert(strcmp((const char*)s[0].seq, "\x1\x1\x1\x1") == 0, "seq(%s)", s[0].seq);
	assert(strcmp((const char*)s[0].qual, "NNNN") == 0, "qual(%s)", s[0].qual);
	assert(strcmp((const char*)s[0].tag, "") == 0, "tag(%s)", s[0].tag);

	assert(s[1].l_name == 5, "l_name(%u)", s[1].l_name);
	assert(s[1].l_seq == 8, "l_seq(%u)", s[1].l_seq);
	assert(s[1].n_tag == 0, "n_tag(%u)", s[1].n_tag);
	assert(strcmp((const char*)s[1].name, "test1") == 0, "name(%s)", s[1].name);
	assert(strcmp((const char*)s[1].seq, "\x1\x8\x1\x8\x2\x4\x2\x4") == 0, "seq(%s)", s[1].seq);
	assert(strcmp((const char*)s[1].qual, "12+3+123") == 0, "qual(%s)", s[1].qual);
	assert(strcmp((const char*)s[1].tag, "") == 0, "tag(%s)", s[1].tag);

	assert(s[2].l_name == 5, "l_name(%u)", s[2].l_name);
	assert(s[2].l_seq == 4, "l_seq(%u)", s[2].l_seq);
	assert(s[2].n_tag == 0, "n_tag(%u)", s[2].n_tag);
	assert(strcmp((const char*)s[2].name, "test2") == 0, "name(%s)", s[2].name);
	assert(strcmp((const char*)s[2].seq, "\x1\x1\x1\x1") == 0, "seq(%s)", s[2].seq);
	assert(strcmp((const char*)s[2].qual, "12@3") == 0, "qual(%s)", s[2].qual);
	assert(strcmp((const char*)s[2].tag, "") == 0, "tag(%s)", s[2].tag);

	assert(s[3].l_name == 5, "l_name(%u)", s[3].l_name);
	assert(s[3].l_seq == 4, "l_seq(%u)", s[3].l_seq);
	assert(s[3].n_tag == 1, "n_tag(%u)", s[3].n_tag);
	assert(strcmp((const char*)s[3].name, "test3") == 0, "name(%s)", s[3].name);
	assert(strcmp((const char*)s[3].seq, "\x1\x2\x4\x8") == 0, "seq(%s)", s[3].seq);
	assert(strcmp((const char*)s[3].qual, "@123") == 0, "qual(%s)", s[3].qual);
	assert(strcmp((const char*)s[3].tag, "COZcomment comment") == 0, "tag(%s)", s[3].tag);

	n_seq = bseq_close(b);
	assert(n_seq == 4);
	remove(filename);
}

unittest( .name = "bseq.fastq.skip" ) {
	char const *filename = "./minialign.unittest.bseq.tmp";
	char const *content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\n12+3\n+123\n"
		"@  test2\n\nAAAA\n+  test2\n\n\n12@3\n\n"
		"@test3  comment comment   \nACGT\n\n+ test3\n@123";

	FILE *fp = fopen(filename, "w");
	assert(fp != NULL);
	fwrite(content, 1, strlen(content), fp);
	fclose(fp);

	uint16_t const tags[1] = { ((uint16_t)'C') | (((uint16_t)'O')<<8) };
	bseq_file_t *b = bseq_open(filename, 256, 0, 0, 1, tags);
	assert(b != NULL);

	uint32_t n_seq = 0;
	void *ptr = NULL;
	uint64_t size = 0;
	bseq_t *s = bseq_read(b, &n_seq, &ptr, &size);

	assert(s != NULL);
	assert(n_seq == 4, "n_seq(%u)", n_seq);
	assert(ptr != NULL);
	assert(size > 0, "size(%lu)", size);

	assert(s[0].l_name == 5, "l_name(%u)", s[0].l_name);
	assert(s[0].l_seq == 4, "l_seq(%u)", s[0].l_seq);
	assert(s[0].n_tag == 0, "n_tag(%u)", s[0].n_tag);
	assert(strcmp((const char*)s[0].name, "test0") == 0, "name(%s)", s[0].name);
	assert(strcmp((const char*)s[0].seq, "\x1\x1\x1\x1") == 0, "seq(%s)", s[0].seq);
	assert(strcmp((const char*)s[0].qual, "") == 0, "qual(%s)", s[0].qual);
	assert(strcmp((const char*)s[0].tag, "") == 0, "tag(%s)", s[0].tag);

	assert(s[1].l_name == 5, "l_name(%u)", s[1].l_name);
	assert(s[1].l_seq == 8, "l_seq(%u)", s[1].l_seq);
	assert(s[1].n_tag == 0, "n_tag(%u)", s[1].n_tag);
	assert(strcmp((const char*)s[1].name, "test1") == 0, "name(%s)", s[1].name);
	assert(strcmp((const char*)s[1].seq, "\x1\x8\x1\x8\x2\x4\x2\x4") == 0, "seq(%s)", s[1].seq);
	assert(strcmp((const char*)s[1].qual, "") == 0, "qual(%s)", s[1].qual);
	assert(strcmp((const char*)s[1].tag, "") == 0, "tag(%s)", s[1].tag);

	assert(s[2].l_name == 5, "l_name(%u)", s[2].l_name);
	assert(s[2].l_seq == 4, "l_seq(%u)", s[2].l_seq);
	assert(s[2].n_tag == 0, "n_tag(%u)", s[2].n_tag);
	assert(strcmp((const char*)s[2].name, "test2") == 0, "name(%s)", s[2].name);
	assert(strcmp((const char*)s[2].seq, "\x1\x1\x1\x1") == 0, "seq(%s)", s[2].seq);
	assert(strcmp((const char*)s[2].qual, "") == 0, "qual(%s)", s[2].qual);
	assert(strcmp((const char*)s[2].tag, "") == 0, "tag(%s)", s[2].tag);

	assert(s[3].l_name == 5, "l_name(%u)", s[3].l_name);
	assert(s[3].l_seq == 4, "l_seq(%u)", s[3].l_seq);
	assert(s[3].n_tag == 18, "n_tag(%u)", s[3].n_tag);
	assert(strcmp((const char*)s[3].name, "test3") == 0, "name(%s)", s[3].name);
	assert(strcmp((const char*)s[3].seq, "\x1\x2\x4\x8") == 0, "seq(%s)", s[3].seq);
	assert(strcmp((const char*)s[3].qual, "") == 0, "qual(%s)", s[3].qual);
	assert(strcmp((const char*)s[3].tag, "COZcomment comment") == 0, "tag(%s)", s[3].tag);

	n_seq = bseq_close(b);
	assert(n_seq == 4);
	remove(filename);
}
/* end of bseq.c */

/* sketch.c */
/**
 * @fn hash64
 */
#if USE_CRC32_HASH != 0
#define hash64(k0, k1, mask)		( (_mm_crc32_u64((k0), (k0)) ^ (k1)) & (mask) )
#else
static _force_inline
uint64_t hash64(uint64_t key, uint64_t unused, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; /* key = (key << 21) - key - 1; */
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; /* key * 265 */
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; /* key * 21 */
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}
#endif

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param p      minimizers; p->a[i].x is the 2k-bit hash value;
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
#if 0
static void mm_sketch(uint8_t const *seq4, uint32_t len, uint32_t w, uint32_t k, uint32_t rid, v4u32_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, k0 = 0, k1 = 0;
	uint32_t i, j, l, buf_pos, min_pos;
	v4u32_t buf[w], min = { .u64 = { UINT64_MAX, UINT64_MAX } };

	assert(len > 0 && w > 0 && k > 0);
	memset(buf, 0xff, w * sizeof(v4u32_t));

	for(i = l = buf_pos = min_pos = 0; i < len; ++i) {
		uint64_t _c = seq4[i];
		v4u32_t info = { .u64 = { UINT64_MAX, UINT64_MAX } };
		if(_c != 0) { // not an ambiguous base
			uint64_t c = 0x03 & ((_c>>1) - (_c>>3));
			k0 = (k0 << 2 | c) & mask;           // forward k-mer
			k1 = (k1 >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if(k0 == k1) continue; // skip "symmetric k-mers" as we don't know it strand
			if(++l >= k) {
				info.u64[0] = hash64(k0 < k1? k0 : k1, mask);
				info.u32[2] = k0 < k1? i : ~i; info.u32[3] = rid;
			}
		} else l = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if(l == w + k - 1) { // special case for the first window - because identical k-mers are not stored yet
			for(j = buf_pos + 1; j < w; ++j)
				if(min.u64[0] == buf[j].u64[0] && buf[j].u64[1] != min.u64[1]) kv_push(v4u32_t, *p, buf[j]);
			for(j = 0; j < buf_pos; ++j)
				if(min.u64[0] == buf[j].u64[0] && buf[j].u64[1] != min.u64[1]) kv_push(v4u32_t, *p, buf[j]);
		}
		if(info.u64[0] <= min.u64[0]) { // a new minimum; then write the old min
			if(l >= w + k) kv_push(v4u32_t, *p, min);
			min = info, min_pos = buf_pos;
		} else if(buf_pos == min_pos) { // old min has moved outside the window
			if(l >= w + k - 1) kv_push(v4u32_t, *p, min);
			for(j = buf_pos + 1, min.u64[0] = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if(min.u64[0] >= buf[j].u64[0]) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for(j = 0; j <= buf_pos; ++j)
				if(min.u64[0] >= buf[j].u64[0]) min = buf[j], min_pos = j;
			if(l >= w + k - 1) { // write identical k-mers
				for(j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if(min.u64[0] == buf[j].u64[0] && min.u64[1] != buf[j].u64[1]) kv_push(v4u32_t, *p, buf[j]);
				for(j = 0; j <= buf_pos; ++j)
					if(min.u64[0] == buf[j].u64[0] && min.u64[1] != buf[j].u64[1]) kv_push(v4u32_t, *p, buf[j]);
			}
		}
		if(++buf_pos == w) buf_pos = 0;
	}
	if(min.u64[0] != UINT64_MAX)
		kv_push(v4u32_t, *p, min);
}
#else
static _force_inline
void mm_sketch(
	uint8_t const *seq4, uint32_t _len,
	uint32_t w, uint32_t k,
	uint32_t rid,
	v4u32_v *p)
{
	uint64_t const kk = k - 1, shift1 = 2*kk, mask = (1ULL<<2*k) - 1;
	uint64_t const len = ((uint64_t)rid<<32) | _len;	/* keep rid in the upper 32bit */
	v4u32_t buf[16];									/* w < 16 must be guaranteed */

	#define _push_kmer(_c) { \
		uint64_t _t = 0x03 & (((_c)>>1) - ((_c)>>3)); \
		k0 = (k0 << 2 | _t) & mask; \
		k1 = (k1 >> 2) | ((3ULL^_t) << shift1); \
	}

	uint64_t i = (uint64_t)rid<<32; seq4 -= (ptrdiff_t)i; i--;
	v4u32_t *q = p->a+p->n, *t = p->a+p->m;
	do {
		/* extend buffer if needed */
		if(q + 64 > t) {
			p->n = q-p->a;
			p->m = MAX2(256, p->m*2);
			p->a = realloc(p->a, p->m*sizeof(v4u32_t));
			q = p->a+p->n;
			t = p->a+p->m;
		}

		uint64_t l = 0, k0 = 0, k1 = 0, min = UINT64_MAX, min_pos = 0;

		/* head */
		while(l < kk) {
			uint64_t c;
			if((c = seq4[++i]) == 0) { goto _loop_tail; }
			_push_kmer(c);
			if(k0 == k1) { continue; }
			buf[l&0x0f] = (v4u32_t){ .u64 = { UINT64_MAX, 0 } };
			l++;
		}
		while(l < kk + w) {
			uint64_t c, h;
			if((c = seq4[++i]) == 0) { goto _loop_tail; }
			_push_kmer(c);
			if(k0 == k1) { continue; }
			buf[l&0x0f] = (v4u32_t){
				.u64 = {
					h = hash64(k0 < k1? k0 : k1, k0 < k1 ? k1 : k0, mask),
					(k0 < k1? 0 : 0xffffffff) ^ i
				}
			};
			if(h <= min) { min = h; min_pos = l & 0x0f; }
			l++;
		}
		for(uint64_t j = kk; j < kk + w; j++) {
			if(buf[j&0x0f].u64[0] == min) { *q++ = buf[j&0x0f]; }
		}
		q--;

		/* body */
		while(1) {
			uint64_t c, h;
			if((c = seq4[++i]) == 0) { goto _loop_tail; }
			_push_kmer(c);
			if(k0 == k1) { continue; }
			buf[l&0x0f] = (v4u32_t){
				.u64 = {
					h = hash64(k0 < k1? k0 : k1, k0 < k1 ? k1 : k0, mask),
					(k0 < k1? 0 : 0xffffffff) ^ i
				}
			};

			if(h <= min) {
				*q++ = buf[min_pos]; min = h; min_pos = l&0x0f;
			} else if(min_pos == ((l-w)&0x0f)) {
				*q++ = buf[min_pos]; min = UINT64_MAX;

				/* note: the following loops degrades performance due to the unpredictable branches */
				for(uint64_t j = l-w+1; j <= l; j++) {
					if(buf[j&0x0f].u64[0] <= min) {
						min = buf[j&0x0f].u64[0]; min_pos = j&0x0f;
					}
				}
				for(uint64_t j = l-w+1; j <= l; j++) {
					if(buf[j&0x0f].u64[0] == min) {
						*q++ = buf[j&0x0f];
					}
				}
				q--;
			}
			l++;

			/* extend buffer if needed */
			if(q+64 <= t) { continue; }
			p->n = q-p->a;
			p->m = MAX2(128, p->m*2);
			p->a = realloc(p->a, p->m*sizeof(v4u32_t));
			q = p->a+p->n;
			t = p->a+p->m;
		}
	_loop_tail:
		/* flush tail */
		if(min != UINT64_MAX) { *q++ = buf[min_pos]; }
	} while(i < len);
	p->n = q-p->a;
	#undef _push_kmer
	return;
}
#endif
/* end of sketch.c */

/* map.c, options */
#define MM_RG			( 0 )				/* Z: read group */
#define MM_CO			( 1 )				/* Z: comment */
#define MM_NH			( 2 )				/* i: #hits */
#define MM_IH 			( 3 )				/* i: index of the record in the hits */
#define MM_AS			( 4 )				/* i: score */
#define MM_XS			( 5 )				/* i: suboptimal score */
#define MM_NM 			( 6 )				/* i: editdist to the reference */
#define MM_SA			( 7 )				/* Z: supplementary records */

#define MM_AVA			( 0x01ULL<<48 )
#define MM_KEEP_QUAL	( 0x02ULL<<48 )
#define MM_CIRCULAR		( 0x04ULL<<48 )
#define MM_OMIT_REP		( 0x08ULL<<48 )		/* omit secondary records */
#define MM_COMP 		( 0x10ULL<<48 )

#define MM_MAF			( 0x01ULL<<56 )
#define MM_BLAST6		( 0x02ULL<<56 )
#define MM_BLASR1		( 0x03ULL<<56 )
#define MM_BLASR4		( 0x04ULL<<56 )
#define MM_PAF			( 0x05ULL<<56 )
#define MM_MHAP 		( 0x06ULL<<56 )
#define MM_FALCON		( 0x07ULL<<56 )


/**
 * @struct mm_mapopt_t
 * @brief parameter container
 * note: implement loggers
 */
typedef struct mm_mapopt_s {
	/* global */
	uint32_t nth;
	uint64_t flag;

	/* sequence parser params */
	uint32_t rmin, qmin;
	uint32_t base_rid, base_qid;

	/* index params */
	uint32_t k, w, b;
	uint32_t n_frq;
	float frq[MAX_FRQ_CNT + 1];

	/* chain and extension */
	int32_t m, x, gi, ge, xdrop, wlen, glen;
	float min_ratio;
	uint32_t min_score, max_cnt;

	/* buffer sizes */
	uint64_t batch_size, outbuf_size;
	
	/* sam optionals */
	uint16_v tags;
	char *rg_line, *rg_id, *arg_line;

	/* logger */
	double inittime;
	uint32_t verbose;
	int (*log)(struct mm_mapopt_s const *, char, char const *, char const *, ...);
	void *fp;
} mm_mapopt_t;

/**
 * @type mm_fprintf_t
 * @brief log printer
 */
typedef int (*mm_log_t)(mm_mapopt_t const *opt, char level, char const *func, char const *fmt, ...);

/**
 * @fn mm_log_printer
 * @brief 0, 1, 2,... for normal message, 8, 9, 10,... for message with timestamp, 16, 17, 18, ... for without header.
 */
static
int mm_log_printer(
	mm_mapopt_t const *opt,	/* option object */
	char level,				/* 'E' and 'W' for error and warning, 0, 1,... for message */
	char const *func,		/* __func__ must be passed */
	char const *fmt,		/* format string */
	...)
{
	if (level < ' ' && (level & 0x07) > opt->verbose) {
		return(0);
	}

	va_list l;
	va_start(l, fmt);

	FILE *fp = (FILE *)opt->fp;
	int r = 0;
	if(level >= ' ' || (level & 0x10) == 0) {
		if(level >= ' ' || (level & 0x08) == 0) {
			r += fprintf(fp, "[%c::%s] ", level < ' '? 'M' : level, func);
		} else {
			r += fprintf(fp, "[%c::%s::%.3f*%.2f] ",
				level < ' '? 'M' : level,					/* 'E' for error */
				func,										/* function name */
				realtime() - opt->inittime,					/* realtime */
				cputime() / (realtime() - opt->inittime));	/* average cpu usage */
		}
	}
	r += vfprintf(fp, fmt, l);								/* body */
	r += fprintf(fp, "\n");
	va_end(l);
	return(r);
}


/**
 * @fn mm_mapopt_destroy
 */
static _force_inline
void mm_mapopt_destroy(mm_mapopt_t *opt)
{
	free(opt->rg_line);
	free(opt->rg_id);
	free(opt->arg_line);
	free(opt);
	return;
}

/**
 * @fn mm_mapopt_init
 */
static _force_inline
mm_mapopt_t *mm_mapopt_init(void)
{
	set_info(0, "[mm_mapopt_init] initialize mapopt object");
	mm_mapopt_t *opt = calloc(1, sizeof(mm_mapopt_t));
	*opt = (mm_mapopt_t){
		/* -V */ .verbose = 1,
		/* -f, -k, -w, -b, -T */ .k = 15, .w = 16, .b = 14, .flag = 0,
		/* -a, -b, -p, -q, -Y */ .m = 1, .x = 1, .gi = 1, .ge = 1, .xdrop = 50,
		/* -s, -m */ .min_score = 50, .min_ratio = 0.3,
		/* -M */ .max_cnt = 0,
		/* -f */ .n_frq = 3, .frq[0] = 0.05, .frq[1] = 0.01, .frq[2] = 0.001,
		/* -t */ .nth = 1,
		/* -R */ .rg_line = NULL, .rg_id = NULL,

		/* -L, -H */.rmin = 0, .qmin = 0,
		/* -W, -G */.wlen = 7000, .glen = 7000,
		// .hlim = 7000, .llim = 7000, .blim = 0, .elim = 200,
		.batch_size = 512 * 1024,
		.outbuf_size = 512 * 1024,
		.base_rid = 0, .base_qid = 0,

		/* initialized time */ .inittime = realtime(),
		/* logger */ .log = (mm_log_t)mm_log_printer, .fp = (void *)stderr
	};
	return(opt);
}

/**
 * @fn mm_mapopt_check
 * @brief check sanity of the option object and print message if not. returns #errors.
 */
static _force_inline
int mm_mapopt_check(mm_mapopt_t *opt)
{
	#define _opt_assert(_cond, _level, ...) { \
		if((_cond)) { \
			ret++; \
			opt->log(opt, _level, __func__, __VA_ARGS__); \
		} \
	}

	int ret = 0;
	/* k-mer length and window size */
	_opt_assert(opt->k >= 32, 								'E', "k must be inside [1,32).");
	_opt_assert(opt->w >= 16,								'E', "w must be inside [1,16).");

	/* score params */
	_opt_assert(opt->m < 1 || opt->m > 5,					'E', "Match award (-a) must be inside [1,5].");
	_opt_assert(opt->x < 1 || opt->x > 5,					'E', "Mismatch penalty (-b) must be inside [1,5].");
	_opt_assert(opt->gi > 5,								'E', "Gap open penalty (-p) must be inside [0,5].");
	_opt_assert(opt->ge < 1 || opt->ge > 5,					'E', "Gap extension penalty (-q) must be inside [1,5].");
	_opt_assert(opt->xdrop < 10 || opt->xdrop > 100,		'E', "Xdrop cutoff must be inside [10,100].");
	_opt_assert(opt->min_score > INT32_MAX, 				'E', "Minimum alignment score must be > 0.");
	_opt_assert(opt->min_ratio < 0.0 || opt->min_ratio > 1.0,
		'E',
		"Minimum alignment score ratio must be inside [0.0,1.0]."
	);

	/* undesired params */
	_opt_assert(opt->gi == 0 && opt->x == 1 && opt->ge == 1,
		'W',
		"(M,X,Gi,Ge) = (1,1,0,1) has positive expected score for two independent random sequences (thus result in false positives)."
		"Please consider using a more stringent score."
	);
	_opt_assert(opt->gi != 0 && opt->x >= (opt->gi + opt->ge),
		'W',
		"Large mismatch penalty with respect to the gap open/extend penalty may cause SEGV or broken CIGAR. [issue #2]"
	);
	_opt_assert(opt->gi != 0 && opt->m + 2*(opt->gi + opt->ge) > 10,
		'W',
		"Large match award or large gap open/extend penalty may cause SEGV or broken CIGAR. [issue #7]"
	);

	/* occurrence filter params */
	_opt_assert(opt->n_frq > MAX_FRQ_CNT,					'E', "Frequency thresholds must be fewer than 16.");
	for(uint64_t i = 0; i < opt->n_frq; i++) {
		_opt_assert(opt->frq[i] < 0.0 || opt->frq[i] > 1.0 || (i != 0 && opt->frq[i-1] < opt->frq[i]),
			'E',
			"Frequency thresholds must be inside [0.0,1.0] and descending."
		);
	}

	_opt_assert(opt->nth < 1 || opt->nth >= MAX_THREADS,
		'E',
		"Thread counts must be inside [1,%u]. For larger values, recompile is needed.", MAX_THREADS
	);
	_opt_assert(opt->batch_size < 64 * 1024,				'E', "Batch size must be > 64k.");
	_opt_assert(opt->outbuf_size < 64 * 1024,				'E', "Output buffer size must be > 64k.");

	_opt_assert(opt->wlen < 100 || opt->wlen >= 100000,		'E', "window edge length must be inside [100,100000).");
	_opt_assert(opt->glen < 100 || opt->glen >= 10000,		'E', "gap chain length must be inside [100,10000).");
	// _opt_assert(opt->sidx >= 16,							'E', "sidx must be inside [0,16).");
	// _opt_assert(opt->eidx >= 16,							'E', "eidx must be inside [0,16).");
	// _opt_assert(opt->hlim < 100 || opt->hlim >= 100000,		'E', "hlim must be inside [100,100000).");
	// _opt_assert(opt->llim < 100 || opt->llim >= 100000,		'E', "llim must be inside [100,100000).");
	// _opt_assert(opt->blim >= 10000,							'E', "blim must be inside [0,10000).");
	// _opt_assert(opt->elim >= 10000,							'E', "elim must be inside [0,10000).");

	#define _dup(x)	({ char *_p = malloc(strlen(x)+1); memcpy(_p, (x), strlen(x)); _p[strlen(x)] = '\0'; _p; })
	if(opt->flag&(0x01ULL<<MM_RG)) {
		if(opt->rg_line == NULL && opt->rg_id == NULL) {
			opt->rg_line = _dup("@RG\tID:1\tSM:default");
			opt->rg_id = _dup("1");
		}
	}
	#undef _dup

	_opt_assert(!(opt->flag&MM_AVA) && (opt->flag&MM_COMP),	'W', "`-A' flag is only effective in all-versus-all mode. ignored.");
	/* opt->flag &= ~MM_COMP; */
	return ret;

	#undef _opt_assert
}

/* end of map.c, options */

/* index.c */

/**
 * @struct mm_idx_seq_t
 * @brief reference-side sequence container
 */
typedef struct {
	uint32_t l_seq, rid, l_name, circular;		/* lengths, id and flag */
	char *name;									/* pointers */
	uint8_t *seq;
} mm_idx_seq_t;
_static_assert(sizeof(mm_idx_seq_t) == 32);
typedef struct { size_t n, m; mm_idx_seq_t *a; } mm_idx_seq_v;

/**
 * @struct mm_idx_bucket_t
 * @brief first-stage data container of double hash table
 */
typedef struct {
	union {
		kh_t h;						/* hash table indexing _p_ and minimizers appearing once */
		v4u32_v a;					/* (minimizer, position) array, for index construction */
	} w;
	uint64_t *p;					/* position array for minimizers appearing >1 times, size in p[0] */
} mm_idx_bucket_t;
_static_assert(sizeof(mm_idx_bucket_t) == 32);

/**
 * @struct mm_idx_t
 * @brief root context container
 */
typedef struct {
	uint64_t mask;					/* must be (1<<b) - 1 */
	mm_idx_bucket_t *bkt;			/* 2nd stage hash table */
	mm_idx_seq_v s;					/* sequence array */
	uint8_t b, w, k, circular;		/* bucket size (in bits), window and k-mer size, circular flag */
	uint32_t base_rid;				/* base sequence id (for all-versus-all) */
	v4u32_v mem;					/* sequence buckets */
} mm_idx_t;

/**
 * @fn mm_idx_init
 */
static _force_inline
mm_idx_t *mm_idx_init(
	uint32_t w, uint32_t k, uint32_t b,
	uint32_t base_rid,
	uint32_t circular)
{
	b = MIN2(k*2, b);				/* clip bucket size */
	mm_idx_t *mi = (mm_idx_t *)calloc(1, sizeof(mm_idx_t));
	*mi = (mm_idx_t){
		.mask = (1<<b) - 1,
		.bkt = (mm_idx_bucket_t *)calloc(1<<b, sizeof(mm_idx_bucket_t)),
		.b = b,
		.w = w<1? 1 : w,
		.k = k,
		.circular = circular,
		.base_rid = base_rid
	};
	return(mi);
}

/**
 * @fn mm_idx_destroy
 */
static _force_inline
void mm_idx_destroy(mm_idx_t *mi)
{
	if(mi == 0) { return; }

	/* buckets */
	for(uint64_t i = 0; i < 1ULL<<mi->b; i++) {
		free(mi->bkt[i].p);
		free(mi->bkt[i].w.h.a);
	}
	free(mi->bkt);

	/* sequence containers */
	for(uint64_t i = 0; i < mi->mem.n; i++) {
		free((void *)mi->mem.a[i].u64[1]);
	}
	free(mi->mem.a);

	/* sequence array */
	free(mi->s.a);
	free(mi);
	return;
}

/**
 * @fn mm_idx_get
 * @brief retrieve element from hash table (hotspot)
 */
static _force_inline
v2u32_t const *mm_idx_get(
	mm_idx_t const *mi,
	uint64_t minier,
	uint32_t *restrict n)
{
	mm_idx_bucket_t *b = &mi->bkt[minier & mi->mask];
	kh_t *h = (kh_t *)&b->w.h;
	uint64_t const *p;

	if(h->a == NULL || (p = kh_get_ptr(h, minier>>mi->b)) == NULL) {
		*n = 0;
		return(NULL);
	}
	if((int64_t)*p >= 0) {
		*n = 1;
		return((v2u32_t const *)p);
	} else {
		*n = (uint32_t)*p;
		return((v2u32_t const *)&b->p[(*p>>32) & 0x7fffffff]);
	}
}

/**
 * @fn mm_idx_cal_max_occ
 * @brief calculate occurrence threshold from frequency
 * (note: redundancy of array construction among multiple calls should be removed)
 */
static _force_inline
uint32_t *mm_idx_cal_max_occ(
	mm_idx_t const *mi,
	uint64_t n_frq,
	float const *frq)
{
	set_info(0, "[mm_idx_cal_max_occ] calculate occurrence thresholds");

	/* calc size of occurrence table */
	uint64_t n = 0;
	for(uint64_t i = 0; i < 1ULL<<mi->b; i++) {
		if(mi->bkt[i].w.h.a != NULL) {
			n += kh_cnt((kh_t *)&mi->bkt[i].w.h);
		}
	}

	/* build occurrence table */
	uint32_t *a = (uint32_t *)malloc(n * sizeof(uint32_t));
	for(uint64_t i = (n = 0); i < 1ULL<<mi->b; i++) {
		kh_t *h = (kh_t *)&mi->bkt[i].w.h;
		if(h->a == NULL) { continue; }
		for(uint64_t k = 0; k < kh_size(h); k++) {
			if(!kh_exist(h, k)) { continue; }
			a[n++] = ((int64_t)kh_val(h, k)) >= 0
				? 1
				: (uint32_t)kh_val(h, k);
		}
	}

	/* calc thresh */
	uint32_t *thres = calloc(1, sizeof(uint32_t) * n_frq);
	for(uint64_t i = 0; i < n_frq; i++) {
		thres[i] = frq[i] <= 0.0
			? UINT32_MAX
			: (ks_ksmall_uint32_t(n, a, (uint32_t)((1.0 - frq[i]) * n)) + 1);
	}
	free(a);
	return(thres);
}

/******************
 * Generate index *
 ******************/

/**
 * @struct mm_idx_pipeline_t
 * @brief multithreaded index construction pipeline context container
 */
typedef struct {
	uint32_t icnt, ocnt;
	bseq_file_t *fp;
	mm_idx_t *mi;
	kvec_t(v4u32_t) hq;
} mm_idx_pipeline_t;

/**
 * @struct mm_idx_step_t
 */
typedef struct {
    uint32_t id, base_rid, n_seq;
	bseq_t const *seq;		/* const!!! */
	void *base;
	uint64_t size;
	v4u32_v a;
} mm_idx_step_t;

/**
 * @fn mm_idx_source
 * @brief sequence reader (source) of the index construction pipeline
 */
static
void *mm_idx_source(uint32_t tid, void *arg)
{
	set_info(tid, "[mm_idx_source] fetch sequence block");
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t*)arg;
	mm_idx_step_t *s = (mm_idx_step_t*)calloc(1, sizeof(mm_idx_step_t));

	/* fetch sequence */
	uint64_t size;
	void *base;
	s->seq = bseq_read(q->fp, &s->n_seq, &base, &size);
	if(s->seq == NULL) {
		free(s); s = NULL;
		return NULL;
	}

	/* assign block id */
	s->id = q->icnt++;

	/* update base_rid and seq array (for all-versus-all) */
	if(q->mi->base_rid == UINT32_MAX) {
		q->mi->base_rid = atoi(s->seq[0].name);		/* assume first seq name is base_rid */
	}
	s->base_rid = q->mi->base_rid + q->mi->s.n;

	/* register fetched block */
	kv_push(v4u32_t, q->mi->mem, ((v4u32_t){ .u64 = { size, (uintptr_t)base }}));

	/* copy sequence pointers */
	if(q->mi->s.n + s->n_seq > q->mi->s.m) {
		q->mi->s.m = MAX2(256, 2*q->mi->s.m);
		q->mi->s.a = realloc(q->mi->s.a, sizeof(mm_idx_seq_t) * q->mi->s.m);
	}
	bseq_t const *src = s->seq;
	mm_idx_seq_t *dst = &q->mi->s.a[q->mi->s.n];
	for(uint64_t i = 0; i < s->n_seq; i++) {
		dst[i] = (mm_idx_seq_t){
			.seq = src[i].seq,   .l_seq = src[i].l_seq,
			.name = src[i].name, .l_name = src[i].l_name,
			.rid = s->base_rid + i,
			.circular = q->mi->circular
		};
	}
	q->mi->s.n += s->n_seq;
	return(s);		/* passed to the minimizer calculation stage */
}

/**
 * @fn mm_idx_worker
 * @brief calculate minimizer
 */
static
void *mm_idx_worker(uint32_t tid, void *arg, void *item)
{
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t *)arg;
	mm_idx_step_t *s = (mm_idx_step_t *)item;

	char buf[128], *p = buf;
	p += _pstr(p, "[mm_idx_worker] bin id "); p += _pnum(uint32_t, p, s->base_rid); p += _pstr(p, ":"); p += _pnum(uint32_t, p, s->base_rid + s->n_seq - 1); *p = '\0';
	set_info(tid, buf);
	for(uint64_t i = 0; i < s->n_seq; i++) {
		mm_sketch(s->seq[i].seq, s->seq[i].l_seq, q->mi->w, q->mi->k, s->base_rid + i, &s->a);
	}
	return(s);
}

/**
 * @fn mm_idx_drain_intl
 */
static _force_inline
void mm_idx_drain_intl(mm_idx_pipeline_t *q, mm_idx_step_t *s)
{
	uint64_t mask = q->mi->mask;
	for(uint64_t i = 0; i < s->a.n; i++) {
		v4u32_v *p = &q->mi->bkt[s->a.a[i].u64[0] & mask].w.a;
		kv_push(v4u32_t, *p, s->a.a[i]);
	}
	free((void *)s->seq);
	free(s->a.a);
	free(s);
	return;
}

/**
 * @fn mm_idx_drain
 */
static
void mm_idx_drain(uint32_t tid, void *arg, void *item)
{
	set_info(tid, "[mm_idx_drain] dump minimizers to pool");
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t *)arg;
	mm_idx_step_t *s = (mm_idx_step_t *)item;

	#if STRICT_STREAM_ORDERING != 0
		/* sorted pipeline (for debugging) */
		kv_hq_push(v4u32_t, incq_comp, q->hq, ((v4u32_t){.u64 = {s->id, (uintptr_t)s}}));
		while(q->hq.n > 1 && q->hq.a[1].u64[0] == q->ocnt) {
			q->ocnt++;
			s = (mm_idx_step_t *)kv_hq_pop(v4u32_t, incq_comp, q->hq).u64[1];
			mm_idx_drain_intl(q, s);
		}
	#else
		/* unsorted pipeline (this may change the results of seed chaining) */
		mm_idx_drain_intl(q, s);
	#endif
	return;
}

/**
 * @struct mm_idx_post_t
 * @brief bucket sorter context
 */
typedef struct {
	mm_idx_t *mi;
	uint32_t from, to;
} mm_idx_post_t;

/**
 * @fn mm_idx_post
 * @brief bucket sorter
 */
static
void *mm_idx_post(uint32_t tid, void *arg, void *item)
{
	set_info(tid, "[mm_idx_post] indexing postprocess");
	mm_idx_post_t *q = (mm_idx_post_t*)item;
	mm_idx_t *mi = q->mi;

	for(uint64_t i = q->from; i < q->to; i++) {
		mm_idx_bucket_t *b = &mi->bkt[i];
		v4u32_t *arr = b->w.a.a;
		uint64_t n_arr = b->w.a.n;

		if(n_arr == 0) {
			_memset_blk_u(b, 0, sizeof(mm_idx_bucket_t));
			continue;
		}

		/* sort by minimizer */
		radix_sort_128x(arr, n_arr);

		/* count keys and preallocate buffer */
		uint64_t n_keys = 0, n_elems = 0;	/* reserve table size at p[0] */
		for(uint64_t j = 1, n = 1; j <= n_arr; j++, n++) {
			if(j != n_arr && arr[j].u64[0] == arr[j - 1].u64[0]) { continue; }
			n_elems += (n > 1)? n : 0;
			n_keys++;
			n = 0;
		}
		kh_init_static(&b->w.h, n_keys / KH_THRESH);
		b->p = (uint64_t *)malloc(sizeof(uint64_t) * (n_elems + 1));
		b->p[0] = n_elems;

		/* create the (2nd-stage) hash table for the bucket */
		for(uint64_t j = 1, n = 1, sp = 1; j <= n_arr; j++, n++) {
			if(j != n_arr && arr[j].u64[0] == arr[j - 1].u64[0]) { continue; }

			v4u32_t *p = &arr[j - n];
			uint64_t key = p->u64[0]>>mi->b, val = p->u64[1];
			if(n != 1) {
				b->p[sp++] = val;
				val = (sp - 1)<<32 | n | 0x01ULL<<63;	/* k = 0 */
				while(n > 1) {
					b->p[sp++] = arr[j - --n].u64[1];
				}
			}
			kh_put(&b->w.h, key, val); n = 0;
		}

		/* deallocate table */
		free(arr);
	}
	return(NULL);
}

/**
 * @fn mm_idx_gen
 * @brief root function of the index construction pipeline
 */
static _force_inline
mm_idx_t *mm_idx_gen(mm_mapopt_t const *opt, bseq_file_t *fp)
{
	mm_idx_pipeline_t pl = {0}, **p;

	set_info(0, "[mm_idx_gen] initialize index object");

	/* init pipeline context */
	pl.icnt = pl.ocnt = 0;
	pl.fp = fp;
	if(pl.fp == NULL) { return(NULL); }		/* fp must be instanciated */

	/* initialize index object */
	pl.mi = mm_idx_init(opt->w, opt->k, opt->b, opt->base_rid, (opt->flag&MM_CIRCULAR) != 0);
	kv_hq_init(pl.hq);

	/* create thread contexts */
	p = (mm_idx_pipeline_t **)calloc(opt->nth, sizeof(mm_idx_pipeline_t*));
	for(uint64_t i = 0; i < opt->nth; i++) { p[i] = &pl; }

	/* read sequence and collect minimizers */
	pt_t *pt = pt_init(opt->nth);
	pt_stream(pt, mm_idx_source, &pl, mm_idx_worker, (void **)p, mm_idx_drain, &pl);
	free(p);
	kv_hq_destroy(pl.hq);
	opt->log(opt, 9, __func__, "collected minimizers");

	/* sort minimizers */
	mm_idx_post_t *q = (mm_idx_post_t *)calloc(opt->nth, sizeof(mm_idx_post_t));
	mm_idx_post_t **qq = (mm_idx_post_t **)calloc(opt->nth, sizeof(mm_idx_post_t*));
	for(uint64_t i = 0; i < opt->nth; i++) {
		q[i].mi = pl.mi;
		q[i].from = (1ULL<<pl.mi->b) * i / opt->nth;
		q[i].to = (1ULL<<pl.mi->b) * (i + 1) / opt->nth;
		qq[i] = &q[i];
	}
	pt_parallel(pt, mm_idx_post, NULL, (void **)qq);
	pt_destroy(pt);
	opt->log(opt, 9, __func__, "sorted minimizers");

	free(q);
	free(qq);
	return(pl.mi);
}

#if 0
/**
 * @fn mm_idx_cmp
 * @brief compare two index objects (for debugging)
 */
static _force_inline
void mm_idx_cmp(mm_mapopt_t const *opt, mm_idx_t const *m1, mm_idx_t const *m2)
{
	for(uint64_t i = 0; i < 1ULL<<opt->b; ++i) {
		mm_idx_bucket_t *bkt1 = &m1->bkt[i], *bkt2 = &m2->bkt[i];
		if(bkt1 == NULL || bkt2 == NULL) {
			if(bkt1 == NULL && bkt2 == NULL) continue;
			if(bkt1) fprintf(stderr, "i(%lu), bkt1 is instanciated but bkt2 is not\n", i);
			if(bkt2) fprintf(stderr, "i(%lu), bkt2 is instanciated but bkt1 is not\n", i);
			continue;
		}

		if(bkt1->n != bkt2->n) fprintf(stderr, "i(%lu), array size differs(%lu, %lu)\n", i, bkt1->n, bkt2->n);
		for(uint64_t j = 0; j < MIN2(bkt1->n, bkt2->n); ++j) {
			if(bkt1->p[j] != bkt2->p[j]) fprintf(stderr, "i(%lu), array differs at j(%lu), (%lu, %lu)\n", i, j, bkt1->p[j], bkt2->p[j]);
		}

		kh_t *h1 = bkt1->h, *h2 = bkt2->h;
		if(h1 == NULL || h2 == NULL) {
			if(h1 == NULL && h2 == NULL) continue;
			if(h1) fprintf(stderr, "i(%lu), h1 is instanciated but h2 is not\n", i);
			if(h2) fprintf(stderr, "i(%lu), h2 is instanciated but h1 is not\n", i);
			continue;
		}
		if(h1->mask != h2->mask) fprintf(stderr, "i(%lu), hash size differs(%u, %u)\n", i, h1->mask, h2->mask);
		if(h1->max != h2->max) fprintf(stderr, "i(%lu), hash max differs(%u, %u)\n", i, h1->max, h2->max);
		if(h1->cnt != h2->cnt) fprintf(stderr, "i(%lu), hash cnt differs(%u, %u)\n", i, h1->cnt, h2->cnt);
		if(h1->ub != h2->ub) fprintf(stderr, "i(%lu), hash ub differs(%u, %u)\n", i, h1->ub, h2->ub);
		for(uint64_t j = 0; j < MIN2(h1->mask, h2->mask)+1; ++j) {
			if(h1->a[j].u64[0] != h2->a[j].u64[0]) fprintf(stderr, "i(%lu), hash key differs at j(%lu), (%lu, %lu)\n", i, j, h1->a[j].u64[0], h2->a[j].u64[0]);
			if(h1->a[j].u64[1] != h2->a[j].u64[1]) fprintf(stderr, "i(%lu), hash val differs at j(%lu), (%lu, %lu)\n", i, j, h1->a[j].u64[1], h2->a[j].u64[1]);
		}
	}
}
#endif

/*************
 * index I/O *
 *************/

#define MM_IDX_MAGIC "MAI\x08"		/* minialign index version 8 */

/**
 * @struct mm_idx_hdr_s
 */
typedef struct mm_idx_hdr_s {
	uint8_t b, w, k, circular;
	uint32_t base_rid;
	uint64_t n_seq, size;
} mm_idx_hdr_t;

/**
 * @fn mm_idx_dump
 * @brief dump index to fp
 */
static _force_inline
void mm_idx_dump(FILE *fp, mm_idx_t const *mi, uint32_t nth)
{
	set_info(0, "[mm_idx_dump] dump index to file (main thread)");

	mm_idx_hdr_t hdr = {
		.b = mi->b, .w = mi->w, .k = mi->k,
		.circular = mi->circular,
		.base_rid = mi->base_rid,
		.n_seq = mi->s.n,
		.size = 0
	};

	/* create compresser threads */
	pg_t *pg = pg_init(fp, nth);

	/* calculate sequence block size */
	for(uint64_t i = 0; i < mi->mem.n; i++) {
		hdr.size += mi->mem.a[i].u64[0];
	}

	/* dump header */
	pgwrite(pg, MM_IDX_MAGIC, strlen(MM_IDX_MAGIC));
	pgwrite(pg, &hdr, sizeof(mm_idx_hdr_t));

	/* dump buckets */
	for(uint64_t i = 0; i < 1ULL<<mi->b; i++) {
		mm_idx_bucket_t *b = &mi->bkt[i];
		uint64_t n = b->p != NULL ? b->p[0] : 0;
		pgwrite(pg, &n, sizeof(uint64_t));					/* value table size */
		debug("i(%llu), n(%llu)", i, n);
		pgwrite(pg, &b->p[1], sizeof(uint64_t) * n);		/* value table content, size in b->p[0] */
		kh_dump(&b->w.h, pg, (khwrite_t)pgwrite);			/* 2nd-stage hash table */
	}

	/* dump sequence contents */
	for(uint64_t i = 0; i < mi->mem.n; i++) {
		pgwrite(pg, (void *)mi->mem.a[i].u64[1], sizeof(char) * mi->mem.a[i].u64[0]);
	}

	/* dump sequence pointer arrays; first convert pointers to offsets */
	for(uint64_t i = 0, j = 0, s = 0; i < mi->s.n; i++) {
		/* calculate offset from head */
		if((uintptr_t)mi->s.a[i].seq < (uintptr_t)mi->mem.a[j].u64[1]
		|| (uintptr_t)mi->s.a[i].seq >= (uintptr_t)mi->mem.a[j].u64[1] + (ptrdiff_t)mi->mem.a[j].u64[0]) {
			s += mi->mem.a[j++].u64[0];
		}

		/* convert pointer to an offset from head */
		mi->s.a[i].name -= (ptrdiff_t)mi->mem.a[j].u64[1]; mi->s.a[i].name += (ptrdiff_t)s;
		mi->s.a[i].seq -= (ptrdiff_t)mi->mem.a[j].u64[1];  mi->s.a[i].seq += (ptrdiff_t)s;
	}
	pgwrite(pg, mi->s.a, sizeof(mm_idx_seq_t) * mi->s.n);

	/* restore pointers */
	for(uint64_t i = 0, j = 0, s = 0; i < mi->s.n; i++) {
		mi->s.a[i].name -= (ptrdiff_t)s; mi->s.a[i].name += (ptrdiff_t)mi->mem.a[j].u64[1];
		mi->s.a[i].seq -= (ptrdiff_t)s;  mi->s.a[i].seq += (ptrdiff_t)mi->mem.a[j].u64[1];
		if((uintptr_t)mi->s.a[i].name > s + mi->mem.a[j].u64[0]) {
			s += mi->mem.a[j++].u64[0];
		}
	}
	pg_destroy(pg);
	return;
}

/**
 * @fn mm_idx_load
 * @brief create index object from file stream
 */
static _force_inline
mm_idx_t *mm_idx_load(FILE *fp, uint32_t nth)
{
	set_info(0, "[mm_idx_load] load index from file (main thread)");

	mm_idx_t *mi = NULL;
	pg_t *pg = pg_init(fp, nth);

	/* read magic */
	char magic[4];
	if(pgread(pg, magic, strlen(MM_IDX_MAGIC)) != strlen(MM_IDX_MAGIC)
	 || strncmp(magic, MM_IDX_MAGIC, strlen(MM_IDX_MAGIC)) != 0) {
		goto _mm_idx_load_fail;
	}

	/* read header */
	mm_idx_hdr_t hdr;
	if(pgread(pg, &hdr, sizeof(mm_idx_hdr_t)) != sizeof(mm_idx_hdr_t)) {
		goto _mm_idx_load_fail;
	}

	/* instanciate index object */
	mi = mm_idx_init(hdr.w, hdr.k, hdr.b, hdr.base_rid, hdr.circular);

	/* load hash tables */
	for(uint64_t i = 0; i < 1ULL<<mi->b; i++) {
		mm_idx_bucket_t *b = &mi->bkt[i];

		/* read value table size */
		uint64_t n;
		if(pgread(pg, &n, sizeof(uint64_t)) != sizeof(uint64_t)) {
			goto _mm_idx_load_fail;
		}

		// debug("i(%llu), n(%llu)", i, n);
		/* read value table content */
		b->p = (uint64_t *)malloc((n + 1) * sizeof(uint64_t));
		b->p[0] = n;
		if(n > 0 && pgread(pg, &b->p[1], sizeof(uint64_t) * n) != sizeof(uint64_t) * n) {
			goto _mm_idx_load_fail;
		}

		/* read hash table */
		kh_load_static(&b->w.h, pg, (khread_t)pgread);
	}

	/* load sequences */
	v4u32_t const mem = {
		.u64 = {
			[0] = hdr.size,
			[1] = (uintptr_t)malloc(sizeof(char) * hdr.size)
		}
	};
	kv_push(v4u32_t, mi->mem, mem);
	if(pgread(pg, (void *)mem.u64[1], sizeof(char) * hdr.size) != sizeof(char) * hdr.size) {
		goto _mm_idx_load_fail;
	}

	/* load sequence pointer array */
	mi->s.n = hdr.n_seq;
	mi->s.a = malloc(sizeof(mm_idx_seq_t) * mi->s.n);
	if(pgread(pg, mi->s.a, sizeof(mm_idx_seq_t) * mi->s.n) != sizeof(mm_idx_seq_t) * mi->s.n) {
		goto _mm_idx_load_fail;
	}

	/* convert offsets to pointers */
	for(uint64_t i = 0; i < mi->s.n; i++) {
		mi->s.a[i].name += (ptrdiff_t)mem.u64[1];
		mi->s.a[i].seq += (ptrdiff_t)mem.u64[1];
	}
	pg_destroy(pg);
	return(mi);
_mm_idx_load_fail:
	pg_destroy(pg);
	mm_idx_destroy(mi);
	return(NULL);
}

/* end of index.c */

/* map.c */

/* forward declaration of mm_align_t */
typedef struct mm_align_s mm_align_t;
typedef struct mm_reg_s mm_reg_t;

/**
 * @struct mm_printer_t
 * @brief alignment result (record) printer interface
 */
typedef struct {
	void (*header)(mm_align_t *);
	void (*unmapped)(mm_align_t *, bseq_t const *);
	void (*mapped)(mm_align_t *, bseq_t const *, mm_reg_t const *);	/* nreg contains #regs in lower 32bit, #uniq in higher 32bit */
	void (*footer)(mm_align_t *);
} mm_printer_t;

/**
 * @struct mm_align_s
 * @brief alignment pipeline context
 */
struct mm_align_s {
	/* output buffer */
	uint8_t *base, *tail, *p;
	uint8_t conv[40];				/* binary -> string conv table */
	mm_printer_t printer;

	/* input */
	bseq_file_t *fp;

	/* worker */
	mm_idx_t const *mi;				/* immutable */
	mm_mapopt_t const *opt;
	uint32_t *occ;
	uint32_t n_occ;
	gaba_t *gaba;

	/* streaming */
	uint32_t icnt, ocnt, base_qid;
	kvec_t(v4u32_t) hq;
	void **t;						/* mm_tbuf_t* array */
	pt_t *pt;
};

/**
 * @struct mm_mini_t
 * @brief minimizer container, alias to v4u32_t and mm_resc_t
 */
typedef struct {
	uint64_t k;
	uint32_t pos, id;
} mm_mini_t;
_static_assert(sizeof(mm_mini_t) == sizeof(v4u32_t));

/**
 * @struct mm_resc_t
 * @brief rescued seed (pos, occ) container, alias to v4u32_t
 */
typedef struct {
	v2u32_t const *p;
	uint32_t qs, n;
} mm_resc_t;
typedef struct { size_t n, m; mm_resc_t *a; } mm_resc_v;
_static_assert(sizeof(mm_resc_t) == sizeof(v4u32_t));

/**
 * @struct mm_seed_t
 * @brief seed position container, alias to v4u32_t
 */
typedef struct {
	uint32_t upos, rid, apos, bpos;
} mm_seed_t;
typedef struct { size_t n, m; mm_seed_t *a; } mm_seed_v;
_static_assert(sizeof(mm_seed_t) == sizeof(v4u32_t));

/**
 * @struct mm_chain_t
 * @brief chain container, alias to v2u32_t
 */
typedef struct {
	uint32_t plen, sid;
} mm_chain_t;
typedef struct { size_t n, m; mm_chain_t *a; } mm_chain_v;
_static_assert(sizeof(mm_chain_t) == sizeof(v2u32_t));

/**
 * @struct mm_res_t
 * @brief result container, alias to v2u32_t and mm_chain_t
 */
typedef struct {
	uint32_t score, bid;
} mm_res_t;
typedef struct { size_t n, m; mm_res_t *a; } mm_res_v;
_static_assert(sizeof(mm_res_t) == sizeof(v2u32_t));

/**
 * @struct mm_map_t
 * @brief mid (kept at every seed pos) -> (chain, ref id) map, alias to v2u32_t
 */
typedef struct {
	uint32_t cid, rid;
} mm_map_t;
typedef struct { size_t n, m; mm_map_t *a; } mm_map_v;
_static_assert(sizeof(mm_map_t) == sizeof(v2u32_t));

/**
 * @struct mm_bin_t
 * @brief result bin, stacked in the bin array
 */
typedef struct {
	uint32_t n_aln, plen;			/* plen holds mapq after postprocess */
	uint32_t lb, ub;				/* (bpos, bpos + blen) */
	gaba_alignment_t const *aln[];
} mm_bin_t;
_static_assert(sizeof(mm_bin_t) % sizeof(void *) == 0);
#define MM_BIN_N					( sizeof(mm_bin_t) / sizeof(void *) )

/**
 * @struct mm_alnset_t
 */
typedef struct {
	// uint32_t score;
	uint32_t rid, mapq;				/* reference id and mapq */
	uint32_t n_aln, b_aln;			/* #alns, base index of reg->bin */
} mm_alnset_t;

typedef struct {
	// uint32_t score;
	uint32_t aid, mapq;				/* alignment id and mapq */
	gaba_alignment_t const a[];
} mm_aln_t;

/**
 * @struct mm_reg_t
 */
struct mm_reg_s {
	uint32_t n_all, n_uniq;
	mm_aln_t const *aln[];
};

/**
 * @struct mm_align_step_t
 * @brief batch object
 */
typedef struct {
	/* query sequences */
	uint32_t id, base_qid, n_seq;
	bseq_t const *seq;				/* const!!! */
	void *base;						/* memory block pointer */
	uint64_t size;					/* memory block size (unused...) */

	/* result buffer */
	lmm_t *lmm;						/* alignment result container */
} mm_align_step_t;

/**
 * @struct mm_tbuf_s
 * @brief per-thread buffer (working context)
 */
typedef struct mm_tbuf_s {
	/* initialized in mm_tbuf_init */
	mm_idx_t mi;					/* (immutable) index object */

	/* frequently-used constants (copied from mm_align_t and mm_mapopt_t; initialized in mm_tbuf_init) */
	/* calculated */
	uint32_t org, thresh;			/* ava filter constants */

	/* copied */
	uint64_t flag;
	uint32_t _pad1, n_occ, occ[MAX_FRQ_CNT + 1];
	int32_t m, x, wlen, glen;		/* chainable window edge length, linkable gap length */
	float min_ratio;
	uint32_t min_score, max_cnt;

	/* initialized in mm_init_query */
	uint32_t adj;					/* qlen + mi->k */
	uint32_t rid, qid;				/* raw qid (used in ava filter) */

	/* query and reference sections (initialized in mm_init_query) */
	uint8_t tail[96];				/* zeros */
	gaba_section_t r[2];			/* [0]: rf, [1]: rr */
	gaba_section_t q[3];			/* [0]: qf, [1]: qr, [2]: qf */
	gaba_section_t t[1];

	/* working buffers */
	mm_resc_v resc;					/* rescued seed array, also used as query minimizer vector */
	mm_resc_t *presc;				/* current head of rescued array */
	mm_seed_v seed;					/* seed array */
	mm_map_v map;					/* (seed -> chain id, seed->rid) map */
	mm_chain_v chain;				/* intervals (chains) on sorted coef, (sid, _ofs(plen)) */
	uint32_t n_res, _pad2;			/* #alignments collected */
	ptr_v bin;						/* gaba_alignment_t* array */
	kh_t pos;						/* alignment dedup hash */
	gaba_dp_t *dp;					/* alignment work */
	gaba_trace_params_t trace;		/* lmm contained */
} mm_tbuf_t;

/**
 * @macro _s, _m
 * @brief extract sign, calc inv
 */
#define _sign(x)		( (x)<0?-1:1 )
#define _smask(x)		( ((int32_t)(x))>>31 )
#define _umask(x)		( 0x7fffffff & (x) )
#define _mabs(x)		( _smask(x)^(x) )
#define _ofs(x)			( (int32_t)0x40000000 - (int32_t)(x) )
// #define _u(x, y)		( (((x)<<1) | 0x40000000) - (y) )
// #define _v(x, y)		( ((x) | 0x40000000) - ((y)<<1) )

#define _u(x, y)		( ((x)<<1) - (y) + ((y) & 0x80000000) + 0x40000000 )
#define _v(x, y)		( (x) - ((y)<<1) + ((y) & 0x80000000) + 0x40000000 )

#define _ud(x, y)		( ((x)<<1) - (y) )
#define _vd(x, y)		( (x) - ((y)<<1) )
#define _p(_ptr)		( ((int32_t const *)(_ptr))[0] + ((int32_t const *)(_ptr))[1] )
#define _q(_ptr)		( ((int32_t const *)(_ptr))[1] - ((int32_t const *)(_ptr))[0] )
#define _inside(_lb, _x, _ub)		( (uint32_t)((_x) - (_lb)) <= (uint32_t)((_ub) - (_lb)) )
#define _key(x, y)		( (x) ^ ((x)>>29) ^ (y) ^ _swap_u64(y) )

/**
 * @fn mm_expand
 * @brief expand minimizer to coef array
 */
static _force_inline
void mm_expand(
	mm_tbuf_t *self,
	uint32_t const n,
	v2u32_t const *r,							/* source array */
	int32_t const qs)							/* query position */
{
	if(n == 0) { return; }

	/* expand array if needed */
	kv_reserve(mm_seed_t, self->seed, self->seed.n + n);

	/* prep mask */
	// uint32_t const qmask = _smask(qs);
	// uint32_t const qs = _mabs(qpos);

	/* iterate over all the collected minimizers */
	for(uint64_t i = 0; i < n; i++) {
		uint32_t rid = r[i].u32[1];
		if(self->org + rid - self->qid < self->thresh) {
			continue;							/* skip if seed is in the lower triangle (all-versus-all) */
		}

		int32_t rs = r[i].u32[0];				/* load reference pos */
		uint32_t rmask = _smask(rs), _rs = rs ^ rmask, _qs = qs ^ rmask;
		self->seed.a[self->seed.n++] = (mm_seed_t){
			.upos = _u(_rs, _qs),				/* coef */
			.rid = rid,							/* rid */
			.apos = _rs, .bpos = _qs			/* rs always positive */
		};
		// debug("pos(%d, %d)", _rs, _qs);
		// qs = rev ? (int32_t)(qf->len - (qs ^ rev) + k - 1) : qs;
		// qs = rev ? (int32_t)(qf->len + qs + k) : qs;
	}
	return;
}

/**
 * @fn mm_collect_seed
 * @brief collect minimizers for the query seq
 * (note: branch misprediction penalty and memory access pattern should be alleviated, but how?)
 */
static _force_inline
void mm_collect_seed(
	mm_tbuf_t *self)
{
	/* clear array and gather minimizers */
	self->resc.n = self->seed.n = 0;
	mm_sketch(self->q[0].base, self->q[0].len, self->mi.w, self->mi.k, 0, (v4u32_v *)&self->resc);
	debug("collected seeds, n(%zu)", self->resc.n);

	uint32_t const max_occ = self->occ[self->n_occ - 1];
	uint32_t const resc_occ = self->occ[0];
	uint64_t n_rescued = 0;

	/* iterate over all the collected minimizers (seeds) on the query */
	for(mm_mini_t *p = (mm_mini_t *)self->resc.a, *t = p + self->resc.n; p < t; p++) {

		/* get minimizer matched on the ref at the current query pos */
		uint32_t n;
		v2u32_t const *r = mm_idx_get(&self->mi, p->k, &n);
		
		/* skip if exceeds repetitive threshold */
		if(n > max_occ) { continue; }

		/* save if less than max but exceeds current threshold */
		int32_t qs = (int32_t)p->pos;
		if(n > resc_occ) {
			self->resc.a[n_rescued++] = (mm_resc_t){
				.p = r, .qs = qs, .n = n
			};
			continue;
		}

		/* append to seed array */
		mm_expand(self, n, r, qs);
	}
	self->resc.n = n_rescued;					/* write back rescued array */
	return;
}

/**
 * @fn mm_seed
 * @brief construct seed array
 */
static _force_inline
uint64_t mm_seed(
	mm_tbuf_t *self,
	uint64_t cnt)								/* iteration count */
{
	debug("seed iteration i(%llu)", cnt);
	if(cnt == 0) {
		/* first collect seeds */
		mm_collect_seed(self);
		self->presc = self->resc.a;				/* init resc pointer */
	} else {
		/* sort rescued array by occurrence */
		debug("sort resc, n(%zu)", self->resc.n);
		if(cnt == 1) { radix_sort_128x((v4u32_t *)self->resc.a, self->resc.n); }

		/* append seeds */
		mm_resc_t *p = self->presc, *t = &self->resc.a[self->resc.n];
		while(p < t && p->n <= self->occ[cnt]) {
			mm_expand(self, p->n, p->p, p->qs);
			p++;
		}
		self->presc = p;						/* write back resc pointer */
	}
	if(self->seed.n == 0) { return(0); }

	/* sort seed array */
	debug("sort seed, n(%zu)", self->seed.n);
	radix_sort_128x((v4u32_t *)self->seed.a, self->seed.n);
	return(self->seed.n);
}

/**
 * @fn mm_chain_core
 * @brief starting from q, link the nearest seed iteratively, open a bin for the chain
 */
static _force_inline
mm_seed_t *mm_chain_core(
	mm_tbuf_t *self,
	mm_seed_t *p)
{
	/* q is the current front, np holds the first unlinked seed */
	mm_seed_t *q = p, *t = self->seed.a + self->seed.n, *np = p + 1;	/* FIXME: update np in the loop for better performance */

	/* save rid, allocate chain bin and map bin */
	uint64_t key = _loadu_u64(&q->upos);
	uint32_t mid = ~self->map.n++;	/* must be NOT (not NEG) to avoid collision at zero */

	/* loop: link seeds which belong to the same (rid, qid) pair */
	do {
		/* load window boundaries */
		uint32_t ubase = _v(q->apos, q->bpos);
		uint32_t uub = q->upos + _ud(self->wlen, 0);
		uint32_t vub = _v(q->apos, q->bpos);
		uint32_t vlb = vub + _vd(0, self->wlen);
		int32_t pd = _p(&q->apos) + 2*self->wlen;

		debug("uub(%d), vlb(%d), vub(%d)", uub, vlb, vub);

		/* loop: pick up the nearest seed */
		mm_seed_t *nq = t;
		for(mm_seed_t *r = q + 1; r < t && r->upos < uub; r++) {
			debug("test seed at i(%llu), rid(%u), (%d, %d), u(%d), v(%d), uub(%d), ubase(%d, %d), np(%ld)",
				r - self->seed.a, r->rid, (int32_t)r->apos, r->bpos, r->upos, _v(r->apos, r->bpos), uub, ubase, ubase + _p(&r->apos), np - self->seed.a);

			// if((int32_t)(r->rid - rid) >= 0) { np = MIN2(np, r);}

			/* terminate if rid changed */
			if(((int64_t)(_loadu_u64(&r->upos) - key)>>30) > 0) {
				debug("rid changed, (%llx, %llx)", _loadu_u64(&r->upos), key);
				break;
			}

			/* test if r is inside the window */
			if(!_inside(vlb, _v(r->apos, r->bpos), vub)) {
				/*np = MIN2(np, r);*/ continue;
			}

			/* update min */
			debug("test pd(%d, %d)", pd, _p(&r->apos));
			if(_p(&r->apos) < pd) {
				nq = r;
				pd = _p(&r->apos);
			}

			/* FIXME: update uub (this seems correct but needs further verified) */
			uub = MIN2(uub, ubase + _p(&r->apos));
			debug("updated uub(%d)", uub);
		}

		/* mark q as chained (is the root for the first itr), update q for the next itr */
		q->rid = mid;							/* store mid, distinguished by its sign */
		if(_unlikely(nq == t)) {
			/* linkable seed not found, open a new bin for the chain */
			uint32_t cid = self->chain.n++;
			self->chain.a[cid] = (mm_chain_t){
				.plen = _ofs(_p(&q->apos) - _p(&p->apos)),
				.sid = ((int32_t)q->apos < 0 ? p : q) - self->seed.a
			};
			debug("linkable seed not found, cid(%u), plen(%d), sid(%u)", cid, _ofs(self->chain.a[cid].plen), self->chain.a[cid].sid);

			/* save rid and cid to map */
			self->map.a[~mid] = (mm_map_t){
				.rid = key>>32, .cid = cid
			};
			return(np);
		}
		debug("chain seed at i(%llu), rid(%u), (%d, %d), u(%d), v(%d), uub(%d)", nq - self->seed.a, nq->rid, nq->apos, nq->bpos, nq->upos, _v(nq->apos, nq->bpos), uub);
		q = nq;
	} while(((_loadu_u64(&q->upos) - key)>>30) == 0);

	/* hit at the middle of the existing chain, shrink chain array by one */
	uint32_t cid = self->map.a[~q->rid].cid;	/* load chain id */
	self->map.a[~mid] = (mm_map_t){
		.rid = key>>32, .cid = cid
	};
	debug("chain recorded, cid(%u), rid(%u)", cid, key>>32);

	/* update plen */
	q = &self->seed.a[self->chain.a[cid].sid];	/* load tail of the joined chain */
	uint32_t plen = _ofs(_p(&q->apos) - _p(&p->apos));
	self->chain.a[cid].plen = MIN2(self->chain.a[cid].plen, plen);
	return(np);
}

/**
 * @fn mm_chain
 * @brief collect chain in self->chain, self->seed must be sorted, self->chain and self->map are not cleared at the head.
 */
static _force_inline
uint64_t mm_chain(
	mm_tbuf_t *self)
{
	/* reserve space for map and chain */
	kv_reserve(mm_map_t, self->map, self->map.n + self->seed.n);
	kv_reserve(mm_chain_t, self->chain, self->chain.n + self->seed.n);

	/* iterate over sorted seed array */
	mm_seed_t *p = self->seed.a, *t = self->seed.a + self->seed.n;
	while(p < t) {
		/* pick up root at p, skip chained */
		if((int32_t)p->rid < 0) { p++; continue; }
		debug("root seed at i(%llu), rid(%u), (%d, %d), u(%d), v(%d)", p - self->seed.a, p->rid, p->apos, p->bpos, p->upos, _v(p->apos, p->bpos));

		p = mm_chain_core(self, p);
	}
	if(self->chain.n == 0) { return(0); }

	/* sort chain by length */
	radix_sort_64x((v2u32_t *)self->chain.a, self->chain.n);
	return(self->chain.n);
}

/**
 * @fn mm_next
 * @brief pick up next seed in the upward region of the current chain
 */
static _force_inline
gaba_pos_pair_t mm_next(
	mm_tbuf_t *self,
	gaba_pos_pair_t cp,
	uint32_t ppos,
	mm_seed_t **pp)
{
	/* load seed array pointers */
	mm_seed_t *p = *pp, *h = self->seed.a;
	uint32_t const cid = self->map.a[~p->rid].cid;
	uint32_t const rid = self->map.a[~p->rid].rid;
	uint32_t const ub = _u(cp.apos, cp.bpos);

	/* skip covered, FIXME: use binary search for better performance */
	while(--p >= h && p->upos > ub && self->map.a[~p->rid].rid == rid) {}

	debug("search start i(%ld)", p - h);

	/* test uncovered */
	while(p >= h && _p(&p->apos) >= ppos && self->map.a[~p->rid].rid == rid) {
		debug("test i(%ld), cid(%u, %u), (%d, %d), u(%d)", p - h, self->map.a[~p->rid].cid, cid, p->apos, p->bpos, p->upos);
		if(self->map.a[~p->rid].cid != cid) { p--; continue; }

		/* hit */
		debug("found next seed at i(%ld), (%d, %d), (%d, %d), u(%d)", p - h, p->apos, p->bpos, p->apos, p->bpos + (_smask(p->bpos) & self->adj), p->upos);
		*pp = p;				/* write back pointer */

		int32_t rev = _smask(p->bpos);
		return((gaba_pos_pair_t){
			.apos = p->apos,
			.bpos = p->bpos + (rev & self->adj)
		});
	}

	/* not found */
	debug("not found");
	return((gaba_pos_pair_t){
		.apos = INT32_MIN,
		.bpos = INT32_MIN
	});
}

/**
 * @fn mm_test_pos
 * @brief test if alignment is already found at the pos,
 * returns score > 0 if found, 0 if evaluated but not meaningful, < 0 if unevaluated.
 */
static _force_inline
int64_t mm_test_pos(
	mm_tbuf_t *self,
	mm_res_t *res,
	gaba_pos_pair_t pos)
{
	v2u32_t *h = (v2u32_t *)kh_get_ptr(&self->pos, _key(_loadu_u64(&pos), _loadu_u64(&self->rid)));
	if(h == NULL) { return(-1); }
	if(h->u64[0] == KH_INIT_VAL) { return(0); }

	/* res_id in lower, aln_id (global) in higher */
	return(((gaba_alignment_t const **)self->bin.a)[h->u32[1]]->score);
}

/**
 * @fn mm_mark_pos
 */
static _force_inline
void mm_mark_pos(
	mm_tbuf_t *self,
	gaba_pos_pair_t pos)
{
	kh_put(&self->pos, _key(_loadu_u64(&pos), _loadu_u64(&self->rid)), KH_INIT_VAL);
	return;
}

/**
 * @fn mm_record
 * @brief record alignment, returns nonzero if new head position found, zero the alignment is duplicated
 */
static _force_inline
uint64_t mm_record(
	mm_tbuf_t *self,
	mm_res_t *res,
	gaba_alignment_t const *a)
{
	/* put head and tail positions to hash */
	uint64_t id = _loadu_u64(&self->rid);
	uint64_t hkey = _key(_loadu_u64(&a->sec->apos), id);
	uint64_t tkey = _key(_loadu_u64(&a->sec->apos) + _loadu_u64(&a->sec->alen), id);
	debug("hash key, h(%llx), t(%llx)", hkey, tkey);

	v2u32_t *h = (v2u32_t *)kh_put_ptr(&self->pos, hkey, 1);
	v2u32_t *t = (v2u32_t *)kh_put_ptr(&self->pos, tkey, 0);		/* noextend */
	uint64_t new = h->u64[0] == KH_INIT_VAL;

	/* open new bin for aln, reuse if the head hit an existing one */
	uint32_t bid = new ? kv_push(void *, self->bin, (void *)a) : h->u32[1];
	gaba_alignment_t const **b = (gaba_alignment_t const **)&self->bin.a[bid];

	debug("id(%u, %u)", a->sec->aid, a->sec->bid);
	debug("record h(%llx, %u), t(%llx, %u), bid(%u)", h[-1].u64[0], h->u32[1], t[-1].u64[0], t->u32[1], bid);

	/* update res */
	res->score -= a->score;

	/* update bin */
	mm_bin_t *bin = (mm_bin_t *)&self->bin.a[res->bid];
	bin->n_aln += h->u64[0] == KH_INIT_VAL;
	bin->plen += _p(&a->sec->alen);
	bin->lb = MIN2(bin->lb, a->sec->bpos);
	bin->ub = MAX2(bin->ub, a->sec->bpos + a->sec->blen);
	debug("update bounds, bid(%u), lb(%u), ub(%u)", res->bid, bin->lb, bin->ub);

	debug("record alignment, (%u, %u) -> (%u, %u), n_aln(%u)", a->sec->apos, a->sec->bpos, a->sec->apos + a->sec->alen, a->sec->bpos + a->sec->blen, bin->n_aln);

	/* update hash */
	if(b[0]->score > a->score) {
		debug("discard, a(%lld), b(%lld)", a->score, b[0]->score);
		*t = (v2u32_t){ .u64[0] = KH_INIT_VAL };	/* re-mark: evaluated but not found */
	} else {
		debug("replace, a(%lld), b(%lld)", a->score, b[0]->score);
		/* replace if the new one is larger than the old one */
		*b = a;
		*h = *t = (v2u32_t){
			.u32 = {
				[0] = res - (mm_res_t *)self->chain.a,
				[1] = bid
			}
		};
	}
	return(new);
}

/**
 * @macro _sec_fw, _sec_rv
 * @brief build forward and reverse section object
 */
#define _sec_fw(_id, _base, _len) ( \
	gaba_build_section((_id)<<1, _base, _len) \
)
#define _sec_rv(_id, _base, _len) ( \
	gaba_build_section( \
		((_id)<<1) + 1, \
		gaba_rev(&(_base)[(_len) - 1], (uint8_t const *)0x800000000000), \
		(_len) \
	) \
)

/**
 * @fn mm_extend_core
 * @brief extension loop, returns fill object with max, never returns NULL.
 */
static _force_inline
gaba_fill_t const *mm_extend_core(
	gaba_dp_t *restrict dp,
	gaba_section_t const *a,
	gaba_section_t const *at,
	gaba_section_t const *b,
	gaba_section_t const *bt,
	gaba_pos_pair_t s)
{
	/* fill root */
	gaba_fill_t const *f = gaba_dp_fill_root(dp, a, s.apos, b, s.bpos);
	if(_unlikely(f == NULL)) { goto _mm_extend_core_abort; }
	debug("fill root, max(%lld), status(%x)", f->max, f->status);

	gaba_fill_t const *m = f;					/* record max */
	uint32_t flag = GABA_STATUS_TERM;
	while((flag & f->status) == 0) {
		/* update section if reached tail */
		if(f->status & GABA_STATUS_UPDATE_A) { a = at; }	/* cmov */
		if(f->status & GABA_STATUS_UPDATE_B) { b = bt; }	/* cmov */

		/* update flag for tail detection */
		flag |= f->status & (GABA_STATUS_UPDATE_A | GABA_STATUS_UPDATE_B);

		/* fill the next section */
		if(_unlikely((f = gaba_dp_fill(dp, f, a, b)) == NULL)) {
			goto _mm_extend_core_abort;
		}
		debug("fill, max(%lld, %lld), status(%x)", f->max, m->max, f->status);
		m = (f->max > m->max)? f : m;
	}
	return(m);									/* never be null */

_mm_extend_core_abort:							/* out-of-memory in libgaba */
	oom_abort(__func__, 0);						/* dump debug print */
	return(NULL);								/* noreturn */
}

/**
 * @fn mm_extend
 */
static _force_inline
uint64_t mm_extend(
	mm_tbuf_t *self)
{
	uint8_t const *lim = (uint8_t const *)0x800000000000;

	/* loop: evaluate chain */
	for(uint64_t k = 0; k < self->chain.n; k++) {
		/* load seed ptr and position */
		mm_seed_t *p = &self->seed.a[self->chain.a[k].sid];		/* load tail */
		int32_t plen = _ofs(self->chain.a[k].plen);
		if(plen * self->m < self->min_score) { continue; }

		int32_t rev = _smask(p->bpos);
		gaba_pos_pair_t cp = {
			.apos = p->apos,
			.bpos = p->bpos + (rev & self->adj)
		};
		// qs = rev ? (int32_t)(qf->len - (qs ^ rev) + k - 1) : qs;
		// qs = rev ? (int32_t)(qf->len + qs + k) : qs;

		/* load ref */
		uint32_t rid = self->map.a[~p->rid].rid;
		mm_idx_seq_t const *ref = &self->mi.s.a[rid - self->mi.base_rid];
		self->rid = rid;
		self->r[0] = _sec_fw(ref->rid, ref->seq, ref->l_seq);
		self->r[1] = _sec_rv(ref->rid, ref->seq, ref->l_seq);

		debug("chain(%llu), cid(%u), rev(%d), plen(%d), (%u, %u), (%d, %d), (%d, %d)", k, self->map.a[~p->rid].cid, rev, plen, self->rid, self->qid, p->apos, p->bpos, cp.apos, cp.bpos);

		/* open result bin */
		uint32_t bid = kv_pushm(void *, self->bin, (void **)&((mm_bin_t){ .lb = UINT32_MAX }), MM_BIN_N);

		/* reuse seed bin for score accumulation */
		mm_res_t *res = (mm_res_t *)&self->chain.a[self->n_res++];
		*res = (mm_res_t){
			.score = _ofs(0),					/* score (negated and offsetted) */
			.bid = bid							/* bin base index */
		};

		/* loop: issue extension until whole chain is covered by alignments */
		for(uint32_t ppos = _p(&cp) - plen, rem = 3, narrow = 0;
			rem > 0 && _p(&cp) >= ppos;
			cp = --rem > 0 ? mm_next(self, cp, ppos, &p) : cp
		) {
			/* reset stack */
			gaba_dp_flush(self->dp, lim, lim);

			/* init tail as zero padded arrays */
			gaba_fill_t const *t = (gaba_fill_t const *)self->tail, *u = t;

			/* downward extension */
			gaba_pos_pair_t tp = cp;
			u = mm_extend_core(&self->dp[narrow], &self->r[0], self->t, &self->q[-rev], self->t, tp);

			/* search pos if extended */
			if(u->max > 0) { tp = gaba_dp_search_max(&self->dp[narrow], u); }
			debug("len(%u, %u), score(%lld), (%u, %u) -> (%u, %u)",
				self->r[0].len, self->q[0].len, u->max, cp.apos, cp.bpos, tp.apos, tp.bpos);

			/* skip if tail is duplicated */
			if(mm_test_pos(self, res, tp) >= 0) {
				narrow = 0;
				debug("duplication detected, try narrower(%u)", narrow);
				continue;			/* try narrower band in the next itr to avoid collision */
			}

			/* upward extension: coordinate reversed here */
			t = mm_extend_core(&self->dp[0], &self->r[1], self->t, &self->q[1+rev], self->t,
				((gaba_pos_pair_t){
					.apos = self->r[0].len - tp.apos - 1,
					.bpos = self->q[0].len - tp.bpos - 1
				})
			);
			if(t->max < self->min_score) {		/* max == 0 indicates alignment was not found */
				mm_mark_pos(self, tp);			/* mark evaluated but not found */
				continue;
			}

			/* generate alignment: coordinates are reversed again, gaps are left-aligned in the resulting path */
			gaba_alignment_t const *a = gaba_dp_trace(&self->dp[0], NULL, t, &self->trace);		/* lmm is contained in self->trace */
			if(a == NULL) {
				debug("failed to generate trace: len(%u, %u), score(%lld), <- (%u, %u)", self->r[0].len, self->q[0].len, t->max, tp.apos, tp.bpos);
				continue;						/* something is wrong... */
			}
			debug("len(%u, %u), score(%lld), (%u, %u) <- (%u, %u)",
				self->r[0].len, self->q[0].len, t->max, a->sec->apos, a->sec->bpos, tp.apos, tp.bpos);

			/* record alignment, update current head position */
			if(mm_record(self, res, a)) { rem = 10; }
			if(_p(&a->sec->apos) < ppos) {
				debug("full length captured, cp(%u, %u), p(%u), ppos(%u)", cp.apos, cp.bpos, _p(&a->sec->apos), ppos);
				break;							/* full length covered, break before calling mm_next */
			}

			/* update itr states */
			cp = *((gaba_pos_pair_t *)&a->sec->apos);
			narrow = 0;
			debug("split detected, try narrower(%u), cp(%u, %u), p(%u), ppos(%u), rem(%u)", narrow, cp.apos, cp.bpos, _p(&cp), ppos, rem);
		}

		/* discard if the score did not exceed the minimum threshold */
		if(res->score > _ofs(self->min_score)) {/* _ofs() inverts the sign */
			self->bin.n = bid;
			self->n_res--;
			debug("remove bin, bid(%llu), n_res(%u)", self->bin.n, self->n_res);
		}
	}
	return(self->n_res);
}

#define MAPQ_DEC	( 4 )
#define MAPQ_COEF	( 1<<MAPQ_DEC )
#define _clip(x)	MAX2(0, MIN2(((uint32_t)(x)), 60 * MAPQ_COEF))
#define _aln(x)		( (gaba_alignment_t const *)(x).u64[1] )

/**
 * @fn mm_prune_regs
 * @brief prune low-scoring alignments, returns #alignments
 */
static _force_inline
uint64_t mm_prune_regs(
	mm_tbuf_t *self)
{
	mm_res_t *res = (mm_res_t *)self->chain.a;	/* alignments, must be sorted */

	uint64_t q = MIN2(self->max_cnt, self->n_res);
	uint32_t min = _ofs((uint32_t)(_ofs(res[0].score) * self->min_ratio));

	while(res[--q].score > min) {
		/* actually nothing to do because all the alignment objects are allocated from lmm */
		mm_bin_t *bin = (mm_bin_t *)&self->bin.a[res[q].bid];
		for(uint64_t i = 0; i < bin->n_aln; i++) {
			gaba_dp_res_free((gaba_alignment_t *)(bin->aln[i] + 1));	/* destroy gaba_alignment_t objects */
		}
	}
	return(q + 1);
}

#if COLLECT_SUPPLEMENTARY != 0
/**
 * @fn mm_collect_supp
 * @brief collect supplementary alignments
 */
static _force_inline
uint64_t mm_collect_supp(
	uint32_t const n_res,
	mm_res_t *restrict res,						/* alignments, must be sorted */
	void **bin)
{
	/* collect supplementaries */
	#define _swap_res(x, y)	{ mm_res_t _tmp = res[x]; res[x] = res[y]; res[y] = _tmp; }
	uint64_t p, q;
	for(p = 1, q = n_res; p < q; p++) {
		uint64_t max = 0;
		for(uint64_t i = p; i < q; i++) {
			/* bin[0] holds a tuple (bpos, blen) */
			mm_bin_t *s = (mm_bin_t *)&bin[res[i].bid];
			int32_t lb = s->lb, ub = s->ub, span = ub - lb;

			debug("bid(%u), s(%u, %u)", res[i].bid, lb, ub);

			for(uint64_t j = 0; j < p; j++) {
				mm_bin_t *t = (mm_bin_t *)&bin[res[j].bid];
				
				/* update boundary */
				if(t->ub < ub) {
					lb = MAX2(lb, t->ub);
				} else {
					ub = MIN2(ub, t->lb);
				}
				debug("bid(%u), t(%u, %u), s(%u, %u)", res[j].bid, t->lb, t->ub, lb, ub);

				/* calculate covered length */
				if(2*(ub - lb) < span) {		/* check if covered by j */
					debug("mark secondary, i(%llu)", i);
					q--; _swap_res(i, q); i--;	/* move to tail */
					goto _loop_tail;
				}
			}

			/* update score */
			max = MAX2(max, ((uint64_t)(2*(ub - lb) - span)<<32) | i);
		_loop_tail:;
		}
		if(max & 0xffffffff) {
			debug("mark supplementary, i(%llu)", max & 0xffffffff);
			_swap_res(p, max & 0xffffffff);		/* move to head, mark supplementary */
		}
	}
	p = MIN2(p, q);
	#undef _swap_res
	return(p);
}

/**
 * @fn mm_post_map
 * @brief mark secondary (repetitive) / supplementary flags
 */
static _force_inline
uint64_t mm_post_map(
	mm_tbuf_t *self)
{
	mm_res_t *res = (mm_res_t *)self->chain.a;	/* alignments, must be sorted */

	/* collect supp */
	uint64_t p = mm_collect_supp(self->n_res, res, self->bin.a);

	/* extract shortest repeat element */
	int32_t usc = 0, lsc = INT32_MAX, tsc = 0;
	for(uint64_t i = p; i < self->n_res; i++) {
		usc = MAX2(usc, res[i].score);
		lsc = MIN2(lsc, res[i].score);
		tsc += res[i].score;
	}
	lsc = (lsc == INT32_MAX)? 0 : lsc;

	/* calc mapq for primary and supplementary alignments */
	double tpc = 1.0;
	for(uint64_t i = 0; i < p; i++) {
		uint32_t score = res[0].score;
		mm_bin_t *bin = (mm_bin_t *)&self->bin.a[res[i].bid];

		/* calc effective length */
		double elen = (double)bin->plen / 2.0;

		/* estimate identity */
		double pid = 1.0 - (double)(elen * self->m - score) / (double)(self->m + self->x) / elen;

		/* estimate unique length */
		double ec = 2.0 / (pid * (double)(self->m + self->x) - (double)self->x);
		double ulen = ec * (score - usc), pe = 1.0 / (ulen * ulen + (double)(self->n_res - p + 1));

		/* estimate mapq */
		bin->plen = _clip(-10.0 * MAPQ_COEF * log10(pe));
		tpc *= 1.0 - pe;
	}

	/* calc mapq for secondary (repetitive) alignments */
	double tpe = MIN2(1.0 - tpc, 1.0);
	for(uint64_t i = p; i < self->n_res; i++) {
		mm_bin_t *bin = (mm_bin_t *)&self->bin.a[res[i].bid];
		bin->plen = _clip(
			  -10.0
			* MAPQ_COEF
			* log10(1.0 - tpe * (double)(res[i].score - lsc + 1) / (double)tsc)
		);
	}
	return(p);		/* #non-repetitive alignments */
}
#else
/**
 * @fn mm_post_map
 */
static _force_inline
uint64_t mm_post_map(
	mm_mapopt_t const *opt,
	uint32_t const n_reg, v4u32_t *reg)
{
	uint32_t i, score = _aln(reg[0])->score, min = score * opt->min_ratio;
	uint32_t ssc = (n_reg>1)? _aln(reg[1])->score : 0;
	uint32_t bsc = (n_reg>1)? _aln(reg[n_reg-1])->score : 0;
	uint32_t tsc = 0;

	double elen = (double)gaba_plen(_aln(reg[0])->sec) / 2.0;
	double pid = 1.0 - (double)(elen * opt->m - score) / (double)(opt->m + opt->x) / elen;
	double ec = 2.0 / (pid * (double)(opt->m + opt->x) - (double)opt->x);
	double ulen = ec * (score - ssc), pe = 1.0 / (ulen * ulen + (double)n_reg);

	for(i = 1; i < n_reg; i++) tsc += _aln(reg[i])->score - bsc + 1;
	reg[0].u32[0] |= _clip(-10.0 * MAPQ_COEF * log10(pe));
	for(i = 1; i < n_reg && (score = _aln(reg[i])->score) >= min; i++) {
		reg[i].u32[0] |=
			  (0x100<<16)
			| _clip(
				  -10.0
				* MAPQ_COEF
				* log10(1.0 - pe * (double)(score - bsc + 1) / (double)tsc)
			  );
	}
	return(1);		/* #non-repetitive alignments */
}
#endif

/**
 * @fn mm_post_ava
 * @brief all-versus-all mapq estimation
 */
static _force_inline
uint64_t mm_post_ava(
	mm_tbuf_t *self)
{
	mm_res_t *res = (mm_res_t *)self->chain.a;	/* alignments, must be sorted */

	uint32_t i, score, min = res[0].score * self->min_ratio;
	for(i = 0; i < self->n_res && (score = res[i].score) >= min; i++) {
		mm_bin_t *bin = (mm_bin_t *)&self->bin.a[res[i].bid];

		/* estimate effective length */
		double elen = (double)bin->plen / 2.0;
		double pid = 1.0 - (double)(elen * self->m - score) / (double)(self->m + self->x) / elen;
		double ec = 2.0 / (pid * (double)(self->m + self->x) - (double)self->x);
		double ulen = ec * score, pe = 1.0 / (ulen + 1);
		bin->plen = _clip(-10.0 * MAPQ_COEF * log10(pe));	// | ((i == 0? 0 : 0x800)<<16);
	}
	return(self->n_res);	/* #non-repetitive alignments */
}
#undef _clip
#undef _aln

/**
 * @fn mm_pack_reg
 * @brief allocate reg array from lmm, copy reg array to it
 */
static _force_inline
mm_reg_t const *mm_pack_reg(
	mm_tbuf_t *self,
	uint32_t n_all,
	uint32_t n_uniq)
{
	/* allocate mem from lmm */
	uint64_t size = sizeof(mm_reg_t) + self->bin.n * sizeof(mm_aln_t *);
	mm_reg_t *reg = lmm_malloc(self->trace.lmm, size);
	mm_aln_t **p = (mm_aln_t **)reg->aln, **b = p;

	/* build head */
	*reg = (mm_reg_t){ 0 };

	/* build reg array (copy from bin) */
	mm_res_t *res = (mm_res_t *)self->chain.a;
	for(uint64_t i = 0; i < n_all; i++) {
		/* copy alignment id and mapq at the head of gaba_alignment_t */
		mm_bin_t *bin = (mm_bin_t *)&self->bin.a[res[i].bid];
		for(uint64_t j = 0; j < bin->n_aln; j++) {
			mm_aln_t *a = (mm_aln_t *)bin->aln[j] - 1;		/* .head_margin = sizeof(mm_aln_t) */
			a->aid = i;
			a->mapq = bin->plen;
			*p++ = a;
/*
			fprintf(stderr, "pack_reg, i(%lu), j(%lu), p(%p), (%u, %u), pos(%u, %u), len(%u, %u), plen(%u)\n",
				i, j, a,
				self->r[0].len, self->q[0].len,
				a->a->sec->apos, a->a->sec->bpos,
				a->a->sec->alen, a->a->sec->blen,
				a->a->path->len);
*/
		}

		/* store #unique alignments */
		if(i == n_uniq - 1) { reg->n_uniq = p - b; }
	}

	/* store #total alignments */
	reg->n_all = p - b;
	return(reg);
}

/**
 * @fn mm_init_query
 */
static _force_inline
void mm_init_query(
	mm_tbuf_t *self,
	uint32_t const l_seq, uint8_t const *seq,	/* query sequence (read) length and pointer (must be 4-bit encoded) */
	uint32_t const qid,							/* query sequence id (used in all-versus-all filter) */
	lmm_t *restrict lmm)
{
	/* clear buffers */
	self->resc.n = 0;
	self->seed.n = 0;
	self->map.n = 0;
	self->chain.n = 0;
	self->bin.n = 0;
	kh_clear(&self->pos);
	self->n_res = 0;

	/* set query seq info */
	self->adj = l_seq + self->mi.k - 1;
	self->qid = qid;
	self->q[0] = _sec_fw(qid, seq, l_seq);
	self->q[1] = _sec_rv(qid, seq, l_seq);
	self->q[2] = _sec_fw(qid, seq, l_seq);

	/* store lmm */
	self->trace.lmm = (void *)lmm;
	return;
}

/**
 * @fn mm_align_seq
 * @brief alignment root function
 */
static _force_inline
mm_reg_t const *mm_align_seq(
	mm_tbuf_t *self,							/* thread-local context */
	uint32_t const l_seq, uint8_t const *seq,	/* query sequence (read) length and pointer (must be 4-bit encoded) */
	uint32_t const qid,							/* query sequence id (used in all-versus-all filter) */
	lmm_t *restrict lmm)						/* memory arena for results */
{
	/* skip unmappable */
	if(l_seq < self->mi.k || l_seq * self->m < self->min_score) {
		return(NULL);
	}

	/* clear buffers */
	mm_init_query(self, l_seq, seq, qid, lmm);

	debug("start: qid(%u)", qid);
	/* seed-chain-extend loop */
	debug("n_occ(%u)", self->n_occ);
	for(uint64_t i = 0; i < self->n_occ; i++) {
		if(mm_seed(self, i) == 0) {
			continue;		/* seed not found */
		}
		if(mm_chain(self) == 0) {
			continue;		/* chain not found */
		}
		if(mm_extend(self) > 0) {
			break;			/* at least one full-length alignment found */
		}
	}

	if(self->n_res == 0) {
		return(NULL);
	}

	/* sort by score in reverse order */
	radix_sort_64x((v2u32_t *)self->chain.a, self->n_res);

	/* prune alignments whose score is less than min_score threshold */
	uint32_t n_all = mm_prune_regs(self);

	/* collect supplementaries (split-read collection) */
	uint32_t n_uniq = ((self->flag & MM_AVA)? mm_post_ava : mm_post_map)(self);
	debug("n_all(%u), n_uniq(%u)", n_all, n_uniq);


	for(uint64_t i = 0; i < n_all; i++) {
		debug("bid(%u)", self->chain.a[i].sid);
		mm_bin_t *bin = (mm_bin_t *)&self->bin.a[self->chain.a[i].sid];

		debug("n_aln(%u), plen(%u), lb(%u), ub(%u)", bin->n_aln, bin->plen, bin->lb, bin->ub);

		for(uint64_t j = 0; j < bin->n_aln; j++) {
			debug("j(%llu), aln(%p)", j, bin->aln[j]);
		}

	}

	/* allocate reg array from memory arena */
	return(mm_pack_reg(self, n_all, n_uniq));
}

/**
 * @struct mm_tmpbuf_t
 q > 0 && * @brief 
 */
typedef struct {
	uint8_t *tail, *p;
	uint8_t base[240];
} mm_tmpbuf_t;

/**
 * @macro _flush
 * @brief flush the buffer if there is no room for(margin + 1) bytes
 */
#define _flush(_buf, _margin) { \
	if(_unlikely((uintptr_t)(_buf)->p + (_margin) >= (uintptr_t)(_buf)->tail)) { \
		fwrite((_buf)->base, sizeof(uint8_t), (_buf)->p - (_buf)->base, stdout); \
		(_buf)->p = (_buf)->base; \
	} \
}

/**
 * @macro _put
 * @brief put a byte to buf
 */
#define _put(_buf, _c) { \
	_flush(_buf, 1); \
	*(_buf)->p++ = (_c); \
}

/**
 * @macro _put
 * @brief format integer n with decimal point after c-th digit
 */
#define _putfi(type, _buf, _n, _c) ({ \
	uint64_t _b = 0; \
	type _m = (type)(_n); \
	int64_t _i = 0; \
	while(_m || _i <= _c) { _b <<= 4; _b += _m % 10, _m /= 10; _i++; } \
	_flush(_buf, _i + 1); \
	for(int64_t _j = _i; _j > (_c); _j--) { *(_buf)->p++ = (_b&0x0f) + '0'; _b>>=4; } \
	*(_buf)->p++ = '.'; \
	for(int64_t _j = (_c); _j > 0; _j--) { *(_buf)->p++ = (_b&0x0f) + '0'; _b>>=4; } \
	_i; \
})

/**
 * @macro _puti, _putn
 * @brief format integer
 */
#define _puti(type, _buf, _n) ({ \
	uint64_t _b = 0; \
	type _m = (type)(_n); \
	int64_t _i = 0; \
	while(_m) { _b <<= 4; _b += _m % 10, _m /= 10; _i++; } \
	_i += (_i==0); \
	_flush(_buf, _i); \
	for(int64_t _j = _i; _j > 0; _j--) { \
		*(_buf)->p++ = (_b&0x0f) + '0'; _b>>=4; \
	} \
	_i; \
})
#define _putn(_buf, _n) _puti(uint32_t, _buf, _n)

/**
 * @macro _putpi
 * @brief format pair of integers into same widths (used in maf formatter for vertical (column) alignment)
 * _buf1 will be flushed, _buf2 is never flushed nor boundary-tested
 */
#define _putpi(type, _buf1, _buf2, _n1, _n2) ({ \
	uint64_t _b1 = 0, _b2 = 0; \
	type _m1 = (type)(_n1), _m2 = (type)(_n2); \
	int64_t _i = 0; \
	while(_m1 | _m2) { \
		_b1 <<= 4; _b1 += _m1 % 10, _m1 /= 10; \
		_b2 <<= 4; _b2 += _m2 % 10, _m2 /= 10; \
		_i++; \
	} \
	_i += (_i==0); \
	_flush(_buf1, _i); \
	for(int64_t _j = _i, _z1 = 0, _z2 = 0; _j > 0; _j--) { \
		uint64_t _c1 = _b1&0x0f, _c2 = _b2&0x0f; \
		_z1 |= _c1 | (_j == 1); _z2 |= _c2 | (_j == 1); \
		*(_buf1)->p++ = _c1 + '0' - (_z1? 0 : 0x10); \
		*(_buf2)->p++ = _c2 + '0' - (_z2? 0 : 0x10); \
		_b1>>=4; _b2>>=4; \
	} \
	_i; \
})
#define _putpn(_buf1, _buf2, _n1, _n2)	_putpi(uint32_t, _buf1, _buf2, _n1, _n2)

/**
 * @macro _puts, _putsk
 * @brief dump string (_puts for '\0'-terminated, _putsk for const)
 */
#define _puts(_buf, _s) { \
	for(uint8_t const *_q = (uint8_t const *)(_s); *_q; _q++) { \
		_put(_buf, *_q); \
	} \
}
#define _putsk(_buf, _s) { \
	uint64_t const _l = strlen(_s); \
	_flush(_buf, _l); \
	for(uint8_t const *_q = (uint8_t const *)(_s), *_t = _q + (_l); _q < _t; _q++) { \
		*(_buf)->p++ = *(_q); \
	} \
}

/**
 * @macro _putsn, _putsnr
 * @brief dump string with simd copy (forward and reverse; for quality string dump)
 */
#define _putsn(_buf, _s, _l) { \
	uint8_t const *_q = (uint8_t const *)(_s); \
	for(uint8_t const *_t = _q + ((_l) & ~0x1fULL); _q < _t; _q += 32) { \
		_flush(_buf, 32); \
		_storeu_v32i8((_buf)->p, _loadu_v32i8(_q)); (_buf)->p += 32; \
	} \
	_flush(_buf, 32); \
	_storeu_v32i8((_buf)->p, _loadu_v32i8(_q)); (_buf)->p += (_l) & 0x1f; \
}
#define _putsnr(_buf, _s, _l) { \
	uint8_t const *_q = (uint8_t const *)(_s) + (_l) - 32; \
	for(; _q >= _s; _q -= 32) { \
		_flush(_buf, 32); \
		_storeu_v32i8((_buf)->p, _swap_v32i8(_loadu_v32i8(_q))); (_buf)->p += 32; \
	} \
	_flush(_buf, 32); \
	_storeu_v32i8((_buf)->p, _swap_v32i8(_loadu_v32i8(_q))); (_buf)->p += (_l) & 0x1f; \
}

/**
 * @macro _putsn32
 * @brief len must be shorter than 32
 */
#define _putsn32(_buf, _s, _len) { \
	_flush(_buf, 32); \
	_storeu_v32i8((_buf)->p, _loadu_v32i8(_s)); \
	(_buf)->p += (_len); \
}

/**
 * @macro _putscn
 * @brief dump constant pattern
 */
#define _putscn(_buf, _c, _l) { \
	v32i8_t const _pattern = _set_v32i8(_c); \
	for(uint64_t _i = 0; _i < ((_l) & ~0x1fULL); _i += 32) { \
		_flush(_buf, 32); \
		_storeu_v32i8((_buf)->p, _pattern); (_buf)->p += 32; \
	} \
	_flush(_buf, 32); \
	_storeu_v32i8((_buf)->p, _pattern); (_buf)->p += (_l) & 0x1f; \
}

/**
 * @macro _putsnt, _putsntr
 * @brief dump string with encoding conversion
 */
#define _putsnt(_buf, _s, _l, _table) { \
	v32i8_t register _conv = _from_v16i8_v32i8(_loadu_v16i8(_table)); \
	v32i8_t _r, _b; \
	uint8_t const *_q = (uint8_t const *)(_s); \
	for(uint8_t const *_t = _q + ((_l) & ~0x1fULL); _q < _t; _q += 32) { \
		_flush(_buf, 32); \
		_r = _loadu_v32i8(_q); \
		_b = _shuf_v32i8(_conv, _r); \
		_storeu_v32i8((_buf)->p, _b); \
		(_buf)->p += 32; \
	} \
	_flush(_buf, 32); \
	_r = _loadu_v32i8(_q); \
	_b = _shuf_v32i8(_conv, _r); \
	_storeu_v32i8((_buf)->p, _b); \
	(_buf)->p += (_l) & 0x1f; \
}
#define _putsntr(_buf, _s, _l, _table) { \
	v32i8_t register _conv = _from_v16i8_v32i8(_loadu_v16i8(_table)); \
	v32i8_t _r, _b; \
	uint8_t const *_q = (uint8_t const *)(_s) + (_l) - 32; \
	for(; _q >= _s; _q -= 32) { \
		_flush(_buf, 32); \
		_r = _loadu_v32i8(_q); \
		_b = _shuf_v32i8(_conv, _swap_v32i8(_r)); \
		_storeu_v32i8((_buf)->p, _b); \
		(_buf)->p += 32; \
	} \
	_flush(_buf, 32); \
	_r = _loadu_v32i8(_q); \
	_b = _shuf_v32i8(_conv, _swap_v32i8(_r)); \
	_storeu_v32i8((_buf)->p, _b); \
	(_buf)->p += (_l) & 0x1f; \
}

/**
 * @macro _putsnt32, _putsntr32
 * @brief dump short string with encoding conversion
 */
#define _putsnt32(_buf, _q, _len, _table) { \
	_flush(_buf, 32); \
	v32i8_t register _conv = _from_v16i8_v32i8(_loadu_v16i8(_table)); \
	_storeu_v32i8((_buf)->p, _shuf_v32i8(_conv, _loadu_v32i8(_q))); \
	(_buf)->p += (_len); \
}
#define _putsntr32(_buf, _q, _len, _table) { \
	_flush(_buf, 32); \
	v32i8_t register _conv = _from_v16i8_v32i8(_loadu_v16i8(_table)); \
	_storeu_v32i8((_buf)->p, _shuf_v32i8(_conv, _swap_v32i8(_loadu_v32i8(_q)))); \
	(_buf)->p += (_len); \
}

/**
 * @macro _putd
 * @brief print direction ('+' for forward, '-' for reverse)
 */
#define _putds(b, _id)		_put(b, ((_id) & 0x01) ? '-' : '+');
#define _putdn(b, _id)		_put(b, ((_id) & 0x01) ? '1' : '0');

/**
 * @macro _t, _c, _sp, _cr
 * @brief print delimiters
 */
#define _t(b)			_put(b, '\t')
#define _c(b)			_put(b, ',')
#define _sp(b)			_put(b, ' ')
#define _cr(b)			_put(b, '\n')

/**
 * @macro _aln, _mapq, _flag
 * @brief type cast utilities
 */
#define _aln(x)		( (gaba_alignment_t const *)(x).u64[1] )
#define _mapq(x)	( (x).u32[0] & 0xffff )
#define _flag(x)	( (x).u32[0]>>16 )

/**
 * @fn mm_load_uint64
 * @brief bit string loader (for path-string parsing)
 */
static inline
uint64_t mm_load_uint64(
	uint64_t const *ptr,
	int64_t pos)
{
	int64_t rem = pos & 63;
	uint64_t a = (ptr[pos>>6]>>rem) | ((ptr[(pos>>6) + 1]<<(63 - rem))<<1);
	return(a);
}
#define _load_uint64(_ptr, _pos) ({ \
	uint64_t _rem = (_pos) & 0x3f; \
	((_ptr)[(_pos)>>6]>>_rem) | (((_ptr)[((_pos)>>6) + 1]<<(63 - _rem))<<1); \
})

/**
 * @fn mm_cigar_printer
 * @brief callback function passed to gaba_dp_print_cigar_*
 */
static
int mm_cigar_printer(void *_b, int64_t len, char c)
{
	mm_align_t *b = (mm_align_t *)_b;
	if(len < 41) {
		_flush(b, 16);
		*b->p++ = (b->conv[len-1]&0x0f) + '0';
		*b->p++ = (b->conv[len-1]>>4) + '0';
		b->p -= len<10;
		*b->p++ = c;
	} else {
		len = _putn(b, (uint32_t)len); _put(b, c);
	}
	return(len+1);
}

/* sam tags */
#define MM_RG			( 0 )		/* Z: read group */
#define MM_CO			( 1 )		/* Z: comment */
#define MM_NH			( 2 )		/* i: #hits */
#define MM_IH 			( 3 )		/* i: index of the record in the hits */
#define MM_AS			( 4 )		/* i: score */
#define MM_XS			( 5 )		/* i: suboptimal score */
#define MM_NM 			( 6 )		/* i: editdist to the reference */
#define MM_SA			( 7 )		/* Z: supplementary records */
#define MM_MD			( 8 )		/* Z: mismatch positions */

/**
 * @fn mm_print_sam_num
 * @brief print number in sam tag format
 */
static _force_inline
uint64_t mm_print_sam_num(
	mm_align_t *b,
	uint8_t type,
	uint8_t const *p)
{
	switch (type) {
		case 'a': _put(b, *p); return(1);
		case 'c': _puti(int8_t, b, (int8_t)*p); return(1);
		case 'C': _puti(uint8_t, b, (uint8_t)*p); return(1);
		case 's': _puti(int16_t, b, *((int16_t*)p)); return(2);
		case 'S': _puti(uint16_t, b, *((uint16_t*)p)); return(2);
		case 'i': _puti(int32_t, b, *((int32_t*)p)); return(4);
		case 'I': _puti(uint32_t, b, *((uint32_t*)p)); return(4);
		case 'f': {
			char buf[32]; uint64_t l = sprintf(buf, "%f", *((float*)p));
			for(uint64_t i = 0; i < l; i++) {
				_put(b, buf[i]);
			}
			return(4);
		}
	}
	return(1);
}

/**
 * @fn mm_restore_sam_tags
 * @brief print saved sam tags
 */
static _force_inline
void mm_restore_sam_tags(
	mm_align_t *b,
	bseq_t const *t)
{
	static uint8_t const tag_size[32] = {
		0, 1, 0xfe, 1,  0, 0, 4, 0,  0xfe, 4, 0, 0,  0, 0, 0, 0,
		0, 0, 0, 2,     0, 0, 0, 0,  0, 0, 0xff, 0,  0, 0, 0, 0,
	};
	// uint8_t const *p = t->tag, *tail = p + t->l_tag;
	uint8_t const *p = t->tag;
	uint64_t n_tag = t->n_tag;
	// while(p < tail) {
	while(n_tag-- > 0) {
		/* print tag label */
		_t(b); _put(b, p[0]); _put(b, p[1]); _put(b, ':'); _put(b, p[2]); _put(b, ':');

		/* print body */
		if(p[2] == 'Z') {
			/* string */
			p += 3;
			while(*p) {
				_put(b, *p); p++;
			}
			p++;
		} else if(tag_size[p[2]&0x1f] == 0xfe) {
			/* array */
			uint8_t type = p[3];
			uint32_t len = *((uint32_t*)&p[4]); p += 8;
			for(uint64_t i = 0; i < len; i++) {
				p += mm_print_sam_num(b, type, p);
				_put(b, ',');
			}
		} else {
			/* number */
			p += mm_print_sam_num(b, p[2], p + 3) + 3;
		}
	}
	return;
}

/**
 * @fn mm_print_sam_header
 * @brief print sam header
 */
static
void mm_print_sam_header(
	mm_align_t *b)
{
	mm_idx_t const *mi = b->mi;
	/* header */
	_putsk(b, "@HD\tVN:1.0\tSO:unsorted\n");

	/* sequences */
	for(uint64_t i = 0; i < mi->s.n; i++) {
		_putsk(b, "@SQ\tSN:");
		_putsn(b, mi->s.a[i].name, mi->s.a[i].l_name);		/* sequence name */
		_putsk(b, "\tLN:");
		_putn(b, mi->s.a[i].l_seq);							/* sequence length */
		_cr(b);
	}

	/* print read group */
	if(b->opt->flag & 0x01ULL<<MM_RG) {
		_puts(b, b->opt->rg_line);
		_cr(b);
	}

	/* program info (note: command line and version should be included) */
	_putsk(b, "@PG\tID:minialign\tPN:minialign\tVN:");
	_puts(b, version());
	_putsk(b, "\tCL:");
	_puts(b, b->opt->arg_line);
	_cr(b);
	return;
}

/**
 * @fn mm_print_sam_unmapped
 * @brief print sam unmapped record
 */
static
void mm_print_sam_unmapped(
	mm_align_t *b,
	bseq_t const *t)
{
	/* print name */
	_putsn(b, t->name, t->l_name);

	/* unmapped fixed string */
	_putsk(b, "\t4\t*\t0\t0\t*\t*\t0\t0\t");

	/* sequence */
	_putsnt(b, t->seq, t->l_seq, "NACMGRSVTWYHKDBN");
	_t(b);

	/* quality if available */
	if(b->opt->flag&MM_KEEP_QUAL && t->qual[0] != '\0') {
		_putsn(b, t->qual, t->l_seq);
	} else {
		_put(b, '*');
	}

	/* tags */
	mm_restore_sam_tags(b, t); _cr(b);
	return;
}

/**
 * @fn mm_print_sam_mapped_core
 */
static _force_inline
void mm_print_sam_mapped_core(
	mm_align_t *b,
	mm_idx_seq_t const *r,
	bseq_t const *t,
	mm_aln_t const *a,
	uint16_t flag)
{
	const gaba_path_section_t *s = &a->a->sec[0];
	uint32_t rs = s->apos;									/* ref start pos */
	uint32_t hl = s->bpos, tl = t->l_seq - s->bpos - s->blen;/* head and tail clip lengths */
	uint32_t qs = (flag & 0x900) ? hl : 0;					/* query start pos */
	uint32_t qe = t->l_seq - ((flag & 0x900) ? tl : 0);		/* query end pos */

	_putsn(b, t->name, t->l_name); _t(b);					/* qname */
	_putn(b, flag | ((s->bid & 0x01)<<4)); _t(b);			/* flag */
	_putsn(b, r->name, r->l_name); _t(b);					/* tname (ref name) */
	_putn(b, rs + 1); _t(b);								/* mapped pos */
	_putn(b, a->mapq>>MAPQ_DEC); _t(b);						/* mapping quality */

	/* print cigar */
	if(hl) {
		_putn(b, hl);
		_put(b, (flag&0x900)? 'H' : 'S');					/* print head clip */
	}

	// fprintf(stderr, "ppos(%u)\n", a->a->sec->ppos);
	gaba_dp_print_cigar_forward(mm_cigar_printer, b, a->a->path->array, 0, a->a->path->len);
	if(tl) {
		_putn(b, tl);
		_put(b, (flag&0x900)? 'H' : 'S');					/* print tail clip */
	}
	_putsk(b, "\t*\t0\t0\t");								/* mate tags, unused */

	/* print sequence */
	if(s->bid & 0x01) {
		_putsntr(b, &t->seq[t->l_seq-qe], qe-qs, "NTGKCYSBAWRDMHVN");
	} else {
		_putsnt(b, &t->seq[qs], qe-qs, "NACMGRSVTWYHKDBN");
	}
	_t(b);

	/* print quality string if available */
	if(b->opt->flag&MM_KEEP_QUAL && t->qual[0] != '\0') {
		if(s->bid & 0x01) {
			_putsnr(b, &t->qual[t->l_seq-qe], qe-qs);
		} else {
			_putsn(b, &t->qual[qs], qe-qs);
		}
	} else {
		_put(b, '*');
	}
	return;
}

/**
 * @fn mm_print_sam_supp
 * @brief print sam supplementary alignment (SA) tag
 */
static _force_inline
void mm_print_sam_supp(
	mm_align_t *b,
	mm_idx_seq_t const *r,
	bseq_t const *t,
	mm_aln_t const *a)
{
	/* rname,pos,strand,CIGAR,mapQ,NM; */
	const gaba_path_section_t *s = &a->a->sec[0];
	uint32_t rs = s->apos;
	uint32_t hl = s->bpos, tl = s->bpos + s->blen;			/* head and tail clips */
	
	_putsn(b, r->name, r->l_name); _c(b);					/* rname */
	_putn(b, rs + 1); _c(b);								/* rpos */
	_put(b, (s->bid & 0x01) ? '+' : '-'); _c(b);			/* direction */

	/* print cigar */
	if(hl != 0) { _putn(b, hl); _put(b, 'H'); }				/* always hard clipped */
	gaba_dp_print_cigar_forward(mm_cigar_printer, b, a->a->path->array, 0, a->a->path->len);
	if(tl != 0) { _putn(b, tl); _put(b, 'H'); }
	_c(b);

	/* mapping quality and editdist */
	_putn(b, a->mapq); _c(b);
	_putn(b, a->a->xcnt + a->a->gecnt); _put(b, ';');
	return;
}

/**
 * @fn mm_print_sam_md
 * @brief print sam MD tag
 */
static _force_inline
void mm_print_sam_md(
	mm_align_t *b,
	mm_idx_seq_t const *r,
	bseq_t const *t,
	mm_aln_t const *a)
{
	/*
	 * print MD tag
	 * note: this operation requires head-to-tail comparison of two sequences,
	 *       which may result in several percent performance degradation.
	 */
	uint64_t const f = b->opt->flag;
	if((f & 0x01ULL<<MM_MD) == 0) { return; }				/* skip if disabled */

	const gaba_path_section_t *s = &a->a->sec[0];
	uint32_t rs = s->apos, qs = s->bpos;

	_puts(b, "\tMD:Z:");
	// uint64_t const *p = (uint64_t const *)a->a->path->array;
	// int64_t pos = 0, lim = a->a->path->len;

	uint64_t const *p = (uint64_t const *)((uint64_t)a->a->path->array & ~(sizeof(uint64_t) - 1));
	uint64_t lim = (((uint64_t)a->a->path->array & sizeof(uint32_t)) ? 32 : 0) + a->a->path->len;
	uint64_t ridx = a->a->path->len;

	uint64_t dir = (s->bid & 0x01)? 1 : 0;
	static uint8_t const comp[16] __attribute__(( aligned(16) )) = {
		0x00, 0x08, 0x04, 0x0c, 0x02, 0x0a, 0x06, 0x0e,
		0x01, 0x09, 0x05, 0x0d, 0x03, 0x0b, 0x07, 0x0f
	};
	static uint8_t const id[16] __attribute__(( aligned(16) )) = {
		0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
		0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f
	};
	v32i8_t const cv = _from_v16i8_v32i8(_load_v16i8(dir ? comp : id));
	uint8_t const *rp = &r->seq[rs], *rb = rp, *qp = &t->seq[dir ? (uint64_t)t->l_seq - qs - 32 : qs];
	while(ridx > 0) {
		/* suppose each indel block is shorter than 32 bases */
		uint64_t arr = _load_uint64(p, lim - ridx);
		uint64_t cnt = tzcnt(~arr) - (arr & 0x01);		/* count #ins */
		ridx -= cnt;
		qp += dir ? -cnt : cnt;

		if((arr & 0x01) == 0) {							/* is_del */
			ridx -= cnt = tzcnt(arr);					/* count #del */
			_putn(b, (int32_t)(rp - rb));
			_put(b, '^');
			_putsnt32(b, rp, cnt, "NACMGRSVTWYHKDBN");
			rp += cnt;
			rb = rp;
		}

		/* match or mismatch */
		uint64_t acnt = 32;
		while(acnt == 32) {
			/* count diagonal */
			arr = _load_uint64(p, lim - ridx);
			acnt = MIN2(tzcnt(arr ^ 0x5555555555555555), ridx)>>1;

			/* load sequence to detect mismatch */
			v32i8_t rv = _shuf_v32i8(cv, _loadu_v32i8(rp)), qv = _loadu_v32i8(qp);
			if(dir) { qv = _swap_v32i8(qv); }

			/* compare and count matches */
			uint64_t mmask = (uint64_t)((v32_masku_t){ .mask = _mask_v32i8(_eq_v32i8(rv, qv)) }).all;
			uint64_t mcnt = MIN2(acnt, tzcnt(~mmask));

			/* adjust pos */
			ridx -= 2*mcnt;
			rp += mcnt;
			qp += dir ? -mcnt : mcnt;

			if(mcnt >= acnt) { continue; }				/* continues longer than 32bp */
			_putn(b, (int32_t)(rp - rb));				/* print match length */
			ridx -= 2*(cnt = MIN2(tzcnt(mmask>>mcnt), acnt - mcnt));
			qp += dir ? -cnt : cnt;
			for(uint64_t i = 0; i < cnt - 1; i++) {
				_put(b, "NACMGRSVTWYHKDBN"[*rp]);		/* print mismatch base */
				_put(b, '0');							/* padding */
				rp++;
			}
			_put(b, "NACMGRSVTWYHKDBN"[*rp]);			/* print mismatch base */
			rp++;
			rb = rp;
		}
	}
	_putn(b, (int32_t)(rp - rb));						/* always print tail even if # == 0 */
	return;
}

/**
 * @fn mm_print_sam_general_tags
 */
static _force_inline
void mm_print_sam_general_tags(
	mm_align_t *b,
	mm_aln_t const *a,
	uint32_t n_reg,
	uint64_t j)
{
	uint64_t const f = b->opt->flag;

	/* read group */
	if(f & 0x01ULL<<MM_RG) {
		_putsk(b, "\tRG:Z:"); _puts(b, b->opt->rg_id);
	}

	/* #hits, including repetitive */
	if(f & 0x01ULL<<MM_NH) {
		_putsk(b, "\tNH:i:"); _putn(b, n_reg);			/* note: (f & MM_OMIT_REP)? n_uniq : n_reg is better? */
	}

	/* record index */
	if(f & 0x01ULL<<MM_IH) {
		_putsk(b, "\tIH:i:"); _putn(b, j);
	}

	/* alignment score */
	if(f & 0x01ULL<<MM_AS) {
		_putsk(b, "\tAS:i:"); _putn(b, a->a->score);
	}

	/* editdist to ref */
	if(f & 0x01ULL<<MM_NM) {
		_putsk(b, "\tNM:i:");
		_putn(b, a->a->xcnt + a->a->gecnt);
	}
	return;
}

/**
 * @fn mm_print_sam_primary_tags
 * @brief print tags only appear in primary alignment, return non-zero if other alignments should be omitted
 */
static _force_inline
uint64_t mm_print_sam_primary_tags(
	mm_align_t *b,
	bseq_t const *t,				/* query sequence */
	mm_reg_t const *reg)			/* alignments, 1..n_uniq is printed in SA */
{
	mm_idx_t const *mi = b->mi;
	uint64_t const f = b->opt->flag;
	uint64_t ret = 0;

	/* second best hit score */
	if((f & 0x01ULL<<MM_XS)) {
		_putsk(b, "\tXS:i:");
		_putn(b, reg->n_all > 1 ? reg->aln[1]->a->score : 0);
	}

	/* supplementary alignment */
	if(f & 0x01ULL<<MM_SA && reg->n_uniq > 1) {
		_putsk(b, "\tSA:Z:");

		/* iterate over the supplementary alignments */
		for(uint64_t i = 1; i < reg->n_uniq; i++) {
			mm_aln_t const *a = reg->aln[i];
			uint32_t rid = a->a->sec->aid>>1;
			mm_idx_seq_t const *r = &mi->s.a[rid - mi->base_rid];

			mm_print_sam_supp(b, r, t, a);
		}
		ret = 1;
	}

	/* saved sam tags */
	mm_restore_sam_tags(b, t);
	return(ret);
}

/**
 * @fn mm_print_sam_mapped
 */
static
void mm_print_sam_mapped(
	mm_align_t *b,
	bseq_t const *t,
	mm_reg_t const *reg)
{
	mm_idx_t const *mi = b->mi;

	/* iterate over alignment sets */
	uint64_t const n = (b->opt->flag & MM_OMIT_REP)? reg->n_uniq : reg->n_all;
	for(uint64_t i = 0, flag = 0; i < n; i++) {
		if(i >= reg->n_uniq) {
			flag = 0x100;				/* mark secondary */
		}

		mm_aln_t const *a = reg->aln[i];
		uint32_t rid = a->a->sec->aid>>1;
		mm_idx_seq_t const *r = &mi->s.a[rid - mi->base_rid];

		/* print body */
		mm_print_sam_mapped_core(b, r, t, a, flag);

		/* print general tags (scores, ...) */
		mm_print_sam_general_tags(b, a, reg->n_all, i);

		/* mismatch position (MD) */
		mm_print_sam_md(b, r, t, a);

		/* primary-specific tags */
		if(i == 0 && mm_print_sam_primary_tags(b, t, reg)) {
			i = n;						/* skip supplementary records when SA tag is enabled */
		}
		_cr(b);

		flag = 0x800;					/* mark supplementary */
	}
	return;
}

/**
 * @fn mm_print_maf_mapped_ref
 * @brief print header line of a maf record
 */
static _force_inline
void mm_print_maf_mapped_ref(
	mm_align_t *b,
	mm_aln_t const *a,
	uint8_t const *rp)
{
	uint64_t const *p = (uint64_t const *)((uint64_t)(a->a->path->array - 1) & ~(sizeof(uint64_t) - 1));
	int64_t pos = a->a->path->len;
	while(pos > 0) {
		uint64_t arr, cnt;

		/* insertion */
		arr = _load_uint64(p, pos);
		pos -= (cnt = lzcnt(~arr));		/* count 1's (insertions) */
		_putsn32(b, "--------------------------------", cnt);

		/* deletion */
		arr = _load_uint64(p, pos);
		pos -= (cnt = lzcnt(arr) - 1);	/* count 0's (deletions) */
		_putsnt32(b, rp, cnt, "NACMGRSVTWYHKDBN");
		rp += cnt;						/* advance pointer */
		do {
			arr = _load_uint64(p, pos);
			cnt = lzcnt(arr^0x5555555555555555)>>1;		/* count diagonals */
			pos -= 2*cnt;
			_putsnt32(b, rp, cnt, "NACMGRSVTWYHKDBN");
			rp += cnt;					/* advance pointer */
		} while(cnt == 32);
	}
	_cr(b);
	return;
}

/**
 * @fn mm_print_maf_query_forward
 */
static _force_inline
void mm_print_maf_query_forward(
	mm_align_t *b,
	mm_aln_t const *a,
	uint8_t const *qp)
{
	uint64_t const *p = (uint64_t const *)((uint64_t)(a->a->path->array - 1) & ~(sizeof(uint64_t) - 1));
	int64_t pos = a->a->path->len;
	while(pos > 0) {
		uint64_t arr, cnt;

		/* insertions */
		arr = _load_uint64(p, pos);
		pos -= (cnt = lzcnt(~arr));
		_putsnt32(b, qp, cnt, "NACMGRSVTWYHKDBN");
		qp += cnt;

		/* deletions */
		arr = _load_uint64(p, pos);
		pos -= (cnt = lzcnt(arr) - 1);
		_putsn32(b, "--------------------------------", cnt);
		
		/* diagonals */
		do {
			arr = _load_uint64(p, pos);
			cnt = lzcnt(arr^0x5555555555555555)>>1;
			pos -= 2*cnt;
			_putsnt32(b, qp, cnt, "NACMGRSVTWYHKDBN");
			qp += cnt;
		} while(cnt == 32);
	}
	return;
}

/**
 * @fn mm_print_maf_query_reverse
 */
static _force_inline
void mm_print_maf_query_reverse(
	mm_align_t *b,
	mm_aln_t const *a,
	uint8_t const *qp)
{
	uint64_t const *p = (uint64_t const *)((uint64_t)(a->a->path->array - 1) & ~(sizeof(uint64_t) - 1));
	int64_t pos = a->a->path->len;
	while(pos > 0) {
		uint64_t arr, cnt;

		/* insertions */
		arr = _load_uint64(p, pos);
		pos -= (cnt = lzcnt(~arr));
		_putsntr32(b, qp, cnt, "NTGKCYSBAWRDMHVN");
		qp -= cnt;

		/* deletions */
		arr = _load_uint64(p, pos);
		pos -= (cnt = lzcnt(arr) - 1);
		_putsn32(b, "--------------------------------", cnt);
		
		/* diagonals */
		do {
			arr = _load_uint64(p, pos);
			cnt = lzcnt(arr^0x5555555555555555)>>1;
			pos -= 2*cnt;
			_putsntr32(b, qp, cnt, "NTGKCYSBAWRDMHVN");
			qp -= cnt;
		} while(cnt == 32);
	}
	return;
}

/**
 * @fn mm_print_maf_mapped_core
 */
static _force_inline
void mm_print_maf_mapped_core(
	mm_align_t *b,
	mm_idx_seq_t const *r,
	bseq_t const *t,
	mm_aln_t const *a)
{
	const gaba_path_section_t *s = &a->a->sec[0];

	/* header */
	_put(b, 'a'); _sp(b); _putsk(b, "score="); _putn(b, a->a->score); _cr(b);

	mm_tmpbuf_t qb;			/* leave buffers uninitialized */
	qb.p = qb.base; qb.tail = qb.base + 240;

	uint32_t const rs = s->apos;
	uint32_t const qs = s->bpos;

	/* heads of sequence record */
	_put(b, 's'); _sp(b);
	_put(&qb, 's'); _sp(&qb);

	/* names */
	uint64_t l = MAX2(r->l_name, t->l_name);
	_putscn(b, ' ', l - r->l_name); _putsn(b, r->name, r->l_name); _sp(b);
	_putscn(&qb, ' ', l - t->l_name); _putsn(&qb, t->name, t->l_name); _sp(&qb);

	/* alignment length */
	_putpn(b, &qb, rs, qs); _sp(b); _sp(&qb);
	_putpn(b, &qb, s->alen, s->blen); _sp(b); _sp(&qb);

	/* directions */
	_put(b, '+'); _sp(b);
	_putds(&qb, s->bid); _sp(&qb);

	/* sequence lengths */
	_putpn(b, &qb, r->l_seq, t->l_seq); _sp(b); _sp(&qb);

	/* reference alignment */
	mm_print_maf_mapped_ref(b, a, &r->seq[rs]);

	/* query */
	_putsn(b, qb.base, qb.p - qb.base);
	if(s->bid & 0x01) {
		mm_print_maf_query_reverse(b, a, &t->seq[(uint64_t)t->l_seq - qs - 32]);
	} else {
		mm_print_maf_query_forward(b, a, &t->seq[qs]);
	}
	_cr(b);
	_cr(b);

}

/**
 * @fn mm_print_maf_mapped
 */
static
void mm_print_maf_mapped(
	mm_align_t *b,
	bseq_t const *t,
	mm_reg_t const *reg)
{
	mm_idx_t const *mi = b->mi;
	uint64_t const n = (b->opt->flag & MM_OMIT_REP)? reg->n_uniq : reg->n_all;
	for(uint64_t i = 0; i < n; i++) {
		mm_aln_t const *a = reg->aln[i];
		uint32_t rid = a->a->sec->aid>>1;
		mm_idx_seq_t const *r = &mi->s.a[rid - mi->base_rid];

		mm_print_maf_mapped_core(b, r, t, a);
	}
	return;
}

/**
 * @fn mm_print_blast6_mapped
 * @brief blast.6 (tabular) format
 * qname rname idt len #x #gi qs qe rs re e-value bitscore
 */
static
void mm_print_blast6_mapped(
	mm_align_t *b,
	bseq_t const *t,
	mm_reg_t const *reg)
{
	mm_idx_t const *mi = b->mi;
	uint64_t const n = (b->opt->flag & MM_OMIT_REP)? reg->n_uniq : reg->n_all;
	for(uint64_t i = 0; i < n; i++) {
		mm_aln_t const *a = reg->aln[i];
		const gaba_path_section_t *s = &a->a->sec[0];
		mm_idx_seq_t const *r = &mi->s.a[(s->aid>>1) - mi->base_rid];

		int32_t dcnt = (a->a->path->len - a->a->gecnt)>>1, slen = dcnt + a->a->gecnt;
		int32_t mid = 100000.0 * (double)(dcnt - a->a->xcnt) / (double)slen;	/* percent identity */
		uint32_t rs = s->bid&0x01? r->l_seq - s->apos : s->apos + 1;
		uint32_t re = s->bid&0x01? s->apos + 1 : r->l_seq - s->apos;
		uint32_t qs = s->bpos + 1, qe = s->bpos + s->blen;

		/* sequence names */
		_putsn(b, t->name, t->l_name); _t(b);
		_putsn(b, r->name, r->l_name); _t(b);

		_putfi(int32_t, b, mid, 3); _t(b);				/* sequence identity */
		_putn(b, slen); _t(b);							/* alignment length */
		_putn(b, a->a->xcnt); _t(b);						/* mismatch count */
		_putn(b, a->a->gicnt); _t(b);					/* gap count */
		
		/* positions */
		_putn(b, qs); _t(b);
		_putn(b, qe); _t(b);
		_putn(b, rs); _t(b);
		_putn(b, re); _t(b);

		/* estimated e-value and bitscore */
		double bit = 1.85 * (double)a->a->score - 0.02;	/* fixme, lambda and k should be calcd from the scoring params */
		int32_t e = 1000.0 * (double)r->l_seq * (double)t->l_seq * pow(2.0, -bit);	/* fixme, lengths are not corrected */
		_putfi(int32_t, b, e, 3); _t(b);
		_putn(b, (int32_t)bit); _cr(b);
	}
	return;
}

/**
 * @fn mm_print_blasr1_mapped
 * @brief blasr.m1 format
 * qname rname qd rd score idt rs re rl qs qe ql #cells
 */
static
void mm_print_blasr1_mapped(
	mm_align_t *b,
	bseq_t const *t,
	mm_reg_t const *reg)
{
	mm_idx_t const *mi = b->mi;
	uint64_t const n = (b->opt->flag & MM_OMIT_REP)? reg->n_uniq : reg->n_all;
	for(uint64_t i = 0; i < n; i++) {
		mm_aln_t const *a = reg->aln[i];
		const gaba_path_section_t *s = &a->a->sec[0];
		mm_idx_seq_t const *r = &mi->s.a[(s->aid>>1) - mi->base_rid];

		int32_t dcnt = (a->a->path->len - a->a->gecnt)>>1, slen = dcnt + a->a->gecnt;
		int32_t mid = 1000000.0 * (double)(dcnt - a->a->xcnt) / (double)slen;	/* percent identity */
		// uint32_t rs = s->bid & 0x01 ? r->l_seq - s->apos - s->alen : s->apos, re = rs + s->alen;
		uint32_t rs = s->apos, re = s->apos + s->alen;
		uint32_t qs = s->bid & 0x01 ? t->l_seq - s->bpos - s->blen : s->bpos, qe = qs + s->blen;

		/* sequence names */
		_putsn(b, t->name, t->l_name); _sp(b);
		_putsn(b, r->name, r->l_name); _sp(b);

		/* directions */
		_put(b, '0'); _sp(b);
		_putdn(b, s->bid); _sp(b);

		/* score in negative, identity */
		_put(b, '-'); _putn(b, a->a->score); _sp(b);
		_putfi(int32_t, b, mid, 4); _sp(b);

		/* positions and lengths */
		_putn(b, rs); _sp(b);
		_putn(b, re); _sp(b);
		_putn(b, r->l_seq); _sp(b);
		_putn(b, qs); _sp(b);
		_putn(b, qe); _sp(b);
		_putn(b, t->l_seq); _sp(b);

		/* #cells, always zero */
		_put(b, '0'); _cr(b);
	}
	return;
}

/**
 * @fn mm_print_blasr4_mapped
 * @brief blasr.m4 format
 * qname rname score idt qd qs qe ql rd rs re rl mapq
 */
static
void mm_print_blasr4_mapped(
	mm_align_t *b,
	bseq_t const *t,
	mm_reg_t const *reg)
{
	mm_idx_t const *mi = b->mi;
	uint64_t const n = (b->opt->flag & MM_OMIT_REP)? reg->n_uniq : reg->n_all;
	for(uint64_t i = 0; i < n; i++) {
		mm_aln_t const *a = reg->aln[i];
		const gaba_path_section_t *s = &a->a->sec[0];
		mm_idx_seq_t const *r = &mi->s.a[(s->aid>>1) - mi->base_rid];

		int32_t dcnt = (a->a->path->len - a->a->gecnt)>>1, slen = dcnt + a->a->gecnt;
		int32_t mid = 1000000.0 * (double)(dcnt - a->a->xcnt) / (double)slen;	/* percent identity */
		// uint32_t rs = s->bid&0x01? r->l_seq-s->apos-s->alen : s->apos, re = rs+s->alen;
		uint32_t rs = s->apos, re = s->apos + s->alen;
		uint32_t qs = s->bid & 0x01 ? t->l_seq - s->bpos - s->blen : s->bpos, qe = qs + s->blen;

		/* sequence names */
		_putsn(b, t->name, t->l_name); _sp(b);
		_putsn(b, r->name, r->l_name); _sp(b);

		/* score in negative, identity */
		_put(b, '-'); _putn(b, a->a->score); _sp(b);
		_putfi(int32_t, b, mid, 4); _sp(b);

		/* query-side */
		_put(b, '0'); _sp(b);
		_putn(b, qs); _sp(b);
		_putn(b, qe); _sp(b);
		_putn(b, t->l_seq); _sp(b);

		/* reference-side */
		_putdn(b, s->bid); _sp(b);
		_putn(b, rs); _sp(b);
		_putn(b, re); _sp(b);
		_putn(b, r->l_seq); _sp(b);

		/* mapping quality */
		_putn(b, (a->mapq>>(MAPQ_DEC - 2)) + ((a->mapq>>(MAPQ_DEC + 2)) & ~0x01)); _cr(b);
	}
	return;
}

/**
 * @fn mm_print_paf_mapped
 * @brief minimap paf format
 * qname ql qs qe qd rname rl rs re #m block_len mapq
 */
static
void mm_print_paf_mapped(
	mm_align_t *b,
	bseq_t const *t,
	mm_reg_t const *reg)
{
	mm_idx_t const *mi = b->mi;
	uint64_t const n = (b->opt->flag & MM_OMIT_REP)? reg->n_uniq : reg->n_all;
	for(uint64_t i = 0; i < n; i++) {
		mm_aln_t const *a = reg->aln[i];
		const gaba_path_section_t *s = &a->a->sec[0];
		mm_idx_seq_t const *r = &mi->s.a[(s->aid>>1) - mi->base_rid];

		uint32_t dcnt = (a->a->path->len - a->a->gecnt)>>1;			/* #matches + #mismatches */
		uint32_t rs = s->apos, re = s->apos + s->alen;
		uint32_t qs = s->bpos, qe = s->bpos + s->blen;

		/* query */
		_putsn(b, t->name, t->l_name); _t(b);
		_putn(b, t->l_seq); _t(b);
		_putn(b, qs); _t(b);
		_putn(b, qe); _t(b);
		_putds(b, s->bid); _t(b);

		/* reference */
		_putsn(b, r->name, r->l_name); _t(b);
		_putn(b, r->l_seq); _t(b);
		_putn(b, rs); _t(b);
		_putn(b, re); _t(b);

		/* #matches, block length, mapping quality */
		_putn(b, dcnt - a->a->xcnt); _t(b);
		_putn(b, dcnt + a->a->gecnt); _t(b);
		_putn(b, a->mapq>>MAPQ_DEC);

		/* print optional tags */
		uint64_t const f = b->opt->flag;
		if(f & 0x01ULL<<MM_AS) { _putsk(b, "\tAS:i:"); _putn(b, a->a->score); }
		if(f & 0x01ULL<<MM_NM) { _putsk(b, "\tNM:i:"); _putn(b, a->a->xcnt + a->a->gecnt); }
		_cr(b);
	}
	return;
}

/**
 * @fn mm_print_mhap_mapped
 * @brief MHAP format
 * qname rname 1-idt score qd qs qe ql rd rs rd rl
 */
static
void mm_print_mhap_mapped(
	mm_align_t *b,
	bseq_t const *t,
	mm_reg_t const *reg)
{
	mm_idx_t const *mi = b->mi;
	uint64_t const n = (b->opt->flag & MM_OMIT_REP)? reg->n_uniq : reg->n_all;
	for(uint64_t i = 0; i < n; i++) {
		mm_aln_t const *a = reg->aln[i];
		const gaba_path_section_t *s = &a->a->sec[0];
		mm_idx_seq_t const *r = &mi->s.a[(s->aid>>1) - mi->base_rid];

		int32_t dcnt = (a->a->path->len - a->a->gecnt)>>1, slen = dcnt + a->a->gecnt;
		int32_t merr = 10000.0 * (1.0 - (double)(dcnt - a->a->xcnt) / (double)slen);	/* error rate */
		// uint32_t rs = s->bid&0x01? r->l_seq-s->apos-s->alen : s->apos, re = rs+s->alen;
		uint32_t rs = s->apos, re = s->apos + s->alen;
		uint32_t qs = s->bid & 0x01 ? t->l_seq - s->bpos - s->blen : s->bpos, qe = qs + s->blen;

		_putsn(b, t->name, t->l_name); _sp(b);			/* qname */
		_putsn(b, r->name, r->l_name); _sp(b);			/* rname */
		_putfi(int32_t, b, merr, 4); _sp(b);			/* identity */
		_putn(b, a->a->score); _sp(b);					/* score */

		/* query-side pos and direction */
		_put(b, '0'); _sp(b);
		_putn(b, qs); _sp(b);
		_putn(b, qe); _sp(b);
		_putn(b, t->l_seq); _sp(b);

		/* reference-side pos and direction */
		_putdn(b, s->bid); _sp(b);
		_putn(b, rs); _sp(b);
		_putn(b, re); _sp(b);
		_putn(b, r->l_seq); _cr(b);
	}
	return;
}

/**
 * @fn mm_print_falcon_mapped
 * @brief falcon_sense input format
 */
static
void mm_print_falcon_mapped(
	mm_align_t *b,
	bseq_t const *t,
	mm_reg_t const *reg)
{
	/* print header line (name seq) */
	_putsn(b, t->name, t->l_name); _sp(b);
	_putsnt(b, t->seq, t->l_seq, "NACMGRSVTWYHKDBN"); _cr(b);
	
	/* print alignment lines */
	mm_idx_t const *mi = b->mi;
	uint64_t const n = (b->opt->flag & MM_OMIT_REP)? reg->n_uniq : reg->n_all;
	for(uint64_t i = 0; i < n; i++) {
		mm_aln_t const *a = reg->aln[i];
		const gaba_path_section_t *s = &a->a->sec[0];
		mm_idx_seq_t const *r = &mi->s.a[(s->aid>>1) - mi->base_rid];

		// uint32_t qs = r->l_seq-s->apos-s->alen, qe = r->l_seq-s->apos;
		uint32_t qs = s->apos, qe = s->apos + s->alen;

		/* name seq */
		_putsn(b, r->name, r->l_name); _sp(b);
		if(s->bid & 0x01) {
			_putsntr(b, &r->seq[r->l_seq-qe], qe-qs, "NTGKCYSBAWRDMHVN");
		} else {
			_putsnt(b, &r->seq[qs], qe-qs, "NACMGRSVTWYHKDBN");
		}
		_cr(b);
	}
	_put(b, '+'); _sp(b); _put(b, '+'); _cr(b);
	return;
}

/**
 * @fn mm_print_falcon_footer
 */
static
void mm_print_falcon_footer(
	mm_align_t *b)
{
	_put(b, '-'); _sp(b); _put(b, '-'); _cr(b);
	return;
}

#undef _d
#undef _t
#undef _c
#undef _cr
#undef _sp
#undef _aln

/**
 * @fn mm_align_source
 * @brief source of the alignment pipeline
 */
static
void *mm_align_source(uint32_t tid, void *arg)
{
	set_info(tid, "[mm_align_source] fetch sequence block");
	mm_align_t *b = (mm_align_t *)arg;
	mm_align_step_t *s = (mm_align_step_t*)calloc(1, sizeof(mm_align_step_t));

	s->lmm = lmm_init(NULL, 512 * 1024);
	s->seq = bseq_read(b->fp, &s->n_seq, &s->base, &s->size);
	if(s->seq == NULL) {
		/* reached the tail (or an error occurred) */
		lmm_clean(s->lmm);
		free(s);
		return(NULL);
	}

	/* assign id */
	s->id = b->icnt++;

	/* update and assign base_qid */
	if(b->base_qid == UINT32_MAX) {
		b->base_qid = atoi(s->seq[0].name);
	}
	s->base_qid = b->base_qid;
	b->base_qid += s->n_seq;
	return(s);
}

/**
 * @fn mm_align_worker
 */
static
void *mm_align_worker(uint32_t tid, void *arg, void *item)
{
	mm_tbuf_t *t = (mm_tbuf_t*)arg;
	mm_align_step_t *s = (mm_align_step_t *)item;

	char buf[128], *p = buf;
	p += _pstr(p, "[mm_align_worker] bin id "); p += _pnum(uint32_t, p, s->base_qid); p += _pstr(p, ":"); p += _pnum(uint32_t, p, s->base_qid + s->n_seq - 1); *p = '\0';
	set_info(tid, buf);

	bseq_t *seq = (bseq_t *)s->seq;
	for(uint64_t i = 0; i < s->n_seq; i++) {
		uint32_t qid = s->base_qid + i;			/* FIXME: parse qid from name with atoi when -M is set */
		seq[i].reg = (void *)mm_align_seq(t, seq[i].l_seq, seq[i].seq, qid, s->lmm);
	}
	return(s);
}

/**
 * @fn mm_align_drainm_intl
 */
static _force_inline
void mm_align_drain_intl(mm_align_t *b, mm_align_step_t *s)
{
	for(uint64_t i = 0; i < s->n_seq; i++) {
		mm_reg_t *r = (mm_reg_t *)s->seq[i].reg;

		/* unmapped */
		if(r == NULL) {
			if(b->printer.unmapped) {
				b->printer.unmapped(b, &s->seq[i]);
			}
			continue;
		}

		/* mapped */
		b->printer.mapped(b, &s->seq[i], r);

		/* FIXME: gaba_alignment_t objects must be freed by gaba_dp_res_free */
		lmm_free(s->lmm, r);					/* actually do nothing */
	}
	free(s->base);
	free((void *)s->seq);
	lmm_clean(s->lmm);
	free(s);
	return;
}

/**
 * @fn mm_align_drain
 * @brief alignment pipeline drain, calls formatter
 */
static
void mm_align_drain(uint32_t tid, void *arg, void *item)
{
	set_info(tid, "[mm_align_drain] dump alignments to sam records");
	mm_align_t *b = (mm_align_t *)arg;
	mm_align_step_t *s = (mm_align_step_t*)item;

	#if STRICT_STREAM_ORDERING != 0
		kv_hq_push(v4u32_t, incq_comp, b->hq, ((v4u32_t){.u64 = {s->id, (uintptr_t)s}}));
		while(b->hq.n > 1 && b->hq.a[1].u64[0] == b->ocnt) {
			b->ocnt++;
			s = (mm_align_step_t*)kv_hq_pop(v4u32_t, incq_comp, b->hq).u64[1];
			mm_align_drain_intl(b, s);
		}
	#else
		mm_align_drain_intl(b, s);
	#endif
	return;
}

/**
 * @fn mm_tbuf_destroy
 * @brief destroy thread-local buffer
 */
static _force_inline
void mm_tbuf_destroy(mm_tbuf_t *t)
{
	if(t == NULL) { return; }
	free(t->resc.a);
	free(t->seed.a);
	free(t->map.a);
	free(t->chain.a);
	free(t->bin.a);
	kh_destroy_static(&t->pos);
	gaba_dp_clean(t->dp);
	free(t);
}

/**
 * @fn mm_tbuf_init
 * @brief create thread-local working buffer
 */
static _force_inline
mm_tbuf_t *mm_tbuf_init(mm_align_t *b)
{
	mm_tbuf_t *t = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	if(t == NULL) {
		return(NULL);
	}

	mm_mapopt_t const *opt = b->opt;
	uint32_t const org = (opt->flag & (MM_AVA | MM_COMP)) == MM_AVA ? 0x40000000 : 0;
	uint32_t const thresh = (opt->flag & MM_AVA? 1 : 0) + org;
	uint8_t const *lim = (uint8_t const *)0x800000000000;

	/* init dp context and copy constants */
	*t = (mm_tbuf_t){
		/* index */
		.mi = *b->mi,

		/* ava filter constants */
		.org = org,
		.thresh = thresh,

		/* copied constants */
		.flag = opt->flag,
		.m = opt->m,
		.x = opt->x,
		.wlen = opt->wlen,
		.glen = opt->glen,
		.min_ratio = opt->min_ratio,
		.min_score = opt->min_score,
		.max_cnt = opt->max_cnt == 0 ? 0xffffffff : opt->max_cnt,
		.n_occ = b->n_occ,

		/* tail section */
		.tail = { 0 },
		.t[0] = _sec_fw(0xfffffffe, t->tail, 32),

		/* dp context */
		.dp = gaba_dp_init(b->gaba, lim, lim)
	};
	if(t->dp == NULL) { goto _fail; }

	/* copy occ */
	memcpy(t->occ, b->occ, t->n_occ * sizeof(uint32_t));
	memset(t->occ + t->n_occ, 0, (MAX_FRQ_CNT + 1 - t->n_occ) * sizeof(uint32_t));

	/* init hash */
	kh_init_static(&t->pos, 0);
	return(t);
_fail:
	mm_tbuf_destroy(t);
	return(NULL);
}

/**
 * @fn mm_align_destroy
 * @brief destroy alignment pipeline context
 */
static _force_inline
void mm_align_destroy(mm_align_t *b)
{
	/* flush output buffer */
	if(b->printer.footer) {
		b->printer.footer(b);
	}
	fwrite(b->base, sizeof(uint8_t), b->p - b->base, stdout);

	/* destroy threads */
	for(uint64_t i = 0; i < b->opt->nth; i++) {
		mm_tbuf_destroy((mm_tbuf_t*)b->t[i]);
	}

	/* destroy contexts */
	free(b->base);
	free(b->occ);
	free(b->t);
	kv_hq_destroy(b->hq);
	gaba_clean(b->gaba);
	pt_destroy(b->pt);
	free(b);
	return;
}

/**
 * @fn mm_align_init
 * @brief create alignment pipeline context
 */
static _force_inline
mm_align_t *mm_align_init(mm_mapopt_t const *opt, mm_idx_t const *mi)
{
	static const mm_printer_t printer[] = {
		[0] = { .header = mm_print_sam_header, .unmapped = mm_print_sam_unmapped, .mapped = mm_print_sam_mapped },
		[MM_MAF>>56] = { .mapped = mm_print_maf_mapped },
		[MM_BLAST6>>56] = { .mapped = mm_print_blast6_mapped },
		[MM_BLASR1>>56] = { .mapped = mm_print_blasr1_mapped },
		[MM_BLASR4>>56] = { .mapped = mm_print_blasr4_mapped },
		[MM_PAF>>56] = { .mapped = mm_print_paf_mapped },
		[MM_MHAP>>56] = { .mapped = mm_print_mhap_mapped },
		[MM_FALCON>>56] = { .mapped = mm_print_falcon_mapped, .footer = mm_print_falcon_footer }
	};

	/* malloc context and output buffer */
	void *p = malloc(sizeof(uint8_t) * opt->outbuf_size);
	mm_align_t *b = calloc(1, sizeof(mm_align_t));
	*b = (mm_align_t){
		/* output buffer */
		.base = p,
		.tail = p + opt->outbuf_size,
		.p = p,

		/* refs */
		.mi = mi,
		.opt = opt,

		/* occurrence table */
		.n_occ = opt->n_frq,
		.occ = mm_idx_cal_max_occ(mi, opt->n_frq, opt->frq),

		/* pipeline contexts */
		.fp = NULL,
		.icnt = 0, .ocnt = 0,
		.base_qid = opt->base_qid,
		.gaba = gaba_init(&((struct gaba_params_s const){
			.m = opt->m,
			.x = opt->x,
			.gi = opt->gi,
			.ge = opt->ge,
			.xdrop = opt->xdrop,
			// .filter_thresh = 7,		/* note: can be disabled? */
			.head_margin = sizeof(mm_aln_t)
		})),
		.hq = { .n = 1, .m = 1, .a = NULL },
		.printer = printer[opt->flag>>56],

		/* threads */
		.t = calloc(opt->nth, sizeof(mm_tbuf_t*)),
		.pt = pt_init(opt->nth)
	};


	/* init output queue, buf and printer */
	if(b->occ == NULL || b->gaba == NULL || b->pt == NULL || b->printer.mapped == NULL) {
		opt->log(opt, 'E', __func__, "internal error. failed to create alignment context. debuginfo: (%p, %p, %p, %p)",
			b->occ, b->gaba, b->pt, b->printer.mapped
		);
		goto _fail;
	}

	/* initialize threads */
	for(uint64_t i = 0; i < b->opt->nth; i++) {
		if((b->t[i] = (void *)mm_tbuf_init(b)) == 0) {
			opt->log(opt, 'E', __func__, "internal error. failed to create worker thread. debuginfo: (%lu)", i);
			goto _fail;
		}
	}

	/* initialize binary->bcd conversion table */
	for(uint64_t i = 0; i < 9; i++) {
		b->conv[i] = (i + 1) % 10;
	}
	for(uint64_t i = 9; i < 40; i++) {
		b->conv[i] = (((i + 1) % 10)<<4) + (i + 1) / 10;
	}

	/* print sam header */
	if(b->printer.header) {
		b->printer.header(b);
	}
	return(b);
_fail:
	mm_align_destroy(b);
	return(NULL);
}

/**
 * @fn mm_align_file
 * @brief multithreaded alignment high-level interface
 */
static _force_inline
int mm_align_file(mm_align_t *b, bseq_file_t *fp)
{
	if(fp == 0) { return(-1); }
	b->fp = fp;
	return(pt_stream(b->pt,
		mm_align_source, b,
		mm_align_worker, (void **)b->t,
		mm_align_drain, b)
	);
}

/* end of map.c */

/* main.c */

/**
 * @fn liftrlimit
 * @brief elevate max virtual memory size
 */
static _force_inline
void liftrlimit()
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}

/**
 * @fn posixly_correct
 * @brief change option parser behavior to posix-compatible
 */
static _force_inline
void posixly_correct()
{
#ifdef __linux__
	setenv("POSIXLY_CORRECT", "1", 0);
#endif
}

/**
 * @fn mm_print_version
 */
static
void mm_print_version(mm_mapopt_t const *opt)
{
	opt->log(opt, 16, __func__, "%s", version());
	return;
}

/**
 * @fn mm_mapopt_load_preset
 * @brief load preset params
 */
static _force_inline
int mm_mapopt_load_preset(mm_mapopt_t *o, char const *arg)
{
	if(strcmp(arg, "pacbio") == 0) {
		o->k = 15; o->w = 10; o->m = 1; o->x = 2; o->gi = 2; o->ge = 1; o->xdrop = 50; o->min_score = 50; o->min_ratio = 0.3;
	} else if(strcmp(arg, "ont") == 0) {
		o->k = 15; o->w = 10; o->m = 1; o->x = 1; o->gi = 1; o->ge = 1; o->xdrop = 50; o->min_score = 50; o->min_ratio = 0.3;
	} else if(strcmp(arg, "ava") == 0) {
		o->k = 14; o->w = 5; o->m = 1; o->x = 2; o->gi = 0; o->ge = 1; o->xdrop = 50; o->min_score = 30; o->min_ratio = 0.05;
		o->flag |= MM_AVA | MM_PAF;
	} else {
		o->log(o, 'W', __func__, "Warning: Unknown preset tag: `%s'.", arg);
		return 1;
	}
	return 0;
}

/**
 * @fn mm_mapopt_parse_threshs
 * @brief parse frequence thresholds
 */
static _force_inline
int mm_mapopt_parse_threshs(
	mm_mapopt_t *o,
	char const *arg)
{
	char const *p = optarg; o->n_frq = 0;
	while(*p != '\0') {
		o->frq[o->n_frq++] = atof(p);
		while(*p != '\0' && *p != ',') { p++; }
		if(*p != '\0' || o->n_frq >= 16) { break; }
		p++;
	}
	for(uint64_t i = 0; i < o->n_frq; i++) {
		if(o->frq[i] < 0.0 || o->frq[i] > 1.0) {
			o->log(o, 'W', __func__, "Warning: Invalid threshold `%f' parsed from `%s'.", o->frq[i], arg);
		}
	}
	return 0;
}

/**
 * @fn mm_mapopt_parse_base_ids
 * @brife parse base id, returns nonzero if failed
 */
static _force_inline
int mm_mapopt_parse_base_ids(
	mm_mapopt_t *o,
	char const *arg)
{
	char const *p = optarg;
	if(!isdigit(*p)) { return(1); }

	/* reference-side id */
	o->base_rid = atoi(p);
	while(*p != '\0' && *p != ',') { p++; }

	/* query-side id */
	if(*p != '\0' && isdigit(p[1])) {
		o->base_qid = atoi(&p[1]);
	}
	return(0);
}

/**
 * @macro _encode
 * @brief map tag string to 16-bit int
 */
#define _encode(p)		( (uint16_t)p[0] | ((uint16_t)p[1]<<8) )

/**
 * @fn mm_mapopt_tag2flag
 */
static _force_inline
uint64_t mm_mapopt_tag2flag(
	mm_mapopt_t *o,
	uint16_t t)
{
	#define _test(_str) ( _encode(_str) == t )

	if(_test("RG")) { return(0x01ULL<<MM_RG); }
	if(_test("CO")) { return(0x01ULL<<MM_CO); }
	if(_test("NH")) { return(0x01ULL<<MM_NH); }
	if(_test("IH")) { return(0x01ULL<<MM_IH); }
	if(_test("AS")) { return(0x01ULL<<MM_AS); }
	if(_test("XS")) { return(0x01ULL<<MM_XS); }
	if(_test("NM")) { return(0x01ULL<<MM_NM); }
	if(_test("SA")) { return(0x01ULL<<MM_SA); }
	if(_test("MD")) { return(0x01ULL<<MM_MD); }

	o->log(o, 'W', __func__, "Unknown tag: `%c%c'.", t & 0xff, t>>8);
	return(0);

	#undef _test
}

/**
 * @fn mm_mapopt_parse_tags
 * @brief parse sam tag identifier and convert to flags
 */
static _force_inline
uint64_t mm_mapopt_parse_tags(
	mm_mapopt_t *o,
	char const *p,
	uint16_v *buf)
{
	uint32_t flag = 0;
	while(*p != '\0') {
		char const *q = p;

		/* split first tag */
		while(*q != '\0' && *q != ',') { q++; }

		/* skip if the length is not 2 */
		if(q - p != 2) {
			o->log(o, 'W', __func__,
				"Unparsable tag: `%.*s', must be two-character string like `RG'.",
				q - p, p
			);
			p = q + 1;
			continue;
		}

		/* save if buffer is available */
		uint16_t t = _encode(p);
		if(buf != NULL) { kv_push(uint16_t, *buf, t); }

		/* set tag */
		uint64_t f = mm_mapopt_tag2flag(o, t);
		if(f == 0 && o->verbose >= 1) {
			o->log(o, 'W', __func__, "Unknown tag: `%.*s'.", 2, p);
		}
		flag |= f;

		/* test tail and advance pointer */
		if(*q == '\0') { break; }
		p = q + 1;
	}
	return(flag);
}

#undef _encode

/**
 * @fn mm_mapopt_parse_rg
 * @brief parse read group string
 */
static _force_inline
void mm_mapopt_parse_rg(
	mm_mapopt_t *o,
	char const *arg)
{
	/* cleanup old ones */
	free(o->rg_line);	o->rg_line = NULL;
	free(o->rg_id);		o->rg_id = NULL;
	o->flag &= ~(0x01ULL<<MM_RG);

	/* first unescape tabs */
	char const *id, *p = arg;
	char b[256], *q = b;
	while(q < b+256 && *p != '\0') {
		if(p[0] == '\\' && p[1] == 't') {
			if(p[1] == '\0') { break; }
			*q++ = '\t';
			p += 2;
		} else {
			*q++ = *p++;
		}
	}
	*q = '\0';

	/* search header and id */
	if(strstr(b, "@RG") != b || (id = strstr(b, "\tID:")) == 0) {
		o->log(o, 'E', __func__, "RG line must start with @RG and contains ID, like `@RG\\tID:1'.");
		return;
	}

	/* copy */
	#define _dup(x)	({ char *_p = malloc(strlen(x)+1); memcpy(_p, (x), strlen(x)); _p[strlen(x)] = '\0'; _p; })
	o->rg_line = _dup(b);
	o->rg_id = _dup(id+4);
	#undef _dup

	/* truncate tail of id line */
	char *t = o->rg_id;
	while(*t != '\t' && *t != '\n' && *t != '\0') { t++; }
	*t = '\0';

	/* enable rg tag */
	o->flag |= 0x01ULL<<MM_RG;
	return;
}

/**
 * @fn mm_mapopt_parse_format
 * @brief convert format string to flag
 */
static _force_inline
uint64_t mm_mapopt_parse_format(
	mm_mapopt_t *o,
	char const *arg)
{
	if(strcmp(arg, "sam") == 0) { return(0); }
	if(strcmp(arg, "maf") == 0) { return(MM_MAF); }
	if(strcmp(arg, "blast6") == 0) { return(MM_BLAST6); }
	if(strcmp(arg, "blasr1") == 0) { return(MM_BLASR1); }
	if(strcmp(arg, "blasr4") == 0) { return(MM_BLASR4); }
	if(strcmp(arg, "paf") == 0) { return(MM_PAF); }
	if(strcmp(arg, "mhap") == 0) { return(MM_MHAP); }
	if(strcmp(arg, "falcon") == 0) { return(MM_FALCON); }
	o->log(o, 'W', __func__, "Unknown output format: `%s'.", arg);
	return 0;
}

/**
 * @fn mm_mapopt_save_args
 */
static _force_inline
void mm_mapopt_save_args(
	mm_mapopt_t *o,
	int argc, char *const *argv)
{
	kvec_t(char) v = { .n = 0, .m = 0, .a = NULL };
	for(uint64_t i = 0; i < (uint64_t)argc; i++) {
		kv_pushm(char, v, argv[i], strlen(argv[i]));
		kv_push(char, v, ' ');
	}
	v.a[v.n - 1] = '\0';
	o->arg_line = v.a;
	return;
}

/**
 * @fn mm_mapopt_parse
 * @brief parse option string and set to mapopt_t object
 */
static _force_inline
int mm_mapopt_parse(
	mm_mapopt_t *o,
	int argc, char *argv[],
	char const **fnr,				/* input index filename */
	char const **fnw,				/* output index filename */
	ptr_v *v)						/* positional argumants */
{
	int ret = 0;
	mm_mapopt_save_args(o, argc, argv);
	while(optind < argc) {
		int ch;

		if((ch = getopt(argc, argv, "t:x:V:c:k:w:f:B:d:l:C:NXAs:m:r:M:a:b:p:q:W:G:L:H:D:Y:O:PQR:T:U:vh")) < 0) {
			/* positional argument */
			kv_push(void *, *v, argv[optind]);
			optind++;
			continue;
		}

		/* option */
		switch (ch) {
			case 't': o->nth = atoi(optarg); break;
			case 'x': mm_mapopt_load_preset(o, optarg); break;
			case 'V':
				if(isdigit(optarg[0])) {
					o->verbose = atoi(optarg);
				} else {
					for(char const *p = optarg; *p == 'V'; p++, o->verbose++) {}
				}
				break;
			case 'k': o->k = atoi(optarg); break;
			case 'w': o->w = atoi(optarg); break;
			case 'f': mm_mapopt_parse_threshs(o, optarg); break;
			case 'B': o->b = atoi(optarg); break;
			case 'd': *fnw = optarg; break;
			case 'l': *fnr = optarg; break;
			case 'C': mm_mapopt_parse_base_ids(o, optarg); break;
			case 'N': o->base_rid = UINT32_MAX; o->base_qid = UINT32_MAX; break;
			case 'X': o->flag |= MM_AVA; break;
			case 'A': o->flag |= MM_COMP; break;
			case 's': o->min_score = atoi(optarg); break;
			case 'm': o->min_ratio = atof(optarg); break;
			case 'M': o->max_cnt = atoi(optarg); break;
			case 'a': o->m = atoi(optarg); break;
			case 'b': o->x = atoi(optarg); break;
			case 'p': o->gi = atoi(optarg); break;
			case 'q': o->ge = atoi(optarg); break;
			case 'W': o->wlen = atoi(optarg); break;
			case 'G': o->glen = atoi(optarg); break;
			case 'L': o->rmin = atoi(optarg); break;
			case 'H': o->qmin = atoi(optarg); break;
			case 'D': o->batch_size = atoi(optarg); break;
			case 'Y': o->xdrop = atoi(optarg); break;
			case 'O': o->flag &= ~(0xffULL<<56); o->flag |= mm_mapopt_parse_format(o, optarg); break;
			case 'P': o->flag |= MM_OMIT_REP; break;
			case 'Q': o->flag |= MM_KEEP_QUAL; break;
			case 'R': mm_mapopt_parse_rg(o, optarg); break;
			case 'T': o->flag |= mm_mapopt_parse_tags(o, optarg, NULL); break;
			case 'U': mm_mapopt_parse_tags(o, optarg, &o->tags); break;
			case 'v': ret = 1; o->verbose = 0; break;
			case 'h': ret = 1; break;
		}
	}
	return(ret);
}


/**
 * @fn mm_print_help
 */
static
void mm_print_help(mm_mapopt_t const *opt)
{
	if(opt->verbose == 0) { return; }

	#define _msg(_level, ...) { \
		opt->log(opt, 16 + _level, __func__, __VA_ARGS__); \
	}

	_msg(1, "\n"
			"  minialign - fast aligner for PacBio and Nanopore long reads\n"
			"\n"
			"minialign is a fast long-read (nucleotide sequence) alignment tool built on\n"
			"the top of minimap long-read overlapper adopting libgaba SIMD-parallelized\n"
			"Smith-Waterman extension algorithm.\n"
			"");
	_msg(1, "Usage:"
			"  first trial:\n"
			"    $ minialign -t4 -xont <ref.fa> <ont2d.{fa,fq,bam}> > mapping.sam\n"
			"\n"
			"  mapping on a prebuilt index (saves ~1min for human genome per run):\n"
			"    $ minialign [indexing options] -d <index.mai> <ref.fa>\n"
			"    $ minialign -l <index.mai> <reads.{fa,fq,bam}> > mapping.sam\n"
			"\n"
			"  all-versus-all alignment in a read set:\n"
			"    $ minialign -X -xava <reads.fa> [<reads.fa> ...] > allvsall.paf\n"
			"");
	_msg(1, "Options:");
	_msg(1, "  Global:");
	_msg(1, "    -x STR       load preset params {pacbio,ont,ava} [ont]");
	_msg(1, "    -t INT       number of threads [%d]", opt->nth);
	_msg(1, "    -X           switch to all-versus-all alignment mode");
	_msg(1, "    -v           show version number [%s]", version());
	_msg(1, "  Indexing:");
	_msg(7, "    -c STR,...   treat specified sequences as circular []");
	_msg(1, "    -k INT       k-mer size [%d]", opt->k);
	_msg(1, "    -w INT       minimizer window size [{-k}*2/3]");
	_msg(1, "    -d FILE      dump index to FILE []");
	_msg(1, "    -l FILE      load index from FILE [] (overriding -k and -w)");
	_msg(2, "    -C INT[,INT] set base rid and qid [%u, %u]", opt->base_rid, opt->base_qid);
	_msg(2, "    -N           treat sequence name as id (seq must be sorted)");
	_msg(2, "    -L INT       sequence length filter: min. ref. length; 0 to disable [%u]", opt->rmin);
	_msg(2, "    -H INT       sequence length filter: min. query length; 0 to disable [%u]", opt->qmin);
	_msg(1, "  Mapping:");
	_msg(2, "    -f FLOAT,... occurrence thresholds [0.05,0.01,0.001]");
	_msg(1, "    -a INT       match award [%d]", opt->m);
	_msg(1, "    -b INT       mismatch penalty [%d]", opt->x);
	_msg(1, "    -p INT       gap open penalty [%d]", opt->gi);
	_msg(1, "    -q INT       gap extension penalty [%d]", opt->ge);
	_msg(2, "    -Y INT       X-drop threshold [%d]", opt->xdrop);
	_msg(1, "    -s INT       minimum alignment score [%d]", opt->min_score);
	_msg(1, "    -m INT       minimum alignment score ratio [%1.2f]", opt->min_ratio);
	_msg(2, "    -M INT       maximum #alignments reported [%d]", opt->max_cnt);
	_msg(2, "    -A           calculate both topright and bottomright triangles (only effective with -X)");
	_msg(1, "  Output:");
	_msg(1, "    -O STR       output format {sam,maf,blast6,blasr1,blasr4,paf,mhap,falcon} [%s]",
		(char const *[]){ "sam", "maf", "blast6", "blasr1", "blasr4", "paf", "mhap", "falcon" }[opt->flag>>56]);
	_msg(2, "    -P           omit secondary (repetitive) alignments");
	_msg(1, "    -Q           include quality string");
	_msg(1, "    -R STR       read group header line, like \"@RG\\tID:1\" [%s]", opt->rg_line? opt->rg_line : "");
	_msg(1, "    -T STR,...   list of optional tags: {RG,AS,XS,NM,NH,IH,SA,MD} []");
	_msg(1, "                   RG is also inferred from -R");
	_msg(1, "                   supp. records are omitted when SA tag is enabled");
	_msg(2, "    -U STR,...   tags to be transferred from the input bam file []");
	_msg(1, "");
	if(opt->verbose < 2) {
		_msg(0, "  Pass -hVV to show all the options.");
		_msg(0, "");
	}

	#undef _msg
	return;
}

/**
 * @fn main
 */
int _export(main)(int argc, char *argv[])
{
	int ret = 1;
	char const *fnr = NULL, *fnw = NULL;
	FILE *fpr = NULL, *fpw = NULL;
	ptr_v query = { 0 };

	/* unittest hook, see unittest.h for the details */
	#if UNITTEST != 0
	if(argc > 1 && strcmp(argv[1], "unittest") == 0) {
		return(unittest_main(argc, argv));
	}
	#endif

	/* enable debug info */
	enable_info(0);
	set_info(0, "[main] parsing arguments");

	/* elevate memory limit */
	liftrlimit();

	/* set option parser behavior to posix-compatible */
	posixly_correct();

	/* init option object (init base time) and parse args */
	mm_mapopt_t *opt = mm_mapopt_init();
	if(mm_mapopt_parse(opt, argc, argv, &fnr, &fnw, &query)) {
		/* set logger to stdout when invoked by option */
		opt->fp = stdout;
		if(opt->verbose == 0) {
			mm_print_version(opt);
		} else {
			mm_print_help(opt);
		}

		/* also exit status is 0 */
		ret = 0;
		goto _final;
	}

	/* if no input is given, print help to stderr and exit with 1 */
	if(!fnr && query.n == 0) {
		mm_print_help(opt);
		goto _final;
	}

	/* if not in index construction mode and query-side file is missing, pour stdin to query */
	if(!fnw && ((fnr && query.n == 0) || (!fnr && query.n == 1 && !(opt->flag&MM_AVA)))) {
		opt->log(opt, 1, __func__, "query-side input redirected to stdin.");
		kv_push(void *, query, "-");
	}

	/* calculate default window size if not specified */
	if(opt->w >= 16) {
		opt->w = (int)(.6666667 * opt->k + .499);
	}

	/* check options */
	if(mm_mapopt_check(opt)) {
		goto _final;
	}

	/* open input index file */
	set_info(0, "[main] open index file");
	if(fnr && (fpr = fopen(fnr, "rb")) == NULL) {
		opt->log(opt, 'E', __func__,
			"failed to open index file `%s'. Please check file path and it exists.",
			fnr
		);
		goto _final;
	}

	/* fnw != NULL switches to index construction mode */
	if(fnw && (fpw = fopen(fnw, "wb")) == NULL) {
		opt->log(opt, 'E', __func__,
			"failed to open index file `%s' in write mode. Please check file path and its permission.",
			fnw
		);
		goto _final;
	}

	/* iterate over index *blocks* */
	for(uint64_t i = 0;
		i < (fpr? UINT64_MAX : (((opt->flag&MM_AVA) || fpw)? query.n : 1));
		i++
	) {
		bseq_file_t *fp = NULL;
		mm_idx_t *mi = NULL;

		/* load/generate index for this index-side iteration */
		if(fpr != NULL) {
			/* fetch the next block */
			mi = mm_idx_load(fpr, opt->nth);
		} else {
			/* fetch the next file */
			fp = bseq_open((char const *)query.a[i], opt->batch_size, 0, opt->rmin, 0, NULL);
			mi = mm_idx_gen(opt, fp);
			opt->base_rid += bseq_close(fp);
		}

		/* check sanity of the index */
		if(mi == NULL) {
			/* when index is loaded from file, mi == NULL indicates the correct tail of the file */
			if(fpr && i > 0) { break; }

			/* otherwise error */
			opt->log(opt, 'E', __func__, "failed to %s `%s'. Please check %s.",
				fpr? "load index file" : "open sequence file",
				fpr? fnr : (const char*)query.a[i],
				fpr? "file path, format and its version" : "file path and format"
			);
			mm_idx_destroy(mi);
			goto _final;
		}
		opt->log(opt, 9, __func__, "loaded/built index for %lu target sequence(s).", mi->s.n);

		/* do the task */
		if(fpw != NULL) {
			mm_idx_dump(fpw, mi, opt->nth);
		} else {
			/* init alignment; create threads and memory arena */
			mm_align_t *aln = mm_align_init(opt, mi);

			/* iterate over all the queries */
			for(uint64_t j = (!fpr && !(opt->flag&MM_AVA)); j < query.n; j++) {
				opt->log(opt, 10, __func__, "fetched `%s'.", (const char *)query.a[j]);

				/* fetch the next query */
				if(!(fp = bseq_open((const char *)query.a[j],
					opt->batch_size,
					(opt->flag&MM_KEEP_QUAL) != 0,
					opt->qmin,
					opt->tags.n, opt->tags.a))
				) {
					/* an error occured in loading query */
					opt->log(opt, 'E', __func__,
						"failed to open sequence file `%s'."
						"Please check file path and format.",
						(const char*)query.a[j]
					);
					mm_align_destroy(aln);
					mm_idx_destroy(mi);
					goto _final;
				}
				mm_align_file(aln, fp);
				bseq_close(fp);

				opt->log(opt, 10, __func__, "finished `%s'.", (const char *)query.a[j]);
			}

			/* finish alignment */
			mm_align_destroy(aln);
		}
		/* cleanup index for the current iteration */
		mm_idx_destroy(mi);
	}

	opt->log(opt, 1, __func__, "Version: %s", MM_VERSION);
	opt->log(opt, 1, __func__, "CMD: %s", opt->arg_line);
	opt->log(opt, 1, __func__, "Real time: %.3f sec; CPU: %.3f sec",
		realtime() - opt->inittime,
		cputime()
	);
	ret = 0;
_final:
	free(query.a);
	mm_mapopt_destroy(opt);
	if(fpr != NULL) { fclose(fpr); }
	if(fpw != NULL) { fclose(fpw); }
	return(ret);
}

/* end of main.c */
