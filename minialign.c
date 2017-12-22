
/**
 * @file minialign.c
 *
 * @brief minialign core implementation
 *
 * @author Hajime Suzuki (original files by Heng Li)
 * @license MIT
 */
#ifndef NDEBUG
#  define DEBUG
#endif
/* configurations */
/**
 * @macro MM_VERSION
 */
#ifndef MM_VERSION
#  define MM_VERSION		"minialign-0.6.0-devel"
#endif

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
#  define UNITTEST 					( 1 )
#endif

/* make sure POSIX APIs are properly activated */
#if defined(__linux__) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE		200112L
#endif
#if defined(__darwin__) && !defined(_DARWIN_C_FULL)
#  define _DARWIN_C_SOURCE		_DARWIN_C_FULL
#endif

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

/* utils.h */
/* max, min, roundup */
#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )
#define _roundup(x, base)			( ((x) + (base) - 1) & ~((base) - 1) )

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
/* utils.h */

#include "sassert.h"
#include "log.h"

/* mm_malloc.c: malloc wrappers */
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
	fprintf(stderr, "[E::%s] Out of memory. (required: %zu MB, maxrss: %ld MB)\n", name, req / 1024, r.ru_maxrss);
	trap();								/* segv trap for debugging; see log.h */
	exit(128);							/* 128 reserved for out of memory */
}

/**
 * @macro mm_malloc
 */
#define mm_malloc(_x) ({ \
	void *_ptr = malloc((size_t)(_x)); \
	if(_unlikely(_ptr == NULL)) { \
		oom_abort(__func__, (_x)); \
	} \
	_ptr; \
})
#define mm_new(type_t, _content) ({ \
	type_t *_ptr = mm_malloc(sizeof(type_t)); \
	*_ptr = (type_t)_content; \
	_ptr; \
})
#define malloc(_x)		mm_malloc(_x)

/**
 * @macro mm_realloc
 */
#define mm_realloc(_x, _y) ({ \
	void *_ptr = realloc((void *)(_x), (size_t)(_y)); \
	if(_unlikely(_ptr == NULL)) { \
		oom_abort(__func__, (_y)); \
	} \
	_ptr; \
})
#define realloc(_x, _y)	mm_realloc(_x, _y)

/**
 * @macro mm_calloc
 */
#define mm_calloc(_x, _y) ({ \
	void *_ptr = calloc((size_t)(_x), (size_t)(_y)); \
	if(_unlikely(_ptr == NULL)) { \
		oom_abort(__func__, (_y)); \
	} \
	_ptr; \
})
#define calloc(_x, _y)	mm_calloc(_x, _y)

/* end of mm_malloc.c */

/* miscellaneous macros and types */
#define UNITTEST_UNIQUE_ID		1
#include "unittest.h"
#include "lmm.h"
#include "gaba_wrap.h"
#include "arch/arch.h"
#include "kvec.h"

unittest_config( .name = "minialign" );
unittest( .name = "mm_malloc" ) {
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

/**
 * sequence encoding table
 * "NTGKCYSBAWRDMHVN"
 * "NACMGRSVTWYHKDBN"
 */
enum alphabet {
	A = 0x00,
	C = 0x01,
	G = 0x02,
	T = 0x03,
	N = 0x80
};
static uint8_t const enc4f[16] = { [1] = A, [2] = C, [4] = G, [8] = T };
static uint8_t const enc4r[16] = { [1] = T, [2] = G, [4] = C, [8] = A };
static uint8_t const encaf[16] = {
	['A' & 0x0f] = A, ['C' & 0x0f] = C,
	['G' & 0x0f] = G, ['T' & 0x0f] = T,
	['U' & 0x0f] = T, ['N' & 0x0f] = N
};
// #define _encaf(_c)			 ( 0x03 & (((_c)>>1) - ((_c)>>3)) )
#define _encaf(_c)			 ( 0x03 & (((_c)>>2) ^ ((_c)>>1)) )
static uint8_t const decaf[16] = { [A] = 'A', [C] = 'C', [G] = 'G', [T] = 'T' };
static uint8_t const decar[16] = { [A] = 'T', [C] = 'G', [G] = 'C', [T] = 'A' };
static uint8_t const idxaf[256] = { ['A'] = 1, ['C'] = 2, ['G'] = 3, ['T'] = 4, ['U'] = 4 };

/**
 * @macro mm_encode_tag
 * @brief map SAM tag string to 16-bit int
 */
#define mm_encode_tag(_p)		( (uint16_t)(_p)[0] | ((uint16_t)(_p)[1]<<8) )

/**
 * @type read_t, write_t
 * @brief abstract I/O function types
 */
typedef uint64_t (*read_t)(void *fp, void *buf, uint64_t size);
typedef uint64_t (*write_t)(void *fp, void const *buf, uint64_t size);

/* time functions */
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
static _force_inline
double realtime(void)
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return(tp.tv_sec + tp.tv_usec * 1e-6);
}

/* random */
static _force_inline
uint64_t mm_rand64(void)
{
	uint64_t bits = 31, n = 0;
	for(uint64_t acc = 0; acc < 64; acc += bits) {
		n <<= bits; n ^= (uint64_t)rand();
	}
	return n;
}

/* string handling */
static _force_inline
char *mm_strndup(char const *p, uint64_t l)
{
	if(!p) { return(NULL); }
	char *q = malloc(sizeof(char) * (l + 1));
	memcpy(q, p, l); q[l] = '\0';
	return(q);
}
static _force_inline
char *mm_strdup(char const *p)
{
	if(!p) { return(NULL); }
	kvec_t(char) b = { 0 }; p--;
	while(*++p) { kv_push(char, b, *p); }
	return(b.a);
}
static _force_inline
char *mm_join(char const *const *p, char c)
{
	kvec_t(char) b = { 0 }; p--;
	while(*++p) {
		char const *q = *p - 1;
		while(*++q) { kv_push(char, b, *q); }
		kv_push(char, b, c);
	}
	b.a[--b.n] = '\0';
	return(b.a);
}
static _force_inline
int mm_startswith(char const *p, char const *prf)
{
	uint64_t l = strlen(p), r = strlen(prf);
	return(l >= r && strncmp(p, prf, r) == 0);
}
static _force_inline
int mm_endswith(char const *p, char const *suf)
{
	uint64_t l = strlen(p), r = strlen(suf);
	return(l >= r && strncmp(p + l - r, suf, r) == 0);
}
static _force_inline
char *mm_append(char *p, char const *suf)
{
	uint64_t l = strlen(p), r = strlen(suf);
	p = realloc(p, sizeof(char) * (l + r + 1));
	memcpy(p + l, suf, r); p[l + r] = '\0';
	return(p);
}
static _force_inline
uint64_t mm_shashn(char const *p, uint64_t l)
{
	uint64_t a = 0x12345678;
	while(*p && l) { a = (a<<5) ^ (a>>3) ^ *p++; l--; }
	return(a);
}
static _force_inline
char const *version(void)
{
	char const *prefix = "minialign-";
	uint64_t spos = mm_startswith(MM_VERSION, prefix) ? 0 : strlen(prefix);	/* match to remove prefix */
	return(&MM_VERSION[spos]);
}
/* end of miscellaneous section */

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

#define kh_ptr(h)		( (h)->a )
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
void kh_dump(kh_t const *h, void *fp, write_t wfp)
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
void kh_load_static(kh_t *h, void *fp, read_t rfp)
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
	if(k0 == k1) { return((kh_bidx_t){ i, 0 }); }	/* duplicated key */

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

	debug("kh_extend called, new_size(%lx), cnt(%u), ub(%u)", size, h->cnt, h->ub);

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
	if(extend != 0 && h->cnt >= h->ub) { kh_extend(h); }
	kh_bidx_t b = kh_allocate(h->a, key, KH_INIT_VAL, h->mask);
	debug("allocated hash bin (%lu, %lu), (%lx, %lx), dst(%ld)", b.idx, b.n, h->a[b.idx].u64[0], h->a[b.idx].u64[1], (b.idx - key) & h->mask);

	h->cnt += b.n;
	// h->ub = ((b.idx - key) & h->mask) > KH_DST_MAX ? 0 : h->ub;
	return(&h->a[b.idx].u64[1]);
}

/**
 * @fn kh_get
 * @brief returns val or UINT64_MAX if not found
 */
static _force_inline
uint64_t kh_get(kh_t const *h, uint64_t key)
{
	uint64_t mask = h->mask, pos = key & mask, k;
	do {
		if((k = h->a[pos].u64[0]) == key) { return(h->a[pos].u64[1]); }
		pos = mask & (pos + 1);
	} while(k + 1 != 0);			/* !is_empty(k) || is_moved(k) */
	return(UINT64_MAX);
}

/**
 * @fn kh_get_ptr
 * @brief returns pointer to val or NULL if not found
 */
static _force_inline
uint64_t *kh_get_ptr(kh_t const *h, uint64_t key)
{
	uint64_t mask = h->mask, pos = key & mask, k;
	do {
		if((k = h->a[pos].u64[0]) == key) { return(&h->a[pos].u64[1]); }
		pos = mask & (pos + 1);
	} while(k + 1 != 0);			/* !is_empty(k) || is_moved(k) */
	return(NULL);
}

/**
 * @struct kh_str_t
 * @brief string -> string hashmap
 */
typedef struct {
	uint8_v s;
	kh_t h;
} kh_str_t;
#define kh_str_ptr(_h)					( kh_ptr(&(_h)->h) )
#define kh_str_size(_h)					( kh_size(&(_h)->h) )
#define kh_str_init_static(_h, _size)	{ (_h)->s = (uint8_v){ 0 }; kh_init_static(&(_h)->h, _size); }
#define kh_str_init(_size)				({ kh_str_t *h = calloc(1, sizeof(kh_str_t)); kh_init_static(h, _size); h; })
#define kh_str_destroy_static(_h)		{ free(((kh_t *)_h)->h.a); free(((kh_t *)_h)->h.s.a); }
#define kh_str_destroy(_h)				{ kh_str_destroy(_h); free(_h); }

/**
 * @fn kh_str_put, kh_str_get
 */
static _force_inline
void kh_str_put(kh_str_t *h, char const *k, uint64_t klen, char const *v, uint64_t vlen)
{
	if(klen == UINT64_MAX) { klen = strlen(k); }
	if(vlen == UINT64_MAX) { vlen = strlen(v); }
	kh_put(&h->h, mm_shashn(k, klen), (kv_size(h->s)<<32) | klen);
	kv_pushm(uint8_t, h->s, k, klen);
	kv_push(uint8_t, h->s, '\0');
	kv_pushm(uint8_t, h->s, v, vlen);
	kv_push(uint8_t, h->s, '\0');
	return;
}
static _force_inline
char const *kh_str_get(kh_str_t const *h, char const *k, uint64_t klen)
{
	if(klen == UINT64_MAX) { klen = strlen(k); }
	uint64_t const *p = kh_get_ptr(&h->h, mm_shashn(k, klen));
	if(p == NULL || (uint32_t)(*p) != klen) { return(0); }		/* not found or length differs */
	char const *q = (char const *)&h->s.a[*p>>32];
	return(strncmp(q, k, klen) == 0 ? q + klen + 1 : NULL);
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
	kh_dump(h, (void *)fp, (write_t)gzwrite);
	gzclose(fp);
	kh_destroy(h);

	/* restore */
	fp = gzopen(filename, "r");
	assert((void *)fp != NULL);
	h = kh_load((void *)fp, (read_t)gzread);
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
#define pt_nth(_pt)	( (_pt)->nth )
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
	if(pt == NULL) { return; }

	/* send termination signal */
	for(uint64_t i = 1; i < pt->nth; i++) {
		pt_enq(pt->c->in, pt->c->tid, PT_EXIT);
	}

	/* wait for the threads terminate */
	for(uint64_t i = 1; i < pt->nth; i++) {
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
	pt_t *pt = calloc(nth, sizeof(pt_t) + sizeof(pt_thread_t));

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
	}
	return(pt);
}

/**
 * @fn pt_set_worker
 * @brief update worker function and argument pointers
 */
static _force_inline
int pt_set_worker(pt_t *pt, void *arg, pt_worker_t wfp)
{
	void *item;

	/* fails when unprocessed object exists in the queue */
	if((item = pt_deq(&pt->in, 0)) != PT_EMPTY) {
		pt_enq(&pt->in, 0, item);
		return(-1);
	}

	/* update pointers */
	for(uint64_t i = 0; i < pt->nth; i++) {
		pt->c[i].wfp = wfp; pt->c[i].warg = arg;
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
int pt_stream(pt_t *pt, void *arg, pt_source_t sfp, pt_worker_t wfp, pt_drain_t dfp)
{
	if(pt_set_worker(pt, arg, wfp)) { return(-1); }

	/* keep balancer between [lb, ub) */
	uint64_t const lb = 2 * pt->nth, ub = 4 * pt->nth;
	uint64_t bal = 0;
	void *it;

	while((it = sfp(0, arg)) != NULL) {
		pt_enq(&pt->in, 0, it);
		if(++bal < ub) { continue; }

		while(bal > lb) {
			/* flush to drain (note: while loop is better?) */
			if((it = pt_deq(&pt->out, 0)) != PT_EMPTY) { bal--; dfp(0, arg, it); }
			/* process one in the master (parent) thread */
			if((it = pt_deq(&pt->in, 0)) != PT_EMPTY) { bal--; dfp(0, arg, wfp(0, arg, it)); }
		}
	}

	/* source depleted, process remainings */
	while((it = pt_deq(&pt->in, 0)) != PT_EMPTY) { pt_enq(&pt->out, 0, wfp(0, arg, it)); }

	/* flush results */
	while(bal > 0) { if((it = pt_deq(&pt->out, 0)) != PT_EMPTY) { bal--; dfp(0, arg, it); } }
	return(0);
}

/**
 * @fn pt_parallel
 */
static _force_inline
int pt_parallel(pt_t *pt, void *arg, pt_worker_t wfp)
{
	/* fails when unprocessed element exists */
	if(pt_set_worker(pt, arg, wfp)) { return(-1); }

	/* push items */
	for(uint64_t i = 1; i < pt->nth; i++) { pt_enq(&pt->in, 0, (void *)i); }
	debug("pushed items");

	/* process the first one in the master (parent) thread */
	wfp(0, arg, (void *)0);
	debug("finished master");

	/* wait for the children done */
	for(uint64_t i = 1; i < pt->nth; i++) {
		while(pt_deq(&pt->out, 0) == PT_EMPTY) {
			struct timespec tv = { .tv_nsec = 512 * 1024 };
			nanosleep(&tv, NULL);	/* wait for a while */
			/* sched_yield(); */
		}
		debug("joined i(%lu)", i);
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
	pt_stream(pt, arr, pt_unittest_source, pt_unittest_worker, pt_unittest_drain);
	assert(icnt == 1024, "icnt(%lu)", icnt);
	assert(ocnt == 512*1025, "ocnt(%lu)", ocnt);

	uint64_t s[4] = { 0 };	//, *sp[1] = { &s[0] };
	pt_parallel(pt, arr, pt_unittest_worker);
	assert(s[0] == 1, "d[0](%lu)", s[0]);
	pt_destroy(pt);
}

unittest( .name = "pt.multi" ) {
	pt_t *pt = pt_init(4);
	assert(pt != NULL);

	uint64_t icnt = 0, ocnt = 0, inc = 1, *arr[4] = { &inc, &inc, &inc, &inc };
	pt_parallel(pt, arr, pt_unittest_worker);
	assert(icnt == 1024, "icnt(%lu)", icnt);
	assert(ocnt == 512*1025, "ocnt(%lu)", ocnt);

	uint64_t s[4] = { 0,1,2,3 };	//, *sp[4] = { &s[0],&s[1],&s[2],&s[3] };
	pt_parallel(pt, arr, pt_unittest_worker);
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
	pg_t *pg = (pg_t *)arg;
	pg_block_t *s = (pg_block_t *)item;
	if(s == NULL || s->len == 0) { return(s); }
	return((s->raw ? pg_deflate : pg_inflate)(s, pg->block_size));
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
pg_t *pg_init(FILE *fp, pt_t *pt)
{
	if(fp == NULL) { return(NULL); }

	/* create context */
	pg_t *pg = calloc(1, sizeof(pg_t) + pt_nth(pt) * sizeof(pg_t *));
	*pg = (pg_t){
		.fp = fp, .pt = pt,
		.lb = pt_nth(pt), .ub = 3 * pt_nth(pt), .bal = 0, .nth = pt_nth(pt),
		.block_size = 1024 * 1024
	};
	kv_hq_init(pg->hq);

	/* init worker args */
	// for(uint64_t i = 0; i < pt_nth(pt); i++) { pg->c[i] = (void *)pg; }
	pt_set_worker(pg->pt, pg, pg_worker);
	return(pg);
}

/**
 * @fn pg_freeze
 * @brief clear ptask queues
 */
static _force_inline
void pg_freeze(pg_t *pg)
{
	/* process current working block */
	pg_block_t *s = pg->s, *t;
	if(s && s->flush == 1 && s->head != 0) {
		s->len = s->head;
		if(pg->nth == 1) {
			pg_write_block(pg, pg_deflate(s, pg->block_size));
		} else {
			pg->bal++; pt_enq(&pg->pt->in, 0, s);
		}
		pg->s = NULL;
	}

	/* process remainings */
	while(pg->bal > 0) {
		if((t = pt_deq(&pg->pt->out, 0)) == PT_EMPTY) {
			/* wait for a while(note: nanosleep is better?) */
			sched_yield(); continue;
		}
		pg->bal--;
		kv_hq_push(v4u32_t, incq_comp, pg->hq, ((v4u32_t){.u64 = {t->id, (uintptr_t)t}}));
	}
}

/**
 * @fn pg_destroy
 */
static _force_inline
void pg_destroy(pg_t *pg)
{
	if(pg == NULL) { return; }
	pg_freeze(pg); free(pg->s);

	/* flush heapqueue */
	while(pg->hq.n > 1) {
		pg->ocnt++;
		pg_block_t *t = (pg_block_t*)kv_hq_pop(v4u32_t, incq_comp, pg->hq).u64[1];
		if(t && t->flush) { pg_write_block(pg, t); } else { free(t); }
	}

	/* write terminator */
	uint64_t z = 0;
	fwrite(&z, sizeof(uint64_t), 1, pg->fp);
	
	/* cleanup contexts */
	kv_hq_destroy(pg->hq);
	free(pg);
	return;
}

/**
 * @fn pgread
 */
static _force_inline
pg_block_t *pg_read_single(pg_t *pg)
{
	/* single-threaded */
	pg_block_t *t;
	if((t = pg_read_block(pg)) == NULL) { pg->eof = 2; return(NULL); }
	return(pg_inflate(t, pg->block_size));
}
static _force_inline
pg_block_t *pg_read_multi(pg_t *pg)
{
	/* multithreaded; read compressed blocks and push them to queue */
	pg_block_t *t;
	while(!pg->eof && pg->bal < pg->ub) {
		if((t = pg_read_block(pg)) == NULL) { pg->eof = 1; break; }
		pg->bal++;
		pt_enq(&pg->pt->in, 0, t);
	}

	/* fetch inflated blocks and push heapqueue to sort */
	while((t = pt_deq(&pg->pt->out, 0)) != PT_EMPTY) {
		pg->bal--;
		kv_hq_push(v4u32_t, incq_comp, pg->hq, ((v4u32_t){.u64 = {t->id, (uintptr_t)t}}));
	}

	/* check if input depleted */
	if(pg->ocnt >= pg->icnt) { pg->eof = 2; return(NULL); }

	/* check if inflated block is available */
	if(pg->hq.n < 2 || pg->hq.a[1].u64[0] > pg->ocnt) {
		sched_yield(); return(NULL);	/* wait for a while */
	}

	/* block is available! */
	pg->ocnt++;
	return((pg_block_t *)kv_hq_pop(v4u32_t, incq_comp, pg->hq).u64[1]);
}
static _force_inline
uint64_t pgread(pg_t *pg, void *dst, uint64_t len)
{
	uint64_t rem = len;
	pg_block_t *s = pg->s;
	if(pg->eof == 2) return(0);
	if(pg->pt->c->wfp != pg_worker && pt_set_worker(pg->pt, pg, pg_worker)) {
		return(0);
	}

	static pg_block_t *(*const fp[2])(pg_t *) = { pg_read_single, pg_read_multi };
	while(rem > 0) {
		/* check and prepare a valid inflated block */
		while(s == NULL || s->head == s->len) {
			free(s); pg->s = s = fp[pg->nth > 1](pg);
			if(pg->eof > 1) { return(len - rem); }
		}

		/* copy to dst buffer */
		uint64_t adv = MIN2(rem, s->len - s->head);
		memcpy(dst + len - rem, s->buf + s->head, adv);
		rem -= adv; s->head += adv;
	}
	return(len);
}

/**
 * @fn pgwrite
 */
static _force_inline
void pg_write_single(pg_t *pg, pg_block_t *s)
{
	/* single-threaded */
	if(s != NULL) { pg_write_block(pg, pg_deflate(s, pg->block_size)); }
	return;
}
static _force_inline
void pg_write_multi(pg_t *pg, pg_block_t *s)
{
	/* push the current block to deflate queue */
	if(s != NULL) { pg->bal++; pt_enq(&pg->pt->in, 0, s); }

	/* fetch copressed block and push to heapqueue to sort */
	pg_block_t *t;
	while(pg->bal > pg->lb) {
		if((t = pt_deq(&pg->pt->out, 0)) == PT_EMPTY) {
			if(pg->bal < pg->ub) { break; }	/* unexpected error */
			sched_yield(); continue;		/* queue full, wait for a while */
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
	return;
}
static _force_inline
uint64_t pgwrite(pg_t *pg, const void *src, uint64_t len)
{
	uint64_t rem = len;
	pg_block_t *s = pg->s;
	if(pg->pt->c->wfp != pg_worker && pt_set_worker(pg->pt, pg, pg_worker)) {
		return(0);
	}

	void (*const fp[2])(pg_t *, pg_block_t *) = { pg_write_single, pg_write_multi };
	while(rem > 0) {
		/* push the current block to queue and prepare an empty one if the current one is full */
		if(s == NULL || s->head == s->len) {
			fp[pg->nth > 1](pg, s);
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

#if 0
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
#endif
/* end of queue.c */
/* bseq.h */
/**
 * @struct bseq_params_t
 */
typedef struct {
	uint64_t batch_size;					/* buffer (block) size */
	uint32_t keep_qual, min_len;			/* 1 to keep quality string, minimum length cutoff (to filter out short seqs) */
	uint32_t n_tag;
	uint16_t const *tag;					/* tags to be preserved (bam), "CO" to comment in fasta */
} bseq_params_t;

/**
 * @struct bseq_seq_t
 * @brief record container
 */
typedef struct {
	uint32_t l_seq;							/* sequence length */
	uint16_t l_name, n_tag;					/* name length, #tags */
	char *name;								/* pointers */
	uint8_t *seq, *qual, *tag;
	uint64_t u64;							/* reserved for alignment result */
} bseq_seq_t;
_static_assert(sizeof(bseq_seq_t) == 48);

/**
 * @struct bseq_t
 */
typedef struct {
	uint64_t u64[3];						/* open for any use */
	uint32_t u32, n_seq;					/* #sequences */
	void *base;								/* sequence memory base pointer */
	uint64_t size;							/* sequence memory size */
	bseq_seq_t seq[];						/* array of sequence objects follows */
} bseq_t;
_static_assert(sizeof(bseq_t) == sizeof(bseq_seq_t));
typedef struct { size_t n, m; bseq_seq_t *a; } bseq_seq_v;
/* bseq.h */

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
	header->text = (char*)calloc(1, header->l_text + 1);
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
	uint64_t n, m;							/* input buffer, aligns to kvec_t */
	uint8_t *a, *p, *t;						/* input buffer, base, current pointer, and tail pointer */
	gzFile fp;
	bam_header_t *bh;
	uint64_t acc;
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
	bseq_params_t *b,
	char const *fn)
{
	/* open stream */
	gzFile f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if(f == NULL) { return(NULL); }

	/* create instance */
	bseq_file_t *fp = (bseq_file_t *)calloc(1, sizeof(bseq_file_t));
	*fp = (bseq_file_t){ .fp = f, .keep_qual = b->keep_qual, .min_len = b->min_len };

	/* determine file type; allow some invalid spaces at the head */
	for(uint64_t i = 0; i < 4; i++) {
		int c;
		if((c = gzgetc(fp->fp)) == 'B') {	/* test bam signature */
			gzungetc(c, fp->fp); fp->bh = bam_read_header(fp->fp); break;
		} else if(c == '>' || c == '@') {	/* test fasta/q delimiter */
			gzungetc(c, fp->fp); fp->delim = c; break;
		}
	}
	if(!fp->bh && !fp->delim) { free(fp); return(NULL); }

	/* init buffer */
	kv_reserve(uint8_t, *fp, b->batch_size);
	fp->n = b->batch_size;

	/* init tags */
	fp->tags = calloc(((fp->l_tags = b->n_tag) + 0x1f) & ~0x1f, sizeof(uint16_t));
	if(b->n_tag && b->tag) { memcpy(fp->tags, b->tag, b->n_tag * sizeof(uint16_t)); }
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
	free(fp->a);
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
	bseq_seq_v *restrict seq,				/* sequence metadata array */
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
	l_tag = ((uint8_t *)fp->t) - stag;							/* FIXME: tag section length */

	/* reserve memory */
	uint64_t req_size = (
		  c->l_qname + 1 + c->l_qseq + 1 						/* name and seq with '\0' */
		+ ((fp->keep_qual && *squal != 0xff)					/* quality string if available */
			? c->l_qseq
			: 0) + 1
		+ (fp->l_tags? l_tag : 0) + 1							/* tag section */
	);
	kv_reserve(uint8_t, *mem, mem->n + req_size);
	bseq_seq_t *s = kv_pushp(bseq_seq_t, *seq);

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
			mem->a[mem->n++] = enc4r[0x0f & (sseq[i>>1]>>((~i&0x01)<<2))];
		}
	} else {
		for(int64_t i = 0; i < c->l_qseq; i++) {
			mem->a[mem->n++] = enc4f[0x0f & (sseq[i>>1]>>((~i&0x01)<<2))];
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
	_p += _len - 32; _q += _len - 32; _q -= _q[-1] == '\r'; _m1>>_len; \
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
 * @fn bseq_read_fasta
 * @brief parse one sequence, returns 0 for success, 1 for buffer starvation, 2 for broken format
 */
static _force_inline
uint64_t bseq_read_fasta(
	bseq_file_t *restrict fp,
	bseq_seq_v *restrict seq,		/* src */
	uint8_v *restrict mem)			/* dst, must have enough space (e.g. 2 * buffer) */
{
	#define _id(x)					(x)
	#define _trans(x)				( _shuf_v32i8(cv, _and_v32i8(fv, x)) )
	#define _forward_state(_state)	fp->state = _state; case _state
	#define _cp()					if(_unlikely(p >= t)) { goto _refill; }

	/* keep them on registers */
	v32i8_t const dv = _set_v32i8(fp->delim == '@'? '+' : fp->delim);
	v32i8_t const sv = _set_v32i8(' '), lv = _set_v32i8('\n'), fv = _set_v32i8(0xf);
	v32i8_t const cv = _from_v16i8_v32i8(_load_v16i8(encaf));

	bseq_seq_t *s = &seq->a[seq->n-1];			/* restore previous states */
	uint8_t *p = fp->p, *q = &mem->a[mem->n];	/* load pointers */
	uint8_t const *t = fp->t;
	uint64_t m;
	debug("enter, state(%u), p(%p), t(%p), eof(%u)", fp->state, p, t, fp->is_eof);
	if(_unlikely(p >= t)) { return(1); }
	switch(fp->state) {							/* dispatcher for the first iteration */
		default:								/* idle or broken */
		_forward_state(1):						/* waiting header */
			if(*p++ != fp->delim) { return(0); }/* broken */
			s = kv_pushp(bseq_seq_t, *seq);		/* create new sequence */
		_forward_state(2):						/* transition to spaces between delim and name */
			_strip(p, t, sv); _cp();
			s->name = (char*)_beg(q, mem->a);
		_forward_state(3):
			m = _readline(p, t, q, sv, _id); _cp();
			p++;								/* skip '\n' or ' ' after sequence name */
			s->l_name = _term(q, mem->a, s->name);
			s->tag = _beg(q, mem->a);			/* prepare room for tag before the third state, to use m before it vanishes */
			s->n_tag = m & fp->keep_comment;	/* set n_tag if comment line found */
			if(s->n_tag == 0) { goto _seq_head; }
		_forward_state(4):						/* spaces between name and comment */
			_strip(p, t, sv); _cp();
			*q++ = 'C'; *q++ = 'O'; *q++ = 'Z';
		_forward_state(5):						/* parsing comment */
			_readline(p, t, q, lv, _id); _cp();	/* refill needed, comment continues */
			p++;								/* skip '\n' after comment */
			while(q[-1] == ' ') { q--; }		/* strip spaces at the tail of the comment */
		_seq_head:
			_term(q, mem->a, s->tag);
			s->seq = _beg(q, mem->a);
		_forward_state(6):						/* parsing seq */
			while(1) {
				m = _readline(p, t, q, dv, _trans);
				/*
				 * EOF, which is not detected by _readline, found. mark by ORing m with the eof status flag (is_eof == 1)
				 * for the sequence array to be correctly capped by the _term macro in the follwing several lines.
				 */
				if(_unlikely(p >= t)) { m |= fp->is_eof; break; }
				if(m & 0x01) { break; }
				p++;							/* skip '\n' */
			}
			if((m & 0x01) == 0) { goto _refill; }/* buffer starved but not yet reached the end of the sequence section */
			s->l_seq = _term(q, mem->a, s->seq);
			s->qual = _beg(q, mem->a);
			if(fp->delim == '>' || _unlikely(p >= t)) { goto _qual_tail; }/* here p >= t only holds when EOF is detected */
		_forward_state(7):
			_skipline(p, t); _cp();
			p++; fp->acc = 0;					/* skip '\n', clear accumulator */
		_forward_state(8):;						/* parsing qual */
			uint64_t acc = fp->acc, lim = s->l_seq;/* load accumulator */
			while(1) {
				if(!fp->keep_qual) { acc += _skipline(p, t); }
				else { uint8_t const *b = q; _readline(p, t, q, lv, _id); acc += q - b; }
				if(_unlikely(p >= t)) { fp->acc = acc; goto _refill; }
				if(acc >= lim) { break; }
				p++;							/* skip '\n' */
			}
		_forward_state(9):
			_cp(); _strip(p, t, lv);
		_qual_tail:
			_term(q, mem->a, s->qual);			/* fall through to go back to idle */
	}
	debug("break, state(%u), eof(%u), name(%s), len(%u)", fp->state, fp->is_eof, mem->a + (uintptr_t)s->name, s->l_seq);
	fp->state = 0;								/* back to idle */
	if((uint32_t)s->l_seq < fp->min_len) { seq->n--; q = mem->a + (uintptr_t)s->name; }		/* squash if seq is short */
_refill:
	debug("return, p(%p), t(%p)", p, t);
	fp->p = p; mem->n = q - mem->a;				/* write back pointers */
	return(p >= t);

	#undef _id
	#undef _trans
	#undef _forward_state
	#undef _cp
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
bseq_t *bseq_read(bseq_file_t *fp)
{
	uint8_v mem = { 0 };
	bseq_seq_v seq = { 0 };
	static uint8_t const margin[BSEQ_MGN] = { 0 };

	/* check the context is valid */
	if(fp->is_eof > 1) { return(NULL); }

	#define _readp(_l)	({ \
		kv_reserve(uint8_t, *fp, _l);		/* dead code when _l == fp->n */ \
		uint64_t _r = gzread(fp->fp, fp->a, _l); \
		fp->p = fp->a; fp->t = fp->a + _r; \
		(_r != _l) ? (fp->is_eof = 1 + (_r == 0)) : 0; \
	})
	#define _reada(type)	({ type _n; if(gzread(fp->fp, &_n, sizeof(type)) != sizeof(type)) { _n = 0; } _n; })

	/* reserve bseq_t space at the head of bseq_seq_t array */
	kv_push(bseq_seq_t, seq, ((bseq_seq_t){ 0 }));			/* (bseq_t){ .base = NULL, .size = 0, .n_seq = 0 } */

	/* reserve mem */
	kv_reserve(uint8_t, mem, 2 * fp->n);
	kv_pushm(uint8_t, mem, margin, BSEQ_MGN);				/* head margin */
	if(fp->bh) {	/* bam */
		while(mem.n < fp->n) {
			uint64_t size = _reada(uint32_t);				/* size might be larger than fp->n */
			if(size == 0 || _readp(size) > 0) { break; }
			bseq_read_bam(fp, &seq, &mem);
		}
	} else {		/* fasta/q */
		while(mem.n < fp->n + BSEQ_MGN) {					/* fetch-and-parse loop */
			while(bseq_read_fasta(fp, &seq, &mem) == 1) {	/* buffer starved */
				if(_readp(fp->n) > 1) { break; }			/* fetch next */
				_storeu_v64i8(fp->t, _zero_v64i8());		/* fill margin of the input buffer */
				kv_reserve(uint8_t, mem, mem.n + 2 * fp->n);/* reserve room for the next parsing unit */
			}
			debug("finished seq, state(%u), is_eof(%u)", fp->state, fp->is_eof);
			if(fp->is_eof > 1) { break; }
			if(fp->state != 0) { goto _bseq_read_fail; }	/* error occurred */
		}
	}
	kv_pushm(uint8_t, mem, margin, BSEQ_MGN);				/* tail margin */

	/* adjust pointers */
	ptrdiff_t b = (ptrdiff_t)mem.a;
	kv_foreach(bseq_seq_t, seq, { p->name += b; p->seq += b; p->qual += b; p->tag += b; });

	/* store info to the bseq_t object, allocated at the head of the bseq_seq_t array */
	bseq_t *s = (bseq_t *)seq.a;
	*s = (bseq_t){ .base = mem.a, .size = mem.n, .n_seq = seq.n - 1, .u32 = 0, .u64 = { 0 } };
	fp->n_seq += seq.n - 1;									/* update sequence counter */
	return(s);

_bseq_read_fail:;
	free(mem.a); fp->is_eof = 2;
	return(NULL);

	#undef _readp
	#undef _reada
}

#if 0
unittest( .name = "bseq.fasta" ) {
	char const *filename = "./minialign.unittest.bseq.tmp";
	char const *content =
		">test0\nAAAA\n"
		"> test1\nATAT\nCGCG\r\n\r\n"
		">  test2\n\nAAAA\n"
		">test3 comment comment  \nACGT\n\n";

	FILE *fp = fopen(filename, "w");
	assert(fp != NULL);
	fwrite(content, 1, strlen(content), fp);
	fclose(fp);

	uint16_t const tags[1] = { ((uint16_t)'C') | (((uint16_t)'O')<<8) };
	bseq_file_t *b = bseq_open(filename, 64, 1, 0, 1, tags);
	assert(b != NULL);
	bseq_t *s = bseq_read(b);

	assert(s != NULL);
	assert(s->n_seq == 4, "n_seq(%u)", s->n_seq);
	assert(s->base != NULL);
	assert(s->size > 0, "size(%lu)", s->size);

	assert(s->seq[0].l_name == 5, "l_name(%u)", s->seq[0].l_name);
	assert(s->seq[0].l_seq == 4, "l_seq(%u)", s->seq[0].l_seq);
	assert(s->seq[0].n_tag == 0, "n_tag(%u)", s->seq[0].n_tag);
	assert(strcmp((const char*)s->seq[0].name, "test0") == 0, "name(%s)", s->seq[0].name);
	assert(strcmp((const char*)s->seq[0].seq, "\x0\x0\x0\x0") == 0, "seq(%s)", s->seq[0].seq);
	assert(strcmp((const char*)s->seq[0].qual, "") == 0, "qual(%s)", s->seq[0].qual);
	assert(strcmp((const char*)s->seq[0].tag, "") == 0, "tag(%s)", s->seq[0].tag);

	assert(s->seq[1].l_name == 5, "l_name(%u)", s->seq[1].l_name);
	assert(s->seq[1].l_seq == 8, "l_seq(%u)", s->seq[1].l_seq);
	assert(s->seq[1].n_tag == 0, "n_tag(%u)", s->seq[1].n_tag);
	assert(strcmp((const char*)s->seq[1].name, "test1") == 0, "name(%s)", s->seq[1].name);
	assert(strcmp((const char*)s->seq[1].seq, "\x0\x3\x0\x3\x1\x2\x1\x2") == 0, "seq(%s)", s->seq[1].seq);
	assert(strcmp((const char*)s->seq[1].qual, "") == 0, "qual(%s)", s->seq[1].qual);
	assert(strcmp((const char*)s->seq[1].tag, "") == 0, "tag(%s)", s->seq[1].tag);

	assert(s->seq[2].l_name == 5, "l_name(%u)", s->seq[2].l_name);
	assert(s->seq[2].l_seq == 4, "l_seq(%u)", s->seq[2].l_seq);
	assert(s->seq[2].n_tag == 0, "n_tag(%u)", s->seq[2].n_tag);
	assert(strcmp((const char*)s->seq[2].name, "test2") == 0, "name(%s)", s->seq[2].name);
	assert(strcmp((const char*)s->seq[2].seq, "\x0\x0\x0\x0") == 0, "seq(%s)", s->seq[2].seq);
	assert(strcmp((const char*)s->seq[2].qual, "") == 0, "qual(%s)", s->seq[2].qual);
	assert(strcmp((const char*)s->seq[2].tag, "") == 0, "tag(%s)", s->seq[2].tag);

	assert(s->seq[3].l_name == 5, "l_name(%u)", s->seq[3].l_name);
	assert(s->seq[3].l_seq == 4, "l_seq(%u)", s->seq[3].l_seq);
	assert(s->seq[3].n_tag == 1, "n_tag(%u)", s->seq[3].n_tag);
	assert(strcmp((const char*)s->seq[3].name, "test3") == 0, "name(%s)", s->seq[3].name);
	assert(strcmp((const char*)s->seq[3].seq, "\x0\x1\x2\x3") == 0, "seq(%s)", s->seq[3].seq);
	assert(strcmp((const char*)s->seq[3].qual, "") == 0, "qual(%s)", s->seq[3].qual);
	assert(strcmp((const char*)s->seq[3].tag, "COZcomment comment") == 0, "tag(%s)", s->seq[3].tag);

	uint64_t n_seq = bseq_close(b);
	assert(n_seq == 4);
	remove(filename);
}

unittest( .name = "bseq.fastq" ) {
	char const *filename = "./minialign.unittest.bseq.tmp";
	char const *content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\n12+3\n+123\r\n"
		"@  test2\n\nAAAA\n+  test2\n\n\n12@3\n\n"
		"@test3  comment comment   \nACGT\n\n+ test3\n@123";

	FILE *fp = fopen(filename, "w");
	assert(fp != NULL);
	fwrite(content, 1, strlen(content), fp);
	fclose(fp);

	uint16_t const tags[1] = { ((uint16_t)'C') | (((uint16_t)'O')<<8) };
	bseq_file_t *b = bseq_open(filename, 256, 1, 0, 1, tags);
	assert(b != NULL);
	bseq_t *s = bseq_read(b);

	assert(s != NULL);
	assert(s->n_seq == 4, "n_seq(%u)", s->n_seq);
	assert(s->base != NULL);
	assert(s->size > 0, "size(%lu)", s->size);

	assert(s->seq[0].l_name == 5, "l_name(%u)", s->seq[0].l_name);
	assert(s->seq[0].l_seq == 4, "l_seq(%u)", s->seq[0].l_seq);
	assert(s->seq[0].n_tag == 0, "n_tag(%u)", s->seq[0].n_tag);
	assert(strcmp((const char*)s->seq[0].name, "test0") == 0, "name(%s)", s->seq[0].name);
	assert(strcmp((const char*)s->seq[0].seq, "\x0\x0\x0\x0") == 0, "seq(%s)", s->seq[0].seq);
	assert(strcmp((const char*)s->seq[0].qual, "NNNN") == 0, "qual(%s)", s->seq[0].qual);
	assert(strcmp((const char*)s->seq[0].tag, "") == 0, "tag(%s)", s->seq[0].tag);

	assert(s->seq[1].l_name == 5, "l_name(%u)", s->seq[1].l_name);
	assert(s->seq[1].l_seq == 8, "l_seq(%u)", s->seq[1].l_seq);
	assert(s->seq[1].n_tag == 0, "n_tag(%u)", s->seq[1].n_tag);
	assert(strcmp((const char*)s->seq[1].name, "test1") == 0, "name(%s)", s->seq[1].name);
	assert(strcmp((const char*)s->seq[1].seq, "\x0\x3\x0\x3\x1\x2\x1\x2") == 0, "seq(%s)", s->seq[1].seq);
	assert(strcmp((const char*)s->seq[1].qual, "12+3+123") == 0, "qual(%s)", s->seq[1].qual);
	assert(strcmp((const char*)s->seq[1].tag, "") == 0, "tag(%s)", s->seq[1].tag);

	assert(s->seq[2].l_name == 5, "l_name(%u)", s->seq[2].l_name);
	assert(s->seq[2].l_seq == 4, "l_seq(%u)", s->seq[2].l_seq);
	assert(s->seq[2].n_tag == 0, "n_tag(%u)", s->seq[2].n_tag);
	assert(strcmp((const char*)s->seq[2].name, "test2") == 0, "name(%s)", s->seq[2].name);
	assert(strcmp((const char*)s->seq[2].seq, "\x0\x0\x0\x0") == 0, "seq(%s)", s->seq[2].seq);
	assert(strcmp((const char*)s->seq[2].qual, "12@3") == 0, "qual(%s)", s->seq[2].qual);
	assert(strcmp((const char*)s->seq[2].tag, "") == 0, "tag(%s)", s->seq[2].tag);

	assert(s->seq[3].l_name == 5, "l_name(%u)", s->seq[3].l_name);
	assert(s->seq[3].l_seq == 4, "l_seq(%u)", s->seq[3].l_seq);
	assert(s->seq[3].n_tag == 1, "n_tag(%u)", s->seq[3].n_tag);
	assert(strcmp((const char*)s->seq[3].name, "test3") == 0, "name(%s)", s->seq[3].name);
	assert(strcmp((const char*)s->seq[3].seq, "\x0\x1\x2\x3") == 0, "seq(%s)", s->seq[3].seq);
	assert(strcmp((const char*)s->seq[3].qual, "@123") == 0, "qual(%s)", s->seq[3].qual);
	assert(strcmp((const char*)s->seq[3].tag, "COZcomment comment") == 0, "tag(%s)", s->seq[3].tag);

	uint64_t n_seq = bseq_close(b);
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
	bseq_t *s = bseq_read(b);

	assert(s != NULL);
	assert(s->n_seq == 4, "n_seq(%u)", s->n_seq);
	assert(s->base != NULL);
	assert(s->size > 0, "size(%lu)", s->size);

	assert(s->seq[0].l_name == 5, "l_name(%u)", s->seq[0].l_name);
	assert(s->seq[0].l_seq == 4, "l_seq(%u)", s->seq[0].l_seq);
	assert(s->seq[0].n_tag == 0, "n_tag(%u)", s->seq[0].n_tag);
	assert(strcmp((const char*)s->seq[0].name, "test0") == 0, "name(%s)", s->seq[0].name);
	assert(strcmp((const char*)s->seq[0].seq, "\x0\x0\x0\x0") == 0, "seq(%s)", s->seq[0].seq);
	assert(strcmp((const char*)s->seq[0].qual, "") == 0, "qual(%s)", s->seq[0].qual);
	assert(strcmp((const char*)s->seq[0].tag, "") == 0, "tag(%s)", s->seq[0].tag);

	assert(s->seq[1].l_name == 5, "l_name(%u)", s->seq[1].l_name);
	assert(s->seq[1].l_seq == 8, "l_seq(%u)", s->seq[1].l_seq);
	assert(s->seq[1].n_tag == 0, "n_tag(%u)", s->seq[1].n_tag);
	assert(strcmp((const char*)s->seq[1].name, "test1") == 0, "name(%s)", s->seq[1].name);
	assert(strcmp((const char*)s->seq[1].seq, "\x0\x3\x0\x3\x1\x2\x1\x2") == 0, "seq(%s)", s->seq[1].seq);
	assert(strcmp((const char*)s->seq[1].qual, "") == 0, "qual(%s)", s->seq[1].qual);
	assert(strcmp((const char*)s->seq[1].tag, "") == 0, "tag(%s)", s->seq[1].tag);

	assert(s->seq[2].l_name == 5, "l_name(%u)", s->seq[2].l_name);
	assert(s->seq[2].l_seq == 4, "l_seq(%u)", s->seq[2].l_seq);
	assert(s->seq[2].n_tag == 0, "n_tag(%u)", s->seq[2].n_tag);
	assert(strcmp((const char*)s->seq[2].name, "test2") == 0, "name(%s)", s->seq[2].name);
	assert(strcmp((const char*)s->seq[2].seq, "\x0\x0\x0\x0") == 0, "seq(%s)", s->seq[2].seq);
	assert(strcmp((const char*)s->seq[2].qual, "") == 0, "qual(%s)", s->seq[2].qual);
	assert(strcmp((const char*)s->seq[2].tag, "") == 0, "tag(%s)", s->seq[2].tag);

	assert(s->seq[3].l_name == 5, "l_name(%u)", s->seq[3].l_name);
	assert(s->seq[3].l_seq == 4, "l_seq(%u)", s->seq[3].l_seq);
	assert(s->seq[3].n_tag == 1, "n_tag(%u)", s->seq[3].n_tag);
	assert(strcmp((const char*)s->seq[3].name, "test3") == 0, "name(%s)", s->seq[3].name);
	assert(strcmp((const char*)s->seq[3].seq, "\x0\x1\x2\x3") == 0, "seq(%s)", s->seq[3].seq);
	assert(strcmp((const char*)s->seq[3].qual, "") == 0, "qual(%s)", s->seq[3].qual);
	assert(strcmp((const char*)s->seq[3].tag, "COZcomment comment") == 0, "tag(%s)", s->seq[3].tag);

	uint64_t n_seq = bseq_close(b);
	assert(n_seq == 4);
	remove(filename);
}
#endif
/* end of bseq.c */

/* sketch.c */
/**
 * @fn hash64
 */
#define hash64(k0, k1, mask)		( (_mm_crc32_u64((k0), (k0)) ^ (k1)) & (mask) )
// #define hash64(k0, k1, mask)		( (k0) & (mask) )

/**
 * @struct mm_sketch_t
 * @brief minimizer window context, can be carried over (connected to) the next segment.
 */
typedef struct {
	uint32_t w, k;
	uint64_v *b;
	uint64_t r[64];					/* UINT64_MAX */
} mm_sketch_t;
typedef struct {
	uint64_t i, u, k0, k1;
} mm_sketch_cap_t;

static _force_inline
void mm_sketch_init(mm_sketch_t *sk, uint32_t w, uint32_t k, uint64_v *b)
{
	sk->w = w; sk->k = k; sk->b = b;
	_memset_blk_u(sk->r, UINT64_MAX, sizeof(uint64_t) * 32);	/* fill UINT64_MAX */
	return;
}

/**
 * @fn mm_sketch
 * @brief calclulate (w,k)-minimizer sketch for seq.
 * minimizers are stored in uint64_t array, where upper 56bit holds actual hash value (minimizer),
 * and lower 8bit contains local index of the length-w bin. direction is encoded in the 7-th bit.
 */
#define _loop_init(_len, _cap) \
	uint64_t const kk = sk->k - 1, shift1 = 2*kk, mask = (1ULL<<2*sk->k) - 1, w = sk->w; \
	kv_reserve(uint64_t, *sk->b, sk->b->n + 4 * (_len) / w);	/* reserve buffer for minimizers, expected number will be 2 * len / w */ \
	uint64_t *q = sk->b->a + sk->b->n; \
	uint8_t const *p = seq, *t = &seq[_len]; \
	uint64_t u = (_cap)->u, k0 = (_cap)->k0, k1 = (_cap)->k1;	/* forward and backward k-mers, min buffer (w < 32 must be guaranteed) */

#define _push_kmer() { \
	uint64_t c = *p++; \
	k0 = (k0 << 2 | c) & mask; k1 = (k1 >> 2) | ((3ULL^c) << shift1); \
}
#define _loop_core() ({ \
	_push_kmer(); \
	uint64_t km = k0 < k1 ? k0 : k1, kx = k0 < k1 ? k1 : k0, m = k0 < k1 ? 0 : 0x80; \
	uint64_t h = hash64(km, kx, mask)<<8 | i | m; f = MIN2(f, h); uint64_t v = MIN2(f, sk->r[i + 1]); \
	if((v == h) | (v - u)) { \
		debug("push, i(%lu), uhfv(%lx, %lx, %lx, %lx), v(%lu, %lu, %lu), p(%p) pos(%lu)", \
			i, u, h, f, v, v>>8, (v>>7)&1, v&0x7f, p, p - seq - 15); \
		*q++ = v; \
	} \
	u = v; h; \
})
#define _push_cap(_i) ({ \
	mm_sketch_cap_t *cap = (mm_sketch_cap_t *)q; \
	*cap = (mm_sketch_cap_t){ .i = ~(_i), .u = u, .k0 = k0, .k1 = k1 }; \
	sk->b->n = (uint64_t *)(cap + 1) - sk->b->a; \
	cap; \
})
// static _force_inline
mm_sketch_cap_t const *mm_sketch(mm_sketch_t *sk, uint8_t const *seq, uint32_t len)
{
	static mm_sketch_cap_t const init = { 0 };
	_loop_init(len, &init);			/* initialize working buffers */
	debug("p(%p), t(%p), len(%u, %lu)", p, t, len, t - p);
	for(uint64_t i = 0; i < kk && p < t; i++) { _push_kmer(); }
	while(t - p >= w) {
		/* calculate hash values for the current block */
		for(uint64_t i = 0, f = UINT64_MAX; i < w; i++) { sk->r[i] = _loop_core(); }
		/* fold backward-min array */
		for(uint64_t i = 0, r = UINT64_MAX; i < w; i++) { r = MIN2(r, sk->r[w-i-1]); sk->r[w-i-1] = r;  }
	}
	uint64_t l = t - p;				/* remainder */
	debug("remainder l(%lu)", l);
	if(l > 0) {
		/* calculate hash values for the tail */
		for(uint64_t i = 0, f = UINT64_MAX; i < l; i++) { sk->r[w+i] = _loop_core() + w; }
		/* fold backward-min array */
		for(uint64_t i = 0, r = UINT64_MAX; i < w; i++) { r = MIN2(r, sk->r[w+l-i-1]); sk->r[w+l-i-1] = r;  }
		/* move to align the head */
		for(uint64_t i = 0; i < w; i++) { sk->r[i] = sk->r[l+i] - l; }
	}
	return(_push_cap(0));
}
static _force_inline
mm_sketch_cap_t const *mm_sketch_cap(mm_sketch_t *sk, mm_sketch_cap_t const *cap, uint8_t const *seq, uint32_t len)
{
	uint64_t ci = ~cap->i, l = MIN2(len, sk->w - ci);
	_loop_init(l, cap); (void)t;	/* t is unused here */
	for(uint64_t i = ci, f = UINT64_MAX; i < w && p < t; i++) { _loop_core(); }
	return(_push_cap(ci + l));
}
#undef _loop_init
#undef _push_kmer
#undef _loop_core
#undef _push_cap
/* end of sketch.c */

/* index.h */
/**
 * @struct mm_idx_params_t
 */
typedef struct {
	uint8_t b, w, k, n_frq;			/* bucket size (in bits), window and k-mer size */
	float frq[MAX_FRQ_CNT];			/* occurrence array */
	kh_str_t circ;					/* circular ref names */
} mm_idx_params_t;

/**
 * @struct mm_idx_seq_t
 * @brief reference-side sequence container
 */
typedef struct {
	uint8_t const *seq;
	char const *name;
	uint32_t l_seq;
	uint16_t l_name, circular;		/* backward offset of name from seq head */
} mm_idx_seq_t;
_static_assert(sizeof(mm_idx_seq_t) == 24);

/**
 * @struct mm_idx_t
 */
typedef struct mm_idx_bkt_s mm_idx_bkt_t;
typedef struct {
	mm_idx_bkt_t *bkt;				/* array of 2nd stage hash table is the 1st stage hash table */
	uint64_t mask;					/* (internal) must be (1<<b) - 1 */
	uint8_t b, w, k, n_occ;			/* bucket size (in bits), window and k-mer size */
	uint32_t occ[MAX_FRQ_CNT];		/* occurrence array */
	uint32_t n_seq, mono;			/* (internal) sequence buckets and monolithic flag */
	mm_idx_seq_t *s;				/* sequence array */
} mm_idx_t;
/* end of index.h */

/* map.h */
/* flags */
#define MM_AVA			( 0x01ULL )
#define MM_OMIT_REP		( 0x08ULL )			/* omit secondary records */
#define MM_COMP 		( 0x10ULL )

/* forward declaration of mm_align_t */
typedef struct mm_align_s mm_align_t;
typedef struct mm_reg_s mm_reg_t;

/**
 * @struct mm_tbuf_params_t
 */
typedef struct {
	mm_idx_t mi;							/* (immutable) index object */
	uint32_t org, thresh;					/* ava filter constants, calculated from opt->base_rid and opt->base_qid */
	uint32_t flag;
	int32_t wlen, glen;						/* chainable window edge length, linkable gap length */
	float min_ratio;
	uint32_t min_score;
	uint32_t base_rid, base_qid;
	gaba_t *ctx;
	gaba_alloc_t alloc;						/* lmm contained */
} mm_tbuf_params_t;

/**
 * @struct mm_align_params_t
 * @brief DP-based alignment calculation parameters
 */
typedef struct {
	uint32_t flag;
	int32_t wlen, glen;						/* chainable window edge length, linkable gap length */
	float min_ratio;
	uint32_t min_score;
	uint32_t base_rid, base_qid;			/* will be updated */
	gaba_params_t p;						/* extension */
} mm_align_params_t;
/* end of map.h */

/* printer.h */
/* tags */
/* sam tags */
#define MM_RG			( 0 )				/* Z: read group */
#define MM_CO			( 1 )				/* Z: comment */
#define MM_NH			( 2 )				/* i: #hits */
#define MM_IH 			( 3 )				/* i: index of the record in the hits */
#define MM_AS			( 4 )				/* i: score */
#define MM_XS			( 5 )				/* i: suboptimal score */
#define MM_NM 			( 6 )				/* i: editdist to the reference */
#define MM_SA			( 7 )				/* Z: supplementary records */
#define MM_MD			( 8 )				/* Z: mismatch positions */

/* formats */
#define MM_SAM			( 0x00ULL )
#define MM_MAF			( 0x01ULL )
#define MM_BLAST6		( 0x02ULL )
#define MM_BLASR1		( 0x03ULL )
#define MM_BLASR4		( 0x04ULL )
#define MM_PAF			( 0x05ULL )
#define MM_MHAP 		( 0x06ULL )
#define MM_FALCON		( 0x07ULL )

/* forward decls */
typedef struct mm_print_s mm_print_t;
typedef void (*mm_print_header_t)(mm_print_t *b, uint32_t n_seq, mm_idx_seq_t const *seq);
typedef void (*mm_print_mapped_t)(mm_print_t *b, mm_idx_seq_t const *ref, bseq_seq_t const *t, mm_reg_t const *reg);
static void mm_print_header(mm_print_t *b, uint32_t n_seq, mm_idx_seq_t const *seq);
static void mm_print_mapped(mm_print_t *b, mm_idx_seq_t const *ref, bseq_seq_t const *t, mm_reg_t const *reg);

/**
 * @struct mm_print_params_t
 */
typedef struct {
	uint64_t outbuf_size;			/* buffer sizes */
	uint32_t flag, format;
	uint32_t n_tag;
	uint16_t const *tag;
	char *arg_line;
	char *rg_line, *rg_id;
} mm_print_params_t;
/* end of printer.h */

/* opt.h */
/**
 * @struct mm_opt_parser_t
 * @brief option parser function prototype
 */
typedef struct mm_opt_s mm_opt_t;
typedef struct {
	uint8_t type;							/* 1 for boolean, 2 for option mandatory */
	void (*fn)(mm_opt_t *o, char const *arg);
} mm_opt_parser_t;
/**
 * @struct mm_opt_t
 * @brief parameter, logger, and parser container
 */
struct mm_opt_s {
	ptr_v parg;
	char *fnw;
	uint32_t nth, help;
	uint16_v tags;
	bseq_params_t b;
	mm_idx_params_t c;						/* index params */
	mm_align_params_t a;					/* chain and align params */
	mm_print_params_t r;

	pt_t *pt;
	double inittime;
	uint32_t verbose, ecnt;					/* error counter */
	int (*log)(struct mm_opt_s const *, char, char const *, char const *, ...);
	void *fp;
	mm_opt_parser_t t[256];					/* parser functions */
};
/**
 * @type mm_fprintf_t
 * @brief log printer
 */
typedef int (*mm_log_t)(mm_opt_t const *o, char level, char const *func, char const *fmt, ...);
/**
 * @fn mm_log_printer
 * @brief 0, 1, 2,... for normal message, 8, 9, 10,... for message with timestamp, 16, 17, 18, ... for without header.
 */
static
int mm_log_printer(
	mm_opt_t const *o,		/* option object */
	char level,				/* 'E' and 'W' for error and warning, 0, 1,... for message */
	char const *func,		/* __func__ must be passed */
	char const *fmt,		/* format string */
	...)
{
	if (level < ' ' && (level & 0x07) > o->verbose) {
		return(0);
	}

	va_list l;
	va_start(l, fmt);

	FILE *fp = (FILE *)o->fp;
	int r = 0;
	if(level >= ' ' || (level & 0x10) == 0) {
		if(level >= ' ' || (level & 0x08) == 0) {
			r += fprintf(fp, "[%c::%s] ", level < ' '? 'M' : level, func);
		} else {
			r += fprintf(fp, "[%c::%s::%.3f*%.2f] ",
				level < ' '? 'M' : level,					/* 'E' for error */
				func,										/* function name */
				realtime() - o->inittime,					/* realtime */
				cputime() / (realtime() - o->inittime));	/* average cpu usage */
		}
	}
	r += vfprintf(fp, fmt, l);								/* body */
	r += fprintf(fp, "\n");
	va_end(l);
	return(r);
}
/* end of opt.h */
/* index.c */
/**
 * @struct mm_mini_t
 * @brief minimizer container, alias to v4u32_t and mm_resc_t
 */
typedef struct {
	uint64_t hrem;
	uint32_t pos, rid;
} mm_mini_t;
_static_assert(sizeof(mm_mini_t) == sizeof(v4u32_t));
typedef struct { size_t n, m; mm_mini_t *a; } mm_mini_v;

/**
 * @struct mm_idx_mem_t
 */
typedef struct { uint64_t size; void *base; } mm_idx_mem_t;

/**
 * @struct mm_idx_intl_t
 * @brief multithreaded index construction pipeline context container
 */
typedef struct {
	mm_idx_t mi;
	uint32_t nth, icnt, ocnt;
	bseq_file_t *fp;
	kh_str_t const *circ;
	uint32_t call, ctest;
	kvec_t(mm_idx_seq_t) svec;
	kvec_t(mm_idx_mem_t) mvec;
	kvec_t(v4u32_t) hq;
	uint32_v *cnt;					/* jagged array, counting minimizer occurrences */
} mm_idx_intl_t;

/**
 * @struct mm_idx_bkt_t
 * @brief first-stage data container of double hash table
 */
struct mm_idx_bkt_s {
	union { kh_t h; mm_mini_v a;} w;						/* 2nd stage hash table, and minimizer array for sorting */
	union { uint64_t *p; struct { uint32_t keys, single; } n; } v;	/* last stage table (size in p[0]), #keys and #single-element keys */
};
_static_assert(sizeof(mm_idx_bkt_t) == 32);

/**
 * @fn mm_idx_destroy
 */
static _force_inline
void mm_idx_destroy(mm_idx_t *mi)
{
	if(mi == NULL || mi->mono == 1) { free(mi); return; }	/* loaded by mm_idx_load */

	/* hash table buckets */
	for(uint64_t i = 0; i < 1ULL<<mi->b; i++) {
		free(mi->bkt[i].w.h.a);
		free(mi->bkt[i].v.p);
	}
	free(mi->bkt);

	/* sequence containers */
	mm_idx_intl_t *mii = (mm_idx_intl_t *)mi;
	for(uint64_t i = 0; i < mii->mvec.n; i++) { free(mii->mvec.a[i].base); }
	free(mii->mvec.a);
	free(mii->svec.a);
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
	mm_idx_bkt_t const *b = &mi->bkt[minier & mi->mask];
	kh_t const *h = &b->w.h;
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
		return((v2u32_t const *)&b->v.p[(*p>>32) & 0x7fffffff]);
	}
}

/******************
 * Generate index *
 ******************/
/**
 * @struct mm_idx_step_t
 */
typedef struct {
	uint64_v a;						/* minimizer array */
	uint32_t id;					/* bin id (to sort packets) */
} mm_idx_step_t;
_static_assert(offsetof(mm_idx_step_t, id) == offsetof(bseq_t, u32));

/**
 * @fn mm_idx_source
 * @brief sequence reader (source) of the index construction pipeline
 */
static 
void *mm_idx_source(uint32_t tid, void *arg)
{
	mm_idx_intl_t *mii = (mm_idx_intl_t *)arg;
	bseq_t *r = bseq_read(mii->fp);			/* fetch sequence block */
	if(r == NULL) { return(NULL); }

	mm_idx_step_t *s = (mm_idx_step_t *)r;	/* create packet, use unused block at the head (see bseq_t) */
	kv_reserve(uint64_t, s->a, 4 * r->size / mii->mi.w);
	s->id = mii->icnt++;					/* assign block id */
	return(s);								/* pass to minimizer calculation stage */
}

/**
 * @fn mm_idx_worker
 * @brief calculate minimizer
 */
static
void *mm_idx_worker(uint32_t tid, void *arg, void *item)
{
	mm_idx_intl_t *mii = (mm_idx_intl_t *)arg;
	mm_idx_step_t *s = (mm_idx_step_t *)item;
	bseq_t *r = (bseq_t *)s;				/* overlapped */

	mm_sketch_t sk;
	for(uint64_t i = 0; i < r->n_seq; i++) {
		mm_sketch_init(&sk, mii->mi.w, mii->mi.k, &s->a);
		mm_sketch_cap_t const *cap = mm_sketch(&sk, r->seq[i].seq, r->seq[i].l_seq);

		/* tail margin for circular sequences */
		uint64_t c = mii->call | (mii->ctest && kh_str_get(mii->circ, r->seq[i].name, r->seq[i].l_name) != NULL);
		if(c) { mm_sketch_cap(&sk, cap, r->seq[i].seq, r->seq[i].l_seq); }			/* nori-shiro */
		r->seq[i].u64 = (s->a.n<<1) | c;	/* (#minimizers: 63, circular:1) */
		debug("c(%lu), n(%lu)", c, s->a.n);
	}
	return(s);
}

/**
 * @fn mm_idx_drain
 * @brief push minimizers to bins
 */
static _force_inline
void mm_idx_drain_intl(mm_idx_intl_t *mii, mm_idx_step_t *s)
{
	debug("drain, n(%zu)", s->a.n);

	/* save memory block pointers */
	bseq_t const *r = (bseq_t const *)s;	/* overlaps with mm_idx_step_t */
	kv_push(mm_idx_mem_t, mii->mvec, ((mm_idx_mem_t){ .size = r->size, .base = r->base }));	/* register fetched block */

	/* copy sequence objects, prepare src and dst */
	bseq_seq_t const *src = r->seq;
	kv_reserve(mm_idx_seq_t, mii->svec, mii->svec.n + r->n_seq);

	mm_idx_bkt_t *bkt = mii->mi.bkt;
	uint64_t mask = mii->mi.mask, w = mii->mi.w, b = mii->mi.b;
	for(uint64_t i = 0, *p = s->a.a; i < r->n_seq; i++) {
		/* push sequence object */
		mii->svec.a[mii->svec.n] = (mm_idx_seq_t){
			.seq = src->seq, .l_seq = src->l_seq,
			.name = src->name, .l_name = src->l_name,
			.circular = src->u64 & 0x01
		};
		/* push minimizers */
		uint64_t base = 0, v = w;
		for(uint64_t *t = s->a.a + (src->u64>>1); p < t; p += 4) {
			for(; (*p>>8) != 0xffffffffffffff; p++) {
				uint64_t u = *p & 0x7f, fr = (*p>>7) & 0x01, h = *p>>8;
				base += u <= v ? w : 0; v = u;
				debug("pos(%ld, %ld)", (-fr) ^ (base + u), ((-fr) ^ (base + u)) - ((-fr) & mii->mi.k));
				kv_push(mm_mini_t, bkt[h & mask].w.a, ((mm_mini_t){
					.hrem = h>>b, .pos = base + u, .rid = (mii->svec.n<<1) + fr
				}));
			}
		}
		src++; mii->svec.n++;				/* update src and dst pointers (update rid) */
	}
	free(s->a.a);
	free(s);
	return;
}
static
void mm_idx_drain(uint32_t tid, void *arg, void *item)
{
	mm_idx_intl_t *mii = (mm_idx_intl_t *)arg;
	mm_idx_step_t *s = (mm_idx_step_t *)item;
	/* sorted pipeline (for debugging) */
	kv_hq_push(v4u32_t, incq_comp, mii->hq, ((v4u32_t){.u64 = {s->id, (uintptr_t)s}}));
	while(mii->hq.n > 1 && mii->hq.a[1].u64[0] == mii->ocnt) {
		mii->ocnt++;
		s = (mm_idx_step_t *)kv_hq_pop(v4u32_t, incq_comp, mii->hq).u64[1];
		mm_idx_drain_intl(mii, s);
	}
	return;
}

/**
 * @fn mm_idx_count_occ
 * @brief sort buckets and count #elements
 */
static
void *mm_idx_count_occ(uint32_t tid, void *arg, void *item)
{
	uint64_t i = (uint64_t)item;
	mm_idx_intl_t *mii = (mm_idx_intl_t *)arg;
	mm_idx_bkt_t *b = &mii->mi.bkt[(1ULL<<mii->mi.b) *  i      / mii->nth] - 1;
	mm_idx_bkt_t *t = &mii->mi.bkt[(1ULL<<mii->mi.b) * (i + 1) / mii->nth];

	uint32_v *cnt = &mii->cnt[tid];
	while(++b < t) {
		uint64_t n_arr = b->w.a.n;
		if(n_arr == 0) { _memset_blk_u(b, 0, sizeof(mm_idx_bkt_t)); continue; }

		/* sort by minimizer */
		mm_mini_t *arr = b->w.a.a;
		radix_sort_128x((v4u32_t *)arr, n_arr);

		kv_reserve(uint32_t, *cnt, cnt->n + n_arr);
		uint32_t *u = &cnt->a[cnt->n];

		/* iterate over minimizers to count #keys */
		uint64_t n_keys = 0, n_single = 0, n = 1, ph = arr[0].hrem;
		for(mm_mini_t *p = &arr[1], *t = &arr[n_arr]; p < t; ph = p++->hrem, n++) {
			if(ph != p->hrem) { n_single += n == 1; *u++ = n; n = 0; n_keys++; }
		}
		b->v.n.single = n_single + (n == 1);
		b->v.n.keys = n_keys + 1;
		// debug("single(%u), keys(%u)", b->v.n.single, b->v.n.keys);
		*u++ = n; cnt->n = u - cnt->a;
	}
	return(NULL);
}

/**
 * @fn mm_idx_build_hash
 */
static
void *mm_idx_build_hash(uint32_t tid, void *arg, void *item)
{
	uint64_t i = (uint64_t)item;
	mm_idx_intl_t *mii = (mm_idx_intl_t *)arg;
	mm_idx_bkt_t *b  = &mii->mi.bkt[(1ULL<<mii->mi.b) *  i      / mii->nth] - 1;
	mm_idx_bkt_t *bt = &mii->mi.bkt[(1ULL<<mii->mi.b) * (i + 1) / mii->nth];

	debug("i(%lu), (%llu, %llu)", i, (1ULL<<mii->mi.b) * i / mii->nth, (1ULL<<mii->mi.b) * (i+1) / mii->nth);
	debug("b(%p), t(%p)", b, bt);
	while(++b < bt) {
		uint64_t n_arr = b->w.a.n;
		if(n_arr == 0) { continue; }

		/* create the (2nd-stage) hash table for the bucket */
		mm_mini_t *arr = b->w.a.a;
		#define _fill_body() { \
			uint64_t key = q->hrem, val = _loadu_u64(&q->pos);	/* uint32_t pos, rid; */ \
			if(++q < p) { \
				r[++sp] = val; val = sp<<32 | 0x01ULL<<63 | 1;	/* swap val with a (base_idx, cnt) tuple */ \
				do { r[++sp] = _loadu_u64(&q->pos); val++; } while(++q < p); \
			} \
			kh_put(&b->w.h, key, val); \
		}
		// debug("i(%lu), p(%p), single(%u), keys(%u)", b - mii->mi.bkt, b->v.p, b->v.n.single, b->v.n.keys);
		kh_init_static(&b->w.h, 1.1 * b->v.n.keys / KH_THRESH);	/* make the hash table occupancy equal to (or slightly smaller than) KH_THRESH */
		uint64_t max_cnt = mii->mi.occ[mii->mi.n_occ - 1], sp = 0, *r = (uint64_t *)arr;	/* reuse minimizer array */
		mm_mini_t *p = arr, *q = p, *t = &arr[n_arr];
		for(uint64_t ph = p++->hrem; p < t; q = p, ph = p++->hrem) {
			if(ph != p->hrem && p - q <= max_cnt) { _fill_body(); }	/* skip if occurs more than the max_cnt threshold */
		}
		if(p - q <= max_cnt) { _fill_body(); }
		#undef _fill_body

		/* shrink table */
		r[0] = n_arr - b->v.n.single;							/* table size saved at p[0] (for use in index serialization) */
		b->v.p = realloc(r, sizeof(uint64_t) * (n_arr - b->v.n.single + 1));	/* shrink array (is this an overhead?) */
		// b->v.p = r;
	}
	return(NULL);
}

/**
 * @fn mm_idx_gen
 * @brief root function of the index construction pipeline
 */
static _force_inline
mm_idx_t *mm_idx_gen(mm_idx_params_t const *o, bseq_file_t *fp, pt_t *pt)
{
	uint8_t b = MIN2(o->k * 2, o->b);		/* clip bucket size */
	mm_idx_intl_t *mmi = calloc(1, sizeof(mm_idx_intl_t));
	*mmi = (mm_idx_intl_t){					/* init pipeline context */
		.mi = (mm_idx_t){
			.bkt = malloc(sizeof(mm_idx_bkt_t) * (1ULL<<b)),
			.mask = (1ULL<<b) - 1, .b = b, .w = o->w, .k = o->k, .n_occ = o->n_frq
		},
		.fp = fp, .nth = pt_nth(pt),
		.cnt   = calloc(pt_nth(pt), sizeof(uint32_v)),
		.circ  = &o->circ,
		.call  = kh_str_ptr(&o->circ) && kh_str_size(&o->circ) == 0,	/* mark all sequences as circular if array is instanciated but no entry found */
		.ctest = kh_str_ptr(&o->circ) && kh_str_size(&o->circ) > 0
	};

	/* read sequence and collect minimizers */
	kv_hq_init(mmi->hq);					/* initialize heapqueue for packet sorting */
	pt_stream(pt, mmi, mm_idx_source, mm_idx_worker, mm_idx_drain);
	kv_hq_destroy(mmi->hq);

	/* sort minimizers then concatenate occurrence arrays */
	pt_parallel(pt, mmi, mm_idx_count_occ);
	for(uint64_t i = 1; i < pt_nth(pt); i++) {
		kv_pushm(uint32_t, mmi->cnt[0], mmi->cnt[i].a, mmi->cnt[i].n);
		free(mmi->cnt[i].a);
	}
	debug("keys(%lu)", mmi->cnt[0].n);

	/* calculate occurrence thresholds */
	for(uint64_t i = 0; i < o->n_frq; i++) {
		mmi->mi.occ[i] = (o->frq[i] <= 0.0
			? UINT32_MAX
			: (ks_ksmall_uint32_t(mmi->cnt[0].n, mmi->cnt[0].a, (uint32_t)((1.0 - o->frq[i]) * mmi->cnt[0].n)) + 1)
		);
	}
	free(mmi->cnt[0].a);

	/* build hash table */
	pt_parallel(pt, mmi, mm_idx_build_hash);

	/* finish */
	mmi->mi.s = mmi->svec.a;
	mmi->mi.n_seq = mmi->svec.n;
	return((mm_idx_t *)mmi);
}

#if 0
/**
 * @fn mm_idx_cmp
 * @brief compare two index objects (for debugging)
 */
static _force_inline
void mm_idx_cmp(mm_opt_t const *opt, mm_idx_t const *m1, mm_idx_t const *m2)
{
	for(uint64_t i = 0; i < 1ULL<<opt->b; ++i) {
		mm_idx_bkt_t *bkt1 = &m1->bkt[i], *bkt2 = &m2->bkt[i];
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

// #define MM_IDX_MAGIC "MAI\x08"		/* minialign index version 8 */
#define MM_IDX_MAGIC	0x0849414d		/* "MAI\x08" in little endian; minialign index version 8 */

/**
 * @fn mm_idx_dump
 * @brief dump index to fp, index is broken after the function call.
 */
#define _up(_x)						( (uintptr_t)(_x) )
#define _inside_ptr(_a, _b, _c)		( (_up(_b) - _up(_a)) < (_up(_c) - _up(_a)) )
static _force_inline
uint64_t mm_idx_dump_calc_size(mm_idx_t const *mi)
{
	uint64_t size = sizeof(mm_idx_t);
	size += sizeof(mm_idx_bkt_t) * (1ULL<<mi->b);
	size += sizeof(mm_idx_seq_t) * mi->n_seq;
	for(uint64_t i = 0; i < 1ULL<<mi->b; i++) {
		mm_idx_bkt_t *b = &mi->bkt[i];
		if(kh_ptr(&b->w.h) == NULL) { continue; }
		size += sizeof(v4u32_t) * kh_size(&b->w.h);
		size += sizeof(uint64_t) * (b->v.p ? *b->v.p + 1 : 0);
	}
	mm_idx_intl_t *mmi = (mm_idx_intl_t *)mi;
	for(uint64_t i = 0; i < mmi->mvec.n; i++) { size += mmi->mvec.a[i].size; }
	return(size);
}
static _force_inline
void mm_idx_dump(mm_idx_t const *mi, void *fp, write_t wfp)
{
	#define _writep(_b, _l)		{ wfp(fp, _b, _l); }
	#define _writea(type, _a)	{ type _t = (_a); _writep(&(_t), sizeof(type)); }

	/* calc size */
	uint64_t size = mm_idx_dump_calc_size(mi);

	/* dump header */
	_writea(uint32_t, MM_IDX_MAGIC); _writea(uint64_t, size);

	/* accumulate offset */
	#define _acc(_bytes)	({ uintptr_t _s = ofs; ofs += (ptrdiff_t)(_bytes); (void *)_s; })
	#define _ofs(_base)		( (ptrdiff_t)((_base) - ofs) )
	uintptr_t ofs = sizeof(mm_idx_t);

	/* dump idx object */
	mm_idx_t mib = *mi;
	mib.bkt = _acc(sizeof(mm_idx_bkt_t) * (1ULL<<mi->b));
	mib.s = _acc(sizeof(mm_idx_seq_t) * mi->n_seq);
	_writea(mm_idx_t, mib);

	/* dump buckets (= first-stage hash table) */
	for(uint64_t i = 0; i < 1ULL<<mi->b; i++) {
		mm_idx_bkt_t b = mi->bkt[i];
		if(kh_ptr(&b.w.h)) {
			b.w.h.a = _acc(sizeof(v4u32_t) * kh_size(&mi->bkt[i].w.h));
			b.v.p = _acc(sizeof(uint64_t) * (*b.v.p + 1));
		}
		_writea(mm_idx_bkt_t, b);
	}

	/* dump sequences */
	mm_idx_intl_t *mmi = (mm_idx_intl_t *)mi;
	mm_idx_mem_t const *q = mmi->mvec.a;
	for(mm_idx_seq_t *p = mi->s, *t = &mi->s[mi->n_seq]; p < t; p++) {
		/* update offset to forward to the next memory block if the sequence does not reside in the current one */
		if(!_inside_ptr(q->base, p->seq, q->base + q->size)) { ofs += q++->size; }
		mm_idx_seq_t s = *p; s.seq -= _ofs(q->base); s.name -= _ofs(q->base);
		_writea(mm_idx_seq_t, s);
	}

	/* dump memory blocks */
	for(mm_idx_bkt_t *b = mi->bkt, *t = &mi->bkt[1ULL<<mi->b]; b < t; b++) {
		if(kh_ptr(&b->w.h) == NULL) { continue; }
		uint64_t const *p = b->v.p ? b->v.p : ((uint64_t const [1]){ 0 });
		_writep(b->w.h.a, sizeof(v4u32_t) * kh_size(&b->w.h));
		_writep(p, sizeof(uint64_t) * (*p + 1));		/* value table, size and content */
	}
	for(mm_idx_mem_t const *p = mmi->mvec.a, *t = &mmi->mvec.a[mmi->mvec.n]; p < t; p++) {
		_writep(p->base, sizeof(uint8_t) * p->size);
	}
	return;
	#undef _writep
	#undef _writea
	#undef _acc
	#undef _ofs
}
#undef _up
#undef _inside_ptr

/**
 * @fn mm_idx_load
 * @brief create index object from file stream
 */
static _force_inline
mm_idx_t *mm_idx_load(void *fp, read_t rfp)
{
	/* read by _l and test if full length is filled, jump to _fail if not */
	#define _readp(_b, _l)	{ if(rfp(fp, _b, _l) != _l) { goto _mm_idx_load_fail; } }
	#define _reada(type)	({ type _n; _readp(&_n, sizeof(type)); _n; })

	mm_idx_t *mi = NULL;
	if(_reada(uint32_t) != MM_IDX_MAGIC) { fprintf(stderr, "magic not correct\n"); goto _mm_idx_load_fail; }
	uint64_t size = _reada(uint64_t);		/* read index size */
	mi = malloc(size); _readp(mi, size);	/* malloc mem and read all FIXME: can be mapped to hugepages to improve performance */
	mi->mono = 1;

	/* restore pointers from offsets */
	#define _rst(_p, _b)	{ (_p) = (void *)((uintptr_t)(_b) + (ptrdiff_t)(_p)); }
	_rst(mi->bkt, mi); _rst(mi->s, mi);		/* keep mi->mem untouched */
	for(mm_idx_bkt_t *b = mi->bkt, *t = &mi->bkt[1ULL<<mi->b]; b < t; b++) {
		if(kh_ptr(&b->w.h) == NULL) { continue; }
		_rst(b->w.h.a, mi); _rst(b->v.p, mi);
	}
	for(uint64_t i = 0; i < mi->n_seq; i++) {
		_rst(mi->s[i].name, mi);
		_rst(mi->s[i].seq, mi);
	}
	#undef _rst
	return(mi);
_mm_idx_load_fail:
	free(mi);
	return(NULL);

	#undef _readp
	#undef _reada
}

/* end of index.c */

/* map.c */
/**
 * @struct mm_resc_t
 * @brief rescued seed (pos, occ) container, alias to v4u32_t
 */
typedef struct {
	uint32_t qs, n;
	v2u32_t const *p;
} mm_resc_t;
typedef struct { size_t n, m; mm_resc_t *a; } mm_resc_v;
_static_assert(sizeof(mm_resc_t) == sizeof(v4u32_t));

/**
 * @struct mm_seed_t, mm_leaf_t
 * @brief seed container / leaf seed container, alias to v4u32_t
 */
typedef struct {
	uint32_t upos, rid, vpos, cid;	/* 2x - y, ref id, x - 2y, previous seed index */
} mm_seed_t;
typedef struct {
	uint32_t rsid, rid, qid, lsid;	/* flag == 0 */
} mm_leaf_t;
typedef struct { size_t n, m; mm_seed_t *a; } mm_seed_v;
_static_assert(sizeof(mm_seed_t) == sizeof(v4u32_t));
_static_assert(sizeof(mm_leaf_t) == sizeof(mm_seed_t));
_static_assert(offsetof(mm_seed_t, cid) == offsetof(mm_leaf_t, lsid));
#define CID_INIT				( -1 )

/**
 * @struct mm_root_t
 * @brief chain head container
 */
typedef struct {
	uint32_t plen, lid;				/* the longest root->leaf path length, root seed index */
} mm_root_t;
typedef struct { size_t n, m; mm_root_t *a; } mm_root_v;
_static_assert(sizeof(mm_root_t) == sizeof(v2u32_t));

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

typedef struct {
	uint32_t aid, mapq;				/* alignment id and mapq */
	gaba_alignment_t const a[];
} mm_aln_t;
struct mm_reg_s {
	uint32_t n_all, n_uniq;
	mm_aln_t const *aln[];
};

#if 0
/**
 * @struct mm_alnset_t
 */
typedef struct {
	uint32_t rid, mapq;				/* reference id and mapq */
	uint32_t n_aln, b_aln;			/* #alns, base index of reg->bin */
} mm_alnset_t;

/* output */
typedef struct mm_aln_s {
	uint32_t nid;
	uint16_t n_aln;
	uint8_t mapq, flag;				/* 0x80 for leaf */
	// uint32_t rid, qid;
	// uint32_t rspos, qspos, repos, qepos;
	struct mm_aln_s const *a;				/* nested structure */
} mm_aln_t;
typedef struct {
	uint32_t nid;
	uint16_t n_aln;
	uint8_t mapq, flag;				/* 0x80 for leaf */
	// uint32_t rid, qid;
	// uint32_t rspos, qspos, repos, qepos;
	gaba_alignment_t const a[1];	/* leaf node contains a gaba_alignment_t object */
} mm_aln_leaf_t;
#endif

/**
 * @struct mm_tbuf_s
 * @brief thread-local context (working buffers)
 */
typedef struct mm_tbuf_s {
	mm_idx_t mi;					/* (immutable) index object */
	int32_t twlen, tglen;			/* leaned chainable window edge length, linkable gap length; converted from wlen and glen */
	float min_ratio;
	uint32_t min_score;
	gaba_dp_t *dp;
	gaba_alloc_t alloc;				/* lmm contained */

	/* initialized in mm_init_query */
	uint32_t rid, qid;				/* raw qid (used in ava filter; ava switching flag is embedded in qid sign bit) */

	/* query and reference sections (initialized in mm_init_query) */
	uint8_t tail[96];				/* zeros or 0x80s which does not match to any bases */
	gaba_section_t r[2];			/* [0]: rf, [1]: rr */
	gaba_section_t q[3];			/* [0]: qf, [1]: qr, [2]: qf */
	gaba_section_t t[1];			/* tail section info */
	gaba_section_t *rtp, *qtp;		/* tail section pointers, &r[0] or &t[0] */

	/* working buffers */
	mm_resc_v resc;					/* rescued seed array */
	mm_resc_t *presc;				/* rescued array pointer */
	mm_seed_v seed;					/* seed array / leaf array */
	uint64_t n_seed;
	mm_root_v root;					/* roots of chain trees */
	uint32_v margin;				/* marginal roots */
	uint32_t n_res;					/* #alignments collected */
	ptr_v bin;						/* gaba_alignment_t* array */
	kh_t pos;						/* alignment dedup hash */
} mm_tbuf_t;

/**
 * @macro _s, _m
 * @brief extract sign, calc inv
 */
#define _sign(_x)		( (_x) < 0 ? -1 : 1 )
#define _smask(_x)		( ((int32_t)(_x))>>31 )
#define _umask(_x)		( 0x7fffffff & (_x) )
#define _mabs(_x)		( _smask(_x)^(_x) )
#define _ofs(_x)		( (int32_t)0x40000000 - (int32_t)(_x) )
// #define _u(x, y)		( (((x)<<1) | 0x40000000) - (y) )
// #define _v(x, y)		( ((x) | 0x40000000) - ((y)<<1) )

#define _u(_x, _y)		( _ud(_x, _y) + _ofs(0) )
#define _v(_x, _y)		( _vd(_x, _y) + _ofs(0) )
#define _ud(_x, _y)		( ((_x)<<1) - (_y) )
#define _vd(_x, _y)		( (_x) - ((_y)<<1) )

#define _bare(_x)		( (_x) - _ofs(0) )
#define _as(_ptr)		( ((_bare((_ptr)->upos)<<1) - _bare((_ptr)->vpos)) / 3 )
#define _bs(_ptr)		( (_bare((_ptr)->upos) - (_bare((_ptr)->vpos)<<1)) / 3 )
#define _ps(_ptr)		( (_ptr)->upos - (_ptr)->vpos )
#define _qs(_ptr)		( ((_ptr)->upos + (_ptr)->vpos) / 3 )

#define _p(_ptr)		( ((int32_t const *)(_ptr))[0] + ((int32_t const *)(_ptr))[1] )
#define _q(_ptr)		( ((int32_t const *)(_ptr))[1] - ((int32_t const *)(_ptr))[0] )
#define _inside(_lb, _x, _ub)		( (uint32_t)((_x) - (_lb)) <= (uint32_t)((_ub) - (_lb)) )
#define _key(_x, _y)	( (_x) ^ ((_x)>>29) ^ (_y) ^ _swap_u64(_y) )
_static_assert(_ud(1000, -5) == -_vd(-5, 1000));	/* invariant condition */
_static_assert(_ud(1000, 10) == -_vd(10, 1000));	/* invariant condition */

/**
 * @fn mm_expand
 * @brief expand minimizer to coef array
 */
static _force_inline
void mm_expand(
	mm_tbuf_t *self,
	uint32_t const n,
	v2u32_t const *r,							/* source array */
	uint32_t const qs)							/* query position */
{
	if(n == 0) { return; }

	/* expand array if needed */
	kv_reserve(mm_seed_t, self->seed, self->seed.n + n);

	/* iterate over all the collected minimizers */
	for(uint64_t i = 0; i < n; i++) {
		uint32_t const rid = r[i].u32[1];
		if(rid < self->qid) { continue; }		/* all-versus-all flag, base_rid, and base_qid are embedded in qid; skip if seed is in the lower triangle (all-versus-all) */
		uint32_t const rs = r[i].u32[0];		/* load reference pos */
		uint32_t const rmask = -(rid & 0x01);
		uint32_t const _rs = rs + (self->mi.k & rmask), _qs = qs ^ rmask;
		self->seed.a[self->seed.n++] = (mm_seed_t){
			.upos = _u(_rs, _qs),				/* coefficients */
			.vpos = _v(_rs, _qs),
			.rid = rid>>1, .cid = CID_INIT		/* first prev node is initialized with zero */
		};
		debug("abpos(%d, %d), pos(%d, %d), uvpos(%d, %d)", rs, qs, _rs, _qs, _bare(_u(_rs, _qs)), _bare(_v(_rs, _qs)));
	}
	return;
}

/**
 * @fn mm_collect_seed
 * @brief collect minimizers for the query seq
 * (note: branch misprediction penalty and memory access pattern should be alleviated, but how?)
 */
// static _force_inline
void mm_collect_seed(
	mm_tbuf_t *self)
{
	/* clear array and gather minimizers */
	self->seed.n = 0;

	mm_sketch_t sk;
	mm_sketch_init(&sk, self->mi.w, self->mi.k, (uint64_v *)&self->root);
	mm_sketch(&sk, self->q[0].base, self->q[0].len);
	debug("collected seeds, n(%zu)", self->root.n);

	/* prepare rescue array */
	kv_reserve(mm_resc_t, self->resc, self->root.n);
	mm_resc_t *s = self->resc.a;

	uint32_t const max_occ = self->mi.occ[self->mi.n_occ - 1];
	uint32_t const resc_occ = self->mi.occ[0];

	/* iterate over all the collected minimizers (seeds) on the query */
	uint64_t w = self->mi.w, pos = 0, v = w;
	for(uint64_t *p = (uint64_t *)self->root.a; (*p>>8) != 0xffffffffffffff; p++) {
		/* first calculate pos for the current one */
		uint64_t u = *p & 0x7f, fr = (*p>>7) & 0x01, h = *p>>8;
		pos += u - v + (u <= v ? w : 0); v = u;
		// debug("fpos(%ld), rpos(%ld), pos(%d), d(%lu)", fpos, rpos, pos, d);

		/* get minimizer matched on the ref at the current query pos */
		uint32_t n;
		v2u32_t const *r = mm_idx_get(&self->mi, h, &n);
		if(n > max_occ) { continue; }			/* skip if exceeds repetitive threshold */
		uint32_t const cpos = (pos + (self->mi.k & -fr)) ^ -fr;
		if(n > resc_occ) {						/* save if less than max but exceeds current threshold */
			*s++ = (mm_resc_t){ .p = r, .qs = cpos, .n = n };
			continue;
		}
		mm_expand(self, n, r, cpos);			/* append to seed array */
	};
	self->resc.n = s - self->resc.a;			/* write back rescued array */
	self->presc = self->resc.a;					/* init resc pointer */
	self->root.n = 0;							/* clear root array */
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
	debug("seed iteration i(%lu)", cnt);
	if(cnt == 0) {
		/* push head sentinel */
		kv_push(mm_seed_t, self->seed, ((mm_seed_t){ .rid = -1, .cid = CID_INIT }));
		mm_collect_seed(self);					/* first collect seeds */

		/* push tail sentinel */
		kv_push(mm_seed_t, self->seed, ((mm_seed_t){ .rid = INT32_MAX, .cid = CID_INIT }));
	} else {
		/* sort rescued array by occurrence */
		debug("sort resc, n(%zu)", self->resc.n);
		if(cnt == 1) { radix_sort_128x((v4u32_t *)self->resc.a, self->resc.n); }

		/* append seeds to seed bin */
		self->seed.n = self->n_seed;			/* clear leaf array */
		mm_resc_t *p = self->presc, *t = &self->resc.a[self->resc.n];
		while(p < t && p->n <= self->mi.occ[cnt]) {
			mm_expand(self, p->n, p->p, p->qs); p++;
		}
		self->presc = p;						/* write back resc pointer */
	}
	self->n_seed = self->seed.n;
	if(self->seed.n == 0) { return(0); }		/* seed not found for this round */

	/* sort seed array */
	debug("sort seed, n(%zu)", self->seed.n);
	radix_sort_128x((v4u32_t *)self->seed.a, self->seed.n);
	for(uint64_t i = 0; i < self->seed.n; i++) {
		debug("rid(%u), upos(%d), vpos(%d)", self->seed.a[i].rid, _bare(self->seed.a[i].upos), _bare(self->seed.a[i].vpos));
	}
	return(self->seed.n);						/* report #seeds found */
}

/**
 * @fn mm_chain_seeds
 */
// static _force_inline
uint64_t mm_chain_seeds(
	mm_tbuf_t *self)
{
	#define _l(_s)			( (mm_leaf_t *)(_s) )
	mm_seed_t *s = self->seed.a;
	mm_root_t *c = self->root.a;
	uint64_t ncid = 0, nlid = self->n_seed, r = 0, nr = 0;
	v4i32_t const tw = _seta_v4i32(-self->twlen, 0, 0, self->twlen);
	v4i32_t const bm = _seta_v4i32(INT32_MIN, INT32_MIN, INT32_MAX, INT32_MIN);
	do {
		#define _ld(_p)			( _shuf_v4i32(_load_v4i32(&s[_p]), _km_v4i32(2, 2, 1, 0)) )
		#define _bw(_x)			( _add_v4i32(_x, tw) )
		#define _vin(_x, _y)	( _mask_v4i32(_gt_v4i32(_x, _y)) )
		#define _vbare(_x)		( _sub_v4i32(_x, _seta_v4i32(_ofs(0), _ofs(0), 0, _ofs(0))) )
		uint64_t cid = ncid, lid = nlid++, rsid = nr; r = nr; nr = UINT64_MAX;
		v4i32_t const base = _min_v4i32(_ld(r), bm);

		/* test branching path */ {
			v4i32_t i = _ld(r);
			_print_v4i32(_vbare(i));
			for(uint64_t q = r, m; ((m = _vin(i, _bw(_ld(q - 1)))) & 0xff) == 0; q--) {/* check if the head can be a branching path */
				_print_v4i32(_vbare(_bw(_ld(q - 1))));
				debug("m(%lx)", m);
				if(m == 0xf000) {
					debug("hit at the middle of chain, cid(%u)", s[q - 1].cid);
					cid = s[q - 1].cid; rsid = _l(s)[c[cid].lid].rsid;
					// _l(s)[lid].rsid = _l(s)[s[q - 1].prev].rsid; s[r].prev = q - 1;
					goto _link_branch_break;
				}
			}

			/* open chain bin, load initial coordinates */
			c[ncid++] = (mm_root_t){ .lid = lid, .plen = _ofs(0) };					/* long enough because used to collect minimizers */
			s[r].cid = ~lid; /*_l(s)[lid].rsid = r;*/ 								/* mark root */

			/* record marginal seed for circularization */
			debug("test margin, as(%d), circular(%u)", _as(&s[r]), self->mi.s[s[r].rid].circular);
			if(_as(&s[r]) < (self->mi.s[s[r].rid].circular ? self->twlen : 0)) {	/* mark marginal if the reference-side seq is circular and the root is in the margin */
				debug("mark marginal");
				kv_push(uint32_t, self->margin, cid);
			}
		}

	_link_branch_break:;
		/* chain loop */
		uint64_t p = r, pi, lsid, m, acc = 0;
		do {
			v4i32_t b = _max_v4i32(base, _bw(_ld(p)));
			_print_v4i32(_vbare(_bw(_ld(p))));
			_print_v4i32(_vbare(b));
			uint64_t q = p + 1;
			do {
				m = _vin(_ld(q), b);
				_print_v4i32(_vbare(_ld(q)));
				_print_v4i32(_gt_v4i32(_ld(q), b));

				acc += 2 * (~(acc + 1) & 0x04);
				debug("m(%lx), acc(%ld, %lx), %s", m, acc>>3, acc&0x03, m == 0xf000 ? "chained" : "not chained");
				if(m == 0xf000) { lsid = q; s[q].cid = cid; acc |= 2; }
				acc |= m != 0xf000; q++;
			} while((m & 0xff) == 0);

			/*
			for(uint64_t q = p + 1; ((m = _vin(_ld(q), b)) & 0xff) == 0; q++) {
				_print_v4i32(_vbare(_ld(q)));
				_print_v4i32(_gt_v4i32(_ld(q), b));
				acc |= m != 0xf000; acc += 2 * (~(acc + 1) & 0x04);
				debug("m(%lx), acc(%ld, %lx), %s", m, acc>>3, acc&0x03, m == 0xf000 ? "chained" : "not chained");
				if(m != 0xf000) { continue; }
				lsid = q; s[q].cid = cid; acc |= 2;
			}
			*/
			pi = s[p + 1].cid != CID_INIT ? 1 : acc>>3; p += pi;
		} while((m & 0xf0) == 0);
		debug("r(%lu), pi(%lu), c(%lu), inc(%lu), nr(%lu, %lu)",
			r, pi, acc>>3, pi == 1 ? acc>>3 : 1, nr, MIN2(nr, r + (pi == 1 ? acc>>3 : 1)));
		nr = MIN2(nr, r + (pi == 1 ? acc>>3 : 1));

		/* open leaf bin */
		// uint32_t rsid = MIN2(r, _l(s)[c[cid].lid].rsid);
		_l(s)[lid] = (mm_leaf_t){				/* record leaf */
			.rsid = rsid, .lsid = lsid, .rid = s[rsid].rid, .qid = self->qid,			/* qid = self->qid */
		};
		debug("record leaf, sid(%u, %u), rid(%u), qid(%u)", _l(s)[lid].rsid, _l(s)[lid].lsid, _l(s)[lid].rid, _l(s)[lid].qid);

		/* update root (chain) */
		uint32_t plen = _ofs(_ps(&s[lsid]) - _ps(&s[rsid]));
		debug("test plen, cid(%lu), plen(%u, %u)", cid, _ofs(c[cid].plen), _ofs(plen));
		if(c[cid].plen > plen) {				/* choose the longer */
			c[cid] = (mm_root_t){ .lid = lid, .plen = plen };
			debug("update plen, lid(%lu), rsid(%u), lsid(%u), rpos(%d, %d), lpos(%d, %d), plen(%u)",
				lid, _l(s)[lid].rsid, _l(s)[lid].lsid,
				_bare(s[_l(s)[lid].rsid].upos), _bare(s[_l(s)[lid].rsid].vpos),
				_bare(s[_l(s)[lid].lsid].upos), _bare(s[_l(s)[lid].lsid].vpos), _ofs(plen));
		}
		#undef _ld
		#undef _bw
		#undef _vin
	} while(s[nr].rid < INT32_MAX);				/* check if the next root is found */
	self->seed.n = nlid;						/* write back tail index */
	self->root.n = ncid;						/* write back root index */
	debug("lid(%lu), cid(%lu)", nlid, ncid);
	return(ncid);								/* #collected chains */
}

/**
 * @fn mm_circularize
 * @brief sweep tails of the chains and link to a head if found
 */
static _force_inline
void mm_circularize(
	mm_tbuf_t *self)
{
	debug("circularize");
	mm_root_t *c = self->root.a;
	mm_seed_t *s = self->seed.a;
	mm_leaf_t const *l = &_l(s)[self->n_seed], *lt = &_l(s)[self->seed.n];
	debug("i(%lu), t(%lu)", l - (mm_leaf_t *)self->seed.a, lt - (mm_leaf_t *)self->seed.a);
	for(uint32_t const *m = self->margin.a, *mt = &self->margin.a[self->margin.n]; m < mt; m++) {
		debug("i(%ld), m(%u), l(%lu)", m - self->margin.a, *m, l - (mm_leaf_t *)self->seed.a);

		mm_leaf_t const *p = &_l(s)[c[*m].lid];
		if(p < l) { continue; }

		uint32_t const rid = s[p->rsid].rid;
		debug("lid(%u), rsid(%u), lsid(%u), rid(%u, %u)", c[*m].lid, p->rsid, p->lsid, s[p->rsid].rid, p->rid);

		uint32_t const uofs = _ud(self->mi.s[rid].l_seq, 0), vofs = _vd(self->mi.s[rid].l_seq, 0);
		uint32_t const ulb = s[p->rsid].upos - self->twlen, uub = s[p->rsid].upos;
		uint32_t const vlb = s[p->rsid].vpos, vub = s[p->rsid].vpos + self->twlen;
		debug("ulb(%d), uub(%d), vlb(%d), vub(%d)", _bare(ulb), _bare(uub), _bare(vlb), _bare(vub));

		/* forward lower boundary */
		while(l < lt && s[l->rsid].rid < rid) {
			debug("i(%lu), rid(%u, %u)", l - (mm_leaf_t *)self->seed.a, s[l->rsid].rid, rid);
			l++;
		}
		while(l < lt && s[l->rsid].vpos - vofs < vlb) {
			debug("i(%lu), vpos(%d, %d)", l - (mm_leaf_t *)self->seed.a, _bare(s[l->rsid].vpos - vofs), _bare(vlb));
			l++;
		}

		/* link */
		for(mm_leaf_t const *q = l; q < lt && s[q->lsid].upos - uofs < s[p->rsid].upos; q++) {
			if(s[q->lsid].upos - uofs < ulb) { continue; }						/* FIXME: squash branches for better performance */
			if(s[q->lsid].upos - uofs > uub) { continue; }
			if(s[q->lsid].vpos - vofs > vub) { continue; }
			debug("link p(%lu) to q(%lu), pos(%d, %d)", p - (mm_leaf_t *)self->seed.a, q - (mm_leaf_t *)self->seed.a, _bare(s[q->lsid].upos), _bare(s[q->lsid].vpos));

			/* circularize */
			uint32_t rcid = s[q->rsid].cid;
			debug("rcid(%d), rsid(%u), lsid(%u)", rcid, q->rsid, q->lsid);
			while((int32_t)rcid >= 0) {
				debug("rcid(%d), rsid(%u), lsid(%u)", s[_l(s)[rcid].rsid].cid, _l(s)[rcid].rsid, _l(s)[rcid].lsid);
				rcid = s[_l(s)[rcid].rsid].cid;
			} rcid = ~rcid;
			debug("found root, rcid(%d)", rcid);

			int32_t pj = _ps(&s[p->rsid]) - _ps(&s[q->lsid]) + self->mi.s[rid].l_seq;
			s[p->rsid].cid = c[rcid].lid;		/* create leaf edge */
			c[rcid].plen -= _ofs(c[*m].plen) + pj;
			c[rcid].lid = c[*m].lid;
			c[*m].plen = _ofs(0);				/* mark invalid */
		}
	}
	return;
}

/**
 * @fn mm_chain
 * @brief collect chain in self->root, self->seed must be sorted, self->root and self->map are not cleared at the head.
 */
static _force_inline
uint64_t mm_chain(
	mm_tbuf_t *self,
	uint64_t cnt)
{
	/* reserve space for map and chain */
	kv_reserve(mm_seed_t, self->seed, self->seed.n + self->seed.n);
	kv_reserve(mm_root_t, self->root, self->root.n + self->seed.n);
	kv_reserve(uint32_t, self->margin, self->margin.n + self->seed.n);

	/* chain */
	if(mm_chain_seeds(self) == 0) { return(0); }
	mm_circularize(self);

	/* sort chain by length */
	radix_sort_64x((v2u32_t *)self->root.a, self->root.n);
	debug("sorted seeds, n(%lu)", self->root.n);
	return(self->root.n);
}

/**
 * @fn mm_next
 * @brief pick up next seed in the upward region of the current chain
 */
static _force_inline
gaba_pos_pair_t mm_next(
	mm_tbuf_t *self,
	gaba_pos_pair_t cp,
	uint32_t lid)					/* leaf id */
{
	mm_seed_t const *s = self->seed.a;
	uint32_t const ulb = _u(cp.apos, cp.bpos) - self->twlen, vub = _v(cp.apos, cp.bpos) + self->twlen;
	do {
		/* lid points at a leaf; FIXME: qid test is incomplete */
		while(_l(s)[lid].qid != 0 && s[_l(s)[lid].lsid].rid != self->rid) {
			lid = s[_l(s)[lid].rsid].cid;
		}
		/* lid points at a seed, vpos overlaps with flag and vpos never be zero for seeds; FIXME: _l(s)[lid].flag == 0 */
		do { lid = s[lid].cid; } while(s[lid].vpos > vub || s[lid].upos < ulb);
	} while((int32_t)_l(s)[lid].rid >= 0);	/* negative for lid in seed_t */

	return((gaba_pos_pair_t){
		.apos = _as(&s[lid]), .bpos = _bs(&s[lid])
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
	uint64_t hkey = _key(_loadu_u64(&a->seg->apos), id);
	uint64_t tkey = _key(_loadu_u64(&a->seg->apos) + _loadu_u64(&a->seg->alen), id);
	debug("hash key, h(%lx), t(%lx)", hkey, tkey);

	v2u32_t *h = (v2u32_t *)kh_put_ptr(&self->pos, hkey, 1);
	v2u32_t *t = (v2u32_t *)kh_put_ptr(&self->pos, tkey, 0);		/* noextend */
	uint64_t new = h->u64[0] == KH_INIT_VAL;

	/* open new bin for aln, reuse if the head hit an existing one */
	uint32_t bid = new ? kv_push(void *, self->bin, (void *)a) : h->u32[1];
	gaba_alignment_t const **b = (gaba_alignment_t const **)&self->bin.a[bid];

	debug("id(%u, %u)", a->seg->aid, a->seg->bid);
	debug("record h(%lx, %u), t(%lx, %u), bid(%u)", h[-1].u64[0], h->u32[1], t[-1].u64[0], t->u32[1], bid);

	/* update res */
	res->score -= a->score;

	/* update bin */
	mm_bin_t *bin = (mm_bin_t *)&self->bin.a[res->bid];
	bin->n_aln += h->u64[0] == KH_INIT_VAL;
	bin->plen += _p(&a->seg->alen);
	bin->lb = MIN2(bin->lb, a->seg->bpos);
	bin->ub = MAX2(bin->ub, a->seg->bpos + a->seg->blen);
	debug("update bounds, bid(%u), lb(%u), ub(%u)", res->bid, bin->lb, bin->ub);

	debug("record alignment, (%u, %u) -> (%u, %u), n_aln(%u)", a->seg->apos, a->seg->bpos, a->seg->apos + a->seg->alen, a->seg->bpos + a->seg->blen, bin->n_aln);

	/* update hash */
	if(b[0]->score > a->score) {
		debug("discard, a(%ld), b(%ld)", a->score, b[0]->score);
		*t = (v2u32_t){ .u64[0] = KH_INIT_VAL };	/* re-mark: evaluated but not found */
	} else {
		debug("replace, a(%ld), b(%ld)", a->score, b[0]->score);
		/* replace if the new one is larger than the old one */
		*b = a;
		*h = *t = (v2u32_t){
			.u32 = {
				[0] = res - (mm_res_t *)self->root.a,
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
	gaba_fill_t const *f = gaba_dp_fill_root(dp, a, s.apos, b, s.bpos, 0);
	if(_unlikely(f == NULL)) { goto _mm_extend_core_abort; }
	debug("fill root, max(%ld), status(%x)", f->max, f->stat);

	gaba_fill_t const *m = f;					/* record max */
	uint32_t flag = GABA_TERM;
	while((flag & f->stat) == 0) {
		/* update section if reached tail */
		if(f->stat & GABA_UPDATE_A) { a = at; }	/* cmov */
		if(f->stat & GABA_UPDATE_B) { b = bt; }	/* cmov */

		/* update flag for tail detection */
		flag |= f->stat & (GABA_UPDATE_A | GABA_UPDATE_B);

		/* fill the next section */
		if(_unlikely((f = gaba_dp_fill(dp, f, a, b, 0)) == NULL)) {
			goto _mm_extend_core_abort;
		}
		debug("fill, max(%ld, %ld), status(%x)", f->max, m->max, f->stat);
		m = f->max > m->max ? f : m;
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
void mm_extend_load_ref(
	mm_tbuf_t *self,
	uint32_t rid)
{
	/* load ref */
	mm_idx_seq_t const *ref = &self->mi.s[rid];
	self->rid = rid;
	self->r[0] = _sec_fw(rid, ref->seq, ref->l_seq);
	self->r[1] = _sec_rv(rid, ref->seq, ref->l_seq);
	self->rtp = ref->circular ? self->r : self->t;
	return;
}
static _force_inline
uint64_t mm_extend(
	mm_tbuf_t *self,
	uint64_t i)
{
	mm_seed_t const *s = self->seed.a;

	/* loop: evaluate chain */
	for(uint64_t k = 0; k < self->root.n; k++) {
		/* load seed ptr and position */
		uint32_t const lid = self->root.a[k].lid;
		int32_t plen = _ofs(self->root.a[k].plen);
		if(plen /* self->m*/ < self->min_score) { continue; }

		/* load seed positions and sequence ids */
		mm_seed_t const *p = &s[_l(s)[lid].lsid];/* load tail */
		gaba_pos_pair_t cp = { .apos = _as(p), .bpos = _bs(p) };
		debug("k(%lu), lid(%u), id(%u, %u), sid(%u, %u), pos(%d, %d), bare(%d, %d), cp(%d, %d)",
			k, lid, _l(s)[lid].rid, _l(s)[lid].qid, _l(s)[lid].rsid, _l(s)[lid].lsid,
			p->upos, p->vpos, _bare(p->upos), _bare(p->vpos),
			cp.apos, cp.bpos);

		int32_t rev = _smask(cp.bpos);
		mm_extend_load_ref(self, _l(s)[lid].rid);
		// debug("chain(%lu), cid(%u), rev(%d), plen(%d), (%u, %u), (%d, %d), (%d, %d)", k, self->map.a[p->mid], rev, plen, self->rid, self->qid, p->apos, p->bpos, cp.apos, cp.bpos);

		/* open result bin */
		uint32_t bid = kv_pushm(void *, self->bin, (void **)&((mm_bin_t){ .lb = UINT32_MAX }), MM_BIN_N);

		/* reuse seed bin for score accumulation */
		mm_res_t *res = (mm_res_t *)&self->root.a[self->n_res++];
		*res = (mm_res_t){
			.score = _ofs(0),					/* score (negated and offsetted) */
			.bid = bid							/* bin base index */
		};

		/* loop: issue extension until whole chain is covered by alignments */
		debug("plen(%d), ppos(%d), p(%u, %d)", plen, _p(&cp) - plen, _p(&cp), _p(&cp) - plen);
		for(uint32_t ppos = _p(&cp) - plen, rem = 3, narrow = 0;
			rem > 0 && _p(&cp) >= ppos;
			cp = --rem > 0 ? mm_next(self, cp, lid) : cp
		) {
			/* reset stack */
			gaba_dp_flush(self->dp);

			/* init tail as zero padded arrays */
			gaba_fill_t const *t = (gaba_fill_t const *)self->tail, *u = t;

			/* downward extension */
			gaba_pos_pair_t tp = cp;
			u = mm_extend_core(&self->dp[narrow], &self->r[0], self->rtp+1, &self->q[-rev], self->qtp-rev, tp);

			/* search pos if extended */
			if(u->max > 0) { tp = gaba_dp_search_max(&self->dp[narrow], u); }
			debug("len(%u, %u), score(%ld), (%u, %u) -> (%u, %u)",
				self->r[0].len, self->q[0].len, u->max, cp.apos, cp.bpos, tp.apos, tp.bpos);

			/* skip if tail is duplicated */
			if(mm_test_pos(self, res, tp) >= 0) {
				narrow = 0;
				debug("duplication detected, try narrower(%u)", narrow);
				continue;			/* try narrower band in the next itr to avoid collision */
			}

			/* upward extension: coordinate reversed here */
			t = mm_extend_core(&self->dp[0], &self->r[1], self->rtp+1, &self->q[1+rev], self->qtp+1+rev,
				((gaba_pos_pair_t){
					.apos = self->r[0].len - tp.apos - 1,
					.bpos = self->q[0].len - tp.bpos - 1
				})
			);
			if(t->max < self->min_score) {	/* max == 0 indicates alignment was not found */
				mm_mark_pos(self, tp);			/* mark evaluated but not found */
				continue;
			}

			/* generate alignment: coordinates are reversed again, gaps are left-aligned in the resulting path */
			gaba_alignment_t const *a = gaba_dp_trace(&self->dp[0], t, &self->alloc);		/* lmm is contained in self->trace */
			if(a == NULL) {
				debug("failed to generate trace: len(%u, %u), score(%ld), <- (%u, %u)", self->r[0].len, self->q[0].len, t->max, tp.apos, tp.bpos);
				continue;						/* something is wrong... */
			}
			debug("len(%u, %u), score(%ld), (%u, %u) <- (%u, %u)",
				self->r[0].len, self->q[0].len, t->max, a->seg->apos, a->seg->bpos, tp.apos, tp.bpos);

			/* record alignment, update current head position */
			if(mm_record(self, res, a)) { rem = 10; }
			if(_p(&a->seg->apos) < ppos) {
				debug("full length captured, cp(%u, %u), p(%u), ppos(%u)", cp.apos, cp.bpos, _p(&a->seg->apos), ppos);
				break;							/* full length covered, break before calling mm_next */
			}

			/* update itr states */
			cp = *((gaba_pos_pair_t *)&a->seg->apos);
			narrow = 0;
			debug("split detected, try narrower(%u), cp(%u, %u), p(%u), ppos(%u), rem(%u)", narrow, cp.apos, cp.bpos, _p(&cp), ppos, rem);
		}

		/* discard if the score did not exceed the minimum threshold */
		if(res->score > _ofs(self->min_score)) {/* _ofs() inverts the sign */
			self->bin.n = bid;
			self->n_res--;
			debug("remove bin, bid(%zu), n_res(%u)", self->bin.n, self->n_res);
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
	mm_res_t *res = (mm_res_t *)self->root.a;	/* alignments, must be sorted */
	uint64_t q = self->n_res;
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
					debug("mark secondary, i(%lu)", i);
					q--; _swap_res(i, q); i--;	/* move to tail */
					goto _loop_tail;
				}
			}

			/* update score */
			max = MAX2(max, ((uint64_t)(2*(ub - lb) - span)<<32) | i);
		_loop_tail:;
		}
		if(max & 0xffffffff) {
			debug("mark supplementary, i(%lu)", max & 0xffffffff);
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
	mm_res_t *res = (mm_res_t *)self->root.a;	/* alignments, must be sorted */

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
		double pid = 1.0 - (double)(elen /* self->m*/ - score) / (double)(2/*self->m + self->x*/) / elen;

		/* estimate unique length */
		double ec = 2.0 / (pid * (double)(2/*self->m + self->x*/) - (double)1/*self->x*/);
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

/**
 * @fn mm_post_ava
 * @brief all-versus-all mapq estimation
 */
static _force_inline
uint64_t mm_post_ava(
	mm_tbuf_t *self)
{
	mm_res_t *res = (mm_res_t *)self->root.a;	/* alignments, must be sorted */

	uint32_t i, score, min = res[0].score * self->min_ratio;
	for(i = 0; i < self->n_res && (score = res[i].score) >= min; i++) {
		mm_bin_t *bin = (mm_bin_t *)&self->bin.a[res[i].bid];

		/* estimate effective length */
		double elen = (double)bin->plen / 2.0;
		double pid = 1.0 - (double)(elen /* self->m*/ - score) / (double)(2/*self->m + self->x*/) / elen;
		double ec = 2.0 / (pid * (double)(2/*self->m + self->x*/) - (double)1/*self->x*/);
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
	mm_reg_t *reg = lmm_malloc((lmm_t *)self->alloc.opaque, size);
	mm_aln_t **p = (mm_aln_t **)reg->aln, **b = p;

	/* build head */
	*reg = (mm_reg_t){ 0 };

	/* build reg array (copy from bin) */
	mm_res_t *res = (mm_res_t *)self->root.a;
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
				a->a->seg->apos, a->a->seg->bpos,
				a->a->seg->alen, a->a->seg->blen,
				a->a->plen);
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
	self->root.n = 0;
	self->margin.n = 0;
	self->bin.n = 0;
	kh_clear(&self->pos);
	self->n_res = 0;

	/* set query seq info */
	self->qid = qid;
	self->q[0] = _sec_fw(qid, seq, l_seq);
	self->q[1] = _sec_rv(qid, seq, l_seq);
	self->q[2] = _sec_fw(qid, seq, l_seq);

	/* store lmm */
	self->alloc.opaque = (void *)lmm;
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
	if(l_seq < self->mi.k || l_seq /* self->m*/ < self->min_score) {
		return(NULL);
	}

	/* clear buffers */
	mm_init_query(self, l_seq, seq, qid, lmm);
	debug("start: qid(%u)", qid);

	/* seed-chain-extend loop */
	debug("n_occ(%u)", self->mi.n_occ);
	for(uint64_t i = 0; i < self->mi.n_occ; i++) {
		if(mm_seed(self, i) == 0) { continue; }	/* seed not found */
		if(mm_chain(self, i) == 0) { continue; }/* chain not found */
		if(mm_extend(self, i) > 0) { break; }	/* at least one full-length alignment found */
	}
	if(self->n_res == 0) { return(NULL); }		/* unmapped */

	/* sort by score in reverse order */
	radix_sort_64x((v2u32_t *)self->root.a, self->n_res);

	/* prune alignments whose score is less than min_score threshold */
	uint32_t n_all = mm_prune_regs(self);

	/* collect supplementaries (split-read collection) */
	uint32_t n_uniq = (1 ? mm_post_ava : mm_post_map)(self);
	debug("n_all(%u), n_uniq(%u)", n_all, n_uniq);


	for(uint64_t i = 0; i < n_all; i++) {
		debug("bid(%u)", self->root.a[i].lid);
		mm_bin_t *bin = (mm_bin_t *)&self->bin.a[self->root.a[i].lid];

		debug("n_aln(%u), plen(%u), lb(%u), ub(%u)", bin->n_aln, bin->plen, bin->lb, bin->ub);

		for(uint64_t j = 0; j < bin->n_aln; j++) {
			debug("j(%lu), aln(%p)", j, bin->aln[j]);
		}

	}

	/* allocate reg array from memory arena */
	return(mm_pack_reg(self, n_all, n_uniq));
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
	free(t->root.a);
	free(t->margin.a);
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
mm_tbuf_t *mm_tbuf_init(mm_tbuf_params_t *u)
{
	mm_tbuf_t *t = (mm_tbuf_t *)calloc(1, sizeof(mm_tbuf_t));
	uint8_t const *lim = (uint8_t const *)0x800000000000;

	/* init dp context and copy constants */
	*t = (mm_tbuf_t){
		.mi = u->mi,
		.twlen = _ud(u->wlen, u->wlen), .tglen = _ud(u->glen, u->glen),
		.min_ratio = u->min_ratio,
		.min_score = u->min_score,
		.dp = gaba_dp_init(u->ctx, lim, lim),
		.alloc = u->alloc,
		.t[0] = _sec_fw(0xfffffffe, t->tail, 32)			/* tail section */
	};
	if(t->dp == NULL) { goto _fail; }

	/* miscellaneous */
	memset(t->tail, N, 96);									/* tail seq array */
	t->qtp = &t->q[1];										/* query tail section info pointer */
	kh_init_static(&t->pos, 128);							/* init hash */
	return(t);
_fail:
	mm_tbuf_destroy(t);
	return(NULL);
}
/* end of map.c */
/* mtmap.c */
/**
 * @struct mm_align_step_t
 * @brief batch object, placed in the unused space at the head of bseq_t
 */
typedef struct {
	uint32_t id, base_qid;
	lmm_t *lmm;						/* alignment result container (local memory allocator) */
} mm_align_step_t;

/**
 * @struct mm_align_s
 * @brief alignment pipeline context
 */
struct mm_align_s {
	bseq_file_t *fp;				/* input, set at the head of mm_align_file */
	mm_tbuf_params_t u;				/* mapper */
	mm_print_t *pr;					/* output */
	/* streaming */
	uint32_t icnt, ocnt;
	kvec_t(v4u32_t) hq;
	pt_t *pt;
	mm_tbuf_t *t[];					/* mm_tbuf_t* array at the tail */
};

/**
 * @fn mm_align_source
 * @brief source of the alignment pipeline
 */
static
void *mm_align_source(uint32_t tid, void *arg)
{
	mm_align_t *b = (mm_align_t *)arg;
	bseq_t *r = bseq_read(b->fp);
	if(r == NULL) { return(NULL); }

	/* update and assign base_qid */
	// if(b->base_qid == UINT32_MAX) { b->base_qid = atoi(r->seq[0].name); }

	/* allocate working buffer */
	mm_align_step_t *s = (mm_align_step_t *)r;		/* overlaps, use unused64[3] */
	*s = (mm_align_step_t){
		.id = b->icnt++,			/* assign id */
		// .base_qid = b->base_qid,
		.lmm = lmm_init(NULL, 512 * 1024)
	};
	// b->base_qid += r->n_seq;		/* update qid */
	return(s);
}

/**
 * @fn mm_align_worker
 */
static
void *mm_align_worker(uint32_t tid, void *arg, void *item)
{
	mm_align_t *b = (mm_align_t *)arg;
	mm_tbuf_t *t = (mm_tbuf_t *)b->t[tid];
	mm_align_step_t *s = (mm_align_step_t *)item;
	bseq_t *r = (bseq_t *)s;
	for(uint64_t i = 0; i < r->n_seq; i++) {
		uint32_t qid = s->base_qid + i;			/* FIXME: parse qid from name with atoi when -M is set */
		r->seq[i].u64 = (uintptr_t)mm_align_seq(t, r->seq[i].l_seq, r->seq[i].seq, qid, s->lmm);
	}
	return(s);
}

/**
 * @fn mm_align_drain_intl
 */
static _force_inline
void mm_align_drain_intl(mm_align_t *b, mm_align_step_t *s)
{
	bseq_t *r = (bseq_t *)s;
	for(uint64_t i = 0; i < r->n_seq; i++) {
		mm_reg_t *reg = (mm_reg_t *)r->seq[i].u64;
		mm_print_mapped(b->pr, b->u.mi.s, &r->seq[i], reg);	/* mapped */
		lmm_free(s->lmm, reg);					/* FIXME: gaba_alignment_t objects must be freed by gaba_dp_res_free */
	}
	free(r->base);
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
	mm_align_t *b = (mm_align_t *)arg;
	mm_align_step_t *s = (mm_align_step_t*)item;

	kv_hq_push(v4u32_t, incq_comp, b->hq, ((v4u32_t){.u64 = {s->id, (uintptr_t)s}}));
	while(b->hq.n > 1 && b->hq.a[1].u64[0] == b->ocnt) {
		b->ocnt++;
		s = (mm_align_step_t*)kv_hq_pop(v4u32_t, incq_comp, b->hq).u64[1];
		mm_align_drain_intl(b, s);
	}
	return;
}

/**
 * @fn mm_align_destroy
 * @brief destroy alignment pipeline context
 */
static _force_inline
void mm_align_destroy(mm_align_t *b)
{
	if(b == NULL) { return; }

	/* destroy threads */
	for(mm_tbuf_t **p = (mm_tbuf_t **)b->t; *p; p++) { mm_tbuf_destroy(*p); }

	/* destroy contexts */
	kv_hq_destroy(b->hq);
	gaba_clean(b->u.ctx);
	free(b);
	return;
}

/**
 * @fn mm_align_init
 * @brief create alignment pipeline context
 */
static _force_inline
mm_align_t *mm_align_init(mm_align_params_t const *a, mm_idx_t const *mi, pt_t *pt)
{
	// uint32_t const org = (a->flag & (MM_AVA | MM_COMP)) == MM_AVA ? 0x40000000 : 0;
	// uint32_t const thresh = (a->flag & MM_AVA ? 1 : 0) + org;

	/* malloc context and output buffer */
	mm_align_t *b = calloc(1, sizeof(mm_align_t) + sizeof(mm_tbuf_t *) * (pt_nth(pt) + 1));
	*b = (mm_align_t){
		/* options */
		#define _cp(_x)		._x = a->_x
		.u = {
			.mi = *mi,	/* .org = org, .thresh = thresh, */
			_cp(flag), _cp(wlen), _cp(glen),
			_cp(min_ratio), _cp(min_score),
			.ctx = gaba_init(&a->p),
			.alloc = {
				.opaque = NULL,
				.lmalloc = (gaba_lmalloc_t)lmm_malloc,
				.lfree = (gaba_lfree_t)lmm_free
			}
		},
		#undef _cp
		/* pipeline contexts */
		.icnt = 0, .ocnt = 0,
		.hq = { .n = 1, .m = 1, .a = NULL },
		/* threads */
		.pt = pt
	};

	/* init output queue, buf and printer */
	if(b->u.ctx == NULL || b->pt == NULL) { goto _fail; }

	/* initialize threads */
	for(uint64_t i = 0; i < pt_nth(pt); i++) {
		if((b->t[i] = (void *)mm_tbuf_init(&b->u)) == 0) { goto _fail; }
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
int mm_align_file(mm_align_t *b, bseq_file_t *fp, mm_print_t *pr)
{
	if(fp == NULL || pr == NULL) { return(-1); }
	b->fp = fp; b->pr = pr;		/* input and output */
	mm_print_header(pr, b->u.mi.n_seq, b->u.mi.s);
	return(pt_stream(b->pt, b,	/* multithreaded mapping */
		mm_align_source, mm_align_worker, mm_align_drain
	));
}
/* end of mtmap.c */

/* printer.c */
/**
 * @struct mm_print_fn_t
 * @brief alignment result (record) printer interface
 */
typedef struct {
	mm_print_header_t header;
	mm_print_mapped_t mapped;		/* nreg contains #regs in lower 32bit, #uniq in higher 32bit */
} mm_print_fn_t;

/**
 * @struct mm_print_t
 */
struct mm_print_s {
	uint8_t *base, *tail, *p;
	uint8_t conv[40];				/* binary -> string conv table */
	mm_print_fn_t fn;
	uint64_t tags;					/* sam optional tags */
	char *arg_line;
	char *rg_line, *rg_id;
};

/**
 * @struct mm_tmpbuf_t
 * @brief temporary string storage for output formatter
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
 * @macro _loadb_u64
 * @brief bit string loader (for path-string parsing)
 */
#define _loadb_u64(_ptr, _pos) ({ \
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
	mm_print_t *b = (mm_print_t *)_b;
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

/**
 * @fn mm_print_sam_num
 * @brief print number in sam tag format
 */
static _force_inline
uint64_t mm_print_sam_num(
	mm_print_t *b,
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
	mm_print_t *b,
	bseq_seq_t const *t)
{
	static uint8_t const tag_size[32] = {
		0, 1, 0xfe, 1,  0, 0, 4, 0,  0xfe, 4, 0, 0,  0, 0, 0, 0,
		0, 0, 0, 2,     0, 0, 0, 0,  0, 0, 0xff, 0,  0, 0, 0, 0,
	};
	uint8_t const *p = t->tag;
	uint64_t n_tag = t->n_tag;
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
void mm_print_sam_header(mm_print_t *b, uint32_t n_seq, mm_idx_seq_t const *ref)
{
	/* header */
	_putsk(b, "@HD\tVN:1.0\tSO:unsorted\n");

	/* sequences */
	for(uint64_t i = 0; i < n_seq; i++) {
		_putsk(b, "@SQ\tSN:");
		_putsn(b, ref[i].name, ref[i].l_name);		/* sequence name */
		_putsk(b, "\tLN:");
		_putn(b, ref[i].l_seq);						/* sequence length */
		_cr(b);
	}

	/* print read group */
	if(b->tags & 0x01ULL<<MM_RG) { _puts(b, b->rg_line); _cr(b); }

	/* program info (note: command line and version should be included) */
	_putsk(b, "@PG\tID:minialign\tPN:minialign\tVN:");
	_puts(b, version());
	_putsk(b, "\tCL:");
	_puts(b, b->arg_line);
	_cr(b);
	return;
}

/**
 * @fn mm_print_sam_unmapped
 * @brief print sam unmapped record
 */
static
void mm_print_sam_unmapped(
	mm_print_t *b,
	bseq_seq_t const *t)
{
	_putsn(b, t->name, t->l_name);					/* print name */
	_putsk(b, "\t4\t*\t0\t0\t*\t*\t0\t0\t");		/* unmapped fixed string */
	_putsnt(b, t->seq, t->l_seq, decaf); _t(b);		/* sequence */
	if(t->qual[0] != '\0') {						/* quality if available */
		_putsn(b, t->qual, t->l_seq);
	} else {
		_put(b, '*');
	}
	mm_restore_sam_tags(b, t); _cr(b);				/* tags */
	return;
}

/**
 * @fn mm_print_sam_mapped_core
 */
static _force_inline
void mm_print_sam_mapped_core(
	mm_print_t *b,
	mm_idx_seq_t const *r,
	bseq_seq_t const *t,
	mm_aln_t const *a,
	uint16_t flag)
{
	const gaba_path_section_t *s = &a->a->seg[0];
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

	// fprintf(stderr, "ppos(%u)\n", a->a->seg->ppos);
	gaba_dp_print_cigar_forward(mm_cigar_printer, b, a->a->path, 0, a->a->plen);
	if(tl) {
		_putn(b, tl);
		_put(b, (flag&0x900)? 'H' : 'S');					/* print tail clip */
	}
	_putsk(b, "\t*\t0\t0\t");								/* mate tags, unused */

	/* print sequence */
	if(s->bid & 0x01) {
		_putsntr(b, &t->seq[t->l_seq-qe], qe-qs, decar);
	} else {
		_putsnt(b, &t->seq[qs], qe-qs, decaf);
	}
	_t(b);

	/* print quality string if available */
	if(t->qual[0] != '\0') {
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
	mm_print_t *b,
	mm_idx_seq_t const *r,
	bseq_seq_t const *t,
	mm_aln_t const *a)
{
	/* rname,pos,strand,CIGAR,mapQ,NM; */
	const gaba_path_section_t *s = &a->a->seg[0];
	uint32_t rs = s->apos;
	uint32_t hl = s->bpos, tl = s->bpos + s->blen;			/* head and tail clips */
	
	_putsn(b, r->name, r->l_name); _c(b);					/* rname */
	_putn(b, rs + 1); _c(b);								/* rpos */
	_put(b, (s->bid & 0x01) ? '+' : '-'); _c(b);			/* direction */

	/* print cigar */
	if(hl != 0) { _putn(b, hl); _put(b, 'H'); }				/* always hard clipped */
	gaba_dp_print_cigar_forward(mm_cigar_printer, b, a->a->path, 0, a->a->plen);
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
	mm_print_t *b,
	mm_idx_seq_t const *r,
	bseq_seq_t const *t,
	mm_aln_t const *a)
{
	/*
	 * print MD tag
	 * note: this operation requires head-to-tail comparison of two sequences,
	 *       which may result in several percent performance degradation.
	 */
	uint64_t const f = b->tags;
	if((f & 0x01ULL<<MM_MD) == 0) { return; }				/* skip if disabled */

	const gaba_path_section_t *s = &a->a->seg[0];
	uint32_t rs = s->apos, qs = s->bpos;

	_puts(b, "\tMD:Z:");
	// uint64_t const *p = (uint64_t const *)a->a->path;
	// int64_t pos = 0, lim = a->a->plen;

	uint64_t const *p = (uint64_t const *)((uint64_t)a->a->path & ~(sizeof(uint64_t) - 1));
	uint64_t lim = (((uint64_t)a->a->path & sizeof(uint32_t)) ? 32 : 0) + a->a->plen;
	uint64_t ridx = a->a->plen;

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
		uint64_t arr = _loadb_u64(p, lim - ridx);
		uint64_t cnt = tzcnt(~arr) - (arr & 0x01);		/* count #ins */
		ridx -= cnt;
		qp += dir ? -cnt : cnt;

		if((arr & 0x01) == 0) {							/* is_del */
			ridx -= cnt = tzcnt(arr);					/* count #del */
			_putn(b, (int32_t)(rp - rb));
			_put(b, '^');
			_putsnt32(b, rp, cnt, decaf);
			rp += cnt;
			rb = rp;
		}

		/* match or mismatch */
		uint64_t acnt = 32;
		while(acnt == 32) {
			/* count diagonal */
			arr = _loadb_u64(p, lim - ridx);
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
				_put(b, decaf[*rp]);		/* print mismatch base */
				_put(b, '0');							/* padding */
				rp++;
			}
			_put(b, decaf[*rp]);			/* print mismatch base */
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
	mm_print_t *b,
	mm_aln_t const *a,
	uint32_t n_reg,
	uint64_t j)
{
	uint64_t const f = b->tags;

	/* read group */
	if(f & 0x01ULL<<MM_RG) {
		_putsk(b, "\tRG:Z:"); _puts(b, b->rg_id);
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
	mm_print_t *b,
	mm_idx_seq_t const *ref,
	bseq_seq_t const *t,				/* query sequence */
	mm_reg_t const *reg)			/* alignments, 1..n_uniq is printed in SA */
{
	uint64_t const f = b->tags;
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
			uint32_t rid = a->a->seg->aid>>1;
			mm_idx_seq_t const *r = &ref[rid];

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
	mm_print_t *b,
	mm_idx_seq_t const *ref,
	bseq_seq_t const *t,
	mm_reg_t const *reg)
{
	if(reg == NULL) { mm_print_sam_unmapped(b, t); return; }

	/* iterate over alignment sets */
	uint64_t const n = (b->tags & MM_OMIT_REP)? reg->n_uniq : reg->n_all;
	for(uint64_t i = 0, flag = 0; i < n; i++) {
		if(i >= reg->n_uniq) {
			flag = 0x100;				/* mark secondary */
		}

		mm_aln_t const *a = reg->aln[i];
		uint32_t rid = a->a->seg->aid>>1;
		mm_idx_seq_t const *r = &ref[rid];

		/* print body */
		mm_print_sam_mapped_core(b, r, t, a, flag);

		/* print general tags (scores, ...) */
		mm_print_sam_general_tags(b, a, reg->n_all, i);

		/* mismatch position (MD) */
		mm_print_sam_md(b, r, t, a);

		/* primary-specific tags */
		if(i == 0 && mm_print_sam_primary_tags(b, r, t, reg)) {
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
	mm_print_t *b,
	mm_aln_t const *a,
	uint8_t const *rp)
{
	uint64_t const *p = (uint64_t const *)((uint64_t)(a->a->path - 1) & ~(sizeof(uint64_t) - 1));
	int64_t pos = a->a->plen;
	while(pos > 0) {
		uint64_t arr, cnt;

		/* insertion */
		arr = _loadb_u64(p, pos);
		pos -= (cnt = lzcnt(~arr));		/* count 1's (insertions) */
		_putsn32(b, "--------------------------------", cnt);

		/* deletion */
		arr = _loadb_u64(p, pos);
		pos -= (cnt = lzcnt(arr) - 1);	/* count 0's (deletions) */
		_putsnt32(b, rp, cnt, decaf);
		rp += cnt;						/* advance pointer */
		do {
			arr = _loadb_u64(p, pos);
			cnt = lzcnt(arr^0x5555555555555555)>>1;		/* count diagonals */
			pos -= 2*cnt;
			_putsnt32(b, rp, cnt, decaf);
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
	mm_print_t *b,
	mm_aln_t const *a,
	uint8_t const *qp)
{
	uint64_t const *p = (uint64_t const *)((uint64_t)(a->a->path - 1) & ~(sizeof(uint64_t) - 1));
	int64_t pos = a->a->plen;
	while(pos > 0) {
		uint64_t arr, cnt;

		/* insertions */
		arr = _loadb_u64(p, pos);
		pos -= (cnt = lzcnt(~arr));
		_putsnt32(b, qp, cnt, decaf);
		qp += cnt;

		/* deletions */
		arr = _loadb_u64(p, pos);
		pos -= (cnt = lzcnt(arr) - 1);
		_putsn32(b, "--------------------------------", cnt);
		
		/* diagonals */
		do {
			arr = _loadb_u64(p, pos);
			cnt = lzcnt(arr^0x5555555555555555)>>1;
			pos -= 2*cnt;
			_putsnt32(b, qp, cnt, decaf);
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
	mm_print_t *b,
	mm_aln_t const *a,
	uint8_t const *qp)
{
	uint64_t const *p = (uint64_t const *)((uint64_t)(a->a->path - 1) & ~(sizeof(uint64_t) - 1));
	int64_t pos = a->a->plen;
	while(pos > 0) {
		uint64_t arr, cnt;

		/* insertions */
		arr = _loadb_u64(p, pos);
		pos -= (cnt = lzcnt(~arr));
		_putsntr32(b, qp, cnt, decar);
		qp -= cnt;

		/* deletions */
		arr = _loadb_u64(p, pos);
		pos -= (cnt = lzcnt(arr) - 1);
		_putsn32(b, "--------------------------------", cnt);
		
		/* diagonals */
		do {
			arr = _loadb_u64(p, pos);
			cnt = lzcnt(arr^0x5555555555555555)>>1;
			pos -= 2*cnt;
			_putsntr32(b, qp, cnt, decar);
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
	mm_print_t *b,
	mm_idx_seq_t const *r,
	bseq_seq_t const *t,
	mm_aln_t const *a)
{
	const gaba_path_section_t *s = &a->a->seg[0];

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
	mm_print_t *b,
	mm_idx_seq_t const *ref,
	bseq_seq_t const *t,
	mm_reg_t const *reg)
{
	uint64_t const n = (b->tags & MM_OMIT_REP)? reg->n_uniq : reg->n_all;
	for(uint64_t i = 0; i < n; i++) {
		mm_aln_t const *a = reg->aln[i];
		uint32_t rid = a->a->seg->aid>>1;
		mm_idx_seq_t const *r = &ref[rid];

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
	mm_print_t *b,
	mm_idx_seq_t const *ref,
	bseq_seq_t const *t,
	mm_reg_t const *reg)
{
	uint64_t const n = (b->tags & MM_OMIT_REP)? reg->n_uniq : reg->n_all;
	for(uint64_t i = 0; i < n; i++) {
		mm_aln_t const *a = reg->aln[i];
		const gaba_path_section_t *s = &a->a->seg[0];
		mm_idx_seq_t const *r = &ref[s->aid>>1];

		int32_t dcnt = (a->a->plen - a->a->gecnt)>>1, slen = dcnt + a->a->gecnt;
		int32_t mid = 100000.0 * (double)(dcnt - a->a->xcnt) / (double)slen;	/* percent identity */
		uint32_t rs = s->bid&0x01? r->l_seq - s->apos : s->apos + 1;
		uint32_t re = s->bid&0x01? s->apos + 1 : r->l_seq - s->apos;
		uint32_t qs = s->bpos + 1, qe = s->bpos + s->blen;

		/* sequence names */
		_putsn(b, t->name, t->l_name); _t(b);
		_putsn(b, r->name, r->l_name); _t(b);

		_putfi(int32_t, b, mid, 3); _t(b);				/* sequence identity */
		_putn(b, slen); _t(b);							/* alignment length */
		_putn(b, a->a->xcnt); _t(b);					/* mismatch count */
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
 * @fn mm_print_paf_mapped
 * @brief minimap paf format
 * qname ql qs qe qd rname rl rs re #m block_len mapq
 */
static
void mm_print_paf_mapped(
	mm_print_t *b,
	mm_idx_seq_t const *ref,
	bseq_seq_t const *t,
	mm_reg_t const *reg)
{
	uint64_t const n = (b->tags & MM_OMIT_REP)? reg->n_uniq : reg->n_all;
	for(uint64_t i = 0; i < n; i++) {
		mm_aln_t const *a = reg->aln[i];
		const gaba_path_section_t *s = &a->a->seg[0];
		mm_idx_seq_t const *r = &ref[s->aid>>1];

		uint32_t dcnt = (a->a->plen - a->a->gecnt)>>1;			/* #matches + #mismatches */
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
		uint64_t const f = b->tags;
		if(f & 0x01ULL<<MM_AS) { _putsk(b, "\tAS:i:"); _putn(b, a->a->score); }
		if(f & 0x01ULL<<MM_NM) { _putsk(b, "\tNM:i:"); _putn(b, a->a->xcnt + a->a->gecnt); }
		_cr(b);
	}
	return;
}
#undef _d
#undef _t
#undef _c
#undef _cr
#undef _sp
#undef _aln

/**
 * @fn mm_print_tag2flag
 */
static _force_inline
uint64_t mm_print_tag2flag(
	uint64_t n_tag,
	uint16_t const *tag)
{
	uint64_t f = 0;
	#define _t(_str)		{ if(mm_encode_tag(#_str) == *p) { f |= 0x01ULL<<MM_##_str; } }
	for(uint16_t const *p = tag, *t = &tag[n_tag]; p < t; p++) {
		_t(RG) _t(CO) _t(NH) _t(IH) _t(AS) _t(XS) _t(NM) _t(SA) _t(MD)
	}
	return(f);
	#undef _t
}

/**
 * @fn mm_print_destroy
 */
static _force_inline
void mm_print_destroy(
	mm_print_t *pr)
{
	fwrite(pr->base, sizeof(uint8_t), pr->p - pr->base, stdout);
	free(pr->base); free(pr);
	return;
}

/**
 * @fn mm_print_init
 */
static _force_inline
mm_print_t *mm_print_init(
	mm_print_params_t const *r)
{
	static const mm_print_fn_t printer[] = {
		[MM_SAM] = { .header = mm_print_sam_header, .mapped = mm_print_sam_mapped },
		[MM_MAF] = { .mapped = mm_print_maf_mapped },
		[MM_PAF] = { .mapped = mm_print_paf_mapped },
		[MM_BLAST6] = { .mapped = mm_print_blast6_mapped }
	};
	mm_print_t *pr = calloc(1, sizeof(mm_print_t));
	void *p = malloc(sizeof(uint8_t) * r->outbuf_size);
	*pr = (mm_print_t){
		.p = p, .base = p, .tail = p + r->outbuf_size,
		.tags = r->flag | mm_print_tag2flag(r->n_tag, r->tag),
		.fn = printer[r->format],
		.arg_line = mm_strdup(r->arg_line),
		.rg_line = mm_strdup(r->rg_line),
		.rg_id = mm_strdup(r->rg_id)
	};
	for(uint64_t i = 0; i < 9; ++i) { pr->conv[i] = (i + 1) % 10; }
	for(uint64_t i = 9; i < 40; ++i) { pr->conv[i] = (((i + 1) % 10)<<4) + (i + 1) / 10; }
	return(pr);
}

/* function dispatchers */
static _force_inline
void mm_print_header(mm_print_t *b, uint32_t n_seq, mm_idx_seq_t const *ref)
{
	if(b->fn.header) { b->fn.header(b, n_seq, ref); }
	return;
}
static _force_inline
void mm_print_mapped(mm_print_t *b, mm_idx_seq_t const *ref, bseq_seq_t const *t, mm_reg_t const *reg)
{
	b->fn.mapped(b, ref, t, reg);
	return;
}
/* end of printer.c */

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
 * @macro mm_split_foreach
 * @brief split string into tokens and pass each to _body.
 * p, len, and i are reserved for pointer, length, and #parsed.
 * break / continue can be used in the _body to terminate / skip the current element.
 */
#define mm_split_foreach(_ptr, _delims, _body) { \
	char const *_q = (_ptr); \
	int64_t i = 0; \
	v16i8_t _dv = _loadu_v16i8(_delims); \
	uint16_t _m, _mask = 0x02<<tzcnt(((v16i8_masku_t){		/* reserve space for '\0' */ \
		.mask = _mask_v16i8(_eq_v16i8(_set_v16i8('\0'), _dv)) \
	}).all); \
	_dv = _bsl_v16i8(_dv, 1); _mask--;						/* push '\0' at the head of the vector */ \
	do { \
		char const *_p = _q; \
		/* test char one by one until dilimiter found */ \
		while(((_m = ((v16i8_masku_t){ .mask = _mask_v16i8(_eq_v16i8(_set_v16i8(*_q), _dv)) }).all) & _mask) == 0) { _q++; } \
		/* delimiter found, pass to _body */ \
		char const *p = _p; \
		uint64_t l = _q++ - _p; \
		if(l > 0) { _body; i++; } \
	} while((_m & 0x01) == 0); \
}

/**
 * @macro oassert
 */
#define oassert(_o, _cond, ...) { \
	if(!(_cond)) { (_o)->log(_o, 'E', __func__, __VA_ARGS__); o->ecnt++; } \
}

/**
 * @fn mm_opt_atoi, mm_opt_atof
 */
static _force_inline
int64_t mm_opt_atoi(mm_opt_t *o, char const *arg, uint64_t len)
{
	if(arg == NULL) { return(0); }
	if(len == UINT32_MAX) { len = strlen(arg); }
	for(char const *p = arg, *t = &arg[len]; p < t && *p != '\0'; p++) {
		if(!isdigit(*p)) { oassert(o, 0, "unparsable number `%.*s'.", len, arg); return(0); }
	}
	return(atoi(arg));
}
static _force_inline
double mm_opt_atof(mm_opt_t *o, char const *arg, uint64_t len)
{
	if(arg == NULL) { return(0); }
	if(len == UINT32_MAX) { len = strlen(arg); }
	static char const allowed[16] __attribute__(( aligned(16) )) = {
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '-', '.', ',', 'e', 'E'
	};
	v16i8_t v = _load_v16i8(allowed);
	for(char const *p = arg, *t = &arg[len]; p < t && *p != '\0'; p++) {
		if(((v16i8_masku_t){ .mask = _mask_v16i8(_eq_v16i8(_set_v16i8(*p), v)) }).all == 0) {
			oassert(o, 0, "unparsable number `%.*s'.", len, arg); return(0);
		}
	}
	return(atof(arg));
}

/**
 * @fn mm_opt_parse_argv
 * @brief parse (argc, argv)-style option arrays. only argv is required here (MUST be NULL-terminated).
 */
static _force_inline
int mm_opt_parse_argv(mm_opt_t *o, char const *const *argv)
{
	char const *const *p = argv - 1;
	char const *q;
	#define _isarg(_q)	( (_q)[0] != '-' || (_q)[1] == '\0' )
	#define _x(x)		( (uint64_t)(x) )
	while((q = *++p)) {											/* p is a jagged array, must be NULL-terminated */
		if(_isarg(q)) {											/* option starts with '-' and longer than 2 letters, such as "-a" and "-ab" */
			kv_push(void *, o->parg, mm_strdup(q)); continue;	/* positional argument */
		}
		while(o->t[_x(*++q)].type == 1) { o->t[_x(*q)].fn(o, NULL); }/* eat boolean options */
		if(!o->t[_x(*q)].fn) { oassert(o, 0, "unknown option `-%c'.", *q); continue; }/* argument option not found */
		char const *r = q[1] ? q+1 : (p[1] && _isarg(p[1]) ? *++p : NULL);/* if the option ends without argument, inspect the next element in the jagged array (originally placed after space(s)) */
		oassert(o, o->t[_x(*q)].type != 2 || r, "missing argument for option `-%c'.", *q);
		if(o->t[_x(*q)].type != 2 || r) { o->t[_x(*q)].fn(o, r); }/* option with argument would be found at the tail */
	}
	kv_push(void *, o->parg, NULL); o->parg.n--;				/* always keep NULL-terminated */
	#undef _isarg
	#undef _x
	return(o->ecnt);
}

/**
 * @fn mm_opt_parse_line
 * @brief parse tab/space-delimited command line arguments
 */
static void mm_opt_parse_line(mm_opt_t *o, char const *arg)
{
	kvec_t(char) str = { 0 };
	kvec_t(char const *) ptr = { 0 };
	mm_split_foreach(arg, " \t\r\n", {
		kv_push(char const *, ptr, (char const *)kv_size(str));
		kv_pushm(char, str, p, l);
		kv_push(char, str, '\0');
	});
	kv_foreach(char const *, ptr, { *p += (ptrdiff_t)str.a; });
	kv_push(char const *, ptr, NULL);
	mm_opt_parse_argv(o, ptr.a);
	free(ptr.a);
	return;
}

/**
 * @fn mm_opt_load_conf
 */
static int mm_opt_load_conf(mm_opt_t *o, char const *arg)
{
	FILE *fp = fopen(arg, "r");
	if(fp == NULL) { oassert(o, 0, "failed to find configuration file `%s'.", arg); return(0); }

	kvec_t(char) str = { 0 };
	kv_reserve(char, str, 1024);
	while(1) {
		if((str.n += fread(str.a, sizeof(char), str.m - str.n, fp)) < str.m) { break; }
		kv_reserve(char, str, 2 * str.n);
	}
	fclose(fp); kv_push(char, str, '\0');
	mm_opt_parse_line(o, str.a);
	return(1);
}

/**
 * @fn mm_opt_preset
 * @brief parse preset line
 */
static void mm_opt_preset(mm_opt_t *o, char const *arg)
{
	struct mm_preset_s {
		char const *key; char const *val;
		struct mm_preset_s const *children[6];
	};
	#define _pre(_k, _v, ...)	&((struct mm_preset_s const){ .key = (_k), .val = (_v), .children = { __VA_ARGS__ } })
	struct mm_preset_s const *presets[] = {
		_pre("pacbio", "-k15 -w10 -a2 -b4 -p4 -q1 -Y50 -s50 -m0.3",
			_pre("clr", "-q1 -r2,2", NULL),
			_pre("ccs", "-b3 -p5", NULL)
		),
		_pre("ont", "-k15 -w10 -a2 -b4 -p4 -q1 -r2,2 -Y50 -s50 -m0.3",
			_pre("r7", "-b2 -p2",
				_pre("1d", "", NULL), _pre("2d", "", NULL)
			),
			_pre("r9", "",
				_pre("1d", "", NULL), _pre("1dsq", "", NULL), _pre("2d", "", NULL)
			),
			_pre("1d", "", NULL), _pre("1dsq", "", NULL), _pre("2d", "", NULL)
		),
		_pre("ava", "-k15 -w5 -a2 -b3 -p0 -q2 -Y50 -s30 -r0.05", NULL),
		NULL
	};
	#undef _pre

	struct mm_preset_s const *const *c = presets - 1;
	mm_split_foreach(arg, ".:", {		/* traverse preset param tree along with parsing */
		while(*++c && (strlen((*c)->key) != l || strncmp(p, (*c)->key, l) != 0)) {}
		if(!*c) {						/* terminate if not matched, not loaded from file */
			oassert(o, mm_opt_load_conf(o, p), "no preset params found for `%.*s'.", l, p);
			break;
		}
		mm_opt_parse_line(o, (*c)->val);/* apply recursively */
		c = (*c)->children - 1;			/* visit child nodes */
	});
	return;
}

/**
 * @fn mm_opt_rg
 * @brief parse RG tag
 */
static void mm_opt_rg(mm_opt_t *o, char const *arg)
{
	/* clear existing tags */
	free(o->r.rg_line); o->r.rg_line = NULL;
	free(o->r.rg_id);   o->r.rg_id = NULL;
	o->r.flag &= ~(0x01ULL<<MM_RG);

	/* parse */
	kvec_t(char) b = { 0 };
	char const *p = arg - 1;
	while(*++p != '\0') {				/* unescape tabs */
		kv_push(char, b, *p == '\\' ? (++p, '\t') : *p);
	}
	mm_split_foreach(b.a, "\t\r\n", {	/* find ID tag */
		if(strstr(p, "ID:") == p) {
			/* RG line sanity confirmed, save line and set flag */
			o->r.rg_line = b.a;
			o->r.rg_id = mm_strndup(p, l);
			o->r.flag |= 0x01ULL<<MM_RG;
			break;
		}
	});
	oassert(o, o->r.rg_id != NULL, "RG line must start with @RG and contains ID, like `@RG\\tID:1'.");
	return;
}

/**
 * @fn mm_opt_tags
 * @brief parse sam tag identifiers and push to tags array.
 */
static void mm_opt_tags(mm_opt_t *o, char const *arg)
{
	mm_split_foreach(arg, ",;:/", {
		oassert(o, l == 2, "unknown tag: `%.*s'.", l, p);
		kv_push(uint16_t, o->tags, mm_encode_tag(p));
	});
	return;
}

/**
 * @fn mm_opt_format
 * @brief convert format string to flag
 */
static void mm_opt_format(mm_opt_t *o, char const *arg)
{
	static struct mm_format_s { char const *k; uint32_t v; } const t[] = {
		{ "sam",    MM_SAM },
		{ "maf",    MM_MAF },
		{ "blast6", MM_BLAST6 },
		{ "paf",    MM_PAF },
		{ NULL, 0xff }
	}, *p = t - 1;
	while((++p)->k && strcmp(p->k, arg) != 0) {}
	o->r.format = p->v;			/* avoid warning */
	oassert(o, o->r.format != 0xff, "unknown output format `%s'.", arg);
	return;
}

/* save index file name */
static void mm_opt_fnw(mm_opt_t *o, char const *arg) { o->fnw = mm_strdup(arg); }

/* flags, global params */
static void mm_opt_keep_qual(mm_opt_t *o, char const *arg) { o->b.keep_qual = 1; }
static void mm_opt_ava(mm_opt_t *o, char const *arg) { o->a.flag |= MM_AVA; }
static void mm_opt_comp(mm_opt_t *o, char const *arg) { o->a.flag |= MM_COMP; }
static void mm_opt_omit_rep(mm_opt_t *o, char const *arg) { o->a.flag |= MM_OMIT_REP; }
static void mm_opt_verbose(mm_opt_t *o, char const *arg) { o->verbose = arg ? (isdigit(*arg) ? mm_opt_atoi(o, arg, UINT32_MAX) : strlen(arg) + 1) : 0; }
static void mm_opt_threads(mm_opt_t *o, char const *arg) {
	o->nth = mm_opt_atoi(o, arg, UINT32_MAX);
	oassert(o, o->nth < MAX_THREADS, "#threads must be less than %d.", MAX_THREADS);
}
static void mm_opt_help(mm_opt_t *o, char const *arg) { o->verbose = 1; o->help++; }

/* input configurations */
static void mm_opt_base_id(mm_opt_t *o, char const *arg) {
	o->a.base_rid = o->a.base_qid = 0x80000000;	/* first clear ids (= infer-from-name state) */
	mm_split_foreach(arg, ",;:/*", {			/* FIXME: '*' should be treated in separete */
		switch(i) {
			case 0: o->a.base_rid = mm_opt_atoi(o, p, l);	/* fall through to set both */
			case 1: o->a.base_qid = mm_opt_atoi(o, p, l);
		}
	});
}
static void mm_opt_circular(mm_opt_t *o, char const *arg) {
	if(kh_str_ptr(&o->c.circ) == NULL) { kh_str_init_static(&o->c.circ, 128); }
	mm_split_foreach(arg, ",;:/*", { kh_str_put(&o->c.circ, p, l, "", 0); });	/* FIXME: '*' should be treated in separete */
}

/* indexing */
static void mm_opt_kmer(mm_opt_t *o, char const *arg) {
	o->c.k = mm_opt_atoi(o, arg, UINT32_MAX);
	oassert(o, o->c.k > 1 && o->c.k < 32, "k must be inside [1,32).");
}
static void mm_opt_window(mm_opt_t *o, char const *arg) {
	o->c.w = mm_opt_atoi(o, arg, UINT32_MAX);
	oassert(o, o->c.w > 1 && o->c.w < 32, "w must be inside [1,32).");
}
static void mm_opt_bin(mm_opt_t *o, char const *arg) {
	o->c.b = mm_opt_atoi(o, arg, UINT32_MAX);
	oassert(o, o->c.b > 1 && o->c.b < 32, "b must be inside [1,32).");
}
static void mm_opt_frq(mm_opt_t *o, char const *arg) {
	o->c.n_frq = 0;			/* clear counter */
	mm_split_foreach(arg, ",;:/", {
		float f = o->c.frq[o->c.n_frq++] = mm_opt_atof(o, p, l);
		oassert(o, i < MAX_FRQ_CNT, "#thresholds must not exceed %d.", MAX_FRQ_CNT);
		oassert(o, f > 0.0 && f < 1.0, "invalid threshold `%f' parsed from `%.*s'.", f, l, p);
		oassert(o, i == 0 || (o->c.frq[i-1] > o->c.frq[i]), "frequency thresholds must be descending.");
	});
}

/* mapping */
static void mm_opt_wlen(mm_opt_t *o, char const *arg) {
	o->a.wlen = mm_opt_atoi(o, arg, UINT32_MAX);
	oassert(o, o->a.wlen > 100 || o->a.wlen < 100000, "window edge length must be inside [100,100000).");
}
static void mm_opt_glen(mm_opt_t *o, char const *arg) {
	o->a.glen = mm_opt_atoi(o, arg, UINT32_MAX);
	oassert(o, o->a.glen > 100 || o->a.glen < 10000, "gap chain length must be inside [100,10000).");
}
static void mm_opt_match(mm_opt_t *o, char const *arg) {
	int32_t m = mm_opt_atoi(o, arg, UINT32_MAX);
	for(uint64_t i = 0; i < 16; i++) {
		if((i&0x03) == (i>>2)) { o->a.p.score_matrix[i] = m; }
	}
	oassert(o, m > 0 && m < 7, "match award (-a) must be inside [1,7].");
}
static void mm_opt_mismatch(mm_opt_t *o, char const *arg) {
	int32_t x = mm_opt_atoi(o, arg, UINT32_MAX);
	for(uint64_t i = 0; i < 16; i++) {
		if((i&0x03) != (i>>2)) { o->a.p.score_matrix[i] = x; }
	}
	oassert(o, x > 0 && x < 7, "mismatch penalty (-b) must be inside [1,7].");
}
static void mm_opt_mod(mm_opt_t *o, char const *arg) {
	mm_split_foreach(arg, ",;:/", {
		oassert(o, idxaf[(uint64_t)p[0]] != 0, "unknown base `%c' in modifier `%.*s'.", p[0], l, p);
		oassert(o, idxaf[(uint64_t)p[1]] != 0, "unknown base `%c' in modifier `%.*s'.", p[1], l, p);
		o->a.p.score_matrix[idxaf[(uint64_t)p[1]] * 4 + idxaf[(uint64_t)p[0]]] += atoi(&p[2]);
	});
}
static void mm_opt_gi(mm_opt_t *o, char const *arg) {
	int32_t gi = mm_opt_atoi(o, arg, UINT32_MAX);
	o->a.p.gi = gi;
	oassert(o, gi < 32, "gap open penalty (-p) must be inside [0,32].");
}
static void mm_opt_ge(mm_opt_t *o, char const *arg) {
	int32_t ge = mm_opt_atoi(o, arg, UINT32_MAX);
	o->a.p.ge = ge;
	oassert(o, ge > 0 && ge < 32, "gap extension penalty (-q) must be inside [1,32].");
}
static void mm_opt_gf(mm_opt_t *o, char const *arg) {
	mm_split_foreach(arg, ",;:/", {
		switch(i) {
			case 0: o->a.p.gfa = mm_opt_atoi(o, p, l);			/* fall through to set both */
			case 1: o->a.p.gfb = mm_opt_atoi(o, p, l);
		}
	});
}
static void mm_opt_xdrop(mm_opt_t *o, char const *arg) {
	int32_t xdrop = mm_opt_atoi(o, arg, UINT32_MAX);
	o->a.p.xdrop = xdrop;
	oassert(o, xdrop > 10 && xdrop < 128, "X-drop cutoff must be inside [10,128].");
}
static void mm_opt_min_len(mm_opt_t *o, char const *arg) {
	o->b.min_len = mm_opt_atoi(o, arg, UINT32_MAX);
	oassert(o, o->b.min_len > 0, "minimum sequence length must be > 0.");
}
static void mm_opt_min_score(mm_opt_t *o, char const *arg) {
	o->a.min_score = mm_opt_atoi(o, arg, UINT32_MAX);
	oassert(o, o->a.min_score > 0, "minimum alignment score must be > 0.");
}
static void mm_opt_min_ratio(mm_opt_t *o, char const *arg) {
	o->a.min_ratio = mm_opt_atof(o, arg, UINT32_MAX);
	oassert(o, o->a.min_ratio > 0.0 && o->a.min_ratio < 1.0, "minimum alignment score ratio must be inside [0.0,1.0].");
}
static void mm_opt_batch(mm_opt_t *o, char const *arg) {
	o->b.batch_size = mm_opt_atoi(o, arg, UINT32_MAX);
	oassert(o, o->b.batch_size > 64 * 1024, "batch size must be > 64k.");
}
static void mm_opt_outbuf(mm_opt_t *o, char const *arg) {
	o->r.outbuf_size = mm_opt_atoi(o, arg, UINT32_MAX);
	oassert(o, o->r.outbuf_size > 64 * 1024, "output buffer size must be > 64k.");
}

/**
 * @fn mm_opt_destroy
 */
static _force_inline
void mm_opt_destroy(mm_opt_t *o)
{
	pt_destroy(o->pt);
	kv_foreach(void *, o->parg, { free(*p); });
	free(o->parg.a);
	free(o->fnw);
	free(o->tags.a);
	free(o->r.rg_line);
	free(o->r.rg_id);
	free(o->r.arg_line);
	free(o);
	return;
}

/**
 * @fn mm_opt_init
 */
static _force_inline
mm_opt_t *mm_opt_init(char const *const *argv)
{
	mm_opt_t *o = calloc(1, sizeof(mm_opt_t));
	*o = (mm_opt_t){
		/* global */
		.nth = 1,
		/* input */
		.b = { .batch_size = 512 * 1024, .min_len = 1, },
		/* indexing params */
		.c = {
			.k = 15, .w = 32, .b = 14,		/* w will be overwritten later */
			.n_frq = 3, .frq = { 0.05, 0.01, 0.001 },
		},
		/* mapping */
		.a = {
			.wlen = 7000, .glen = 7000,
			.min_score = 50, .min_ratio = 0.3,
			.p = {
				.score_matrix = { 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1 },
				.gi = 1, .ge = 1, .gfa = 0, .gfb = 0, .xdrop = 50
			},
		},
		/* output */
		.r = { .outbuf_size = 512 * 1024, .arg_line = mm_join(argv, ' '), },
		/* initialized time and loggers */
		.inittime = realtime(),
		.verbose = 1, .fp = (void *)stderr, .log = (mm_log_t)mm_log_printer,
		/* parsers: mapping from option character to functionality */
		.t = {
			['x'] = { 2, mm_opt_preset },
			['R'] = { 2, mm_opt_rg },
			['T'] = { 2, mm_opt_tags },
			['O'] = { 2, mm_opt_format },
			['d'] = { 2, mm_opt_fnw },

			['X'] = { 1, mm_opt_ava },
			['A'] = { 1, mm_opt_comp },
			['P'] = { 1, mm_opt_omit_rep },
			['Q'] = { 1, mm_opt_keep_qual },
			['v'] = { 3, mm_opt_verbose },
			['h'] = { 1, mm_opt_help },
			['t'] = { 2, mm_opt_threads },

			['k'] = { 2, mm_opt_kmer },
			['w'] = { 2, mm_opt_window },
			['c'] = { 2, mm_opt_circular },
			['f'] = { 2, mm_opt_frq },
			['B'] = { 2, mm_opt_bin },
			['C'] = { 3, mm_opt_base_id },
			['L'] = { 2, mm_opt_min_len },

			['W'] = { 2, mm_opt_wlen },
			['G'] = { 2, mm_opt_glen },
			['a'] = { 2, mm_opt_match },
			['b'] = { 2, mm_opt_mismatch },
			['e'] = { 2, mm_opt_mod },
			['p'] = { 2, mm_opt_gi },
			['q'] = { 2, mm_opt_ge },
			['r'] = { 2, mm_opt_gf },
			['Y'] = { 2, mm_opt_xdrop },
			['s'] = { 2, mm_opt_min_score },
			['m'] = { 2, mm_opt_min_ratio },
			['1'] = { 2, mm_opt_batch },
			['2'] = { 2, mm_opt_outbuf }
		}
	};
	if(mm_opt_parse_argv(o, ++argv)) { free(o); return(NULL); }
	if(o->c.w >= 32) { o->c.w = (int)(2.0/3.0 * o->c.k + .499); }		/* calc. default window size (proportional to kmer length) if not specified */
	if((o->pt = pt_init(o->nth)) == NULL) {
		mm_opt_destroy(o); return(NULL);
	}
	return(o);
}

/**
 * @fn mm_print_help
 */
static
void mm_print_help(mm_opt_t const *o)
{
	if(o->verbose == 0) { return; }

	#define _msg(_level, ...) { \
		o->log(o, 16 + _level, __func__, __VA_ARGS__); \
	}

	_msg(1, "\n"
			"  minialign - fast and accurate alignment tool for long reads\n"
			"");
	_msg(1, "Usage:\n"
			"  first trial:\n"
			"    $ minialign -t4 -xont.r9.1d ref.fa ont_r9.4_1d.fq > mapping.sam\n"
			"");
	_msg(1, "  mapping on a prebuilt index (saves ~1min for human genome per run):\n"
			"    $ minialign [indexing options] -d index.mai ref.fa\n"
			"    $ minialign index.mai reads.fq > mapping.sam\n"
			"");
	_msg(2, "  all-versus-all alignment in a read set:\n"
			"    $ minialign -X -xava reads.fa [reads.fa ...] > allvsall.paf\n"
			"");
	_msg(1, "Options:");
	_msg(1, "  Global:");
	_msg(1, "    -x STR       load preset params [ont]");
	_msg(1, "                   {pacbio.{clr,ccs},ont.{r7,r9}.{1d,1dsq,2d},ava}");
	_msg(1, "    -t INT       number of threads [%d]", o->nth);
	_msg(1, "    -d FILE      index construction mode, dump index to FILE");
	_msg(2, "    -X           all-versus-all mode, `-XA' to calculate both A~B and B~A.");
	_msg(1, "    -v [INT]     show version number / set verbose level");
	_msg(1, "  Indexing:");
	_msg(1, "    -k INT       k-mer size [%d]", o->c.k);
	_msg(1, "    -w INT       minimizer window size [{-k}*2/3]");
	_msg(1, "    -c STR,...   circular reference name, `*' to mark all as circular []");
	_msg(2, "    -B INT       1st stage hash table size base [%u]", o->c.b)
	_msg(2, "    -C INT[,INT] set base rid and qid, `*' to infer from seq. name [%u, %u]", o->a.base_rid, o->a.base_qid);
	_msg(2, "    -L INT       min seq length; 0 to disable [%u]", o->b.min_len);
	_msg(1, "  Mapping:");
	_msg(2, "    -f FLOAT,... occurrence thresholds [%rf]", o->c.n_frq, o->c.frq);
	_msg(1, "    -a INT       match award [%d]", o->a.p.score_matrix[0]);
	_msg(1, "    -b INT       mismatch penalty [%d]", o->a.p.score_matrix[1]);
	_msg(1, "    -e STR,...   score matrix modifier, `GA+3' adds 3 to (r,q)=(G,A) pair");
	_msg(1, "    -p INT       gap open penalty offset for large indels [%d]", o->a.p.gi);
	_msg(1, "    -q INT       per-base penalty for large indels [%d]", o->a.p.ge);
	_msg(1, "    -r INT[,INT] per-base penalty for small ins[,del] (0 to disable) [%d,%d]", o->a.p.gfa, o->a.p.gfb);
	_msg(2, "    -Y INT       X-drop threshold [%d]", o->a.p.xdrop);
	_msg(1, "    -s INT       minimum score [%d]", o->a.min_score);
	_msg(1, "    -m INT       minimum score ratio to max [%1.2f]", o->a.min_ratio);
	_msg(1, "  Output:");
	_msg(1, "    -O STR       output format {sam,maf,blast6,paf} [%s]",
		(char const *[]){ "sam", "maf", "blast6", "blasr1", "blasr4", "paf", "mhap", "falcon" }[o->r.format]);
	_msg(2, "    -P           omit secondary (repetitive) alignments");
	_msg(1, "    -Q           include quality string");
	_msg(2, "    -R STR       read group header line, such as `@RG\\tID:1' [%s]", o->r.rg_line ? o->r.rg_line : "");
	_msg(2, "    -T STR,...   optional tags: {RG,AS,XS,NM,NH,IH,SA,MD} []");
	_msg(2, "                   RG is also inferred from `-R'");
	_msg(2, "                   supp. records are omitted when SA is enabled");
	_msg(2, "                   tags in the input BAM file will also transferred");
	_msg(1, "");
	if(o->verbose < 2) {
		_msg(0, "  Pass -hvv to show all the options.");
		_msg(0, "");
	}

	#undef _msg
	return;
}

/**
 * @fn main_index
 */
static _force_inline
int main_index(mm_opt_t *o)
{
	char const *fn;
	if(!mm_endswith(o->fnw, ".mai")) {
		o->log(o, 'W', __func__, "index filename does not ends with `.mai' (added).");
		o->fnw = mm_append(o->fnw, ".mai");
	}
	pg_t *pg = pg_init(fopen(fn = o->fnw, "wb"), o->pt);
	if(!pg) { goto _main_index_fail; }

	/* iterate over index *blocks* */
	bseq_params_t br = o->b;			/* copy to local stack */
	br.keep_qual = 0; br.n_tag = 0;		/* overwrite */
	kv_foreach(void *, o->parg, {
		bseq_file_t *fp = bseq_open(&br, *p);
		mm_idx_t *mi = mm_idx_gen(&o->c, fp, o->pt);
		o->a.base_rid += bseq_close(fp);

		/* check sanity of the index */
		if(mi == NULL) { fn = *p; goto _main_index_fail; }

		/* dump index */
		o->log(o, 9, __func__, "built index for %lu target sequence(s).", mi->n_seq);
		mm_idx_dump(mi, pg, (write_t)pgwrite);
		mm_idx_destroy(mi);
	});
	pg_destroy(pg);
	return(0);

_main_index_fail:;
	uint64_t is_idx = fn == o->fnw;
	o->log(o, 'E', __func__, "failed to open %s file `%s'%s. Please check file path and its %s.",
		is_idx ? "index" : "sequence",
		fn,
		is_idx ? " in write mode" : "",
		is_idx ? "permission" : "format"
	);
	return(1);
}

/**
 * @fn main_align
 */
static _force_inline
void main_align_error(mm_opt_t *o, int stat, char const *fn, char const *file)
{
	switch(stat) {
	case 1: o->log(o, 'E', fn, "failed to instanciate alignment context."); break;
	case 2: o->log(o, 'E', fn, "failed to open index file `%s'. Please check file path and its version.", file); break;
	case 3: o->log(o, 'E', fn, "failed to open sequence file `%s'. Please check file path and its format.", file); break;
	}
	return;
}
static _force_inline
int main_align(mm_opt_t *o)
{
	pg_t *pg = NULL;
	mm_idx_t *mi = NULL;
	mm_align_t *aln = NULL;
	mm_print_t *pr = NULL;

	/* first test if prebuilt index is available, then instanciate pg reader. pg != NULL indicates prebuilt index is available for this batch */
	if(mm_endswith(*o->parg.a, ".mai") && (pg = pg_init(fopen(*o->parg.a, "rb"), o->pt)) == NULL) {
		goto _main_align_fail;
	}
	uint64_t rt = 1, qh = 1;					/* tail of reference-side arguments, head of query-side arguments */
	if((o->a.flag & MM_AVA) && pg == NULL) {	/* all-versus-all mode without prebuilt index is a special case */
		rt = o->parg.n; qh = 0;					/* calc all-versus-all between the arguments */
	}
	if(qh == o->parg.n) {						/* if query-side file is missing, pour stdin to query */
		o->log(o, 1, __func__, "query-side input redirected to stdin.");
		kv_push(void *, o->parg, mm_strdup("-"));
		kv_push(void *, o->parg, NULL); o->parg.n--;
	}

	/* file I/O and error handlings */
	#define _bseq_open_wrap(_b, _fn) ({ \
		bseq_file_t *_fp = bseq_open(_b, _fn); \
		if(_fp == NULL) { main_align_error(o, 2, __func__, _fn); goto _main_align_fail; } \
		_fp; \
	})
	#define _mm_idx_load_wrap(_pg, _r) ({ \
		mm_idx_t *_mi = NULL; \
		if(_pg) { \
			_mi = mm_idx_load(_pg, (read_t)pgread); \
			pg_freeze(_pg);		/* release thread worker */ \
		} else if(*(_r)) { \
			debug("ref(%s)", *(_r)); \
			bseq_file_t *_fp = _bseq_open_wrap(&br, *(_r)); \
			_mi = mm_idx_gen(&o->c, _fp, o->pt); \
			o->a.base_rid += bseq_close(_fp); \
			if(_mi == NULL) { main_align_error(o, 3, __func__, *(_r)); goto _main_align_fail; } \
			(_r)++;		/* increment r when in the on-the-fly mode and the index is correctly built */ \
		} \
		_mi; \
	})

	bseq_params_t br = o->b, bq = o->b;
	br.keep_qual = 0; br.n_tag = 0;
	pr = mm_print_init(&o->r);

	/* iterate over index *blocks* */
	char const *const *r = (char const *const *)o->parg.a;
	char const *const *t = (char const *const *)&o->parg.a[rt];
	while(r < t && (mi = _mm_idx_load_wrap(pg, r))) {
		o->log(o, 9, __func__, "loaded/built index for %lu target sequence(s).", mi->n_seq);
		/* initialize alignment context for this batch */
		if((aln = mm_align_init(&o->a, mi, o->pt)) == NULL) {
			main_align_error(o, 1, __func__, NULL);
			goto _main_align_fail;
		}
		/* iterate over queries */
		for(char const *const *q = (char const *const *)&o->parg.a[qh]; *q; q++) {
			debug("query(%s)", *q);
			bseq_file_t *fp = _bseq_open_wrap(&bq, *q);
			mm_align_file(aln, fp, pr);
			bseq_close(fp);
			o->log(o, 9, __func__, "finished mapping `%s' onto `%s'.", *q, pg ? *o->parg.a : r[-1]);
		}
		mm_align_destroy(aln);
		mm_idx_destroy(mi);
	}
	return(0);

_main_align_fail:;
	pg_destroy(pg);
	mm_idx_destroy(mi);
	mm_print_destroy(pr);
	mm_align_destroy(aln);
	return(1);
}

/**
 * @fn main
 */
int _export(main)(int argc, char *argv[])
{
	int ret = 1;

	/* unittest hook, see unittest.h for the details */
	#if UNITTEST != 0
	if(argc > 1 && strcmp(argv[1], "unittest") == 0) {
		return(unittest_main(argc, argv));
	}
	#endif

	/* elevate memory limit */
	liftrlimit();

	/* init option object (init base time) and parse args */
	mm_opt_t *o = mm_opt_init((char const *const *)argv);
	if(o == NULL) { return(1); }
	o->log(o, 0, __func__, "Version: %s, Build: %s", MM_VERSION, MM_ARCH);		/* always print version */
	if(o->help || o->parg.n == 0) {								/* when -h is passed or no input file is given, print help message */
		if(o->help) { o->fp = stdout; ret = 0; }				/* set logger to stdout when invoked by -h option, also exit status is 0 (not an error) */
		mm_print_help(o);
		goto _main_final;
	}
	if((ret = (o->fnw ? main_index : main_align)(o)) == 0) {	/* dispatch tasks and get return code */
		o->log(o, 1, __func__, "Command: %s", o->r.arg_line);	/* print log when succeeded */
		o->log(o, 1, __func__, "Real time: %.3f sec; CPU: %.3f sec", realtime() - o->inittime, cputime());
	}
_main_final:;
	mm_opt_destroy(o);
	return(ret);
}

/* end of main.c */
