
/* make sure posix apis are properly activated */
#if defined(__linux__) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE		200112L
#endif
#if defined(__darwin__) && !defined(_DARWIN_C_FULL)
#  define _DARWIN_C_SOURCE	_DARWIN_C_FULL
#endif

/* set non-zero value to ensure the order of output records */
#define STRICT_STREAM_ORDERING		( 1 )

/* collect suppementary alignments, set zero for compatibility with 0.4.x */
#define COLLECT_SUPPLEMENTARY		( 1 )

/* use crc32 for hash64, set zero for compatibility with 0.4.x */
#define USE_CRC32_HASH				( 1 )

/* max #threads, set larger value if needed */
#define MAX_THREADS					( 64 )

/* set zero to disable unittests */
#define UNITTEST 					( 0 )

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

/* global consts */
#ifndef MM_VERSION
#  define MM_VERSION		"minialign-0.5.3-unknown"
#endif

/* max, min */
#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )

/* _likely, _unlikely */
#define _likely(x)		__builtin_expect(!!(x), 1)
#define _unlikely(x)	__builtin_expect(!!(x), 0)

/* timer */
double mm_realtime0;

static double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

static double realtime()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

/* malloc wrappers */
typedef struct {
	uint64_t enabled;
	const char *msg;
	void *_pad[14];
} mm_info_t;
_static_assert(sizeof(mm_info_t) == 128);	// fits in a cache line
mm_info_t info[MAX_THREADS+1] __attribute__(( aligned(64) ));

#define enable_info(t)		{ info[t].enabled = 1; }
#define disable_info(t)		{ info[t].enabled = 0; }
#define set_info(t, x)		{ info[t].msg = (const char *)(x); }
static void oom_abort(const char *name, size_t req)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	fprintf(stderr, "[E::%s] Out of memory. (required: %zu B, maxrss: %ld MB)\n", name, req, r.ru_maxrss);
	for (uint64_t i = 0; i < MAX_THREADS+1; ++i) {
		if (info[i].enabled)
			fprintf(stderr, "[E::%s]  thread %" PRIu64 ": %s\n", name, i, info[i].msg? info[i].msg : "No information available.");
	}
	exit(128);	// 128 reserved for out of memory
}

#define mm_malloc(x) ({ \
	void *_ptr = malloc((size_t)(x)); \
	if (_ptr == NULL) { oom_abort(__func__, (x)); } \
	_ptr; \
})
#define malloc(x)		mm_malloc(x)

#define mm_realloc(x, y) ({ \
	void *_ptr = realloc((void *)(x), (size_t)(y)); \
	if (_ptr == NULL) { oom_abort(__func__, (y)); } \
	_ptr; \
})
#define realloc(x, y)	mm_realloc(x, y)

#define mm_calloc(x, y) ({ \
	void *_ptr = calloc((size_t)(x), (size_t)(y)); \
	if (_ptr == NULL) { oom_abort(__func__, (y)); } \
	_ptr; \
})
#define calloc(x, y)	mm_calloc(x, y)

#define _pnum(type, _buf, _n) ({ \
	uint8_t _b[16] = {0}; \
	type _m = (type)(_n); int64_t _i = 0; \
	while (_m) _b[_i++] = _m % 10, _m /= 10; \
	uint64_t _len = _i + (_i == 0); \
	for (int64_t _j = 0; _j < _len; _j++) { _buf[_j] = _b[_len - _j - 1] + '0'; } \
	_len; \
})

#define _pstr(_buf, _str) ({ \
	const char *_p = (_str); \
	char *_q = (_buf); \
	while ((*_q++ = *_p++)) {} \
	_q - _buf - 1; \
})

/* 64bit random number generator (for unittest) */
uint64_t mm_rand64(void)
{
	uint64_t bits = 31, n = 0;
	for (uint64_t acc = 0; acc < 64; acc += bits) {
		n <<= bits; n ^= (uint64_t)rand();
	}
	return n;
}

/* end of misc.c */

#include "kvec.h"
#include "gaba.h"
#include "lmm.h"
#include "arch/arch.h"

#define UNITTEST_UNIQUE_ID		1
#include "unittest.h"

unittest_config( .name = "minialign" );

/* test for timers and mallocs */
unittest( .name = "misc.timer" ) {
	double cpu = cputime(), real = realtime();	// make sure they do not raise segv
	assert(isnan(cpu) == 0);
	assert(isnan(real) == 0);
}

unittest( .name = "misc.malloc" ) {
	uint64_t const size = 1024 * 1024 * 1024;
	uint8_t *p = malloc(size);
	assert(p != NULL);

	memset(p, 0, size);		// make sure we can touch this area

	p = realloc(p, 2*size);
	assert(p != NULL);

	p = realloc(p, size/2);
	assert(p != NULL);
	free(p);
}
/* end of unittest */

/* minimap.h */
typedef struct { uint32_t x[2]; } v2u32_t;
typedef struct { uint32_t x[4]; } v4u32_t;
typedef struct { uint64_t x[2]; } v2u64_t;
_static_assert(sizeof(v2u32_t) == 8);
_static_assert(sizeof(v4u32_t) == 16);
_static_assert(sizeof(v2u64_t) == 16);

typedef union {
	uint64_t u64[2];
	uint32_t u32[4];
} mm128_t;
_static_assert(sizeof(mm128_t) == 16);

typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; v2u32_t *a; } v2u32_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; uint32_t *a; } uint32_v;
typedef struct { size_t n, m; uint16_t *a; } uint16_v;
typedef struct { size_t n, m; uint8_t *a; } uint8_v;
typedef struct { size_t n, m; void **a; } ptr_v;

#include "ksort.h"
#define sort_key_64x(a) ((a).x[0])
KRADIX_SORT_INIT(64x, v2u32_t, sort_key_64x, 4)
#define sort_key_128x(a) ((a).u64[0])
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8) 
KSORT_INIT_GENERIC(uint32_t)


/* hash.c */
typedef struct kh_s {
	uint32_t mask, max, cnt, ub;		// size of the array, element count, upper bound
	mm128_t *a;
} kh_t;
#define KH_SIZE			( 256 )
#define KH_THRESH		( 0.4 )
#define kh_size(h)		( (h)->mask + 1 )
#define kh_cnt(h)		( (h)->cnt )
#define kh_exist(h, i)	( (h)->a[i].u64[0] + 2 >= 2 )
#define kh_key(h, i)	( (h)->a[i].u64[0] )
#define kh_val(h, i)	( (h)->a[i].u64[1] )

static kh_t *kh_init(uint64_t size)
{
	kh_t *h = calloc(1, sizeof(kh_t));
	// roundup
	size = size < KH_SIZE? KH_SIZE : (0x8000000000000000>>(lzcnt(size - 1) - 1));
	h->mask = size - 1; h->max = size; h->cnt = 0; h->ub = size * KH_THRESH;
	h->a = malloc(sizeof(mm128_t) * size);
	for (uint64_t i = 0; i < size; ++i) h->a[i].u64[0] = UINT64_MAX, h->a[i].u64[1] = 0;
	return h;
}

static void kh_destroy(kh_t *h)
{
	if (h == NULL) { return; }
	free(h->a); free(h);
	return;
}

typedef uint64_t (*khwrite_t)(void *, void *, uint64_t);
static void kh_dump(kh_t *h, void *fp, khwrite_t wfp)
{
	uint32_t x[2] = { 0 };
	if (h == 0) { wfp(fp, x, sizeof(uint32_t) * 2); return; }
	x[0] = h->mask + 1; x[1] = h->cnt;
	wfp(fp, x, sizeof(uint32_t) * 2);
	wfp(fp, h->a, sizeof(mm128_t) * x[0]);
	return;
}

typedef uint64_t (*khread_t)(void *, void *, uint64_t);
static kh_t *kh_load(void *fp, khread_t rfp)
{
	uint32_t x[2] = {0};
	if ((rfp(fp, x, sizeof(uint32_t) * 2)) != sizeof(uint32_t) * 2 || x[0] == 0) return NULL;
	kh_t *h = calloc(1, sizeof(kh_t));
	h->mask = x[0] - 1; h->max = x[0]; h->cnt = x[1]; h->ub = h->max * KH_THRESH;
	h->a = malloc(sizeof(mm128_t) * x[0]);
	if ((rfp(fp, h->a, sizeof(mm128_t) * x[0])) != sizeof(mm128_t) * x[0]) { free(h->a); free(h); return NULL; }
	return h;
}

static void kh_clear(kh_t *h)
{
	if (h == 0) { return; }
	h->mask = KH_SIZE - 1; h->cnt = 0; h->ub = KH_SIZE * KH_THRESH;
	for (uint64_t i = 0; i < KH_SIZE; ++i) h->a[i].u64[0] = UINT64_MAX, h->a[i].u64[1] = 0;
	return;
}

static uint64_t kh_put_intl(mm128_t *a, uint64_t k, uint64_t v, uint64_t mask)
{
	uint64_t kc = k, pos = k & mask;
	while (1) {
		uint64_t kt = a[pos].u64[0];
		if (kc <= kt) {				// k < kt || is_empty(kt) || is_moved(kt)
			uint64_t vt = a[pos].u64[1];
			a[pos] = (mm128_t){ .u64 = { k, v } };
			if (kt + 2 < 2) return 1;// is_empty(kt) || is_moved(kt)
			if (k == kt) return 0;	// k == kt, replaced
			k = kc = kt; v = vt;	// robinhood swap
		}
		kc -= (mask + 1) & (pos + 1); pos = mask & (pos + 1);
	}
	return 0;
}

static void kh_extend(kh_t *h)
{
	uint64_t prev_size = h->mask + 1, size = 2 * prev_size, mask = size - 1;
	h->mask = mask; h->ub = size * KH_THRESH;
	if (size > h->max) { h->a = realloc(h->a, sizeof(mm128_t) * size); h->max = size; }
	for (uint64_t i = 0; i < prev_size; ++i) h->a[i + prev_size].u64[0] = UINT64_MAX, h->a[i + prev_size].u64[1] = 0;
	for (uint64_t i = 0; i < size; ++i) {
		uint64_t k = h->a[i].u64[0];	// key
		if (k + 2 < 2 || (k & mask) == i) continue;
		uint64_t v = h->a[i].u64[1]; h->a[i].u64[0] = UINT64_MAX-1;	// moved
		kh_put_intl(h->a, k, v, mask);
	}
	return;
}

static void kh_put(kh_t *h, uint64_t key, uint64_t val)
{
	if (h->cnt >= h->ub) kh_extend(h);
	h->cnt += kh_put_intl(h->a, key, val, h->mask);
	return;
}

static uint64_t kh_get(kh_t *h, uint64_t key)
{
	uint64_t mask = h->mask, pos = key & mask, k;
	do {
		if ((k = h->a[pos].u64[0]) == key) return h->a[pos].u64[1];
		pos = mask & (pos + 1);
	} while (k + 1 != 0);			// !is_empty(k) || is_moved(k)
	return UINT64_MAX;
}

static const uint64_t *kh_get_ptr(kh_t *h, uint64_t key)
{
	uint64_t mask = h->mask, pos = key & mask, k;
	do {
		if ((k = h->a[pos].u64[0]) == key) return &h->a[pos].u64[1];
		pos = mask & (pos + 1);
	} while (k + 1 != 0);			// !is_empty(k) || is_moved(k)
	return NULL;
}

unittest( .name = "kh.base" ) {
	kh_t *h = kh_init(0);
	assert(h != NULL);

	const uint64_t kmask = mm_rand64(), vmask = mm_rand64(), cnt = 1024 * 1024;
	const uint64_t *p;

	assert(kh_cnt(h) == 0, "cnt(%lu)", kh_cnt(h));
	for (uint64_t i = 0; i < 2*cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p == NULL);
		assert(kh_get(h, i^kmask) == UINT64_MAX);
	}

	// put key-val pairs
	for (uint64_t i = 0; i < cnt; i++) kh_put(h, i^kmask, i^vmask);
	assert(kh_cnt(h) == cnt, "cnt(%lu)", kh_cnt(h));
	for (uint64_t i = 0; i < cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p != NULL, "i(%lu)", i);
		assert(*p == (i^vmask), "i(%lu), val(%lu, %lu)", i, *p, (i^vmask));
		assert(kh_get(h, i^kmask) == (i^vmask), "i(%lu), val(%lu, %lu)", i, *p, (i^vmask));
	}
	for (uint64_t i = cnt; i < 2*cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p == NULL);
		assert(kh_get(h, i^kmask) == UINT64_MAX);
	}

	// clear
	kh_clear(h);
	assert(kh_cnt(h) == 0, "cnt(%lu)", kh_cnt(h));
	for (uint64_t i = 0; i < cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p == NULL);
		assert(kh_get(h, i^kmask) == UINT64_MAX);
	}

	// add another set of key-val pairs
	for (uint64_t i = cnt; i < 2*cnt; i++) kh_put(h, i^kmask, i^vmask);
	assert(kh_cnt(h) == cnt, "cnt(%lu)", kh_cnt(h));
	for (uint64_t i = 0; i < cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p == NULL);
		assert(kh_get(h, i^kmask) == UINT64_MAX);
	}
	for (uint64_t i = cnt; i < 2*cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p != NULL, "i(%lu)", i);
		assert(*p == (i^vmask), "i(%lu), val(%lu, %lu)", i, *p, (i^vmask));
		assert(kh_get(h, i^kmask) == (i^vmask), "i(%lu), val(%lu, %lu)", i, *p, (i^vmask));
	}

	kh_destroy(h);
}

unittest( .name = "kh.io" ) {
	kh_t *h = kh_init(0);
	const uint64_t kmask = mm_rand64(), vmask = mm_rand64(), cnt = 1024 * 1024;
	for (uint64_t i = 0; i < cnt; i++) kh_put(h, i^kmask, i^vmask);

	// dump
	const char *filename = "./minialign.unittest.kh.tmp";
	gzFile fp = gzopen(filename, "w");
	assert((void*)fp != NULL);
	kh_dump(h, (void*)fp, (khwrite_t)gzwrite);
	gzclose(fp);
	kh_destroy(h);

	// restore
	fp = gzopen(filename, "r");
	assert((void*)fp != NULL);
	h = kh_load((void*)fp, (khread_t)gzread);
	assert(h != NULL);
	gzclose(fp);

	const uint64_t *p;
	assert(kh_cnt(h) == cnt, "cnt(%lu)", kh_cnt(h));
	for (uint64_t i = 0; i < cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p != NULL, "i(%lu)", i);
		assert(*p == (i^vmask), "i(%lu), val(%lu, %lu)", i, *p, (i^vmask));
		assert(kh_get(h, i^kmask) == (i^vmask), "i(%lu), val(%lu, %lu)", i, *p, (i^vmask));
	}
	for (uint64_t i = cnt; i < 2*cnt; i++) {
		p = kh_get_ptr(h, i^kmask);
		assert(p == NULL);
		assert(kh_get(h, i^kmask) == UINT64_MAX);
	}
	kh_destroy(h);
	remove(filename);
}
/* end of hash.c */

/* queue.c */
typedef void *(*pt_source_t)(uint32_t tid, void *arg);
typedef void *(*pt_worker_t)(uint32_t tid, void *arg, void *item);
typedef void (*pt_drain_t)(uint32_t tid, void *arg, void *item);

typedef struct pt_q_s {
	uint64_t lock, head, tail, size;
	void **elems;
	uint64_t _pad1[3];
} pt_q_t;
_static_assert(sizeof(pt_q_t) == 64);

typedef struct pt_thread_s {
	pthread_t th;
	uint64_t tid;
	pt_q_t *in, *out;
	volatile pt_worker_t wfp;
	volatile void *warg;
} pt_thread_t;

typedef struct pt_s {
	pt_q_t in, out;
	uint32_t nth;
	pt_thread_t c[];	// [0] is reserved for master
} pt_t;

#define PT_EMPTY	( (void *)(UINT64_MAX) )
#define PT_EXIT		( (void *)(UINT64_MAX-1) )

uint64_t pt_enq(pt_q_t *q, uint64_t tid, void *elem)
{
	uint64_t z, ret = UINT64_MAX;
	do { z = UINT32_MAX; } while (!cas(&q->lock, &z, tid));
	uint64_t head = q->head, tail = q->tail, size = q->size;
	if (((head + 1) % size) != tail) {
		q->elems[head] = elem; q->head = (head + 1) % size; ret = 0;
	}
	do { z = tid; } while (!cas(&q->lock, &z, UINT32_MAX));
	return ret;
}

void *pt_deq(pt_q_t *q, uint64_t tid)
{
	void *elem = PT_EMPTY;
	uint64_t z;
	do { z = UINT32_MAX; } while (!cas(&q->lock, &z, tid));
	uint64_t head = q->head, tail = q->tail, size = q->size;
	if (head != tail) {
		elem = q->elems[tail]; q->tail = (tail + 1) % size;
	}
	do { z = tid; } while (!cas(&q->lock, &z, UINT32_MAX));
	return elem;
}

static void *pt_dispatch(void *s)
{
	pt_thread_t *c = (pt_thread_t *)s;
	struct timespec tv = { .tv_nsec = 512 * 1024 };
	void *ping = PT_EMPTY, *pong = PT_EMPTY;
	while (1) {
		ping = pt_deq(c->in, c->tid);
		if (ping == PT_EMPTY && pong == PT_EMPTY) { nanosleep(&tv, NULL); }
		if (pong != PT_EMPTY) pt_enq(c->out, c->tid, c->wfp(c->tid, (void*)c->warg, pong));
		if (ping == PT_EXIT) break;
		pong = pt_deq(c->in, c->tid);
		if (ping == PT_EMPTY && pong == PT_EMPTY) { nanosleep(&tv, NULL); }
		if (ping != PT_EMPTY) pt_enq(c->out, c->tid, c->wfp(c->tid, (void*)c->warg, ping));
		if (pong == PT_EXIT) break;
	}
	return NULL;
}

static void pt_destroy(pt_t *pt)
{
	void *status;
	for (uint64_t i = 1; i < pt->nth; ++i) pt_enq(pt->c->in, pt->c->tid, PT_EXIT);
	for (uint64_t i = 1; i < pt->nth; ++i) { disable_info(i); pthread_join(pt->c[i].th, &status); }
	while (pt_deq(pt->c->in, pt->c->tid) != PT_EMPTY) {}
	while (pt_deq(pt->c->out, pt->c->tid) != PT_EMPTY) {}
	free(pt->in.elems); free(pt->out.elems); free(pt);
}

static pt_t *pt_init(uint32_t nth)
{
	nth = (nth == 0)? 1 : nth;
	pt_t *pt = calloc(1, sizeof(pt_t) + sizeof(pt_thread_t) * nth);
	pt->in = (pt_q_t){ .lock = UINT32_MAX, .elems = calloc(8 * nth, sizeof(void*)), .size = 8 * nth };
	pt->out = (pt_q_t){ .lock = UINT32_MAX, .elems = calloc(8 * nth, sizeof(void*)), .size = 8 * nth };

	pt->nth = nth; pt->c[0].tid = 0; pt->c[0].in = &pt->in; pt->c[0].out = &pt->out;
	for (uint64_t i = 1; i < nth; ++i) {
		pt->c[i].tid = i; pt->c[i].in = &pt->in; pt->c[i].out = &pt->out;
		pthread_create(&pt->c[i].th, NULL, pt_dispatch, (void *)&pt->c[i]);
		enable_info(i);
	}
	return pt;
}

static int pt_set_worker(pt_t *pt, pt_worker_t wfp, void **warg)
{
	void *item;
	if ((item = pt_deq(pt->c->in, pt->c->tid)) != PT_EMPTY) { pt_enq(pt->c->in, pt->c->tid, item); return -1; }
	for (uint64_t i = 0; i < pt->nth; ++i) pt->c[i].wfp = wfp, pt->c[i].warg = warg? warg[i] : NULL;
	fence();
	return 0;
}

// sfp and dfp are called in the same thread
static int pt_stream(pt_t *pt, pt_source_t sfp, void *sarg, pt_worker_t wfp, void **warg, pt_drain_t dfp, void *darg)
{
	if (pt_set_worker(pt, wfp, warg)) return -1;
	uint64_t bal = 0, lb = 2 * pt->nth, ub = 4 * pt->nth;
	void *item;
	while ((item = sfp(pt->c->tid, sarg)) != NULL) {
		pt_enq(pt->c->in, pt->c->tid, item);
		if (++bal < ub) continue;
		while (bal > lb) {
			if ((item = pt_deq(pt->c->out, pt->c->tid)) != PT_EMPTY) dfp(pt->c->tid, darg, item), bal--;
			if ((item = pt_deq(pt->c->in, pt->c->tid)) != PT_EMPTY) dfp(pt->c->tid, darg, wfp(pt->c->tid, warg? *warg : NULL, item)), bal--;
		}
	}
	while ((item = pt_deq(pt->c->in, pt->c->tid)) != PT_EMPTY) pt_enq(pt->c->out, pt->c->tid, wfp(pt->c->tid, warg? *warg : NULL, item));
	while (bal > 0) if ((item = pt_deq(pt->c->out, pt->c->tid)) != PT_EMPTY) dfp(pt->c->tid, darg, item), bal--;
	return 0;
}

static int pt_parallel(pt_t *pt, pt_worker_t wfp, void **warg, void **item)
{
	if (pt_set_worker(pt, wfp, warg)) return -1;
	for (uint64_t i = 1; i < pt->nth; ++i) pt_enq(pt->c->in, pt->c->tid, item? item[i] : NULL);
	wfp(pt->c->tid, warg? warg[0] : NULL, item? item[0] : NULL);
	for (uint64_t i = 1; i < pt->nth; ++i) {
		while (pt_deq(pt->c->out, pt->c->tid) == PT_EMPTY) sched_yield();
	}
	return 0;
}

static void *pt_unittest_source(uint32_t tid, void *arg)
{
	uint64_t *s = (uint64_t*)arg;
	if (*s >= 1024) return NULL;
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
	pt_stream(pt, pt_unittest_source, (void *)&icnt, pt_unittest_worker, (void**)arr, pt_unittest_drain, (void*)&ocnt);
	assert(icnt == 1024, "icnt(%lu)", icnt);
	assert(ocnt == 512*1025, "ocnt(%lu)", ocnt);

	uint64_t s[4] = { 0 }, *sp[1] = { &s[0] };
	pt_parallel(pt, pt_unittest_worker, (void**)arr, (void**)sp);
	assert(s[0] == 1, "d[0](%lu)", s[0]);
	pt_destroy(pt);
}

unittest( .name = "pt.multi" ) {
	pt_t *pt = pt_init(4);
	assert(pt != NULL);

	uint64_t icnt = 0, ocnt = 0, inc = 1, *arr[4] = { &inc, &inc, &inc, &inc };
	pt_stream(pt, pt_unittest_source, (void *)&icnt, pt_unittest_worker, (void**)arr, pt_unittest_drain, (void*)&ocnt);
	assert(icnt == 1024, "icnt(%lu)", icnt);
	assert(ocnt == 512*1025, "ocnt(%lu)", ocnt);

	uint64_t s[4] = { 0,1,2,3 }, *sp[4] = { &s[0],&s[1],&s[2],&s[3] };
	pt_parallel(pt, pt_unittest_worker, (void**)arr, (void**)sp);
	for (uint64_t i = 0; i < 4; i++) {
		assert(s[i] == (i+1), "i(%lu), d[i](%lu)", i, s[i]);
	}
	pt_destroy(pt);
}

/* stdio stream with multithreaded compression / decompression */
typedef struct {
	uint64_t head, len;
	uint32_t id, raw, flush, pad;
	uint8_t buf[];
} pg_block_t;

typedef struct {
	FILE *fp;
	pt_t *pt;
	pg_block_t *s;
	uint32_t ub, lb, bal, icnt, ocnt, eof, nth;
	uint64_t block_size;
	kvec_t(mm128_t) hq;
	void *c[];
} pg_t;
#define incq_comp(a, b)		( (int64_t)(a).u64[0] - (int64_t)(b).u64[0] )

static pg_block_t *pg_deflate(pg_block_t *in, uint64_t block_size)
{
	uint64_t buf_size = block_size * 1.2;
	pg_block_t *out = malloc(sizeof(pg_block_t) + buf_size);
	z_stream zs = {
		.next_in = in->buf, .avail_in = in->len,
		.next_out = out->buf, .avail_out = buf_size
	};
	deflateInit2(&zs, 1, Z_DEFLATED, 15, 8, Z_DEFAULT_STRATEGY);
	deflate(&zs, Z_FINISH);
	deflateEnd(&zs);
	out->head = 0; out->len = buf_size - zs.avail_out;
	out->id = in->id; out->raw = 0; out->flush = 1;
	free(in);
	return out;
}

static pg_block_t *pg_inflate(pg_block_t *in, uint64_t block_size)
{
	uint64_t buf_size = block_size * 1.2;
	pg_block_t *out = malloc(sizeof(pg_block_t) + buf_size);
	z_stream zs = {
		.next_in = in->buf, .avail_in = in->len,
		.next_out = out->buf, .avail_out = buf_size
	};
	inflateInit2(&zs, 15);
	inflate(&zs, Z_FINISH);
	inflateEnd(&zs);
	out->head = 0; out->len = buf_size - zs.avail_out;
	out->id = in->id; out->raw = 1; out->flush = 0;
	free(in);
	return out;
}

static void *pg_worker(uint32_t tid, void *arg, void *item)
{
	pg_t *pg = (pg_t*)arg;
	pg_block_t *s = (pg_block_t*)item;
	if (s == NULL || s->len == 0) return s;

	char buf[128], *p = buf;
	p += _pstr(p, "[pg_worker] bucket id: "); p += _pnum(uint32_t, p, s->id); *p = '\0';
	set_info(tid, buf);
	return (s->raw? pg_deflate : pg_inflate)(s, pg->block_size);
}

static pg_block_t *pg_read_block(pg_t *pg)
{
	pg_block_t *s = malloc(sizeof(pg_block_t) + pg->block_size);
	if (fread(&s->len, sizeof(uint64_t), 1, pg->fp) != 1 || s->len == 0) { goto _fail; }
	if (fread(s->buf, sizeof(uint8_t), s->len, pg->fp) != s->len) { goto _fail; }
	s->head = 0; s->id = pg->icnt++; s->raw = 0; s->flush = 0;
	return s;
_fail:
	free(s);
	return NULL;
}

static void pg_write_block(pg_t *pg, pg_block_t *s)
{
	if (s->len == 0) return;
	fwrite(&s->len, sizeof(uint64_t), 1, pg->fp);
	fwrite(s->buf, sizeof(uint8_t), s->len, pg->fp);
	free(s);
	return;
}

static pg_t *pg_init(FILE *fp, uint32_t nth)
{
	pg_t *pg = calloc(1, sizeof(pg_t) + nth * sizeof(pg_t*));
	*pg = (pg_t){
		.fp = fp,
		.pt = pt_init(nth),
		.lb = nth, .ub = 3 * nth, .bal = 0, .nth = nth,
		.block_size = 1024 * 1024
	};
	kv_hq_init(pg->hq);
	for (uint64_t i = 0; i < nth; ++i) pg->c[i] = (void*)pg;
	pt_set_worker(pg->pt, pg_worker, pg->c);
	return pg;
}

static void pg_destroy(pg_t *pg)
{
	pg_block_t *s = pg->s, *t;
	if (s && s->flush == 1 && s->head != 0) {
		s->len = s->head;
		if (pg->nth == 1) { pg_write_block(pg, pg_deflate(s, pg->block_size)); }
		else { pt_enq(pg->pt->c->in, pg->pt->c->tid, s); pg->bal++; }
	} else free(s);
	while (pg->bal > 0) {
		if ((t = pt_deq(pg->pt->c->out, pg->pt->c->tid)) == PT_EMPTY) { sched_yield(); continue; }
		pg->bal--;
		if (t && !t->flush) { free(t); continue; }
		kv_hq_push(mm128_t, incq_comp, pg->hq, ((mm128_t){.u64 = {t->id, (uintptr_t)t}}));
	}
	while (pg->hq.n > 1) {
		pg->ocnt++;
		t = (pg_block_t*)kv_hq_pop(mm128_t, incq_comp, pg->hq).u64[1];
		pg_write_block(pg, t);
	}
	uint64_t z = 0;
	fwrite(&z, sizeof(uint64_t), 1, pg->fp);	/* terminator */
	pt_destroy(pg->pt); kv_hq_destroy(pg->hq); free(pg);
	return;
}

static uint64_t pgread(pg_t *pg, void *dst, uint64_t len)
{
	uint64_t rem = len;
	pg_block_t *s = pg->s, *t;
	if (pg->eof == 2) return 0;
	while (rem > 0) {
		while (!s || s->head == s->len) {
			free(s); s = NULL;
			if (pg->nth == 1) {
				if ((t = pg_read_block(pg)) == NULL) { pg->eof = 2; return len-rem; }
				pg->s = s = pg_inflate(t, pg->block_size);
			} else {
				while (!pg->eof && pg->bal < pg->ub) {
					if ((t = pg_read_block(pg)) == NULL) { pg->eof = 1; break; }
					pg->bal++;
					pt_enq(pg->pt->c->in, pg->pt->c->tid, t);
				}
				while ((t = pt_deq(pg->pt->c->out, pg->pt->c->tid)) != PT_EMPTY) {
					pg->bal--;
					kv_hq_push(mm128_t, incq_comp, pg->hq, ((mm128_t){.u64 = {t->id, (uintptr_t)t}}));
				}
				if (pg->ocnt >= pg->icnt) { pg->eof = 2; return len-rem; }
				if (pg->hq.n < 2 || pg->hq.a[1].u64[0] > pg->ocnt) { sched_yield(); continue; }
				pg->ocnt++;
				pg->s = s = (pg_block_t*)kv_hq_pop(mm128_t, incq_comp, pg->hq).u64[1];
			}
		}
		uint64_t adv = MIN2(rem, s->len-s->head);
		memcpy(dst + len - rem, s->buf + s->head, adv);
		rem -= adv; s->head += adv;
	}
	return len;
}

static uint64_t pgwrite(pg_t *pg, const void *src, uint64_t len)
{
	uint64_t rem = len;
	pg_block_t *s = pg->s, *t;
	while (rem > 0) {
		if (!s || s->head == s->len) {
			if (pg->nth == 1) {
				if (s) pg_write_block(pg, pg_deflate(s, pg->block_size));
			} else {
				while (pg->bal > pg->lb) {
					if ((t = pt_deq(pg->pt->c->out, pg->pt->c->tid)) == PT_EMPTY) {
						if (pg->bal >= pg->ub) { sched_yield(); continue; } else break;
					}
					pg->bal--;
					kv_hq_push(mm128_t, incq_comp, pg->hq, ((mm128_t){.u64 = {t->id, (uintptr_t)t}}));
				}
				while (pg->hq.n > 1 && pg->hq.a[1].u64[0] <= pg->ocnt) {
					pg->ocnt++;
					t = (pg_block_t*)kv_hq_pop(mm128_t, incq_comp, pg->hq).u64[1];
					pg_write_block(pg, t);
				}
				if (s) {
					pg->bal++;
					pt_enq(pg->pt->c->in, pg->pt->c->tid, s);
				}
			}
			s = malloc(sizeof(pg_block_t) + pg->block_size);
			s->head = 0; s->len = pg->block_size; s->id = pg->icnt++; s->raw = 1; s->flush = 1;
		}
		uint64_t adv = MIN2(rem, s->len-s->head);
		memcpy(s->buf + s->head, src + len - rem, adv);
		rem -= adv; s->head += adv;
	}
	pg->s = s;
	return len;
}

unittest( .name = "pg.single" ) {
	const uint64_t size = 1024 * 1024 * 1024;
	uint8_t *p = malloc(size), *q = malloc(size);
	for (uint64_t i = 0; i < size; i++) p[i] = i % 253;

	const char *filename = "./minialign.unittest.pg.tmp";
	FILE *fp = fopen(filename, "w");
	assert(fp != NULL);
	pg_t *pg = pg_init(fp, 1);
	assert(pg != NULL);

	for (uint64_t i = 0; i < 3; i++) {
		uint64_t l = pgwrite(pg, p, size);
		assert(l == size, "l(%lu), size(%lu)", l, size);
	}
	pg_destroy(pg);
	fclose(fp);


	fp = fopen(filename, "r");
	assert(fp != NULL);
	pg = pg_init(fp, 1);
	assert(pg != NULL);

	for (uint64_t i = 0; i < 3; i++) {
		uint64_t l = pgread(pg, q, size);
		assert(l == size, "l(%lu), size(%lu)", l, size);
		assert(memcmp(p, q, size) == 0);
	}
	pg_destroy(pg);
	fclose(fp);
	free(p); free(q); remove(filename);
}

unittest( .name = "pg.multi" ) {
	const uint64_t size = 1024 * 1024 * 1024;
	uint8_t *p = malloc(size), *q = malloc(size);
	for (uint64_t i = 0; i < size; i++) p[i] = i % 253;

	const char *filename = "./minialign.unittest.pg.tmp";
	FILE *fp = fopen(filename, "w");
	assert(fp != NULL);
	pg_t *pg = pg_init(fp, 4);
	assert(pg != NULL);

	for (uint64_t i = 0; i < 3; i++) {
		uint64_t l = pgwrite(pg, p, size);
		assert(l == size, "l(%lu), size(%lu)", l, size);
	}
	pg_destroy(pg);
	fclose(fp);


	fp = fopen(filename, "r");
	assert(fp != NULL);
	pg = pg_init(fp, 4);
	assert(pg != NULL);

	for (uint64_t i = 0; i < 3; i++) {
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
typedef struct {
	int32_t n_targets;
	char **target_name;
	uint32_t *target_len;
	size_t l_text, n_text;
	char *text;
} bam_header_t;

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

static void bam_header_destroy(bam_header_t *header)
{
	if (header == 0) return;
	if (header->target_name) {
		for (uint64_t i = 0; i < (uint64_t)header->n_targets; ++i)
			if (header->target_name[i]) free(header->target_name[i]);
		if (header->target_len) free(header->target_len);
		free(header->target_name);
	}
	if (header->text) free(header->text);
	free(header);
}

static bam_header_t *bam_read_header(gzFile fp)
{
	bam_header_t *header;
	char buf[4];
	int magic_len;
	// read "BAM1"
	magic_len = gzread(fp, buf, 4);
	if (magic_len != 4 || strncmp(buf, "BAM\001", 4) != 0) return NULL;
	header = calloc(1, sizeof(bam_header_t));
	// read plain text and the number of reference sequences
	if (gzread(fp, &header->l_text, 4) != 4) goto _bam_read_header_fail;
	header->text = (char*)calloc(header->l_text + 1, 1);
	if ((size_t)gzread(fp, header->text, header->l_text) != header->l_text) goto _bam_read_header_fail;
	if (gzread(fp, &header->n_targets, 4) != 4) goto _bam_read_header_fail;
	// read reference sequence names and lengths
	header->target_name = (char**)calloc(header->n_targets, sizeof(char*));
	header->target_len = (uint32_t*)calloc(header->n_targets, 4);
	for (uint64_t i = 0; i != (uint64_t)header->n_targets; ++i) {
		int32_t name_len;
		if (gzread(fp, &name_len, 4) != 4) goto _bam_read_header_fail;
		header->target_name[i] = (char*)calloc(name_len, 1);
		if (gzread(fp, header->target_name[i], name_len) != name_len) {
			goto _bam_read_header_fail;
		}
		if (gzread(fp, &header->target_len[i], 4) != 4) goto _bam_read_header_fail;
	}
	return header;
_bam_read_header_fail:
	bam_header_destroy(header);
	return NULL;
}

/* end of bamlite.c */

/* bseq.c */
#define BSEQ_MGN			( 64 )
typedef struct {
	gzFile fp;
	bam_header_t *bh;
	uint8_t *p, *base, *tail;
	uint64_t size, acc;	// equal to batch_size
	uint16_t *tags;
	uint32_t l_tags, n_seq, min_len;
	uint8_t is_eof, delim, keep_qual, keep_comment, state;
} bseq_file_t;

// sam optional tags
typedef struct {
	uint32_t size;	// size of data array
	char tag[2], type;
	uint8_t flag;	// 1 if only contained in primary
	uint8_t data[];
} bseq_tag_t;

typedef struct {
	// lengths
	uint32_t l_seq, l_name, l_tag, _pad;
	// pointers
	char *name;
	uint8_t *seq, *qual, *tag;
} bseq_t;
_static_assert(sizeof(bseq_t) == 48);
typedef struct { size_t n, m; bseq_t *a; } bseq_v;

static uint64_t bseq_search_tag(uint32_t l_tags, const uint16_t *tags, uint16_t t1, uint16_t t2)
{
	v32i16_t tv = _set_v32i16((t2<<8) | t1);
	for (uint64_t i = 0; i < ((l_tags + 0x1f) & ~0x1f); i+=0x20) {
		if (((v32_masku_t){ .mask = _mask_v32i16(_eq_v32i16(_loadu_v32i16(&tags[i]), tv)) }).all != 0) return 1;
	}
	return 0;
}

static bseq_file_t *bseq_open(const char *fn, uint64_t batch_size, uint32_t keep_qual, uint32_t min_len, uint32_t l_tags, const uint16_t *tags)
{
	int c;
	bseq_file_t *fp;
	gzFile f;

	set_info(0, "[bseq_open] initialize bseq object");
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (f == NULL) return NULL;
	fp = (bseq_file_t*)calloc(1, sizeof(bseq_file_t));
	fp->n_seq = 0; fp->keep_qual = keep_qual; fp->min_len = min_len;
	fp->fp = f; fp->bh = NULL; fp->delim = 0;
	for (uint64_t i = 0; i < 4; i++) {	// allow some invalid spaces at the head
		if ((c = gzgetc(fp->fp)) == 'B') {	// test bam signature
			gzungetc(c, fp->fp); fp->bh = bam_read_header(fp->fp); break;
		} else if (c == '>' || c == '@') {
			gzungetc(c, fp->fp); fp->delim = c; break;
		}
	}
	if (!fp->bh && !fp->delim) { free(fp); return NULL; }

	// init buffer
	fp->tail = fp->p = (fp->base = malloc((fp->size = batch_size) + 2*BSEQ_MGN)) + batch_size;
	memset(fp->base, 0, BSEQ_MGN);
	memset(fp->base + batch_size + BSEQ_MGN, 0, BSEQ_MGN);

	fp->tags = calloc(((fp->l_tags = l_tags) + 0x1f) & ~0x1f, sizeof(uint16_t));
	if (l_tags && tags) memcpy(fp->tags, tags, l_tags * sizeof(uint16_t));
	fp->keep_comment = bseq_search_tag(fp->l_tags, fp->tags, 'C', 'O');
	return fp;
}

static uint32_t bseq_close(bseq_file_t *fp)
{
	if (fp == NULL) return 0;
	uint32_t n_seq = fp->n_seq;
	gzclose(fp->fp); bam_header_destroy(fp->bh); free(fp->base); free(fp->tags); free(fp);
	return n_seq;
}

static uint64_t bseq_save_tags(uint32_t l_tags, const uint16_t *tags, uint32_t l_arr, uint8_t *arr, uint8_t *q)
{
	static const uint8_t tag_size[32] = {
		0, 1, 0xfe, 1,  0, 0, 4, 0,  0xfe, 4, 0, 0,  0, 0, 0, 0,
		0, 0, 0, 2,     0, 0, 0, 0,  0, 0, 0xff, 0,  0, 0, 0, 0,
	};
	const uint8_t *p = arr, *tail = arr + l_arr;
	while (p < tail) {
		uint64_t keep = bseq_search_tag(l_tags, tags, p[0], p[1]), size;
		if ((size = tag_size[p[2]&0x1f]) == 0xff) {
			if (keep) while ((*q++ = *p++)) {}
			else while (*p++) {}
		} else {
			uint64_t len =  (size == 0xfe)? 8 + tag_size[p[3]&0x1f] * *((uint32_t*)&p[4]) : 3 + size;
			if (keep) { memcpy(q, p, len); q += len; }
			p += len;
		}
	}
	return p - arr;
}

static uint64_t bseq_read_bam(bseq_file_t *fp, uint64_t size, bseq_v *seq, uint8_v *mem)
{
	// parse bam
	bseq_t *s;
	uint8_t *sname, *sseq, *squal, *stag;
	uint32_t l_tag;
	bam_core_t *c;

	if ((c = (bam_core_t*)fp->p)->flag&0x900) return 0;	// skip supp/secondary
	if ((uint32_t)c->l_qseq < fp->min_len) return 0;	// skip short reads
	sname = fp->p + sizeof(bam_core_t);
	sseq = sname + c->l_qname + sizeof(uint32_t)*c->n_cigar;
	squal = sseq + (c->l_qseq+1) / 2;
	stag = squal + c->l_qseq; l_tag = ((uint8_t*)fp->p + size) - stag;

	kv_pushp(bseq_t, *seq, &s);
	kv_reserve(uint8_t, *mem, mem->n + c->l_qname + c->l_qseq + ((fp->keep_qual && *squal != 0xff)? c->l_qseq : 0) + (fp->l_tags? l_tag : 0) + 3);
	s->l_seq = c->l_qseq;
	s->l_name = c->l_qname - 1;		// remove tail '\0'
	s->l_tag = fp->l_tags? l_tag : 0;

	s->name = (char *)mem->n;
	memcpy(mem->a + mem->n, sname, c->l_qname); mem->n += c->l_qname;

	s->seq = (uint8_t *)mem->n;
	if (c->flag&0x10) {
		for (int64_t i = c->l_qseq-1; i >= 0; --i) mem->a[mem->n++] = "\x0\x8\x4\x0\x2\x0\x0\x0\x1"[0x0f & (sseq[i>>1]>>((~i&0x01)<<2))];
	} else {
		for (int64_t i = 0; i < c->l_qseq; ++i) mem->a[mem->n++] = 0x0f & (sseq[i>>1]>>((~i&0x01)<<2));
	}
	mem->a[mem->n++] = '\0';

	s->qual = (uint8_t *)mem->n;
	if (fp->keep_qual && *squal != 0xff) {
		if (c->flag&0x10) { for (int64_t i = c->l_qseq-1; i >= 0; --i) mem->a[mem->n++] = squal[i] + 33; }
		else { for (int64_t i = 0; i < c->l_qseq; ++i) mem->a[mem->n++] = squal[i] + 33; }
	}
	mem->a[mem->n++] = '\0';

	s->tag = (uint8_t *)mem->n;
	if (fp->l_tags && l_tag) mem->n += bseq_save_tags(fp->l_tags, fp->tags, l_tag, stag, mem->a + mem->n);
	mem->a[mem->n++] = '\0';
	return 0;
}

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
	const v32i8_t _lv = _set_v32i8('\n'); \
	do { \
		v32i8_t _r = _loadu_v32i8(_p), _s = _op(_r); _storeu_v32i8(_q, _s); \
		_m1 = _match(_r, _dv); _m2 = _match(_r, _lv); \
		ZCNT_RESULT uint64_t _l = MIN2(tzcnt(_m1 | _m2), _t - _p); \
		_len = _l; _p += 32; _q += 32; \
	} while (_len >= 32); \
	_p += _len - 32; _q += _len - 32; _q -= _q[-1] == 0x0f; _m1>>_len; \
})
#define _skipline(_p, _t) ({ \
	uint64_t _m; \
	uint64_t _len; \
	const uint8_t *_b = _p; \
	const v32i8_t _lv = _set_v32i8('\n'); \
	do { \
		v32i8_t _r = _loadu_v32i8(_p); \
		_m = _match(_r, _lv); \
		ZCNT_RESULT uint64_t _l = MIN2(tzcnt(_m), _t - _p); _len = _l; _p += 32; \
	} while (_len >= 32); \
	_p += _len - 32; _p - _b; \
})
#define _beg(_q, _b)		( (uint8_t *)(_q - _b) )
#define _term(_q, _b, _ofs) ({ \
	uint64_t _len = (uint64_t)(_q - &(_b)[(uint64_t)(_ofs)]); \
	*_q++ = '\0'; _len; \
})

static uint64_t bseq_read_sam(bseq_file_t *fp, bseq_v *seq, uint8_v *mem)
{
	return 0;
}

// returns 0 when correctly finished, 1 when buffer starved, >=2 when broken
static uint64_t bseq_read_fasta(bseq_file_t *fp, bseq_v *seq, uint8_v *mem)
{
	// mem must have enough space (e.g. 2 * buffer)
	// keep them on registers
	#define _id(x)			(x)
	#define _trans(x)		( _shuf_v32i8(cv, _and_v32i8(fv, x)) )
	const v32i8_t dv = _set_v32i8(fp->delim == '@'? '+' : fp->delim);
	const v32i8_t sv = _set_v32i8(' '), lv = _set_v32i8('\n'), fv = _set_v32i8(0xf);
	const v32i8_t cv = _from_v16i8_v32i8(_seta_v16i8(0,0,15,0,0,0,0,0,4,0,0,8,2,0,1,0));

	bseq_t *s = &seq->a[seq->n-1];
	uint8_t *p = fp->p, *q = &mem->a[mem->n];
	const uint8_t *t = fp->tail;
	uint64_t m, acc, lim, ret = 1;
	if (p >= t) return 1;	// refill needed
	switch (fp->state) {
		case 0:				// idle
			if (*p++ != fp->delim) return 2;	// broken
			kv_pushp(bseq_t, *seq, &s);			// create new sequence
			fp->state = 1;						// transition to spaces between delim and name
		case 1:
			_strip(p, t, sv);
			if (p >= t) goto _refill;
			s->name = (char*)_beg(q, mem->a);
			fp->state = 2;
		case 2:
			m = _readline(p, t, q, sv, _id);
			if (p >= t) goto _refill;
			p++;								// skip '\n' or ' '
			s->l_name = _term(q, mem->a, s->name);
			s->tag = _beg(q, mem->a);
			if ((m & 0x01) == 0) goto _seq_head;
			fp->state = 3;
		case 3:
			_strip(p, t, sv);
			if (p >= t) goto _refill;
			*q++ = 'C'; *q++ = 'O'; *q++ = 'Z';
			fp->state = 4;
		case 4:									// parsing comment
			_readline(p, t, q, lv, _id);
			if (p >= t) goto _refill;			// refill needed, comment continues
			p++;								// skip '\n'
			while (*--q == ' ') {} q++;			// strip spaces
			if (!fp->keep_comment) q = mem->a + (uint64_t)s->tag;
		_seq_head:
			s->l_tag = _term(q, mem->a, s->tag);
			s->seq = _beg(q, mem->a);
			fp->state = 5;
		case 5:									// parsing seq
			while (1) {
				m = _readline(p, t, q, dv, _trans);
				if (p >= t) { m |= fp->is_eof; break; }
				if (m & 0x01) break;
				p++;							// skip '\n'
			}
			if ((m & 0x01) == 0) goto _refill;
			s->l_seq = _term(q, mem->a, s->seq);
			if (p >= t) break;
			s->qual = _beg(q, mem->a);
			if (fp->delim == '>') goto _qual_tail;
			fp->state = 6;
		case 6:
			_skipline(p, t);
			if (p >= t) goto _refill;
			p++;
			fp->state = 7; fp->acc = 0;
		case 7:									// parsing qual
			acc = fp->acc, lim = s->l_seq;
			if (fp->keep_qual) {
				while (1) {
					const uint8_t *b = q; _readline(p, t, q, lv, _id); acc += q - b;
					if (p >= t) { fp->acc = acc; goto _refill; }
					if (acc >= lim) break;
					p++;						// skip '\n'
				}
			} else {
				while (1) {
					acc += _skipline(p, t);
					if (p >= t) { fp->acc = acc; goto _refill; }
					if (acc >= lim) break;
					p++;						// skip '\n'
				}
			}
			fp->state = 8;
		case 8:
			if (p >= t) goto _refill;
			_strip(p, t, lv);
		_qual_tail:
			_term(q, mem->a, s->qual);
			break;
		default:			// invalid state
			return 2;		// broken
	}
	ret = fp->is_eof; fp->state = 0;	// back to idle
	if ((uint32_t)s->l_seq < fp->min_len) seq->n--, q = mem->a + (uint64_t)s->name;
_refill:
	fp->p = p; mem->n = q - mem->a;
	return ret;

	#undef _id
	#undef _trans
}
#undef _match
#undef _strip
#undef _readline
#undef _skipline
#undef _beg
#undef _term

static bseq_t *bseq_read(bseq_file_t *fp, uint32_t *n, void **base, uint64_t *size)
{
	uint8_v mem = {0};
	bseq_v seq = {0};
	static const uint8_t margin[64] = {0};

	set_info(0, "[bseq_read] read sequence block from file");
	if (fp->is_eof && fp->p >= fp->tail) return NULL;
	kv_reserve(uint8_t, mem, 2*fp->size);
	kv_pushm(uint8_t, mem, margin, 64);
	if (fp->bh) {
		while (mem.n < fp->size) {
			uint32_t size;
			if (gzread(fp->fp, &size, 4) != 4) { fp->is_eof = 1; break; }
			if (fp->size < size) fp->tail = (fp->base = fp->p = realloc(fp->base, fp->size = size)) + size;
			if (gzread(fp->fp, fp->p, size) != size) { seq.n = 0; fp->is_eof = 2; break; }
			bseq_read_bam(fp, size, &seq, &mem);
		}
	} else {
		while (mem.n < fp->size + 64) {
			uint64_t ret;
			while ((ret = bseq_read_fasta(fp, &seq, &mem)) == 1) {
				if (fp->is_eof) goto _tail;							// correctly reached the end
				fp->p = fp->base + BSEQ_MGN;
				fp->tail = fp->p + gzread(fp->fp, fp->p, fp->size);
				memset(fp->tail, 0, BSEQ_MGN);
				if (fp->tail < fp->base + BSEQ_MGN + fp->size) fp->is_eof = 1;
				if (mem.n + 2*fp->size > mem.m) mem.a = realloc(mem.a, mem.m *= 2);
			}
			if (ret > 1) { seq.n = 0; fp->is_eof = 2; goto _tail; }	// error occurred
		}
	_tail:;
	}

	kv_pushm(uint8_t, mem, margin, 64);
	for (uint64_t i = 0; i < seq.n; i++) {
		seq.a[i].name += (ptrdiff_t)mem.a;
		seq.a[i].seq += (ptrdiff_t)mem.a;
		seq.a[i].qual += (ptrdiff_t)mem.a;
		seq.a[i].tag += (ptrdiff_t)mem.a;
	}
	if (seq.n == 0) free(mem.a), mem.a = 0, fp->is_eof = 1;
	fp->n_seq += seq.n;
	*n = seq.n;
	*base = (void*)mem.a;
	*size = mem.n;
	return seq.a;
}

unittest( .name = "bseq.fasta" ) {
	const char *filename = "./minialign.unittest.bseq.tmp";
	const char *content =
		">test0\nAAAA\n"
		"> test1\nATAT\nCGCG\n\n"
		">  test2\n\nAAAA\n"
		">test3 comment comment  \nACGT\n\n";

	FILE *fp = fopen(filename, "w");
	assert(fp != NULL);
	fwrite(content, 1, strlen(content), fp);
	fclose(fp);

	const uint16_t tags[1] = { ((uint16_t)'C') | (((uint16_t)'O')<<8) };
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
	assert(s[0].l_tag == 0, "l_tag(%u)", s[0].l_tag);
	assert(strcmp((const char*)s[0].name, "test0") == 0, "name(%s)", s[0].name);
	assert(strcmp((const char*)s[0].seq, "\x1\x1\x1\x1") == 0, "seq(%s)", s[0].seq);
	assert(strcmp((const char*)s[0].qual, "") == 0, "qual(%s)", s[0].qual);
	assert(strcmp((const char*)s[0].tag, "") == 0, "tag(%s)", s[0].tag);

	assert(s[1].l_name == 5, "l_name(%u)", s[1].l_name);
	assert(s[1].l_seq == 8, "l_seq(%u)", s[1].l_seq);
	assert(s[1].l_tag == 0, "l_tag(%u)", s[1].l_tag);
	assert(strcmp((const char*)s[1].name, "test1") == 0, "name(%s)", s[1].name);
	assert(strcmp((const char*)s[1].seq, "\x1\x8\x1\x8\x2\x4\x2\x4") == 0, "seq(%s)", s[1].seq);
	assert(strcmp((const char*)s[1].qual, "") == 0, "qual(%s)", s[1].qual);
	assert(strcmp((const char*)s[1].tag, "") == 0, "tag(%s)", s[1].tag);

	assert(s[2].l_name == 5, "l_name(%u)", s[2].l_name);
	assert(s[2].l_seq == 4, "l_seq(%u)", s[2].l_seq);
	assert(s[2].l_tag == 0, "l_tag(%u)", s[2].l_tag);
	assert(strcmp((const char*)s[2].name, "test2") == 0, "name(%s)", s[2].name);
	assert(strcmp((const char*)s[2].seq, "\x1\x1\x1\x1") == 0, "seq(%s)", s[2].seq);
	assert(strcmp((const char*)s[2].qual, "") == 0, "qual(%s)", s[2].qual);
	assert(strcmp((const char*)s[2].tag, "") == 0, "tag(%s)", s[2].tag);

	assert(s[3].l_name == 5, "l_name(%u)", s[3].l_name);
	assert(s[3].l_seq == 4, "l_seq(%u)", s[3].l_seq);
	assert(s[3].l_tag == 18, "l_tag(%u)", s[3].l_tag);
	assert(strcmp((const char*)s[3].name, "test3") == 0, "name(%s)", s[3].name);
	assert(strcmp((const char*)s[3].seq, "\x1\x2\x4\x8") == 0, "seq(%s)", s[3].seq);
	assert(strcmp((const char*)s[3].qual, "") == 0, "qual(%s)", s[3].qual);
	assert(strcmp((const char*)s[3].tag, "COZcomment comment") == 0, "tag(%s)", s[3].tag);

	n_seq = bseq_close(b);
	assert(n_seq == 4);
	remove(filename);
}

unittest( .name = "bseq.fastq" ) {
	const char *filename = "./minialign.unittest.bseq.tmp";
	const char *content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\n12+3\n+123\n"
		"@  test2\n\nAAAA\n+  test2\n\n\n12@3\n\n"
		"@test3  comment comment   \nACGT\n\n+ test3\n@123";

	FILE *fp = fopen(filename, "w");
	assert(fp != NULL);
	fwrite(content, 1, strlen(content), fp);
	fclose(fp);

	const uint16_t tags[1] = { ((uint16_t)'C') | (((uint16_t)'O')<<8) };
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
	assert(s[0].l_tag == 0, "l_tag(%u)", s[0].l_tag);
	assert(strcmp((const char*)s[0].name, "test0") == 0, "name(%s)", s[0].name);
	assert(strcmp((const char*)s[0].seq, "\x1\x1\x1\x1") == 0, "seq(%s)", s[0].seq);
	assert(strcmp((const char*)s[0].qual, "NNNN") == 0, "qual(%s)", s[0].qual);
	assert(strcmp((const char*)s[0].tag, "") == 0, "tag(%s)", s[0].tag);

	assert(s[1].l_name == 5, "l_name(%u)", s[1].l_name);
	assert(s[1].l_seq == 8, "l_seq(%u)", s[1].l_seq);
	assert(s[1].l_tag == 0, "l_tag(%u)", s[1].l_tag);
	assert(strcmp((const char*)s[1].name, "test1") == 0, "name(%s)", s[1].name);
	assert(strcmp((const char*)s[1].seq, "\x1\x8\x1\x8\x2\x4\x2\x4") == 0, "seq(%s)", s[1].seq);
	assert(strcmp((const char*)s[1].qual, "12+3+123") == 0, "qual(%s)", s[1].qual);
	assert(strcmp((const char*)s[1].tag, "") == 0, "tag(%s)", s[1].tag);

	assert(s[2].l_name == 5, "l_name(%u)", s[2].l_name);
	assert(s[2].l_seq == 4, "l_seq(%u)", s[2].l_seq);
	assert(s[2].l_tag == 0, "l_tag(%u)", s[2].l_tag);
	assert(strcmp((const char*)s[2].name, "test2") == 0, "name(%s)", s[2].name);
	assert(strcmp((const char*)s[2].seq, "\x1\x1\x1\x1") == 0, "seq(%s)", s[2].seq);
	assert(strcmp((const char*)s[2].qual, "12@3") == 0, "qual(%s)", s[2].qual);
	assert(strcmp((const char*)s[2].tag, "") == 0, "tag(%s)", s[2].tag);

	assert(s[3].l_name == 5, "l_name(%u)", s[3].l_name);
	assert(s[3].l_seq == 4, "l_seq(%u)", s[3].l_seq);
	assert(s[3].l_tag == 18, "l_tag(%u)", s[3].l_tag);
	assert(strcmp((const char*)s[3].name, "test3") == 0, "name(%s)", s[3].name);
	assert(strcmp((const char*)s[3].seq, "\x1\x2\x4\x8") == 0, "seq(%s)", s[3].seq);
	assert(strcmp((const char*)s[3].qual, "@123") == 0, "qual(%s)", s[3].qual);
	assert(strcmp((const char*)s[3].tag, "COZcomment comment") == 0, "tag(%s)", s[3].tag);

	n_seq = bseq_close(b);
	assert(n_seq == 4);
	remove(filename);
}

unittest( .name = "bseq.fastq.skip" ) {
	const char *filename = "./minialign.unittest.bseq.tmp";
	const char *content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\n12+3\n+123\n"
		"@  test2\n\nAAAA\n+  test2\n\n\n12@3\n\n"
		"@test3  comment comment   \nACGT\n\n+ test3\n@123";

	FILE *fp = fopen(filename, "w");
	assert(fp != NULL);
	fwrite(content, 1, strlen(content), fp);
	fclose(fp);

	const uint16_t tags[1] = { ((uint16_t)'C') | (((uint16_t)'O')<<8) };
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
	assert(s[0].l_tag == 0, "l_tag(%u)", s[0].l_tag);
	assert(strcmp((const char*)s[0].name, "test0") == 0, "name(%s)", s[0].name);
	assert(strcmp((const char*)s[0].seq, "\x1\x1\x1\x1") == 0, "seq(%s)", s[0].seq);
	assert(strcmp((const char*)s[0].qual, "") == 0, "qual(%s)", s[0].qual);
	assert(strcmp((const char*)s[0].tag, "") == 0, "tag(%s)", s[0].tag);

	assert(s[1].l_name == 5, "l_name(%u)", s[1].l_name);
	assert(s[1].l_seq == 8, "l_seq(%u)", s[1].l_seq);
	assert(s[1].l_tag == 0, "l_tag(%u)", s[1].l_tag);
	assert(strcmp((const char*)s[1].name, "test1") == 0, "name(%s)", s[1].name);
	assert(strcmp((const char*)s[1].seq, "\x1\x8\x1\x8\x2\x4\x2\x4") == 0, "seq(%s)", s[1].seq);
	assert(strcmp((const char*)s[1].qual, "") == 0, "qual(%s)", s[1].qual);
	assert(strcmp((const char*)s[1].tag, "") == 0, "tag(%s)", s[1].tag);

	assert(s[2].l_name == 5, "l_name(%u)", s[2].l_name);
	assert(s[2].l_seq == 4, "l_seq(%u)", s[2].l_seq);
	assert(s[2].l_tag == 0, "l_tag(%u)", s[2].l_tag);
	assert(strcmp((const char*)s[2].name, "test2") == 0, "name(%s)", s[2].name);
	assert(strcmp((const char*)s[2].seq, "\x1\x1\x1\x1") == 0, "seq(%s)", s[2].seq);
	assert(strcmp((const char*)s[2].qual, "") == 0, "qual(%s)", s[2].qual);
	assert(strcmp((const char*)s[2].tag, "") == 0, "tag(%s)", s[2].tag);

	assert(s[3].l_name == 5, "l_name(%u)", s[3].l_name);
	assert(s[3].l_seq == 4, "l_seq(%u)", s[3].l_seq);
	assert(s[3].l_tag == 18, "l_tag(%u)", s[3].l_tag);
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
#if USE_CRC32_HASH != 0
#define hash64(k0, k1, mask)		( (_mm_crc32_u64((k0), (k0)) ^ (k1)) & (mask) )
#else
static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
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
static void mm_sketch(const uint8_t *seq4, uint32_t len, uint32_t w, uint32_t k, uint32_t rid, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, k0 = 0, k1 = 0;
	uint32_t i, j, l, buf_pos, min_pos;
	mm128_t buf[w], min = { .u64 = { UINT64_MAX, UINT64_MAX } };

	assert(len > 0 && w > 0 && k > 0);
	memset(buf, 0xff, w * sizeof(mm128_t));

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		uint64_t _c = seq4[i];
		mm128_t info = { .u64 = { UINT64_MAX, UINT64_MAX } };
		if (_c != 0) { // not an ambiguous base
			uint64_t c = 0x03 & ((_c>>1) - (_c>>3));
			k0 = (k0 << 2 | c) & mask;           // forward k-mer
			k1 = (k1 >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (k0 == k1) continue; // skip "symmetric k-mers" as we don't know it strand
			if (++l >= k) {
				info.u64[0] = hash64(k0 < k1? k0 : k1, mask);
				info.u32[2] = k0 < k1? i : ~i; info.u32[3] = rid;
			}
		} else l = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.u64[0] == buf[j].u64[0] && buf[j].u64[1] != min.u64[1]) kv_push(mm128_t, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.u64[0] == buf[j].u64[0] && buf[j].u64[1] != min.u64[1]) kv_push(mm128_t, *p, buf[j]);
		}
		if (info.u64[0] <= min.u64[0]) { // a new minimum; then write the old min
			if (l >= w + k) kv_push(mm128_t, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1) kv_push(mm128_t, *p, min);
			for (j = buf_pos + 1, min.u64[0] = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.u64[0] >= buf[j].u64[0]) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.u64[0] >= buf[j].u64[0]) min = buf[j], min_pos = j;
			if (l >= w + k - 1) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.u64[0] == buf[j].u64[0] && min.u64[1] != buf[j].u64[1]) kv_push(mm128_t, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.u64[0] == buf[j].u64[0] && min.u64[1] != buf[j].u64[1]) kv_push(mm128_t, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.u64[0] != UINT64_MAX)
		kv_push(mm128_t, *p, min);
}
#else
static void mm_sketch(const uint8_t *seq4, uint32_t _len, uint32_t w, uint32_t k, uint32_t rid, mm128_v *p)
{
	const uint64_t kk = k - 1, shift1 = 2*kk, mask = (1ULL<<2*k) - 1, len = ((uint64_t)rid<<32) | _len;
	mm128_t buf[16];	// w < 16 must be guaranteed
	#define _push_kmer(_c) { \
		uint64_t _t = 0x03 & (((_c)>>1) - ((_c)>>3)); \
		k0 = (k0 << 2 | _t) & mask; \
		k1 = (k1 >> 2) | ((3ULL^_t) << shift1); \
	}

	uint64_t i = (uint64_t)rid<<32; seq4 -= (ptrdiff_t)i; i--;
	mm128_t *q = p->a+p->n, *t = p->a+p->m;
	do {
		if (q+64 > t) { p->n = q-p->a; p->m = MAX2(256, p->m*2); p->a = realloc(p->a, p->m*sizeof(mm128_t)); q = p->a+p->n; t = p->a+p->m; }
		uint64_t l = 0, k0 = 0, k1 = 0, min = UINT64_MAX, min_pos = 0;
		while (l < kk) {
			uint64_t c;
			if ((c = seq4[++i]) == 0) goto _loop_tail;
			_push_kmer(c); if (k0 == k1) continue;
			buf[l&0x0f] = (mm128_t){ .u64 = { UINT64_MAX, 0 } };
			l++;
		}
		while (l < kk + w) {
			uint64_t c, h;
			if ((c = seq4[++i]) == 0) goto _loop_tail;
			_push_kmer(c); if (k0 == k1) continue;
			buf[l&0x0f] = (mm128_t){ .u64 = { h = hash64(k0 < k1? k0 : k1, k0 < k1 ? k1 : k0, mask), (k0 < k1? 0 : 0xffffffff) ^ i } };
			if (h <= min) min = h, min_pos = l & 0x0f;
			l++;
		}
		for (uint64_t j = kk; j < kk + w; ++j) if (buf[j&0x0f].u64[0] == min) *q++ = buf[j&0x0f];
		q--;
		while (1) {
			uint64_t c, h;
			if ((c = seq4[++i]) == 0) goto _loop_tail;
			_push_kmer(c); if (k0 == k1) continue;
			buf[l&0x0f] = (mm128_t){ .u64 = { h = hash64(k0 < k1? k0 : k1, k0 < k1 ? k1 : k0, mask), (k0 < k1? 0 : 0xffffffff) ^ i } };

			if (h <= min) {
				*q++ = buf[min_pos], min = h, min_pos = l&0x0f;
			} else if (min_pos == ((l-w)&0x0f)) {
				*q++ = buf[min_pos]; min = UINT64_MAX;
				for (uint64_t j = l-w+1; j <= l; ++j) if (buf[j&0x0f].u64[0] <= min) min = buf[j&0x0f].u64[0], min_pos = j&0x0f;
				for (uint64_t j = l-w+1; j <= l; ++j) if (buf[j&0x0f].u64[0] == min) *q++ = buf[j&0x0f];
				q--;
			}
			l++;
			if (q+64 <= t) continue;
			p->n = q-p->a; p->m = MAX2(128, p->m*2); p->a = realloc(p->a, p->m*sizeof(mm128_t)); q = p->a+p->n; t = p->a+p->m;
		}
	_loop_tail:
		if (min != UINT64_MAX) *q++ = buf[min_pos];
	} while (i < len);
	p->n = q-p->a;
	#undef _push_kmer
	return;
}
#endif
/* end of sketch.c */

/* map.c, options */
#define MM_RG			( 0 )		// Z: read group
#define MM_CO			( 1 )		// Z: comment
#define MM_NH			( 2 )		// i: #hits
#define MM_IH 			( 3 )		// i: index of the record in the hits
#define MM_AS			( 4 )		// i: score
#define MM_XS			( 5 )		// i: suboptimal score
#define MM_NM 			( 6 )		// i: editdist to the reference
#define MM_SA			( 7 )		// Z: supplementary records

#define MM_AVA			( 0x01ULL<<48 )
#define MM_KEEP_QUAL	( 0x02ULL<<48 )
#define MM_CIRCULAR		( 0x04ULL<<48 )
#define MM_OMIT_REP		( 0x08ULL<<48 )		// omit secondary records
#define MM_COMP 		( 0x10ULL<<48 )

#define MM_MAF			( 0x01ULL<<56 )
#define MM_BLAST6		( 0x02ULL<<56 )
#define MM_BLASR1		( 0x03ULL<<56 )
#define MM_BLASR4		( 0x04ULL<<56 )
#define MM_PAF			( 0x05ULL<<56 )
#define MM_MHAP 		( 0x06ULL<<56 )
#define MM_FALCON		( 0x07ULL<<56 )

typedef struct {
	uint32_t sidx, eidx, nth, min, k, w, b, verbose;
	uint64_t flag;
	float min_ratio;
	int32_t m, x, gi, ge, xdrop, rmin, qmin, llim, hlim, elim, blim, max_cnt;
	uint32_t n_frq;
	float frq[16];
	uint16_v tags;
	uint64_t batch_size, outbuf_size;
	char *rg_line, *rg_id;
	uint32_t base_rid, base_qid;
} mm_mapopt_t;

static void mm_mapopt_destroy(mm_mapopt_t *opt)
{
	free(opt->rg_line); free(opt->rg_id); free(opt);
	return;
}

static mm_mapopt_t *mm_mapopt_init(void)
{
	set_info(0, "[mm_mapopt_init] initialize mapopt object");
	mm_mapopt_t *opt = calloc(1, sizeof(mm_mapopt_t));
	*opt = (mm_mapopt_t){
		/* -V */ .verbose = 1,
		/* -f, -k, -w, -b, -T */ .k = 15, .w = 16, .b = 14, .flag = 0,
		/* -a, -b, -p, -q, -Y */ .m = 1, .x = 1, .gi = 1, .ge = 1, .xdrop = 50,
		/* -s, -m */ .min = 50, .min_ratio = 0.3,
		/* -M */ .max_cnt = 0,
		/* -f */ .n_frq = 3, .frq[0] = 0.05, .frq[1] = 0.01, .frq[2] = 0.001,
		/* -t */ .nth = 1,
		/* -R */ .rg_line = NULL, .rg_id = NULL,

		/* -S, -E, -L, -H */.sidx = 0, .eidx = 3, .rmin = 0, .qmin = 0,
		.hlim = 7000, .llim = 7000, .blim = 0, .elim = 200,
		.batch_size = 512 * 1024,
		.outbuf_size = 512 * 1024,
		.base_rid = 0, .base_qid = 0
	};
	return opt;
}

static int mm_mapopt_check(mm_mapopt_t *opt, int (*_fprintf)(FILE*,const char*,...), FILE *_fp)
{
	int ret = 0;
	if (opt->k >= 32) _fprintf(_fp, "[E::%s] k must be inside [1,32).\n", __func__), ret = 1;
	if (ret) return(ret);
	if (opt->w >= 16) _fprintf(_fp, "[E::%s] w must be inside [1,16).\n", __func__), ret = 1;
	if (ret) return(ret);

	if (opt->m < 1 || opt->m > 5) _fprintf(_fp, "[E::%s] Match award must be inside [1,5].\n", __func__), ret = 1;
	if (opt->x < 1 || opt->x > 5) _fprintf(_fp, "[E::%s] Mismatch penalty must be inside [1,5].\n", __func__), ret = 1;
	if (opt->gi > 5) _fprintf(_fp, "[E::%s] Gap open penalty must be inside [0,5].\n", __func__), ret = 1;
	if (opt->ge < 1 || opt->ge > 5) _fprintf(_fp, "[E::%s] Gap extension penalty must be inside [1,5].\n", __func__), ret = 1;
	if (ret) return(ret);

	if (opt->gi == 0 && opt->x == 1 && opt->ge == 1)
		_fprintf(_fp, "[I::%s] (M,X,Gi,Ge) = (1,1,0,1) has positive expected score for two independent random sequences (thus result in false positives). Please consider using a more stringent score.\n", __func__);
	if (opt->gi != 0 && opt->x >= (opt->gi + opt->ge))
		_fprintf(_fp, "[I::%s] Large mismatch penalty with respect to the gap open/extend penalty may cause SEGV or broken CIGAR. [issue #2]\n", __func__);
	if (opt->gi != 0 && opt->m + 2*(opt->gi + opt->ge) > 10)
		_fprintf(_fp, "[I::%s] Large match award or large gap open/extend penalty may cause SEGV or broken CIGAR. [issue #7]\n", __func__);
	if (opt->xdrop < 10 || opt->xdrop > 100) _fprintf(_fp, "[E::%s] Xdrop cutoff must be inside [10,100].\n", __func__), ret = 1;
	if (opt->min > INT32_MAX) _fprintf(_fp, "[E::%s] Minimum alignment score must be > 0.\n", __func__), ret = 1;
	if (opt->min_ratio < 0.0 || opt->min_ratio > 1.0) _fprintf(_fp, "[E::%s] Minimum alignment score ratio must be inside [0.0,1.0].\n", __func__), ret = 1;
	if (opt->n_frq >= 16) _fprintf(_fp, "[E::%s] Frequency thresholds must be fewer than 16.\n", __func__), ret = 1;
	for (uint64_t i = 0; i < opt->n_frq; ++i)
		if (opt->frq[i] < 0.0 || opt->frq[i] > 1.0 || (i != 0 && opt->frq[i-1] < opt->frq[i]))
			_fprintf(_fp, "[E::%s] Frequency thresholds must be inside [0.0,1.0] and descending.\n", __func__), ret = 1;
	if (ret) return(ret);

	if (opt->nth < 1 || opt->nth >= MAX_THREADS) _fprintf(_fp, "[E::%s] Thread counts must be inside [1,%u]. For larger values, recompile is needed.\n", __func__, MAX_THREADS), ret = 1;
	if (opt->batch_size < 64 * 1024) _fprintf(_fp, "[E::%s] Batch size must be > 64k.\n", __func__), ret = 1;
	if (opt->outbuf_size < 64 * 1024) _fprintf(_fp, "[E::%s] Output buffer size must be > 64k.\n", __func__), ret = 1;
	if (ret) return(ret);

	if (opt->sidx >= 16) _fprintf(_fp, "[E::%s] sidx must be inside [0,16).\n", __func__), ret = 1;
	if (opt->eidx >= 16) _fprintf(_fp, "[E::%s] eidx must be inside [0,16).\n", __func__), ret = 1;
	if (opt->hlim < 100 || opt->hlim >= 100000) _fprintf(_fp, "[E::%s] hlim must be inside [100,100000).\n", __func__), ret = 1;
	if (opt->llim < 100 || opt->llim >= 100000) _fprintf(_fp, "[E::%s] llim must be inside [100,100000).\n", __func__), ret = 1;
	if (opt->blim >= 10000) _fprintf(_fp, "[E::%s] blim must be inside [0,10000).\n", __func__), ret = 1;
	if (opt->elim >= 10000) _fprintf(_fp, "[E::%s] elim must be inside [0,10000).\n", __func__), ret = 1;
	if (ret) return(ret);

	#define _dup(x)	({ char *_p = malloc(strlen(x)+1); memcpy(_p, (x), strlen(x)); _p[strlen(x)] = '\0'; _p; })
	if (opt->flag&(0x01ULL<<MM_RG)) {
		if (opt->rg_line == NULL && opt->rg_id == NULL) {
			opt->rg_line = _dup("@RG\tID:1\tSM:default"); opt->rg_id = _dup("1");
		}
	}
	#undef _dup

	if (!(opt->flag&MM_AVA) && (opt->flag&MM_COMP)) _fprintf(_fp, "[W::%s] `-A' flag is only effective in all-versus-all mode. ignored.\n", __func__), opt->flag &= ~MM_COMP;
	return ret;
}

/* end of map.c, options */

/* index.c */
typedef struct {
	// lengths and flag
	uint32_t l_seq, rid, l_name, circular;
	// pointers
	char *name;
	uint8_t *seq;
} mm_idx_seq_t;
_static_assert(sizeof(mm_idx_seq_t) == 32);
typedef struct { size_t n, m; mm_idx_seq_t *a; } mm_idx_seq_v;

typedef struct {
	mm128_v a;   // (minimizer, position) array
	uint64_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	kh_t *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

typedef struct {
	uint64_t mask;
	mm_idx_bucket_t *bkt;
	mm_idx_seq_v s;
	uint8_t b, w, k, circular;
	uint32_t base_rid;

	// work
	mm128_v a;
	uint64_v size;
	ptr_v base;
} mm_idx_t;

static mm_idx_t *mm_idx_init(uint32_t w, uint32_t k, uint32_t b, uint32_t base_rid, uint32_t circular)
{
	mm_idx_t *mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w<1? 1 : w; mi->k = k; mi->b = MIN2(k*2, b); mi->mask = (1<<mi->b) - 1;
	mi->circular = circular; mi->base_rid = base_rid;
	mi->bkt = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	return mi;
}

static void mm_idx_destroy(mm_idx_t *mi)
{
	if (mi == 0) return;
	for (uint64_t i = 0; i < 1ULL<<mi->b; ++i) {
		free(mi->bkt[i].p);
		free(mi->bkt[i].a.a);
		kh_destroy((kh_t*)mi->bkt[i].h);
	}
	free(mi->bkt);
	for (uint64_t i = 0; i < mi->base.n; ++i) free(mi->base.a[i]);
	free(mi->base.a); free(mi->size.a);
	free(mi->s.a); free(mi->a.a); free(mi);
}

static const v2u32_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, uint64_t *n)
{
	mm_idx_bucket_t *b = &mi->bkt[minier & mi->mask];
	kh_t *h = (kh_t*)b->h;
	const uint64_t *p;
	if (h == NULL || (p = kh_get_ptr(h, minier>>mi->b)) == NULL) {
		*n = 0; return NULL;
	}
	if ((int64_t)*p >= 0) {
		*n = 1;
		return (const v2u32_t*)p;
	} else {
		*n = (uint32_t)*p;
		return (const v2u32_t*)&b->p[(*p>>32) & 0x7fffffff];
	}
}

static uint32_t mm_idx_cal_max_occ(const mm_idx_t *mi, double f)
{
	uint64_t n = 0;
	uint32_t thres;

	set_info(0, "[mm_idx_cal_max_occ] calculate occurrence thresholds");
	if (f <= 0.) return UINT32_MAX;
	for (uint64_t i = (n = 0); i < 1ULL<<mi->b; ++i)
		if (mi->bkt[i].h) n += kh_cnt((kh_t*)mi->bkt[i].h);
	
	uint32_t *a = (uint32_t*)malloc(n * sizeof(uint32_t));
	for (uint64_t i = (n = 0); i < 1ULL<<mi->b; ++i) {
		kh_t *h = (kh_t*)mi->bkt[i].h;
		if (h == 0) continue;
		for (uint64_t k = 0; k < kh_size(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = ((int64_t)kh_val(h, k))>=0? 1 : (uint32_t)kh_val(h, k);
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

/******************
 * Generate index *
 ******************/

typedef struct {
	uint32_t icnt, ocnt;
	bseq_file_t *fp;
	mm_idx_t *mi;
	kvec_t(mm128_t) hq;
} mm_idx_pipeline_t;

typedef struct {
    uint32_t id, base_rid, n_seq;
	const bseq_t *seq;	// const!
	void *base;
	uint64_t size;
	mm128_v a;
} mm_idx_step_t;

static void *mm_idx_source(uint32_t tid, void *arg)
{
	set_info(tid, "[mm_idx_source] fetch sequence block");
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t*)arg;
	mm_idx_step_t *s = (mm_idx_step_t*)calloc(1, sizeof(mm_idx_step_t));
	uint64_t size;
	void *base;

	s->seq = bseq_read(q->fp, &s->n_seq, &base, &size);
	if (s->seq == 0) { free(s), s = NULL; return NULL; }
	s->id = q->icnt++;

	// register fetched block
	kv_push(void*, q->mi->base, base);
	kv_push(uint64_t, q->mi->size, size);

	// update base_rid and seq array
	if (q->mi->base_rid == UINT32_MAX) q->mi->base_rid = atoi(s->seq[0].name);	// assume first seq name is base_rid
	s->base_rid = q->mi->base_rid + q->mi->s.n;
	if (q->mi->s.n + s->n_seq > q->mi->s.m) {
		q->mi->s.m = MAX2(256, 2*q->mi->s.m);
		q->mi->s.a = realloc(q->mi->s.a, sizeof(mm_idx_seq_t) * q->mi->s.m);
	}

	const bseq_t *src = s->seq;
	mm_idx_seq_t *dst = &q->mi->s.a[q->mi->s.n];
	for (uint64_t i = 0; i < s->n_seq; ++i) {
		dst[i] = (mm_idx_seq_t){
			.l_seq = src[i].l_seq, .l_name = src[i].l_name,
			.rid = s->base_rid + i, .circular = q->mi->circular,
			.name = src[i].name, .seq = src[i].seq
		};
	}
	q->mi->s.n += s->n_seq;
	return s;
}

static void *mm_idx_worker(uint32_t tid, void *arg, void *item)
{
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t*)arg;
	mm_idx_step_t *s = (mm_idx_step_t*)item;

	char buf[128], *p = buf;
	p += _pstr(p, "[mm_idx_worker] bin id "); p += _pnum(uint32_t, p, s->base_rid); p += _pstr(p, ":"); p += _pnum(uint32_t, p, s->base_rid + s->n_seq - 1); *p = '\0';
	set_info(tid, buf);
	for (uint64_t i = 0; i < s->n_seq; ++i)
		mm_sketch(s->seq[i].seq, s->seq[i].l_seq, q->mi->w, q->mi->k, s->base_rid + i, &s->a);
	return s;
}

static void mm_idx_drain_intl(mm_idx_pipeline_t *q, mm_idx_step_t *s)
{
	uint64_t mask = q->mi->mask;
	for (uint64_t i = 0; i < s->a.n; ++i) {
		mm128_v *p = &q->mi->bkt[s->a.a[i].u64[0]&mask].a;
		kv_push(mm128_t, *p, s->a.a[i]);
	}
	free((void*)s->seq); free(s->a.a); free(s);
	return;
}

static void mm_idx_drain(uint32_t tid, void *arg, void *item)
{
	set_info(tid, "[mm_idx_drain] dump minimizers to pool");
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t*)arg;
	mm_idx_step_t *s = (mm_idx_step_t*)item;

	#if STRICT_STREAM_ORDERING != 0
		kv_hq_push(mm128_t, incq_comp, q->hq, ((mm128_t){.u64 = {s->id, (uintptr_t)s}}));
		while (q->hq.n > 1 && q->hq.a[1].u64[0] == q->ocnt) {
			s = (mm_idx_step_t*)kv_hq_pop(mm128_t, incq_comp, q->hq).u64[1]; q->ocnt++;
			mm_idx_drain_intl(q, s);
		}
	#else
		mm_idx_drain_intl(q, s);
	#endif
	return;
}

typedef struct {
	mm_idx_t *mi;
	uint32_t from, to;
} mm_idx_post_t;

static void *mm_idx_post(uint32_t tid, void *arg, void *item)
{
	set_info(tid, "[mm_idx_post] indexing postprocess");
	mm_idx_post_t *q = (mm_idx_post_t*)item;
	mm_idx_t *mi = q->mi;

	for (uint64_t i = q->from; i < q->to; ++i) {
		mm_idx_bucket_t *b = &mi->bkt[i];
		if (b->a.n == 0) continue;

		// sort by minimizer
		radix_sort_128x(b->a.a, b->a.a + b->a.n);

		// count and preallocate
		uint64_t n_keys = 0; b->n = 0;
		for (uint64_t j = 1, n = 1; j <= b->a.n; ++j, ++n) {
			if (j != b->a.n && b->a.a[j].u64[0] == b->a.a[j-1].u64[0]) continue;
			b->n += (n > 1)? n : 0; ++n_keys; n = 0;
		}
		kh_t *h = kh_init(n_keys / KH_THRESH);
		b->p = (uint64_t*)malloc(sizeof(uint64_t) * b->n);

		// create the hash table
		for (uint64_t j = 1, n = 1, sp = 0; j <= b->a.n; ++j, ++n) {
			if (j != b->a.n && b->a.a[j].u64[0] == b->a.a[j-1].u64[0]) continue;

			mm128_t *p = &b->a.a[j-n];
			uint64_t key = p->u64[0]>>mi->b, val = p->u64[1];
			if (n != 1) {
				b->p[sp++] = val; val = (sp-1)<<32 | n | 0x01ULL<<63;	// k = 0
				while (n > 1) b->p[sp++] = b->a.a[j - --n].u64[1];
			}
			kh_put(h, key, val); n = 0;
		}
		b->h = h;

		// deallocate and clear b->a
		free(b->a.a);
		b->a.n = b->a.m = 0, b->a.a = 0;
	}
	return NULL;
}

static mm_idx_t *mm_idx_gen(const mm_mapopt_t *opt, bseq_file_t *fp)
{
	mm_idx_pipeline_t pl = {0}, **p;

	set_info(0, "[mm_idx_gen] initialize index object");
	pl.icnt = pl.ocnt = 0;
	pl.fp = fp;
	if (pl.fp == NULL) return NULL;
	pl.mi = mm_idx_init(opt->w, opt->k, opt->b, opt->base_rid, (opt->flag&MM_CIRCULAR) != 0);
	kv_hq_init(pl.hq);

	p = (mm_idx_pipeline_t**)calloc(opt->nth, sizeof(mm_idx_pipeline_t*));
	for (uint64_t i = 0; i < opt->nth; ++i) p[i] = &pl;

	pt_t *pt = pt_init(opt->nth);
	pt_stream(pt, mm_idx_source, &pl, mm_idx_worker, (void**)p, mm_idx_drain, &pl);
	free(p); kv_hq_destroy(pl.hq);
	if (opt->verbose >= 1)
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	mm_idx_post_t *q = (mm_idx_post_t*)calloc(opt->nth, sizeof(mm_idx_post_t));
	mm_idx_post_t **qq = (mm_idx_post_t**)calloc(opt->nth, sizeof(mm_idx_post_t*));
	for (uint64_t i = 0; i < opt->nth; ++i) {
		q[i].mi = pl.mi;
		q[i].from = (1ULL<<pl.mi->b)*i/opt->nth;
		q[i].to = (1ULL<<pl.mi->b)*(i+1)/opt->nth;
		qq[i] = &q[i];
	}
	pt_parallel(pt, mm_idx_post, NULL, (void**)qq);
	pt_destroy(pt);

	free(q); free(qq);
	if (opt->verbose >= 1)
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
	return pl.mi;
}

#if 0
static void mm_idx_cmp(const mm_mapopt_t *opt, const mm_idx_t *m1, const mm_idx_t *m2)
{
	for (uint64_t i = 0; i < 1ULL<<opt->b; ++i) {
		mm_idx_bucket_t *bkt1 = &m1->bkt[i], *bkt2 = &m2->bkt[i];
		if (bkt1 == NULL || bkt2 == NULL) {
			if (bkt1 == NULL && bkt2 == NULL) continue;
			if (bkt1) fprintf(stderr, "i(%lu), bkt1 is instanciated but bkt2 is not\n", i);
			if (bkt2) fprintf(stderr, "i(%lu), bkt2 is instanciated but bkt1 is not\n", i);
			continue;
		}

		if (bkt1->n != bkt2->n) fprintf(stderr, "i(%lu), array size differs(%lu, %lu)\n", i, bkt1->n, bkt2->n);
		for (uint64_t j = 0; j < MIN2(bkt1->n, bkt2->n); ++j) {
			if (bkt1->p[j] != bkt2->p[j]) fprintf(stderr, "i(%lu), array differs at j(%lu), (%lu, %lu)\n", i, j, bkt1->p[j], bkt2->p[j]);
		}

		kh_t *h1 = bkt1->h, *h2 = bkt2->h;
		if (h1 == NULL || h2 == NULL) {
			if (h1 == NULL && h2 == NULL) continue;
			if (h1) fprintf(stderr, "i(%lu), h1 is instanciated but h2 is not\n", i);
			if (h2) fprintf(stderr, "i(%lu), h2 is instanciated but h1 is not\n", i);
			continue;
		}
		if (h1->mask != h2->mask) fprintf(stderr, "i(%lu), hash size differs(%u, %u)\n", i, h1->mask, h2->mask);
		if (h1->max != h2->max) fprintf(stderr, "i(%lu), hash max differs(%u, %u)\n", i, h1->max, h2->max);
		if (h1->cnt != h2->cnt) fprintf(stderr, "i(%lu), hash cnt differs(%u, %u)\n", i, h1->cnt, h2->cnt);
		if (h1->ub != h2->ub) fprintf(stderr, "i(%lu), hash ub differs(%u, %u)\n", i, h1->ub, h2->ub);
		for (uint64_t j = 0; j < MIN2(h1->mask, h2->mask)+1; ++j) {
			if (h1->a[j].u64[0] != h2->a[j].u64[0]) fprintf(stderr, "i(%lu), hash key differs at j(%lu), (%lu, %lu)\n", i, j, h1->a[j].u64[0], h2->a[j].u64[0]);
			if (h1->a[j].u64[1] != h2->a[j].u64[1]) fprintf(stderr, "i(%lu), hash val differs at j(%lu), (%lu, %lu)\n", i, j, h1->a[j].u64[1], h2->a[j].u64[1]);
		}
	}
}
#endif

/*************
 * index I/O *
 *************/

#define MM_IDX_MAGIC "MAI\7"		/* minialign index version 7 */

static void mm_idx_dump(FILE *fp, const mm_idx_t *mi, uint32_t nth)
{
	uint64_t y[2] = { mi->s.n, 0 };
	set_info(0, "[mm_idx_dump] dump index to file (main thread)");

	for (uint64_t i = 0; i < mi->size.n; ++i) y[1] += mi->size.a[i];

	pg_t *pg = pg_init(fp, nth);
	pgwrite(pg, MM_IDX_MAGIC, strlen(MM_IDX_MAGIC));
	pgwrite(pg, &mi->b, sizeof(uint8_t) * 4);	// b, w, k, circular
	pgwrite(pg, &mi->base_rid, sizeof(uint32_t));
	pgwrite(pg, y, sizeof(uint64_t) * 2);
	for (uint64_t i = 0; i < 1ULL<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->bkt[i];
		pgwrite(pg, &b->n, sizeof(uint32_t) * 1);
		pgwrite(pg, b->p, sizeof(uint64_t) * b->n);
		kh_dump((kh_t*)b->h, pg, (khwrite_t)pgwrite);
	}
	for (uint64_t i = 0; i < mi->base.n; ++i)
		pgwrite(pg, mi->base.a[i], sizeof(char) * mi->size.a[i]);
	for (uint64_t i = 0, j = 0, s = 0; i < mi->s.n; ++i) {
		if ((uintptr_t)mi->s.a[i].seq < (uintptr_t)mi->base.a[j]
		|| (uintptr_t)mi->s.a[i].seq >= (uintptr_t)mi->base.a[j] + (ptrdiff_t)mi->size.a[j]) {
			s += mi->size.a[j++];
		}
		mi->s.a[i].name -= (ptrdiff_t)mi->base.a[j], mi->s.a[i].seq -= (ptrdiff_t)mi->base.a[j];
		mi->s.a[i].name += (ptrdiff_t)s, mi->s.a[i].seq += (ptrdiff_t)s;
	}
	pgwrite(pg, mi->s.a, sizeof(mm_idx_seq_t) * mi->s.n);
	// restore pointers
	for (uint64_t i = 0, j = 0, s = 0; i < mi->s.n; ++i) {
		mi->s.a[i].name -= (ptrdiff_t)s, mi->s.a[i].seq -= (ptrdiff_t)s;
		mi->s.a[i].name += (ptrdiff_t)mi->base.a[j], mi->s.a[i].seq += (ptrdiff_t)mi->base.a[j];
		if ((uintptr_t)mi->s.a[i].name > s + mi->size.a[j]) s += mi->size.a[j++];
	}
	pg_destroy(pg);
	return;
}

static mm_idx_t *mm_idx_load(FILE *fp, uint32_t nth)
{
	char magic[4];
	uint64_t bsize, y[2];
	mm_idx_t *mi, b;
	set_info(0, "[mm_idx_load] load index from file (main thread)");

	pg_t *pg = pg_init(fp, nth);
	if (pgread(pg, magic, strlen(MM_IDX_MAGIC)) != strlen(MM_IDX_MAGIC)) return NULL;
	if (strncmp(magic, MM_IDX_MAGIC, strlen(MM_IDX_MAGIC)) != 0) return NULL;
	if (pgread(pg, &b.b, sizeof(uint8_t) * 4) != sizeof(uint8_t) * 4) return NULL;	// b, w, k, circular
	if (pgread(pg, &b.base_rid, sizeof(uint32_t)) != sizeof(uint32_t)) return NULL;
	if (pgread(pg, y, sizeof(uint64_t) * 2) != sizeof(uint64_t) * 2) return NULL;
	mi = mm_idx_init(b.w, b.k, b.b, b.base_rid, b.circular); mi->s.n = y[0]; bsize = y[1];
	for (uint64_t i = 0; i < 1ULL<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->bkt[i];
		if (pgread(pg, &b->n, sizeof(uint32_t) * 1) != sizeof(uint32_t) * 1) goto _mm_idx_load_fail;
		b->p = (uint64_t*)malloc(b->n * sizeof(uint64_t));
		if (pgread(pg, b->p, sizeof(uint64_t) * b->n) != sizeof(uint64_t) * (size_t)b->n) goto _mm_idx_load_fail;
		b->h = kh_load(pg, (khread_t)pgread);
	}
	mi->base.n = mi->size.n = 1;
	mi->base.a = malloc(sizeof(void*) * mi->base.n);
	mi->size.a = malloc(sizeof(uint64_t) * mi->size.n);
	mi->base.a[0] = malloc(sizeof(char) * bsize);
	mi->size.a[0] = bsize;

	const uint64_t chunk_size = 1024 * 1024 * 1024;	// 1G to fit in signed int
	for (uint64_t b = 0; b < mi->size.a[0]; b += chunk_size) {
		uint64_t size = MIN2(chunk_size, mi->size.a[0] - b);
		if (pgread(pg, mi->base.a[0] + b, sizeof(char) * size) != sizeof(char) * size) goto _mm_idx_load_fail;
	}
	mi->s.a = malloc(sizeof(mm_idx_seq_t) * mi->s.n);
	if (pgread(pg, mi->s.a, sizeof(mm_idx_seq_t) * mi->s.n) != sizeof(mm_idx_seq_t) * mi->s.n) goto _mm_idx_load_fail;
	for (uint64_t i = 0; i < mi->s.n; ++i) mi->s.a[i].name += (ptrdiff_t)mi->base.a[0], mi->s.a[i].seq += (ptrdiff_t)mi->base.a[0];
	pg_destroy(pg);
	return mi;
_mm_idx_load_fail:
	pg_destroy(pg);
	mm_idx_destroy(mi);
	return NULL;
}

/* end of index.c */

/* map.c */
typedef struct mm_align_s mm_align_t;
typedef struct {
	void (*header)(mm_align_t*);
	void (*unmapped)(mm_align_t*, const bseq_t*);
	void (*mapped)(mm_align_t*, const bseq_t*, uint64_t, const mm128_t*);	// nreg contains #regs in lower 32bit, #uniq in higher 32bit
	void (*footer)(mm_align_t*);
} mm_printer_t;

struct mm_align_s {
	uint8_t *base, *tail, *p;
	uint8_t conv[40];	// binary -> string conv table
	const mm_idx_t *mi;
	const mm_mapopt_t *opt;
	uint64_t *occ;
	uint32_t n_occ;
	uint32_t icnt, ocnt, base_qid;
	bseq_file_t *fp;
	kvec_t(mm128_t) hq;
	mm_printer_t printer;
	gaba_t *gaba;
	void **t;	// mm_tbuf_t **
	pt_t *pt;
};

typedef struct {
	// query sequences
	uint32_t id, base_qid, n_seq;
	const bseq_t *seq;	// const!
	void *base;
	uint64_t size;
	// results
	lmm_t *lmm;		// alignment result container
	mm128_v reg;
} mm_align_step_t;

typedef struct mm_tbuf_s { // per-thread buffer
	const mm_align_t *b;
	gaba_section_t qf, qr, t;
	mm128_v resc; // also used as query minimizer vector
	mm128_v coef; // Hough transform coefficient
	v2u32_v intv; // intervals on sorted coef
	kh_t *pos;
	gaba_dp_t *dp;	// alignment work
} mm_tbuf_t;

#define _s(x)		( (x)<0?-1:1)
#define _m(x)		( (((int32_t)(x))>>31)^(x) )
static void mm_expand(uint32_t n, const v2u32_t *r, uint32_t qid, int32_t qs, uint32_t org, uint32_t thresh, mm128_v *coef)
{
	const int32_t ofs = 0x40000000;
	if (n == 0) return;
	kv_reserve(mm128_t, *coef, coef->n + n);
	for (uint64_t i = 0; i < n; ++i) {	// iterate over all the collected minimizers
		if (org + r[i].x[1] - qid < thresh) continue;
		int32_t rs = (int32_t)r[i].x[0], _qs = (rs>>31) ^ qs, _rs = (rs>>31) ^ rs;
		mm128_t *p = &coef->a[coef->n++];
		p->u32[0] = ofs + _rs - (_qs>>1); p->u32[1] = r[i].x[1]; p->u32[2] = _rs; p->u32[3] = _qs;
	}
	return;
}

static void mm_collect(const mm_idx_t *mi, uint32_t l_seq, const uint8_t *seq, uint32_t qid, uint32_t max_occ, uint32_t resc_occ, uint32_t org, uint32_t thresh, mm128_v *mini, mm128_v *coef)
{
	uint64_t n_resc = 0;
	mini->n = coef->n = 0;
	mm_sketch(seq, l_seq, mi->w, mi->k, 0, mini);

	for (uint64_t i = 0; i < mini->n; ++i) {
		uint64_t n;
		const v2u32_t *r;
		int32_t qs = (int32_t)mini->a[i].u32[2];
		r = mm_idx_get(mi, mini->a[i].u64[0], &n);	// get minimizer at current pos
		if (n > max_occ) continue;	// skip if exceeds repetitive threshold
		if (n > resc_occ) {
			mm128_t *q = &mini->a[n_resc++];
			q->u32[0] = qs; q->u32[1] = n; q->u64[1] = (uintptr_t)r;
			continue;
		}
		mm_expand(n, r, qid, qs, org, thresh, coef);
	}
	mini->n = n_resc;	// write back rescued array
	return;
}

static void mm_chain(uint64_t l_coef, mm128_t *coef, uint32_t llim, uint32_t hlim, uint32_t min, v2u32_v *intv)
{
	const int32_t ofs = 0x40000000;
	const uint32_t chained = 0x80000000, mask = 0x7fffffff;
	for (uint64_t i = intv->n = 0, j, k, n; i < l_coef; i = MIN2(j, k)) {
		uint32_t rid = coef[i].u32[1];
		int32_t rs = coef[i].u32[2], qs = coef[i].u32[3], re = rs, qe = qs;
		int32_t l, h, lub = coef[i].u32[0] + llim, hub = rs - (qs<<1) + ofs, hlb = hub - hlim;
		uint32_t len = 0;
		for (j = i+1, k = UINT64_MAX, n = i; j < l_coef && (l = (int32_t)coef[j].u32[0]) < lub; j++) {
			if ((int32_t)coef[j].u32[2] < 0) continue;
			re = coef[j].u32[2], qe = coef[j].u32[3];	/* tagged seeds are skipped, no need to mask */
			if (rid != coef[j].u32[1] || (qs^qe)&chained || (h = ofs + re - (qe<<1)) < hlb || h > hub) { k = MIN2(j, k); continue; }	// out of range, skip
			lub = l + llim; hub = h; hlb = h - hlim;
			coef[j].u32[2] |= chained; n = j;
		}
		re = coef[n].u32[2] & mask; qe = coef[n].u32[3]; qs = _m(qs); qe = _m(qe);
		len = re-rs+_s(qe-qs)*(qe-qs);
		if (len < min) continue;
		v2u32_t *p;
		kv_pushp(v2u32_t, *intv, &p);
		p->x[0] = (uint32_t)ofs - len; p->x[1] = i;
	}
	return;
}

static uint64_t mm_short_chain(uint64_t l_coef, mm128_t *coef, uint32_t llim, uint32_t hlim, uint32_t rid, uint32_t eidx)
{
	const int32_t ofs = 0x40000000;
	const uint32_t chained = 0x80000000, mask = 0x7fffffff;
	int32_t rs = coef->u32[2] & mask, qs = coef->u32[3], re = rs, qe = qs;
	int32_t l, h, lub = coef->u32[0] + llim, hub = rs - (qs<<1) + ofs, hlb = hub - hlim;
	uint32_t len = 0;
	if (rid != coef->u32[1]) return 1;
	rid = coef->u32[1];
	uint64_t j, k;
	for (j = 1, k = UINT64_MAX; j < l_coef && (l = (int32_t)coef[j].u32[0]) < lub; j++) {
		re = coef[j].u32[2] & mask, qe = coef[j].u32[3];
		if (rid != coef[j].u32[1] || (qs^qe)&chained || (h = ofs + re - (qe<<1)) < hlb || h > hub) { k = MIN2(j, k); continue; }
		lub = l + llim; hub = h; hlb = h - hlim;	// n = j;
		if (++len > eidx) break;
	}
	return MIN2(j, k);
}

static uint64_t mm_rescue(const mm_idx_seq_t *ref, uint64_t l_coef, mm128_t *coef, mm128_t *r, uint32_t llim, uint32_t hlim, uint32_t elim, uint32_t blim, uint32_t eidx)
{
	const int32_t ofs = 0x40000000, mask = 0x7fffffff;
	int32_t re = r->u32[2], qe = r->u32[3], l, h;
	int32_t lt = ofs + re - (qe>>1), ht = ofs + re - (qe<<1), lb = lt - blim, le = lt + elim, hb = ht + blim, he = ht - elim, lub = lt + llim;
	for (uint64_t i = 0; i < l_coef && (int32_t)coef[i].u32[0] < lb; ++i) {}
	for (uint64_t i = 0; i < l_coef && (l = (int32_t)coef[i].u32[0]) < lub; ++i) {
		uint32_t rs = coef[i].u32[2] & mask, qs = coef[i].u32[3];
		if ((h = ofs + rs - (qs<<1)) > hb || (l < le && h > he && ((l > lt) ^ (h < ht)))) continue;
		return mm_short_chain(l_coef - i, coef + i, llim, hlim, ref->rid, eidx);
	}
	return 1;
}

static const gaba_alignment_t *mm_extend(
	const mm_idx_seq_t *ref, uint32_t l_coef, mm128_t *coef, uint32_t k, uint32_t min, uint32_t sidx, uint32_t eidx,
	gaba_dp_t *dp, gaba_section_t *qf, gaba_section_t *qr, gaba_section_t *t, kh_t *pos, lmm_t *lmm)
{
	const uint32_t mask = 0x7fffffff;
	const uint8_t *lim = (const uint8_t*)0x800000000000;
	gaba_section_t *qu, *qd, *r, *q;
	gaba_section_t rf = { .id = ref->rid<<1, .len = ref->l_seq, .base = (const uint8_t*)ref->seq };
	gaba_section_t rr = { .id = (ref->rid<<1)+1, .len = ref->l_seq, .base = gaba_rev((const uint8_t*)ref->seq+ref->l_seq-1, lim) };
	gaba_pos_pair_t p = {0};
	gaba_dp_flush(dp, lim, lim);
	uint64_t key = (uint64_t)ref->rid<<32;
	gaba_alignment_t *a = NULL;
	for (uint64_t i = sidx; i < eidx && i < l_coef; ++i) {
		if (i != 0 && (int32_t)coef[i].u32[2] >= 0) continue;	// skip head
		const gaba_stack_t *stack = gaba_dp_save_stack(dp);
		int32_t rs = coef[i].u32[2] & mask, qs = coef[i].u32[3];
		uint64_t rev = qs<0; qu = rev? qf : qr; qd = rev? qr : qf;
		qs = _m(qs); qs = rev? (int32_t)(qf->len-qs+k-1) : qs;
		// upward extension
		gaba_fill_t *f = gaba_dp_fill_root(dp, r = &rr, ref->l_seq-rs-1, q = qu, qu->len-qs-1), *m = f;
		if (f == NULL) goto _abort;
		uint32_t flag = GABA_STATUS_TERM;
		do {
			if (f->status & GABA_STATUS_UPDATE_A) r = t;
			if (f->status & GABA_STATUS_UPDATE_B) q = t;
			flag |= f->status & (GABA_STATUS_UPDATE_A | GABA_STATUS_UPDATE_B);
			if ((f = gaba_dp_fill(dp, f, r, q)) == NULL) goto _abort;
			m = (f->max > m->max)? f : m;
		} while (!(flag & f->status));
		// find max
		p = gaba_dp_search_max(dp, m);
		// check duplicate
		key |= (uint32_t)(p.apos - (p.bpos>>1));
		if (kh_get_ptr(pos, key) != NULL) return 0;	// already evaluated
		kh_put(pos, key, (uintptr_t)NULL);	// mark evaluated
		// downward extension from max
		gaba_dp_flush_stack(dp, stack);
		if ((m = f = gaba_dp_fill_root(dp, r = &rf, ref->l_seq-p.apos-1, q = qd, qd->len-p.bpos-1)) == NULL) goto _abort;
		flag = GABA_STATUS_TERM;
		do {
			if (f->status & GABA_STATUS_UPDATE_A) r = t;
			if (f->status & GABA_STATUS_UPDATE_B) q = t;
			flag |= f->status & (GABA_STATUS_UPDATE_A | GABA_STATUS_UPDATE_B);
			if ((f = gaba_dp_fill(dp, f, r, q)) == NULL) goto _abort;
			m = (f->max > m->max)? f : m;
		} while (!(flag & f->status));
		if (m->max < min) { key &= 0xffffffff00000000; continue; }
		// convert alignment to cigar
		a = gaba_dp_trace(dp, NULL, m, GABA_TRACE_PARAMS( .lmm = lmm ));	// might be NULL
		break;
	}
	// record head
	if (a) kh_put(pos, key, (uintptr_t)a);
	return a;
_abort:;
	oom_abort(__func__, 0);
	return 0;
}

#define MAPQ_DEC	( 4 )
#define MAPQ_COEF	( 1<<MAPQ_DEC )
#define _clip(x)	MAX2(0, MIN2(((uint32_t)(x)), 60 * MAPQ_COEF))
#define _aln(x)		( (const gaba_alignment_t*)(x).u64[1] )
static uint64_t mm_prune_regs(const mm_mapopt_t *opt, void *lmm, uint32_t n_reg, mm128_t *reg)
{
	uint64_t q = opt->max_cnt? MIN2(opt->max_cnt, n_reg) : n_reg;
	uint32_t min = (uint32_t)(_aln(reg[0])->score * opt->min_ratio);
	while (_aln(reg[q-1])->score < min) q--;
	for (uint64_t i = q; i < n_reg; ++i) { lmm_free(lmm, (void*)reg[i].u64[1]); }		// destroy gaba_alignment_t objects
	return q;
}

#if COLLECT_SUPPLEMENTARY != 0
static uint64_t mm_post_map(const mm_mapopt_t *opt, uint32_t n_reg, mm128_t *reg)
{
	// collect supplementaries
	#define _swap_128(x, y)	{ mm128_t _tmp = reg[x]; reg[x] = reg[y]; reg[y] = _tmp; }
	uint64_t p, q;
	for (p = 1, q = n_reg; p < q; ++p) {
		uint64_t max = 0;
		for (uint64_t i = p; i < q; ++i) {
			const gaba_path_section_t *s = &_aln(reg[i])->sec[0];
			int32_t lb = s->bpos, ub = s->bpos + s->blen, span = ub - lb;

			for (uint64_t j = 0; j < p; ++j) {
				const gaba_path_section_t *t = &_aln(reg[j])->sec[0];
				if (t->bpos + t->blen < ub) lb = MAX2(lb, t->bpos + t->blen);
				else ub = MIN2(ub, t->bpos);
				if (2*(ub - lb) < span) {	// covered by j
					q--; _swap_128(i, q); i--; reg[q].u32[0] |= 0x100<<16;
					goto _loop_tail;
				}
			}
			max = MAX2(max, ((uint64_t)(2*(ub - lb) - span)<<32) | i);
		_loop_tail:;
		}
		if (max&0xffffffff) { _swap_128(p, max&0xffffffff); reg[p].u32[0] |= 0x800<<16; }	// move to head, mark supplementary
	}
	p = MIN2(p, q);
	#undef _swap_128

	// extract shortest repeat element
	int32_t usc = 0, lsc = INT32_MAX, tsc = 0;
	for (uint64_t i = p; i < n_reg; ++i) {
		usc = MAX2(usc, _aln(reg[i])->score);
		lsc = MIN2(lsc, _aln(reg[i])->score);
		tsc += _aln(reg[i])->score;
	}
	lsc = (lsc==INT32_MAX)? 0 : lsc;

	// calc mapq for primary and supplementary alignments
	double tpc = 1.0;
	for (uint64_t i = 0; i < p; ++i) {
		uint32_t score = _aln(reg[0])->score;
		double elen = (double)gaba_plen(_aln(reg[i])->sec) / 2.0, pid = 1.0 - (double)(elen * opt->m - score) / (double)(opt->m + opt->x) / elen;
		double ec = 2.0 / (pid * (double)(opt->m + opt->x) - (double)opt->x), ulen = ec * (score - usc), pe = 1.0 / (ulen * ulen + (double)(n_reg-p+1));
		reg[i].u32[0] |= _clip(-10.0 * MAPQ_COEF * log10(pe)); tpc *= 1.0 - pe;
	}

	// calc mapq for secondary (repetitive) alignments
	double tpe = MIN2(1.0 - tpc, 1.0);
	for (uint64_t i = p; i < n_reg; ++i) {
		reg[i].u32[0] |= _clip(-10.0 * MAPQ_COEF * log10(1.0 - tpe * (double)(_aln(reg[i])->score - lsc + 1) / (double)tsc));
	}
	return p;
}
#else
static uint64_t mm_post_map(const mm_mapopt_t *opt, uint32_t n_reg, mm128_t *reg)
{
	uint32_t i, score = _aln(reg[0])->score, min = score * opt->min_ratio;
	uint32_t ssc = (n_reg>1)? _aln(reg[1])->score : 0, bsc = (n_reg>1)? _aln(reg[n_reg-1])->score : 0, tsc = 0;
	double elen = (double)gaba_plen(_aln(reg[0])->sec) / 2.0, pid = 1.0 - (double)(elen * opt->m - score) / (double)(opt->m + opt->x) / elen;
	double ec = 2.0 / (pid * (double)(opt->m + opt->x) - (double)opt->x), ulen = ec * (score - ssc), pe = 1.0 / (ulen * ulen + (double)n_reg);

	for (i = 1; i < n_reg; ++i) tsc += _aln(reg[i])->score - bsc + 1;
	reg[0].u32[0] |= _clip(-10.0 * MAPQ_COEF * log10(pe));
	for (i = 1; i < n_reg && (score = _aln(reg[i])->score) >= min; ++i)
		reg[i].u32[0] |= (0x100<<16) | _clip(-10.0 * MAPQ_COEF * log10(1.0 - pe * (double)(score - bsc + 1) / (double)tsc));
	return 1;
}
#endif
static uint64_t mm_post_ava(const mm_mapopt_t *opt, uint32_t n_reg, mm128_t *reg)
{
	uint32_t i, score, min = _aln(reg[0])->score * opt->min_ratio;
	for (i = 0; i < n_reg && (score = _aln(reg[i])->score) >= min; ++i) {
		double elen = (double)gaba_plen(_aln(reg[i])->sec) / 2.0, pid = 1.0 - (double)(elen * opt->m - score) / (double)(opt->m + opt->x) / elen;
		double ec = 2.0 / (pid * (double)(opt->m + opt->x) - (double)opt->x), ulen = ec * score, pe = 1.0 / (ulen + 1);
		reg[i].u32[0] |= _clip(-10.0 * MAPQ_COEF * log10(pe)) | ((i == 0? 0 : 0x800)<<16);
	}
	return n_reg;
}
#undef _clip
#undef _aln

static const mm128_t *mm_align_seq(
	mm_tbuf_t *b, const mm_mapopt_t *opt, const mm_idx_t *mi, uint32_t l_seq, const uint8_t *seq, uint32_t qid,
	uint32_t n_occ, uint64_t *occ, lmm_t *lmm, uint64_t *n_reg)
{
	int32_t min = opt->min;
	const int32_t ofs = 0x40000000;
	const uint32_t mask = 0x7fffffff;
	const uint32_t org = (opt->flag&(MM_AVA|MM_COMP))==MM_AVA? ofs : 0, thresh = (opt->flag&MM_AVA? 1 : 0) + org;

	mm128_v reg = { 0 };
	// prepare section info for alignment
	uint8_t tail[96];
	const uint8_t *lim = (const uint8_t*)0x800000000000;
	gaba_section_t qf, qr, mg;
	if (l_seq == 0) return NULL;

	memset(tail, 0, 96);
	qf.id = 0; qf.len = l_seq; qf.base = (const uint8_t*)seq;
	qr.id = 1; qr.len = l_seq; qr.base = gaba_rev((const uint8_t*)seq+l_seq-1, lim);
	mg.id = 4; mg.len = 32; mg.base = tail+32;
	b->coef.n = b->resc.n = b->intv.n = 0; kh_clear(b->pos);
	for (uint64_t i = 0, j = 0; reg.n == 0 && i < n_occ; ++i) {
		if (i == 0) {
			mm_collect(mi, l_seq, seq, qid, occ[n_occ-1], occ[0], org, thresh, &b->resc, &b->coef);
		} else {
			if (i == 1) radix_sort_128x(b->resc.a, b->resc.a + b->resc.n);
			for (uint64_t k = 0; k < b->coef.n; ++k) b->coef.a[k].u32[2] &= mask;
			for (; j < b->resc.n && b->resc.a[j].u32[1] <= occ[i]; ++j)
				mm_expand(b->resc.a[j].u32[1], (const v2u32_t*)b->resc.a[j].u64[1], qid, b->resc.a[j].u32[0], org, thresh, &b->coef);
		}
		radix_sort_128x(b->coef.a, b->coef.a + b->coef.n);

		mm_chain(b->coef.n, b->coef.a, opt->llim, opt->hlim, min>>2, &b->intv);
		radix_sort_64x(b->intv.a, b->intv.a + b->intv.n);
		if (b->intv.a == 0) continue;
		for (uint64_t k = 0; k < b->intv.n && ofs - (int32_t)b->intv.a[k].x[0] > min>>2; ++k) {
			const gaba_alignment_t *a = 0;
			mm128_t *t = &b->coef.a[b->coef.n], r;
			uint64_t p = b->coef.n-b->intv.a[k].x[1], q = p, l = 0, m = ofs - b->intv.a[k].x[0];
			mm_idx_seq_t *ref = &mi->s.a[(t-p)->u32[1]-mi->base_rid];
			do {
				a = mm_extend(ref, q, t-q, mi->k, min, opt->sidx, opt->eidx, b->dp, &qf, &qr, &mg, b->pos, lmm);
				if (p == q && a == 0) break;	// skip chain if first extension did not result in meaningful alignment
				r.u32[2] = (t-p)->u32[2] & mask, r.u32[3] = (t-p)->u32[3];
				if (a != 0) {
					mm128_t *s;
					kv_pushp(mm128_t, reg, &s);
					s->u32[0] = (a->sec->bid&0x01? 0 : 0x10)<<16;
					s->u32[1] = ofs - a->score;
					s->u64[1] = (uintptr_t)a;
					r.u32[2] = ref->l_seq-a->sec->apos; r.u32[3] = (a->sec->bid&0x01)? l_seq-a->sec->bpos : -a->sec->bpos+mi->k+1;
					min = MAX2(min, opt->min_ratio*a->score);
					l += a->sec->alen+a->sec->blen;
				}
				q = p;
			} while (l < m && (p -= mm_rescue(ref, p, t-p, &r, opt->llim, opt->hlim, opt->elim, opt->blim, opt->eidx)) < q - 1);
		}
	}
	*n_reg = reg.n;
	if (reg.n == 0) return NULL;
	#if 0
	radix_sort_128x(reg.a, reg.a + reg.n);
	*n_reg = ((opt->flag & MM_AVA)? mm_post_ava : mm_post_map)(opt, reg.n, reg.a);
	for (uint64_t i = *n_reg; i < reg.n; ++i) { lmm_free(lmm, (void*)reg.a[i].u64[1]); }
	#else
	radix_sort_128x(reg.a, reg.a + reg.n);	// sort by score in reverse order
	*n_reg = reg.n = mm_prune_regs(opt, lmm, reg.n, reg.a);	// prune alignment whose score is less than min_score threshold
	*n_reg |= ((opt->flag & MM_AVA)? mm_post_ava : mm_post_map)(opt, reg.n, reg.a)<<32;
	#endif
	return reg.a;
}

#undef _s
#undef _m

// flush the buffer if there is no room for (margin + 1) bytes
#define _flush(_buf, _margin) { \
	if (_unlikely((uintptr_t)(_buf)->p + (_margin) >= (uintptr_t)(_buf)->tail)) { \
		fwrite((_buf)->base, sizeof(uint8_t), (_buf)->p - (_buf)->base, stdout); \
		(_buf)->p = (_buf)->base; \
	} \
}
#define _put(_buf, _c) { \
	_flush(_buf, 1); \
	*(_buf)->p++ = (_c); \
}
#define _putfi(type, _buf, _n, _c) ({ \
	uint64_t _b = 0; \
	type _m = (type)(_n); int64_t _i = 0; \
	while (_m) { _b <<= 4; _b += _m % 10, _m /= 10; _i++; } \
	_i += (_i==0); _i += (_c - _i + 1 > 0)? _c - _i + 1 : 0; \
	_flush(_buf, _i + 1); \
	for (int64_t _j = _i; _j > (_c); _j--) { *(_buf)->p++ = (_b&0x0f) + '0'; _b>>=4; } \
	*(_buf)->p++ = '.'; \
	for (int64_t _j = (_c); _j > 0; _j--) { *(_buf)->p++ = (_b&0x0f) + '0'; _b>>=4; } \
	_i; \
})
#define _puti(type, _buf, _n) ({ \
	uint64_t _b = 0; \
	type _m = (type)(_n); int64_t _i = 0; \
	while (_m) { _b <<= 4; _b += _m % 10, _m /= 10; _i++; } \
	_i += (_i==0); \
	_flush(_buf, _i); \
	for (int64_t _j = _i; _j > 0; _j--) { *(_buf)->p++ = (_b&0x0f) + '0'; _b>>=4; } \
	_i; \
})
#define _putn(_buf, _n) _puti(uint32_t, _buf, _n)
#define _puts(_buf, _s) { \
	for (const uint8_t *_q = (const uint8_t*)(_s); *_q; _q++) { _put(_buf, *_q); } \
}
#define _putsk(_buf, _s) { \
	const uint64_t _l = strlen(_s); \
	_flush(_buf, _l); \
	for (const uint8_t *_q = (const uint8_t*)(_s), *_t = _q + (_l); _q < _t; _q++) *(_buf)->p++ = *(_q); \
}
#define _putsn(_buf, _s, _l) { \
	const uint8_t *_q = (const uint8_t *)(_s); \
	for (const uint8_t *_t = _q + ((_l) & ~0x1fULL); _q < _t; _q += 32) { \
		_flush(_buf, 32); \
		_storeu_v32i8((_buf)->p, _loadu_v32i8(_q)); (_buf)->p += 32; \
	} \
	_flush(_buf, 32); \
	_storeu_v32i8((_buf)->p, _loadu_v32i8(_q)); (_buf)->p += (_l) & 0x1f; \
}
#define _putsnr(_buf, _s, _l) { \
	const uint8_t *_q = (const uint8_t *)(_s) + (_l) - 32; \
	for (; _q >= _s; _q -= 32) { \
		_flush(_buf, 32); \
		_storeu_v32i8((_buf)->p, _swap_v32i8(_loadu_v32i8(_q))); (_buf)->p += 32; \
	} \
	_flush(_buf, 32); \
	_storeu_v32i8((_buf)->p, _swap_v32i8(_loadu_v32i8(_q))); (_buf)->p += (_l) & 0x1f; \
}
#define _putsnt(_buf, _s, _l, _table) { \
	v32i8_t register _conv = _from_v16i8_v32i8(_loadu_v16i8(_table)); \
	v32i8_t _r, _b; \
	const uint8_t *_q = (const uint8_t *)(_s); \
	for (const uint8_t *_t = _q + ((_l) & ~0x1fULL); _q < _t; _q += 32) { \
		_flush(_buf, 32); \
		_r = _loadu_v32i8(_q); _b = _shuf_v32i8(_conv, _r); \
		_storeu_v32i8((_buf)->p, _b); (_buf)->p += 32; \
	} \
	_flush(_buf, 32); \
	_r = _loadu_v32i8(_q); _b = _shuf_v32i8(_conv, _r); \
	_storeu_v32i8((_buf)->p, _b); (_buf)->p += (_l) & 0x1f; \
}
#define _putsntr(_buf, _s, _l, _table) { \
	v32i8_t register _conv = _from_v16i8_v32i8(_loadu_v16i8(_table)); \
	v32i8_t _r, _b; \
	const uint8_t *_q = (const uint8_t *)(_s) + (_l) - 32; \
	for (; _q >= _s; _q -= 32) { \
		_flush(_buf, 32); \
		_r = _loadu_v32i8(_q); _b = _shuf_v32i8(_conv, _swap_v32i8(_r)); \
		_storeu_v32i8((_buf)->p, _b); (_buf)->p += 32; \
	} \
	_flush(_buf, 32); \
	_r = _loadu_v32i8(_q); _b = _shuf_v32i8(_conv, _swap_v32i8(_r)); \
	_storeu_v32i8((_buf)->p, _b); (_buf)->p += (_l) & 0x1f; \
}
#define _putd(b, _id)		_put(b, ((_id)&0x01)? '+' : '-');

#define MM_RG			( 0 )		// Z: read group
#define MM_CO			( 1 )		// Z: comment
#define MM_NH			( 2 )		// i: #hits
#define MM_IH 			( 3 )		// i: index of the record in the hits
#define MM_AS			( 4 )		// i: score
#define MM_XS			( 5 )		// i: suboptimal score
#define MM_NM 			( 6 )		// i: editdist to the reference
#define MM_SA			( 7 )		// Z: supplementary records
#define MM_MD			( 8 )		// Z: mismatch positions

static uint64_t mm_print_num(mm_align_t *b, uint8_t type, const uint8_t *p)
{
	switch (type) {
		case 'a': _put(b, *p); return 1;
		case 'c': _puti(int8_t, b, (int8_t)*p); return 1;
		case 'C': _puti(uint8_t, b, (uint8_t)*p); return 1;
		case 's': _puti(int16_t, b, *((int16_t*)p)); return 2;
		case 'S': _puti(uint16_t, b, *((uint16_t*)p)); return 2;
		case 'i': _puti(int32_t, b, *((int32_t*)p)); return 4;
		case 'I': _puti(uint32_t, b, *((uint32_t*)p)); return 4;
		case 'f': {
			char buf[32]; uint64_t l = sprintf(buf, "%f", *((float*)p));
			for (uint64_t i = 0; i < l; ++i) { _put(b, buf[i]); }
			return 4;
		}
	}
	return 1;
}

#define _t(b)			_put(b, '\t')
#define _c(b)			_put(b, ',')
#define _sp(b)			_put(b, ' ')
#define _cr(b)			_put(b, '\n')
static void mm_restore_tags(mm_align_t *b, const bseq_t *t)
{
	static const uint8_t tag_size[32] = {
		0, 1, 0xfe, 1,  0, 0, 4, 0,  0xfe, 4, 0, 0,  0, 0, 0, 0,
		0, 0, 0, 2,     0, 0, 0, 0,  0, 0, 0xff, 0,  0, 0, 0, 0,
	};
	const uint8_t *p = t->tag, *tail = p + t->l_tag;
	while (p < tail) {
		_t(b); _put(b, p[0]); _put(b, p[1]); _put(b, ':'); _put(b, p[2]); _put(b, ':');
		if (p[2] == 'Z') {
			p += 3; while(*p) { _put(b, *p); p++; } p++;
		} else if (tag_size[p[2]&0x1f] == 0xfe) {
			uint8_t type = p[3];
			uint32_t len = *((uint32_t*)&p[4]); p += 8;
			for (uint64_t i = 0; i < len; ++i) { p += mm_print_num(b, type, p); _put(b, ','); }
		} else {
			p += mm_print_num(b, p[2], p + 3) + 3;
		}
	}
	return;
}

static void mm_print_header_sam(mm_align_t *b)
{
	const mm_idx_t *mi = b->mi;
	_putsk(b, "@HD\tVN:1.0\tSO:unsorted\n");
	for (uint64_t i = 0; i < mi->s.n; ++i) {
		_putsk(b, "@SQ\tSN:"); _putsn(b, mi->s.a[i].name, mi->s.a[i].l_name);
		_putsk(b, "\tLN:"); _putn(b, mi->s.a[i].l_seq);
		_cr(b);
	}
	if (b->opt->flag & 0x01ULL<<MM_RG) { _puts(b, b->opt->rg_line); _cr(b); }
	_putsk(b, "@PG\tID:minialign\tPN:minialign\n");
	return;
}

static void mm_print_unmapped_sam(mm_align_t *b, const bseq_t *t)
{
	_putsn(b, t->name, t->l_name); _putsk(b, "\t4\t*\t0\t0\t*\t*\t0\t0\t");
	_putsnt(b, t->seq, t->l_seq, "NACMGRSVTWYHKDBN");
	_t(b);
	if (b->opt->flag&MM_KEEP_QUAL && t->qual[0] != '\0') { _putsn(b, t->qual, t->l_seq); }
	else { _put(b, '*'); }
	mm_restore_tags(b, t); _cr(b);
	return;
}

static int mm_cigar_printer(void *_b, int64_t len, char c)
{
	mm_align_t *b = (mm_align_t*)_b;
	if (len < 41) {
		_flush(b, 16);
		*b->p++ = (b->conv[len-1]&0x0f) + '0';
		*b->p++ = (b->conv[len-1]>>4) + '0';
		b->p -= len<10;
		*b->p++ = c;
	} else {
		len = _putn(b, (uint32_t)len); _put(b, c);
	}
	return len+1;
}

#define _aln(x)		( (const gaba_alignment_t*)(x).u64[1] )
#define _mapq(x)	( (x).u32[0] & 0xffff )
#define _flag(x)	( (x).u32[0]>>16 )
static void mm_print_mapped_sam_core(mm_align_t *b, const mm_idx_seq_t *r, const bseq_t *t, const gaba_alignment_t *a, uint16_t mapq, uint16_t flag)
{
	const gaba_path_section_t *s = &a->sec[0];
	uint32_t rs = r->l_seq-s->apos-s->alen, hl = t->l_seq-s->bpos-s->blen, tl = s->bpos;
	uint32_t qs = (flag&0x900)? hl : 0, qe = t->l_seq - ((flag&0x900)? tl : 0);

	_putsn(b, t->name, t->l_name); _t(b); _putn(b, flag); _t(b); _putsn(b, r->name, r->l_name); _t(b);
	_putn(b, rs+1); _t(b); _putn(b, mapq>>MAPQ_DEC); _t(b);
	if (hl) { _putn(b, hl); _put(b, (flag&0x900)? 'H' : 'S'); }
	gaba_dp_print_cigar_reverse(mm_cigar_printer, b, a->path->array, 0, a->path->len);
	if (tl) { _putn(b, tl); _put(b, (flag&0x900)? 'H' : 'S'); }
	_putsk(b, "\t*\t0\t0\t");

	if (flag&0x10) { _putsntr(b, &t->seq[t->l_seq-qe], qe-qs, "NTGKCYSBAWRDMHVN"); }
	else { _putsnt(b, &t->seq[qs], qe-qs, "NACMGRSVTWYHKDBN"); }
	_t(b);
	if (b->opt->flag&MM_KEEP_QUAL && t->qual[0] != '\0') {
		if (flag&0x10) { _putsnr(b, &t->qual[t->l_seq-qe], qe-qs); }
		else { _putsn(b, &t->qual[qs], qe-qs); }
	} else {
		_put(b, '*');
	}
	return;
}

static void mm_print_sam_supp(mm_align_t *b, const mm_idx_seq_t *r, const bseq_t *t, const gaba_alignment_t *a, uint16_t mapq, uint16_t flag)
{
	// rname,pos,strand,CIGAR,mapQ,NM;
	const gaba_path_section_t *s = &a->sec[0];
	uint32_t rs = r->l_seq-s->apos-s->alen, hl = t->l_seq-s->bpos-s->blen, tl = s->bpos;
	_putsn(b, r->name, r->l_name); _c(b); _putn(b, rs+1); _c(b); _put(b, (flag&0x04)? '-' : '+'); _c(b);
	if (hl) { _putn(b, hl); _put(b, (flag&0x900)? 'H' : 'S'); }
	gaba_dp_print_cigar_reverse(mm_cigar_printer, b, a->path->array, 0, a->path->len);
	if (tl) { _putn(b, tl); _put(b, (flag&0x900)? 'H' : 'S'); } _c(b);
	_putn(b, mapq); _c(b); _putn(b, a->xcnt + a->gecnt); _put(b, ';');
}

static inline
uint64_t parse_load_uint64(
	uint64_t const *ptr,
	int64_t pos)
{
	int64_t rem = pos & 63;
	uint64_t a = (ptr[pos>>6]>>rem) | ((ptr[(pos>>6) + 1]<<(63 - rem))<<1);
	return(a);
}

void mm_print_sam_md(mm_align_t *b, const mm_idx_seq_t *r, const bseq_t *t, const gaba_alignment_t *a, uint16_t flag)
{
	/*
	 * print MD tag
	 * note: this operation requires head-to-tail comparison of two sequences,
	 *       which may result in several percent performance degradation.
	 */
	const gaba_path_section_t *s = &a->sec[0];
	uint32_t rs = r->l_seq-s->apos-s->alen, qs = t->l_seq-s->bpos-s->blen;

	#define _load(_ptr, _pos) ({ \
		uint64_t _rem = (_pos) & 0x3f; \
		((_ptr)[(_pos)>>6]>>_rem) | (((_ptr)[((_pos)>>6) + 1]<<(63 - _rem))<<1); \
	})

	_puts(b, "\tMD:Z:");
	uint64_t const *p = (uint64_t const *)((uint64_t)(a->path->array - 1) & ~(sizeof(uint64_t) - 1));
	int64_t pos = a->path->len;

	uint64_t dir = (flag&0x10)? 1 : 0;
	static uint8_t const comp[16] __attribute__(( aligned(16) )) = {
		0x00, 0x08, 0x04, 0x0c, 0x02, 0x0a, 0x06, 0x0e,
		0x01, 0x09, 0x05, 0x0d, 0x03, 0x0b, 0x07, 0x0f
	};
	static uint8_t const id[16] __attribute__(( aligned(16) )) = {
		0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
		0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f
	};
	const v32i8_t cv = _from_v16i8_v32i8(_load_v16i8(dir? comp : id));
	const uint8_t *rp = &r->seq[rs], *rb = rp, *qp = &t->seq[dir? (uint64_t)t->l_seq-qs-32 : qs];
	while (pos > 0) {
		// suppose each indel block is shorter than 32 bases
		uint64_t arr = _load(p, pos), cnt = lzcnt(~arr);	// count #ins
		pos -= cnt; qp += dir? -cnt : cnt;
		if (((((int64_t)arr)>>62)^0x01) > 0) {				// is_del
			pos -= (cnt = lzcnt(arr) - 1);
			_putn(b, (int32_t)(rp-rb)); _put(b, '^'); _putsnt(b, rp, cnt, "NACMGRSVTWYHKDBN");
			rp += cnt; rb = rp;
		}
		// match or mismatch
		uint64_t acnt = 32;
		while (acnt == 32) {
			arr = _load(p, pos); acnt = lzcnt(arr^0x5555555555555555)>>1;
			v32i8_t rv = _shuf_v32i8(cv, _loadu_v32i8(rp)), qv = _loadu_v32i8(qp);
			if (dir) qv = _swap_v32i8(qv);
			uint64_t mmask = (uint64_t)((v32_masku_t){ .mask = _mask_v32i8(_eq_v32i8(rv, qv)) }).all;
			uint64_t mcnt = MIN2(pos>>1, tzcnt(~mmask));
			pos -= 2*mcnt; rp += mcnt; qp += dir? -mcnt : mcnt;
			if (mcnt >= acnt) continue;
			_putn(b, (int32_t)(rp-rb));
			pos -= 2*(cnt = MIN2(tzcnt(mmask>>mcnt), acnt - mcnt)); qp += dir? -cnt : cnt;
			for (uint64_t i = 0; i < cnt-1; ++i) { _put(b, "NACMGRSVTWYHKDBN"[*rp]); _put(b, '0'); rp++; }
			_put(b, "NACMGRSVTWYHKDBN"[*rp]); rp++;
			rb = rp;
		}
	}
	if (rp-rb) { _putn(b, (int32_t)(rp-rb)); }
	#undef _load
	return;
}

static void mm_print_mapped_sam(mm_align_t *b, const bseq_t *t, uint64_t n_reg, const mm128_t *reg)
{
	const uint64_t f = b->opt->flag;
	const mm_idx_t *mi = b->mi;
	const uint32_t n_uniq = n_reg>>32;
	n_reg &= 0xffffffff;
	for (uint64_t j = 0; j < ((f & MM_OMIT_REP)? n_uniq : n_reg); ++j) {
		uint16_t mapq = _mapq(reg[j]), flag = _flag(reg[j]);
		const gaba_alignment_t *a = _aln(reg[j]);
		const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];

		// print body
		mm_print_mapped_sam_core(b, r, t, a, mapq, flag);

		// print tags
		if (f & 0x01ULL<<MM_RG) { _putsk(b, "\tRG:Z:"); _puts(b, b->opt->rg_id); }
		if (f & 0x01ULL<<MM_NH) { _putsk(b, "\tNH:i:"); _putn(b, n_reg); }	// (f & MM_OMIT_REP)? n_uniq : n_reg
		if (f & 0x01ULL<<MM_IH) { _putsk(b, "\tIH:i:"); _putn(b, j); }
		if (f & 0x01ULL<<MM_AS) { _putsk(b, "\tAS:i:"); _putn(b, a->score); }
		if ((f & 0x01ULL<<MM_XS) && (flag & 0x900) == 0) { _putsk(b, "\tXS:i:"); _putn(b, n_reg>1? _aln(reg[1])->score : 0); }
		if (f & 0x01ULL<<MM_NM) { _putsk(b, "\tNM:i:"); _putn(b, a->xcnt + a->gecnt); }
		if (f & 0x01ULL<<MM_SA && (flag & 0x900) == 0 && n_uniq > 1) {
			_putsk(b, "\tSA:Z:");
			for (uint64_t k = 1; k < n_uniq; ++k)
				mm_print_sam_supp(b, &mi->s.a[(_aln(reg[k])->sec->aid>>1) - mi->base_rid], t, _aln(reg[k]), _mapq(reg[k]), _flag(reg[k]));
			j = n_uniq;		// skip supplementary records when SA tag is enabled
		}
		if (f & 0x01ULL<<MM_MD) mm_print_sam_md(b, r, t, a, flag);
		if ((flag & 0x900) == 0) mm_restore_tags(b, t);
		_cr(b);
	}
	return;
}

static void mm_print_mapped_maf_core(mm_align_t *b, const mm_idx_seq_t *r, const bseq_t *t, const gaba_alignment_t *a, uint16_t mapq, uint16_t flag)
{
	#define _load(_ptr, _pos) ({ \
		uint64_t _rem = (_pos) & 0x3f; \
		((_ptr)[(_pos)>>6]>>_rem) | (((_ptr)[((_pos)>>6) + 1]<<(63 - _rem))<<1); \
	})
	#define _putsnf32(_buf, _q, _len) { \
		_flush(_buf, 32); \
		v32i8_t register _conv = _from_v16i8_v32i8(_loadu_v16i8("NACMGRSVTWYHKDBN")); \
		_storeu_v32i8((_buf)->p, _shuf_v32i8(_conv, _loadu_v32i8(_q))); (_buf)->p += (_len); \
	}
	#define _putsnr32(_buf, _q, _len) { \
		_flush(_buf, 32); \
		v32i8_t register _conv = _from_v16i8_v32i8(_loadu_v16i8("NTGKCYSBAWRDMHVN")); \
		_storeu_v32i8((_buf)->p, _shuf_v32i8(_conv, _swap_v32i8(_loadu_v32i8(_q)))); (_buf)->p += (_len); \
	}
	#define _putgn32(_buf, _len) { \
		_flush(_buf, 32); \
		_storeu_v32i8((_buf)->p, _loadu_v32i8("--------------------------------")); (_buf)->p += (_len); \
	}
	const gaba_path_section_t *s = &a->sec[0];
	uint32_t rs = r->l_seq-s->apos-s->alen, qs = t->l_seq-s->bpos-s->blen;
	uint64_t arr, cnt;
	_put(b, 'a'); _sp(b); _putsk(b, "score="); _putn(b, a->score); _cr(b);	// header

	// ref
	_put(b, 's'); _sp(b); _putsn(b, r->name, r->l_name); _sp(b); _putn(b, rs); _sp(b);
	_putn(b, s->alen); _sp(b); _put(b, '+'); _sp(b); _putn(b, r->l_seq); _sp(b);

	uint64_t const *p = (uint64_t const *)((uint64_t)(a->path->array - 1) & ~(sizeof(uint64_t) - 1));
	int64_t pos = a->path->len;
	const uint8_t *rp = &r->seq[rs];
	while (pos > 0) {
		arr = _load(p, pos); pos -= (cnt = lzcnt(~arr)); _putgn32(b, cnt);
		arr = _load(p, pos); pos -= (cnt = lzcnt(arr) - 1); _putsnf32(b, rp, cnt); rp += cnt;
		do {
			arr = _load(p, pos); cnt = lzcnt(arr^0x5555555555555555)>>1; pos -= 2*cnt;
			_putsnf32(b, rp, cnt); rp += cnt;
		} while (cnt == 32);
	}
	_cr(b);

	// query
	_put(b, 's'); _sp(b); _putsn(b, t->name, t->l_name); _sp(b); _putn(b, qs); _sp(b);
	_putn(b, s->blen); _sp(b); _putd(b, s->bid); _sp(b); _putn(b, t->l_seq); _sp(b);

	p = (uint64_t const *)((uint64_t)(a->path->array - 1) & ~(sizeof(uint64_t) - 1));
	pos = a->path->len;
	if (flag&0x10) {
		const uint8_t *qp = &t->seq[(uint64_t)t->l_seq-qs-32];
		while (pos > 0) {
			arr = _load(p, pos); pos -= (cnt = lzcnt(~arr)); _putsnr32(b, qp, cnt); qp -= cnt;
			arr = _load(p, pos); pos -= (cnt = lzcnt(arr) - 1); _putgn32(b, cnt);
			do {
				arr = _load(p, pos); cnt = lzcnt(arr^0x5555555555555555)>>1; pos -= 2*cnt;
				_putsnr32(b, qp, cnt); qp -= cnt;
			} while (cnt == 32);
		}
	} else {
		const uint8_t *qp = &t->seq[qs];
		while (pos > 0) {
			arr = _load(p, pos); pos -= (cnt = lzcnt(~arr)); _putsnf32(b, qp, cnt); qp += cnt;
			arr = _load(p, pos); pos -= (cnt = lzcnt(arr) - 1); _putgn32(b, cnt);
			do {
				arr = _load(p, pos); cnt = lzcnt(arr^0x5555555555555555)>>1; pos -= 2*cnt;
				_putsnf32(b, qp, cnt); qp += cnt;
			} while (cnt == 32);
		}
	}
	_cr(b); _cr(b);
	#undef _load
	#undef _putsnf32
	#undef _putsnr32
	#undef _putgn32
	return;
}

static void mm_print_mapped_maf(mm_align_t *b, const bseq_t *t, uint64_t n_reg, const mm128_t *reg)
{
	const uint64_t f = b->opt->flag;
	const mm_idx_t *mi = b->mi;
	const uint32_t n_uniq = n_reg>>32;
	n_reg &= 0xffffffff;
	for (uint64_t j = 0; j < ((f & MM_OMIT_REP)? n_uniq : n_reg); ++j) {
		uint16_t mapq = _mapq(reg[j]), flag = _flag(reg[j]);
		const gaba_alignment_t *a = _aln(reg[j]);
		const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];
		mm_print_mapped_maf_core(b, r, t, a, mapq, flag);
	}
}

// qname rname idt len #x #gi qs qe rs re e-value bitscore
static void mm_print_mapped_blast6(mm_align_t *b, const bseq_t *t, uint64_t n_reg, const mm128_t *reg)
{
	const mm_idx_t *mi = b->mi;
	for (uint64_t j = 0; j < (uint32_t)n_reg; ++j) {
		const gaba_alignment_t *a = _aln(reg[j]);
		const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];
		const gaba_path_section_t *s = &a->sec[0];
		int32_t dcnt = (a->path->len - a->gecnt)>>1, slen = dcnt + a->gecnt;
		int32_t mid = 100000.0 * (double)(dcnt - a->xcnt) / (double)slen;	// percent identity
		uint32_t rs = s->bid&0x01? r->l_seq-s->apos-s->alen+1 : r->l_seq-s->apos;
		uint32_t re = s->bid&0x01? r->l_seq-s->apos : r->l_seq-s->apos-s->alen+1;
		uint32_t qs = t->l_seq-s->bpos-s->blen+1, qe = t->l_seq-s->bpos;

		_putsn(b, t->name, t->l_name); _t(b); _putsn(b, r->name, r->l_name); _t(b);
		_putfi(int32_t, b, mid, 3); _t(b); _putn(b, slen); _t(b); _putn(b, a->xcnt); _t(b); _putn(b, a->gicnt); _t(b);
		_putn(b, qs); _t(b); _putn(b, qe); _t(b); _putn(b, rs); _t(b); _putn(b, re); _t(b);

		double bit = 1.85 * (double)a->score - 0.02;	// fixme, lambda and k should be calcd from the scoring params
		int32_t e = 1000.0 * (double)r->l_seq * (double)t->l_seq * pow(2.0, -bit);	// fixme, lengths are not corrected
		_putfi(int32_t, b, e, 3); _t(b); _putn(b, (int32_t)bit); _cr(b);
	}
	return;
}

// qname rname qd rd score idt rs re rl qs qe ql #cells
static void mm_print_mapped_blasr1(mm_align_t *b, const bseq_t *t, uint64_t n_reg, const mm128_t *reg)
{
	const mm_idx_t *mi = b->mi;
	for (uint64_t j = 0; j < (uint32_t)n_reg; ++j) {
		const gaba_alignment_t *a = _aln(reg[j]);
		const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];
		const gaba_path_section_t *s = &a->sec[0];
		int32_t dcnt = (a->path->len - a->gecnt)>>1, slen = dcnt + a->gecnt;
		int32_t mid = 1000000.0 * (double)(dcnt - a->xcnt) / (double)slen;	// percent identity
		uint32_t rs = s->bid&0x01? r->l_seq-s->apos-s->alen : s->apos, re = rs+s->alen;
		uint32_t qs = s->bid&0x01? t->l_seq-s->bpos-s->blen : s->bpos, qe = qs+s->blen;

		_putsn(b, t->name, t->l_name); _sp(b); _putsn(b, r->name, r->l_name); _sp(b);
		_put(b, '0'); _sp(b); _put(b, s->bid&0x01? '0' : '1'); _sp(b);
		_put(b, '-'); _putn(b, a->score); _sp(b);	// score in negative
		_putfi(int32_t, b, mid, 4); _sp(b);
		_putn(b, rs); _sp(b); _putn(b, re); _sp(b); _putn(b, r->l_seq); _sp(b);
		_putn(b, qs); _sp(b); _putn(b, qe); _sp(b); _putn(b, t->l_seq); _sp(b);
		_put(b, '0'); _cr(b);	// #cells is always 0
	}
	return;
}

// qname rname score idt qd qs qe ql rd rs re rl mapq
static void mm_print_mapped_blasr4(mm_align_t *b, const bseq_t *t, uint64_t n_reg, const mm128_t *reg)
{
	const mm_idx_t *mi = b->mi;
	for (uint64_t j = 0; j < (uint32_t)n_reg; ++j) {
		uint16_t mapq = _mapq(reg[j]);
		const gaba_alignment_t *a = _aln(reg[j]);
		const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];
		const gaba_path_section_t *s = &a->sec[0];
		int32_t dcnt = (a->path->len - a->gecnt)>>1, slen = dcnt + a->gecnt;
		int32_t mid = 1000000.0 * (double)(dcnt - a->xcnt) / (double)slen;	// percent identity
		uint32_t rs = s->bid&0x01? r->l_seq-s->apos-s->alen : s->apos, re = rs+s->alen;
		uint32_t qs = s->bid&0x01? t->l_seq-s->bpos-s->blen : s->bpos, qe = qs+s->blen;

		_putsn(b, t->name, t->l_name); _sp(b); _putsn(b, r->name, r->l_name); _sp(b);
		_put(b, '-'); _putn(b, a->score); _sp(b); _putfi(int32_t, b, mid, 4); _sp(b);	// add '-' before the number to negate score
		_put(b, '0'); _sp(b); _putn(b, qs); _sp(b); _putn(b, qe); _sp(b); _putn(b, t->l_seq); _sp(b);
		_put(b, s->bid&0x01? '0' : '1'); _sp(b); _putn(b, rs); _sp(b); _putn(b, re); _sp(b); _putn(b, r->l_seq); _sp(b);
		_putn(b, (mapq>>(MAPQ_DEC - 2)) + ((mapq>>(MAPQ_DEC + 2)) & ~0x01)); _cr(b);
	}
	return;
}

// qname ql qs qe qd rname rl rs re #m block_len mapq
static void mm_print_mapped_paf(mm_align_t *b, const bseq_t *t, uint64_t n_reg, const mm128_t *reg)
{
	const mm_idx_t *mi = b->mi;
	for (uint64_t j = 0; j < (uint32_t)n_reg; ++j) {
		uint16_t mapq = _mapq(reg[j]);
		const gaba_alignment_t *a = _aln(reg[j]);
		const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];
		const gaba_path_section_t *s = &a->sec[0];
		uint32_t dcnt = (a->path->len - a->gecnt)>>1;
		uint32_t rs = r->l_seq-s->apos-s->alen, re = r->l_seq-s->apos, qs = t->l_seq-s->bpos-s->blen, qe = t->l_seq-s->bpos;
		_putsn(b, t->name, t->l_name); _t(b); _putn(b, t->l_seq); _t(b); _putn(b, qs); _t(b); _putn(b, qe); _t(b); _putd(b, s->bid); _t(b);
		_putsn(b, r->name, r->l_name); _t(b); _putn(b, r->l_seq); _t(b); _putn(b, rs); _t(b); _putn(b, re); _t(b);
		_putn(b, dcnt - a->xcnt); _t(b); _putn(b, dcnt + a->gecnt); _t(b); _putn(b, mapq>>MAPQ_DEC);

		const uint64_t f = b->opt->flag;
		if (f & 0x01ULL<<MM_AS) { _putsk(b, "\tAS:i:"); _putn(b, a->score); }
		if (f & 0x01ULL<<MM_NM) { _putsk(b, "\tNM:i:"); _putn(b, a->xcnt + a->gecnt); }
		_cr(b);
	}
	return;
}

// qname rname 1-idt score qd qs qe ql rd rs rd rl
static void mm_print_mapped_mhap(mm_align_t *b, const bseq_t *t, uint64_t n_reg, const mm128_t *reg)
{
	const mm_idx_t *mi = b->mi;
	for (uint64_t j = 0; j < (uint32_t)n_reg; ++j) {
		const gaba_alignment_t *a = _aln(reg[j]);
		const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];
		const gaba_path_section_t *s = &a->sec[0];
		int32_t dcnt = (a->path->len - a->gecnt)>>1, slen = dcnt + a->gecnt;
		int32_t mid = 1000000.0 * (double)(dcnt - a->xcnt) / (double)slen;	// percent identity
		uint32_t rs = s->bid&0x01? r->l_seq-s->apos-s->alen : s->apos, re = rs+s->alen;
		uint32_t qs = s->bid&0x01? t->l_seq-s->bpos-s->blen : s->bpos, qe = qs+s->blen;

		_putsn(b, t->name, t->l_name); _sp(b); _putsn(b, r->name, r->l_name); _sp(b);
		_putfi(int32_t, b, 1.0-mid, 4); _sp(b); _putn(b, a->score); _sp(b);
		_put(b, '0'); _sp(b); _putn(b, qs); _sp(b); _putn(b, qe); _sp(b); _putn(b, t->l_seq); _sp(b);
		_put(b, s->bid&0x01? '0' : '1'); _sp(b); _putn(b, rs); _sp(b); _putn(b, re); _sp(b); _putn(b, r->l_seq); _cr(b);
	}
	return;
}

static void mm_print_mapped_falcon(mm_align_t *b, const bseq_t *t, uint64_t n_reg, const mm128_t *reg)
{
	// print header line
	_putsn(b, t->name, t->l_name); _sp(b); _putsnt(b, t->seq, t->l_seq, "NACMGRSVTWYHKDBN"); _cr(b);
	// print alignment lines
	const mm_idx_t *mi = b->mi;
	for (uint64_t j = 0; j < (uint32_t)n_reg; ++j) {
		const gaba_alignment_t *a = _aln(reg[j]);
		const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];
		const gaba_path_section_t *s = &a->sec[0];
		uint32_t qs = r->l_seq-s->apos-s->alen, qe = r->l_seq-s->apos;
		_putsn(b, r->name, r->l_name); _sp(b);
		if (s->bid&0x01) { _putsnt(b, &r->seq[qs], qe-qs, "NACMGRSVTWYHKDBN"); }
		else { _putsntr(b, &r->seq[r->l_seq-qe], qe-qs, "NTGKCYSBAWRDMHVN"); }
		_cr(b);
	}
	_put(b, '+'); _sp(b); _put(b, '+'); _cr(b);
	return;
}
static void mm_print_footer_falcon(mm_align_t *b)
{
	_put(b, '-'); _sp(b); _put(b, '-'); _cr(b);
	return;
}

#undef _d
#undef _t
#undef _c
#undef _cr
#undef _put
#undef _putn
#undef _puts
#undef _putsn
#undef _putsk
#undef _aln

static void *mm_align_source(uint32_t tid, void *arg)
{
	set_info(tid, "[mm_align_source] fetch sequence block");
	mm_align_t *b = (mm_align_t*)arg;
	mm_align_step_t *s = (mm_align_step_t*)calloc(1, sizeof(mm_align_step_t));

	s->lmm = lmm_init(NULL, 512 * 1024);
	if (s->lmm == NULL) return NULL;
	s->seq = bseq_read(b->fp, &s->n_seq, &s->base, &s->size);
	if (s->seq == 0) { lmm_clean(s->lmm), free(s), s = NULL; return NULL; }
	if (b->base_qid == UINT32_MAX) b->base_qid = atoi(s->seq[0].name);
	s->id = b->icnt++;
	s->base_qid = b->base_qid; b->base_qid += s->n_seq;
	return s;
}

static void *mm_align_worker(uint32_t tid, void *arg, void *item)
{
	mm_tbuf_t *t = (mm_tbuf_t*)arg;
	mm_align_step_t *s = (mm_align_step_t*)item;

	char buf[128], *p = buf;
	p += _pstr(p, "[mm_align_worker] bin id "); p += _pnum(uint32_t, p, s->base_qid); p += _pstr(p, ":"); p += _pnum(uint32_t, p, s->base_qid + s->n_seq - 1); *p = '\0';
	set_info(tid, buf);

	for (uint64_t i = 0; i < s->n_seq; ++i) {
		mm128_t *r;
		kv_pushp(mm128_t, s->reg, &r);
		uint32_t qid = s->base_qid + i;			// FIXME: parse qid from name with atoi when -M is set
		r->u64[0] = (uintptr_t)mm_align_seq(t, t->b->opt, t->b->mi, s->seq[i].l_seq, s->seq[i].seq, qid, t->b->n_occ, t->b->occ, s->lmm, &r->u64[1]);
	}
	return s;
}

static void mm_align_drain_intl(mm_align_t *b, mm_align_step_t *s)
{
	for (uint64_t i = 0; i < s->n_seq; ++i) {
		mm128_t *r = (mm128_t*)s->reg.a[i].u64[0];
		uint64_t n_reg = s->reg.a[i].u64[1];	// n_reg in lower 32bit, n_uniq in higher 32bit
		if (r == 0) {
			if (b->printer.unmapped) b->printer.unmapped(b, &s->seq[i]);
			continue;
		}
		b->printer.mapped(b, &s->seq[i], n_reg, r);
		for (uint64_t j = 0; j < (uint32_t)n_reg; ++j) { lmm_free(s->lmm, (void*)r[j].u64[1]); }
		free(r);
	}
	free(s->reg.a); free(s->base); free((void*)s->seq); lmm_clean(s->lmm); free(s);
	return;
}

static void mm_align_drain(uint32_t tid, void *arg, void *item)
{
	set_info(tid, "[mm_align_drain] dump alignments to sam records");
	mm_align_t *b = (mm_align_t*)arg;
	mm_align_step_t *s = (mm_align_step_t*)item;

	#if STRICT_STREAM_ORDERING != 0
		kv_hq_push(mm128_t, incq_comp, b->hq, ((mm128_t){.u64 = {s->id, (uintptr_t)s}}));
		while (b->hq.n > 1 && b->hq.a[1].u64[0] == b->ocnt) {
			s = (mm_align_step_t*)kv_hq_pop(mm128_t, incq_comp, b->hq).u64[1]; b->ocnt++;
			mm_align_drain_intl(b, s);
		}
	#else
		mm_align_drain_intl(b, s);
	#endif
	return;
}

static void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == NULL) return;
	free(b->resc.a); free(b->coef.a); free(b->intv.a);
	kh_destroy(b->pos); gaba_dp_clean(b->dp);
	free(b);
}

static mm_tbuf_t *mm_tbuf_init(mm_align_t *b)
{
	mm_tbuf_t *t = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	if (t == 0) return NULL;
	const uint8_t *lim = (const uint8_t *)0x800000000000;
	t->b = b;
	if ((t->pos = kh_init(0)) == 0) goto _fail;
	if ((t->dp = gaba_dp_init(b->gaba, lim, lim)) == 0) goto _fail;
	return t;
_fail:
	mm_tbuf_destroy(t);
	return NULL;
}

static void mm_align_destroy(mm_align_t *b)
{
	if (b->printer.footer) b->printer.footer(b);
	fwrite(b->base, sizeof(uint8_t), b->p - b->base, stdout);
	for (uint64_t i = 0; i < b->opt->nth; ++i) mm_tbuf_destroy((mm_tbuf_t*)b->t[i]);
	free(b->base); kv_hq_destroy(b->hq); free(b->occ); free(b->t); gaba_clean(b->gaba); pt_destroy(b->pt);
	free(b);
	return;
}

static mm_align_t *mm_align_init(const mm_mapopt_t *opt, const mm_idx_t *mi)
{
	static const mm_printer_t printer[] = {
		[0] = { .header = mm_print_header_sam, .unmapped = mm_print_unmapped_sam, .mapped = mm_print_mapped_sam },
		[MM_MAF>>56] = { .mapped = mm_print_mapped_maf },
		[MM_BLAST6>>56] = { .mapped = mm_print_mapped_blast6 },
		[MM_BLASR1>>56] = { .mapped = mm_print_mapped_blasr1 },
		[MM_BLASR4>>56] = { .mapped = mm_print_mapped_blasr4 },
		[MM_PAF>>56] = { .mapped = mm_print_mapped_paf },
		[MM_MHAP>>56] = { .mapped = mm_print_mapped_mhap },
		[MM_FALCON>>56] = { .mapped = mm_print_mapped_falcon, .footer = mm_print_footer_falcon }
	};

	mm_align_t *b = calloc(1, sizeof(mm_align_t));
	if (b == NULL) return NULL;

	// init output queue, buf and printer
	b->tail = (b->base = b->p = malloc(sizeof(uint8_t) * opt->outbuf_size)) + opt->outbuf_size;
	if (b->base == 0) goto _fail;
	b->icnt = b->ocnt = 0; b->base_qid = opt->base_qid;
	kv_hq_init(b->hq);
	b->printer = printer[opt->flag>>56];
	if (!b->printer.mapped) goto _fail;
	for (uint64_t i = 0; i < 9; ++i) b->conv[i] = (i + 1) % 10;
	for (uint64_t i = 9; i < 40; ++i) b->conv[i] = (((i + 1) % 10)<<4) + (i + 1) / 10;

	// set consts
	b->mi = mi, b->opt = opt;

	// calc occ
	if ((b->occ = calloc(b->opt->n_frq, sizeof(uint64_t))) == 0) goto _fail;
	b->n_occ = b->opt->n_frq;
	for (uint64_t i = 0; i < b->n_occ; ++i) b->occ[i] = mm_idx_cal_max_occ(mi, b->opt->frq[i]);

	// initialize alignment context
	struct gaba_params_s p = {
		.m = opt->m, .x = opt->x, .gi = opt->gi, .ge = opt->ge,
		.xdrop = opt->xdrop, .filter_thresh = 7
	};
	if ((b->gaba = gaba_init(&p)) == 0) goto _fail;

	// initialize threads
	if ((b->t = (void**)calloc(b->opt->nth, sizeof(mm_tbuf_t*))) == 0) goto _fail;
	for (uint64_t i = 0; i < b->opt->nth; ++i) { if ((b->t[i] = (void*)mm_tbuf_init(b)) == 0) goto _fail; }
	if ((b->pt = pt_init(b->opt->nth)) == 0) goto _fail;

	// print sam header
	if (b->printer.header) b->printer.header(b);
	return b;
_fail:
	mm_align_destroy(b);
	return NULL;
}

static int mm_align_file(mm_align_t *b, bseq_file_t *fp)
{
	if (fp == 0) return -1;
	b->fp = fp;
	return pt_stream(b->pt, mm_align_source, b, mm_align_worker, (void**)b->t, mm_align_drain, b);
}

/* end of map.c */

/* main.c */

static void liftrlimit()
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}

static void posixly_correct()
{
#ifdef __linux__
	setenv("POSIXLY_CORRECT", "1", 0);
#endif
}

static const char *mm_get_version(void)
{
	const char *prefix = "minialign-";
	uint64_t spos = strncmp(MM_VERSION, prefix, MIN2(strlen(MM_VERSION), strlen(prefix)))? 0 : strlen(prefix);
	return &MM_VERSION[spos];
}

static void mm_print_version(void)
{
	puts(mm_get_version());
	return;
}

static void mm_print_help(const mm_mapopt_t *opt)
{
	if (opt->verbose == 0) return;
	fprintf(stderr, "\n"
					"  minialign - fast aligner for PacBio and Nanopore long reads\n"
					"\n"
					"minialign is a fast long-read (nucleotide sequence) alignment tool built on\n"
					"the top of minimap long-read overlapper adopting libgaba SIMD-parallelized\n"
					"Smith-Waterman extension algorithm.\n"
					"\n");
	fprintf(stderr, "Usage:\n"
					"  first trial:\n"
					"    $ minialign -t4 -xont <ref.fa> <ont2d.{fa,fq,bam}> > mapping.sam\n"
					"\n"
					"  mapping on a prebuilt index (saves ~1min for human genome per run):\n"
					"    $ minialign [indexing options] -d <index.mai> <ref.fa>\n"
					"    $ minialign -l <index.mai> <reads.{fa,fq,bam}> > mapping.sam\n"
					"\n"
					"  all-versus-all alignment in a read set:\n"
					"    $ minialign -X -xava <reads.fa> [<reads.fa> ...] > allvsall.paf\n"
					"\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  Global:\n");
	fprintf(stderr, "    -x STR       load preset params {pacbio,ont,ava} [ont]\n");
	fprintf(stderr, "    -t INT       number of threads [%d]\n", opt->nth);
	fprintf(stderr, "    -X           switch to all-versus-all alignment mode\n");
	fprintf(stderr, "    -v           show version number [%s]\n", mm_get_version());
	fprintf(stderr, "  Indexing:\n");
	// fprintf(stderr, "    -c STR,...   treat specified sequences as circular []\n");
	fprintf(stderr, "    -k INT       k-mer size [%d]\n", opt->k);
	fprintf(stderr, "    -w INT       minimizer window size [{-k}*2/3]\n");
	fprintf(stderr, "    -d FILE      dump index to FILE []\n");
	fprintf(stderr, "    -l FILE      load index from FILE [] (overriding -k and -w)\n");
	if (opt->verbose >= 2) {
		fprintf(stderr, "    -C INT[,INT] set base rid and qid [%u, %u]\n", opt->base_rid, opt->base_qid);
		fprintf(stderr, "    -N           treat sequence name as id (seq must be sorted)\n");
		fprintf(stderr, "    -L INT       sequence length filter: min. ref. length; 0 to disable [%u]\n", opt->rmin);
		fprintf(stderr, "    -H INT       sequence length filter: min. query length; 0 to disable [%u]\n", opt->qmin);
	}
	fprintf(stderr, "  Mapping:\n");
	if (opt->verbose >= 2)
		fprintf(stderr, "    -f FLOAT,... occurrence thresholds [0.05,0.01,0.001]\n");
	fprintf(stderr, "    -a INT       match award [%d]\n", opt->m);
	fprintf(stderr, "    -b INT       mismatch penalty [%d]\n", opt->x);
	fprintf(stderr, "    -p INT       gap open penalty [%d]\n", opt->gi);
	fprintf(stderr, "    -q INT       gap extension penalty [%d]\n", opt->ge);
	if (opt->verbose >= 2)
		fprintf(stderr, "    -Y INT       X-drop threshold [%d]\n", opt->xdrop);
	fprintf(stderr, "    -s INT       minimum alignment score [%d]\n", opt->min);
	fprintf(stderr, "    -m INT       minimum alignment score ratio [%1.2f]\n", opt->min_ratio);
	if (opt->verbose >= 2) {
		fprintf(stderr, "    -M INT       maximum #alignments reported [%d]\n", opt->max_cnt);
		fprintf(stderr, "    -A           calculate both topright and bottomright triangles (only effective with -X)\n");
	}
	fprintf(stderr, "  Output:\n");
	fprintf(stderr, "    -O STR       output format {sam,maf,blast6,blasr1,blasr4,paf,mhap,falcon} [%s]\n",
		(const char *[]){ "sam", "maf", "blast6", "blasr1", "blasr4", "paf", "mhap", "falcon" }[opt->flag>>56]);
	if (opt->verbose >= 2)
		fprintf(stderr, "    -P           omit secondary (repetitive) alignments\n");
	fprintf(stderr, "    -Q           include quality string\n");
	fprintf(stderr, "    -R STR       read group header line, like \"@RG\\tID:1\" [%s]\n", opt->rg_line? opt->rg_line : "");
	fprintf(stderr, "    -T STR,...   list of optional tags: {RG,AS,XS,NM,NH,IH,SA,MD} []\n");
	fprintf(stderr, "                   RG is also inferred from -R\n");
	fprintf(stderr, "                   supp. records are omitted when SA tag is enabled\n");
	if (opt->verbose >= 2)
		fprintf(stderr, "    -U STR,...   tags to be transferred from the input bam file []\n");
	fprintf(stderr, "\n");
	if (opt->verbose < 2) {
		fprintf(stderr, "  Pass -hVV to show all the options.\n");
		fprintf(stderr, "\n");
	}
	return;
}

static int mm_mapopt_load_preset(mm_mapopt_t *o, const char *arg)
{
	if (strcmp(arg, "pacbio") == 0) {
		o->k = 15; o->w = 10; o->m = 1; o->x = 2; o->gi = 2; o->ge = 1; o->xdrop = 50; o->min = 50; o->min_ratio = 0.3;
	} else if (strcmp(arg, "ont") == 0) {
		o->k = 15; o->w = 10; o->m = 1; o->x = 1; o->gi = 1; o->ge = 1; o->xdrop = 50; o->min = 50; o->min_ratio = 0.3;
	} else if (strcmp(arg, "ava") == 0) {
		o->k = 14; o->w = 5; o->m = 1; o->x = 2; o->gi = 0; o->ge = 1; o->xdrop = 50; o->min = 30; o->min_ratio = 0.05;
		o->flag |= MM_AVA | MM_PAF;
	} else {
		if (o->verbose >= 1) fprintf(stderr, "[M::%s] Warning: Unknown preset tag: `%s'.\n", __func__, arg);
		return 1;
	}
	return 0;
}

static int mm_mapopt_parse_threshs(mm_mapopt_t *o, const char *arg)
{
	const char *p = optarg; o->n_frq = 0;
	while (*p) {
		o->frq[o->n_frq++] = atof(p);
		while (*p && *p != ',') p++;
		if (!*p || o->n_frq >= 16) break;
		p++;
	}
	for (uint64_t i = 0; i < o->n_frq; ++i) {
		if (o->frq[i] < 0.0 || o->frq[i] > 1.0)
			fprintf(stderr, "[M::%s] Warning: Invalid threshold `%f' parsed from `%s'.\n", __func__, o->frq[i], arg);
	}
	return 0;
}

static int mm_mapopt_parse_base_ids(mm_mapopt_t *o, const char *arg)
{
	const char *p = optarg;
	if (!isdigit(*p)) return 1;
	o->base_rid = atoi(p);
	while (*p && *p != ',') { p++; }
	if (*p && isdigit(p[1])) o->base_qid = atoi(&p[1]);
	return 0;
}

static uint64_t mm_mapopt_parse_tags(mm_mapopt_t *o, const char *p, uint16_v *buf)
{
	uint32_t flag = 0;
	while (*p != '\0') {
		const char *q = p;
		while (*q != ',' && *q != '\0') q++;
		if (q - p == 2) {
			if (buf) kv_push(uint16_t, *buf, (uint16_t)p[0] | ((uint16_t)p[1]<<8));
			if (strncmp(p, "RG", 2) == 0) flag |= 0x01ULL<<MM_RG;
			else if (strncmp(p, "CO", 2) == 0) flag |= 0x01ULL<<MM_CO;
			else if (strncmp(p, "NH", 2) == 0) flag |= 0x01ULL<<MM_NH;
			else if (strncmp(p, "IH", 2) == 0) flag |= 0x01ULL<<MM_IH;
			else if (strncmp(p, "AS", 2) == 0) flag |= 0x01ULL<<MM_AS;
			else if (strncmp(p, "XS", 2) == 0) flag |= 0x01ULL<<MM_XS;
			else if (strncmp(p, "NM", 2) == 0) flag |= 0x01ULL<<MM_NM;
			else if (strncmp(p, "SA", 2) == 0) flag |= 0x01ULL<<MM_SA;
			else if (strncmp(p, "MD", 2) == 0) flag |= 0x01ULL<<MM_MD;
			else if (o->verbose >= 1) fprintf(stderr, "[M::%s] Unknown tag: `%c%c'.\n", __func__, p[0], p[1]);
		}
		if (*q == '\0') break;
		p = q + 1;
	}
	return flag;
}

static void mm_mapopt_parse_rg(mm_mapopt_t *o, const char *arg)
{
	free(o->rg_line); free(o->rg_id);
	o->rg_line = NULL; o->rg_id = NULL;
	o->flag &= ~(0x01ULL<<MM_RG);

	const char *id, *p = arg;
	char b[256], *q = b;
	while (q < b+256 && *p != '\0') {
		if (p[0] == '\\' && p[1] == 't') { if (p[1] == '\0') { break; } *q++ = '\t'; p += 2; }
		else { *q++ = *p++; }
	}
	*q = '\0';
	if (strstr(b, "@RG") != b) return;
	if ((id = strstr(b, "\tID:")) == 0) return;

	#define _dup(x)	({ char *_p = malloc(strlen(x)+1); memcpy(_p, (x), strlen(x)); _p[strlen(x)] = '\0'; _p; })
	o->rg_line = _dup(b); o->rg_id = _dup(id+4);
	char *t = o->rg_id;
	while (*t != '\t' && *t != '\n' && *t != '\0') { t++; }
	*t = '\0';
	o->flag |= 0x01ULL<<MM_RG;
	#undef _dup
	return;
}

static uint64_t mm_mapopt_parse_format(mm_mapopt_t *o, const char *arg)
{
	if (strcmp(arg, "sam") == 0) return 0;
	else if (strcmp(arg, "maf") == 0) return MM_MAF;
	else if (strcmp(arg, "blast6") == 0) return MM_BLAST6;
	else if (strcmp(arg, "blasr1") == 0) return MM_BLASR1;
	else if (strcmp(arg, "blasr4") == 0) return MM_BLASR4;
	else if (strcmp(arg, "paf") == 0) return MM_PAF;
	else if (strcmp(arg, "mhap") == 0) return MM_MHAP;
	else if (strcmp(arg, "falcon") == 0) return MM_FALCON;
	else if (o->verbose >= 1) fprintf(stderr, "[M::%s] Unknown output format: `%s'.\n", __func__, arg);
	return 0;
}

static int mm_mapopt_parse(mm_mapopt_t *o, int argc, char *argv[], const char **fnr, const char **fnw, ptr_v *v)
{
	int ret = 0;
	while (optind < argc) {
		int ch;
		if ((ch = getopt(argc, argv, "t:x:V:c:k:w:f:B:d:l:C:NXAs:m:r:M:a:b:p:q:L:H:I:J:S:E:D:Y:O:PQR:T:U:vh")) < 0) {
			kv_push(void*, *v, argv[optind]); optind++; continue;
		}
		switch (ch) {
			case 't': o->nth = atoi(optarg); break;
			case 'x': mm_mapopt_load_preset(o, optarg); break;
			case 'V':
				if (isdigit(optarg[0])) { o->verbose = atoi(optarg); }
				else { for (const char *p = optarg; *p; p++, o->verbose++) {} }
				break;
			case 'k': o->k = atoi(optarg); break;
			case 'w': o->w = atoi(optarg); break;
			case 'f': mm_mapopt_parse_threshs(o, optarg); break;
			case 'B': o->b = atoi(optarg); break;
			case 'd': *fnw = optarg; break;
			case 'l': *fnr = optarg; break;
			case 'C': mm_mapopt_parse_base_ids(o, optarg); break;
			case 'N': o->base_rid = UINT32_MAX, o->base_qid = UINT32_MAX; break;
			case 'X': o->flag |= MM_AVA; break;
			case 'A': o->flag |= MM_COMP; break;
			case 's': o->min = atoi(optarg); break;
			case 'm': o->min_ratio = atof(optarg); break;
			case 'M': o->max_cnt = atoi(optarg); break;
			case 'a': o->m = atoi(optarg); break;
			case 'b': o->x = atoi(optarg); break;
			case 'p': o->gi = atoi(optarg); break;
			case 'q': o->ge = atoi(optarg); break;
			case 'F': o->llim = atoi(optarg); break;
			case 'G': o->hlim = atoi(optarg); break;
			case 'L': o->rmin = atoi(optarg); break;
			case 'H': o->qmin = atoi(optarg); break;
			case 'I': o->blim = atoi(optarg); break;
			case 'J': o->elim = atoi(optarg); break;
			case 'S': o->sidx = atoi(optarg); break;
			case 'E': o->eidx = atoi(optarg); break;
			case 'D': o->batch_size = atoi(optarg); break;
			case 'Y': o->xdrop = atoi(optarg); break;
			case 'O': o->flag &= ~(0xffULL<<56), o->flag |= mm_mapopt_parse_format(o, optarg); break;
			case 'P': o->flag |= MM_OMIT_REP; break;
			case 'Q': o->flag |= MM_KEEP_QUAL; break;
			case 'R': mm_mapopt_parse_rg(o, optarg); break;
			case 'T': o->flag |= mm_mapopt_parse_tags(o, optarg, NULL); break;
			case 'U': mm_mapopt_parse_tags(o, optarg, &o->tags); break;
			case 'v': ret = 1; break;
			case 'h': ret = 2; break;
		}
	}
	return ret;
}

int main(int argc, char *argv[])
{
	int ret = 1;
	const char *fnr = NULL, *fnw = NULL;
	bseq_file_t *fp = NULL;
	ptr_v v = { 0 };
	FILE *fpr = NULL, *fpw = NULL;

	// unittest hook
	#if UNITTEST != 0
	if (argc > 1 && strcmp(argv[1], "unittest") == 0) return unittest_main(argc, argv);
	#endif

	enable_info(0);
	set_info(0, "[main] parsing arguments");
	liftrlimit(); posixly_correct();
	mm_realtime0 = realtime();
	mm_mapopt_t *opt = mm_mapopt_init();
	switch (mm_mapopt_parse(opt, argc, argv, &fnr, &fnw, &v)) {
		case 1: mm_print_version(); ret = 0; goto _final;
		case 2: mm_print_help(opt); ret = 0; goto _final;
	}
	if (!fnr && v.n == 0) { mm_print_help(opt); goto _final; }
	if (!fnw && ((fnr && v.n == 0) || (!fnr && v.n == 1 && !(opt->flag&MM_AVA)))) {
		if (opt->verbose >= 1) fprintf(stderr, "[M::%s] query-side input redirected to stdin.\n", __func__);
		kv_push(void*, v, "-");
	}
	if (opt->w >= 16) opt->w = (int)(.6666667 * opt->k + .499);
	if (mm_mapopt_check(opt, fprintf, stderr)) goto _final;

	set_info(0, "[main] open index file");
	if (fnr && (fpr = fopen(fnr, "rb")) == NULL) {
		fprintf(stderr, "[E::%s] failed to open index file `%s'. Please check file path and it exists.\n", __func__, fnr);
		goto _final;
	}
	if (fnw && (fpw = fopen(fnw, "wb")) == NULL) {
		fprintf(stderr, "[E::%s] failed to open index file `%s' in write mode. Please check file path and its permission.\n", __func__, fnw);
		goto _final;
	}
	for (uint64_t i = 0; i < (fpr? UINT64_MAX : (((opt->flag&MM_AVA) || fpw)? v.n : 1)); ++i) {
		mm_idx_t *mi = NULL;
		mm_align_t *aln = NULL;

		// load/generate index for this index-side iteration
		if (fpr) {
			mi = mm_idx_load(fpr, opt->nth);
		} else {
			fp = bseq_open((const char *)v.a[i], opt->batch_size, 0, opt->rmin, 0, NULL);
			mi = mm_idx_gen(opt, fp);
			opt->base_rid += bseq_close(fp);
		}

		// check sanity of the index
		if (mi == 0) {
			if (fpr && i > 0) break;
			fprintf(stderr, "[E::%s] failed to %s `%s'. Please check %s.\n", __func__,
				fpr? "load index file" : "open sequence file",
				fpr? fnr : (const char*)v.a[i], fpr? "file path, format and its version" : "file path and format");
			goto _final;
		}
		if (opt->verbose >= 1) fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built index for %lu target sequence(s)\n", __func__,
			realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->s.n);

		// do the task
		if (fpw) {
			mm_idx_dump(fpw, mi, opt->nth);
		} else {
			aln = mm_align_init(opt, mi);
			for (uint64_t j = (!fpr && !(opt->flag&MM_AVA)); j < v.n; ++j) {
				if (!(fp = bseq_open((const char*)v.a[j], opt->batch_size, (opt->flag&MM_KEEP_QUAL)!=0, opt->qmin, opt->tags.n, opt->tags.a))) {
					fprintf(stderr, "[E::%s] failed to open sequence file `%s'. Please check file path and format.\n", __func__, (const char*)v.a[j]);
					goto _final;
				}
				mm_align_file(aln, fp);
				bseq_close(fp);
			}
			mm_align_destroy(aln);
		}
		mm_idx_destroy(mi);
	}
	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (uint64_t i = 0; i < (uint64_t)argc; ++i) fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
	ret = 0;
_final:
	free(v.a); mm_mapopt_destroy(opt);
	if (fpr) fclose(fpr);
	if (fpw) fclose(fpw);
	return ret;
}

/* end of main.c */
