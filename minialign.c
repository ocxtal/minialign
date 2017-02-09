#include <getopt.h>
#include <unistd.h>
#include <zlib.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/types.h>


#define MM_VERSION "0.4.5-rc1"

#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )

int mm_verbose = 3;
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

_Thread_local const char *info;	// thread-local comment on the current tasks
#define set_info(x)		( info = (const char *)(x) )
static void oom_abort(const char *name)
{
	fprintf(stderr, "[M::%s] ERROR: Out of memory. (%s)\n", name, info == NULL? "No additional information available" : info);
	exit(128);	// 128 reserved for out of memory
}

#define mm_malloc(x) ({ \
	void *_ptr = malloc((size_t)(x)); \
	if (_ptr == NULL) { oom_abort(__func__); } \
	_ptr; \
})
#define malloc(x)		mm_malloc(x)

#define mm_realloc(x, y) ({ \
	void *_ptr = realloc((void *)(x), (size_t)(y)); \
	if (_ptr == NULL) { oom_abort(__func__); } \
	_ptr; \
})
#define realloc(x, y)	mm_realloc(x, y)

#define mm_calloc(x, y) ({ \
	void *_ptr = calloc((size_t)(x), (size_t)(y)); \
	if (_ptr == NULL) { oom_abort(__func__); } \
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

/* end of misc.c */

#include "kvec.h"
#include "gaba.h"
#include "lmm.h"
#include "sassert.h"
#include "arch/arch.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

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
#define KH_SIZE			( 16 )
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
	size |= size>>1; size |= size>>2; size |= size>>4;
	size |= size>>8; size |= size>>16; size |= size>>32;

	h->mask = KH_SIZE - 1; h->max = KH_SIZE; h->cnt = 0; h->ub = KH_SIZE * KH_THRESH;
	h->a = malloc(sizeof(mm128_t) * KH_SIZE);
	for (uint64_t i = 0; i < KH_SIZE; ++i) h->a[i].u64[0] = (int64_t)-2, h->a[i].u64[1] = 0;
	return h;
}

static void kh_destroy(kh_t *h)
{
	if (h == NULL) { return; }
	free(h->a); free(h);
	return;
}

static void kh_dump(gzFile fp, kh_t *h)
{
	uint32_t x[2] = {0};
	if (h == 0) { gzwrite(fp, x, sizeof(uint32_t) * 2); return; }
	x[0] = h->mask + 1; x[1] = h->cnt;
	gzwrite(fp, x, sizeof(uint32_t) * 2);
	gzwrite(fp, h->a, sizeof(mm128_t) * x[0]);
	return;
}

static kh_t *kh_load(gzFile fp)
{
	uint32_t x[2] = {0};
	if ((gzread(fp, x, sizeof(uint32_t) * 2)) != sizeof(uint32_t) * 2 || x[0] == 0) return NULL;
	kh_t *h = calloc(1, sizeof(kh_t));
	h->mask = x[0] - 1; h->max = x[0]; h->cnt = x[1]; h->ub = h->max * KH_THRESH;
	h->a = malloc(sizeof(mm128_t) * x[0]);
	if ((gzread(fp, h->a, sizeof(mm128_t) * x[0])) != sizeof(mm128_t) * x[0]) { free(h->a); free(h); return NULL; }
	return h;
}

static void kh_clear(kh_t *h)
{
	if (h == 0) { return; }
	h->mask = KH_SIZE - 1; h->cnt = 0; h->ub = KH_SIZE * KH_THRESH;
	for (uint64_t i = 0; i < KH_SIZE; ++i) h->a[i].u64[0] = (int64_t)-2, h->a[i].u64[1] = 0;
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
	for (uint64_t i = 0; i < prev_size; ++i) h->a[i + prev_size].u64[0] = (int64_t)-2, h->a[i + prev_size].u64[1] = 0;
	for (uint64_t i = 0; i < size; ++i) {
		uint64_t k = h->a[i].u64[0];	// key
		if (k + 2 < 2 || (k & mask) == i) continue;
		uint64_t v = h->a[i].u64[1]; h->a[i].u64[0] = (int64_t)-1;	// moved
		kh_put_intl(h->a, k, v, mask);
	}
	return;
}

static void kh_put(kh_t *h, uint64_t key, uint64_t val)
{
	uint64_t mask = h->mask;
	if (h->cnt >= h->ub) kh_extend(h);
	h->cnt += kh_put_intl(h->a, key, val, mask);
	return;
}

static uint64_t kh_get(kh_t *h, uint64_t key)
{
	uint64_t mask = h->mask, pos = key & mask, k;
	do {
		if ((k = h->a[pos].u64[0]) == key) return h->a[pos].u64[1];
		pos = mask & (pos + 1);
	} while (k + 1 < (int64_t)-1);	// !is_empty(k) || is_moved(k)
	return (int64_t)-1;
}

static uint64_t *kh_get_ptr(kh_t *h, uint64_t key)
{
	uint64_t mask = h->mask, pos = key & mask, k;
	do {
		if ((k = h->a[pos].u64[0]) == key) return &h->a[pos].u64[1];
		pos = mask & (pos + 1);
	} while (k + 1 < (int64_t)-1);	// !is_empty(k) || is_moved(k)
	return NULL;
}
/* end of hash.c */

/* queue.c */
#define _cas(ptr, cmp, val) __atomic_compare_exchange_n(ptr, cmp, val, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)
typedef void *(*pt_source_t)(void *arg);
typedef void *(*pt_worker_t)(void *arg, void *item);
typedef void (*pt_drain_t)(void *arg, void *item);

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
	pt_worker_t wfp;
	void *warg;
} pt_thread_t;

typedef struct pt_s {
	pt_q_t in, out;
	uint32_t nth;
	pt_thread_t c[];	// [0] is reserved for master
} pt_t;

#define PT_EMPTY	( (void *)((int64_t)-1) )
#define PT_EXIT		( (void *)((int64_t)-2) )

uint64_t pt_enq(pt_q_t *q, uint64_t tid, void *elem)
{
	uint64_t z, ret = (int64_t)-1;
	do { z = 0xffffffff; } while (!_cas(&q->lock, &z, tid));
	uint64_t head = q->head, tail = q->tail, mask = q->size - 1;
	if (((head + 1) & mask) != tail) {
		q->elems[head] = elem; q->head = (head + 1) & mask; ret = 0;
	}
	do { z = tid; } while (!_cas(&q->lock, &z, 0xffffffff));
	return ret;
}

void *pt_deq(pt_q_t *q, uint64_t tid)
{
	void *elem = PT_EMPTY;
	uint64_t z;
	do { z = 0xffffffff; } while (!_cas(&q->lock, &z, tid));
	uint64_t head = q->head, tail = q->tail, mask = q->size - 1;
	if (head != tail) {
		elem = q->elems[tail]; q->tail = (tail + 1) & mask;
	}
	do { z = tid; } while (!_cas(&q->lock, &z, 0xffffffff));
	return elem;
}

static void *pt_dispatch(void *s)
{
	pt_thread_t *c = (pt_thread_t *)s;
	void *ping = PT_EMPTY, *pong = PT_EMPTY;
	while(1) {
		ping = pt_deq(c->in, c->tid);
		if (ping == PT_EMPTY && pong == PT_EMPTY) sched_yield();
		if (pong != PT_EMPTY) pt_enq(c->out, c->tid, c->wfp(c->warg, pong));
		if (ping == PT_EXIT) break;
		pong = pt_deq(c->in, c->tid);
		if (ping == PT_EMPTY && pong == PT_EMPTY) sched_yield();
		if (ping != PT_EMPTY) pt_enq(c->out, c->tid, c->wfp(c->warg, ping));
		if (pong == PT_EXIT) break;
	}
	return NULL;
}

static void pt_destoroy(pt_t *pt)
{
	void *status;
	for (uint64_t i = 1; i < pt->nth; ++i) pt_enq(pt->c->in, pt->c->tid, PT_EXIT);
	for (uint64_t i = 1; i < pt->nth; ++i) pthread_join(pt->c[i].th, &status);
	while (pt_deq(pt->c->in, pt->c->tid) != PT_EMPTY) {}
	while (pt_deq(pt->c->out, pt->c->tid) != PT_EMPTY) {}
	free(pt->in.elems); free(pt->out.elems); free(pt);
}

static pt_t *pt_init(uint32_t nth)
{
	nth = (nth == 0)? 1 : nth;
	pt_t *pt = calloc(1, sizeof(pt_t) + sizeof(pt_thread_t) * nth);
	pt->in = (pt_q_t){ .lock = 0xffffffff, .elems = calloc(8 * nth, sizeof(void*)), .size = 8 * nth };
	pt->out = (pt_q_t){ .lock = 0xffffffff, .elems = calloc(8 * nth, sizeof(void*)), .size = 8 * nth };

	pt->nth = nth; pt->c[0].tid = 0; pt->c[0].in = &pt->in; pt->c[0].out = &pt->out;
	for (uint64_t i = 1; i < nth; ++i) {
		pt->c[i].tid = i; pt->c[i].in = &pt->in; pt->c[i].out = &pt->out;
		pthread_create(&pt->c[i].th, NULL, pt_dispatch, (void *)&pt->c[i]);
	}
	return pt;
}

static int pt_set_worker(pt_t *pt, pt_worker_t wfp, void **warg)
{
	void *item;
	if ((item = pt_deq(pt->c->in, pt->c->tid)) != PT_EMPTY) { pt_enq(pt->c->in, pt->c->tid, item); return -1; }
	for (uint64_t i = 0; i < pt->nth; ++i) pt->c[i].wfp = wfp, pt->c[i].warg = warg[i];
	return 0;
}

static int pt_stream(pt_t *pt, pt_source_t sfp, void *sarg, pt_worker_t wfp, void **warg, pt_drain_t dfp, void *darg)
{
	if (pt_set_worker(pt, wfp, warg)) return -1;
	uint64_t bal = 0, lb = 2 * pt->nth, ub = 4 * pt->nth;
	void *item;
	while ((item = sfp(sarg)) != NULL) {
		pt_enq(pt->c->in, pt->c->tid, item);
		if (++bal < ub) continue;
		while (bal > lb) {
			if ((item = pt_deq(pt->c->out, pt->c->tid)) != PT_EMPTY) dfp(darg, item), bal--;
			if ((item = pt_deq(pt->c->in, pt->c->tid)) != PT_EMPTY) dfp(darg, wfp(*warg, item)), bal--;
		}
	}
	while ((item = pt_deq(pt->c->in, pt->c->tid)) != PT_EMPTY) pt_enq(pt->c->out, pt->c->tid, wfp(*warg, item));
	while (bal > 0) if ((item = pt_deq(pt->c->out, pt->c->tid)) != PT_EMPTY) dfp(darg, item), bal--;
	return 0;
}

static int pt_parallel(pt_t *pt, pt_worker_t wfp, void **warg, void **src, void **dst)
{
	if (pt_set_worker(pt, wfp, warg)) return -1;
	for (uint64_t i = 1; i < pt->nth; ++i) pt_enq(pt->c->in, pt->c->tid, src? src[i] : NULL);
	void *res = wfp(*warg, src? src[0] : NULL); if (dst && dst[0]) dst[0] = res;
	for (uint64_t i = 1; i < pt->nth; ++i) {
		while ((res = pt_deq(pt->c->out, pt->c->tid)) == PT_EMPTY) {}
		if (dst && dst[i]) dst[i] = res;
	}
	return 0;
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
	if (magic_len != 4 || strncmp(buf, "BAM\001", 4) != 0) {
		// fprintf(stderr, "[bam_read_header] invalid BAM binary header (this is not a BAM file).\n");
		return NULL;
	}
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
typedef struct {
	gzFile fp;
	kseq_t *ks;
	bam_header_t *bh;
	uint8_v buf;
	uint16_t *tags;
	uint32_t l_tags, base_rid;
	uint8_t is_eof, keep_qual, keep_comment;
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
	uint32_t l_seq, rid, l_name, l_tag;
	// pointers
	char *name;
	uint8_t *seq, *qual, *tag;
} bseq_t;
_static_assert(sizeof(bseq_t) == 48);
typedef struct { size_t n, m; bseq_t *a; } bseq_v;

// ascii to 4bit conversion
static unsigned char seq_nt4_table_4bit[32] = {
	0, 1, 0, 2,  0, 0, 0, 4,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  8, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
};

static uint64_t bseq_search_tag(uint32_t l_tags, const uint16_t *tags, uint16_t t1, uint16_t t2)
{
	v32i16_t tv = _set_v32i16((t2<<8) | t1);
	for (uint64_t i = 0; i < ((l_tags + 0x1f) & ~0x1f); i+=0x20) {
		if (((v32_masku_t){ .mask = _mask_v32i16(_eq_v32i16(_load_v32i16(&tags[i]), tv)) }).all != 0) return 1;
	}
	return 0;
}

static bseq_file_t *bseq_open(const char *fn, uint32_t base_rid, uint32_t keep_qual, uint32_t l_tags, const uint16_t *tags)
{
	int c;
	bseq_file_t *fp;
	gzFile f;

	set_info("[bseq_open] initialize bseq object");
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (f == 0) return 0;
	fp = (bseq_file_t*)calloc(1, sizeof(bseq_file_t));
	fp->base_rid = base_rid; fp->keep_qual = keep_qual; fp->fp = f;
	for (uint64_t i = 0; i < 4; i++) {	// allow some invalid spaces at the head
		if ((c = gzgetc(fp->fp)) == 'B') {	// test bam signature
			gzungetc(c, fp->fp), fp->bh = bam_read_header(fp->fp); break;
		} else if (c == '>' || c == '@') {
			gzungetc(c, fp->fp), fp->ks = kseq_init(fp->fp); break;
		}
	}
	if (!fp->bh && !fp->ks) { free(fp); return 0; }

	fp->tags = calloc(((fp->l_tags = l_tags) + 0x1f) & ~0x1f, sizeof(uint16_t));
	if (l_tags && tags) memcpy(fp->tags, tags, l_tags * sizeof(uint16_t));
	fp->keep_comment = bseq_search_tag(fp->l_tags, fp->tags, 'C', 'O');
	return fp;
}

static uint32_t bseq_close(bseq_file_t *fp)
{
	if (fp == 0) return 0;
	uint32_t base_rid = fp->base_rid;
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp->tags); free(fp);
	return base_rid;
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

static void bseq_read_bam(bseq_file_t *fp, uint64_t size, bseq_v *seq, uint8_v *mem)
{
	// parse bam
	bseq_t *s;
	uint8_t *sname, *sseq, *squal, *stag;
	uint32_t l_tag;
	bam_core_t *c;

	if (fp->buf.n < size) kv_reserve(uint8_t, fp->buf, size);
	gzread(fp->fp, fp->buf.a, size);
	if ((c = (bam_core_t*)fp->buf.a)->flag&0x900) return;	// skip supp/secondary
	sname = fp->buf.a + sizeof(bam_core_t);
	sseq = sname + c->l_qname + sizeof(uint32_t)*c->n_cigar;
	squal = sseq + (c->l_qseq+1) / 2;
	stag = squal + c->l_qseq; l_tag = ((uint8_t*)fp->buf.a + size) - stag;

	kv_pushp(bseq_t, *seq, &s);
	kv_reserve(uint8_t, *mem, mem->n + c->l_qname + c->l_qseq + ((fp->keep_qual && *squal != 0xff)? c->l_qseq : 0) + (fp->l_tags? l_tag : 0) + 3);
	s->l_seq = c->l_qseq;
	s->rid = seq->n + fp->base_rid - 1;
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
	// if (fp->l_tags && l_tag) { memcpy(mem->a + mem->n, stag, l_tag); mem->n += l_tag; }
	if (fp->l_tags && l_tag) mem->n += bseq_save_tags(fp->l_tags, fp->tags, l_tag, stag, mem->a + mem->n);
	mem->a[mem->n++] = '\0';
	return;
}

static void bseq_read_fasta(bseq_file_t *fp, bseq_v *seq, uint8_v *mem)
{
	bseq_t *s;
	kseq_t *ks = fp->ks;

	uint64_t l_tag = (fp->keep_comment && ks->comment.l)? ks->comment.l+3 : 0;
	kv_pushp(bseq_t, *seq, &s);
	kv_reserve(uint8_t, *mem, mem->n + ks->name.l + ks->seq.l + ((fp->keep_qual && ks->qual.l)? ks->seq.l : 0) + l_tag + 4);
	s->l_seq = ks->seq.l;
	s->rid = seq->n + fp->base_rid - 1;
	s->l_name = ks->name.l;
	s->l_tag = l_tag;	// comment is transferred to CO:Z: tag if bseq_open is initialized with "CO" tag option

	s->name = (char *)mem->n;
	memcpy(mem->a + mem->n, ks->name.s, ks->name.l); mem->n += ks->name.l;
	mem->a[mem->n++] = '\0';

	s->tag = (uint8_t *)mem->n;
	if (fp->keep_comment && ks->comment.l) {
		mem->a[mem->n++] = 'C'; mem->a[mem->n++] = 'O'; mem->a[mem->n++] = 'Z';
		memcpy(mem->a + mem->n, ks->comment.s, ks->comment.l); mem->n += ks->comment.l;
	}
	mem->a[mem->n++] = '\0';

	s->seq = (uint8_t *)mem->n;
	for (uint64_t i = 0; i < ks->seq.l; ++i) mem->a[mem->n++] = seq_nt4_table_4bit[0x1f&ks->seq.s[i]];
	mem->a[mem->n++] = '\0';

	s->qual = (uint8_t *)mem->n;
	if (fp->keep_qual && ks->qual.l) { memcpy(mem->a + mem->n, ks->qual.s, ks->qual.l); mem->n += ks->qual.l; }
	mem->a[mem->n++] = '\0';
	return;
}

static bseq_t *bseq_read(bseq_file_t *fp, uint64_t chunk_size, uint32_t *n, void **base, uint64_t *size)
{
	uint8_v mem = {0};
	bseq_v seq = {0};
	static const uint8_t margin[64] = {0};

	set_info("[bseq_read] read sequence block from file");
	kv_reserve(uint8_t, mem, chunk_size + 128);
	kv_pushm(uint8_t, mem, margin, 64);
	if (fp->bh) {
		uint32_t size;
		while (gzread(fp->fp, &size, 4) == 4) {
			bseq_read_bam(fp, size, &seq, &mem);
			if (mem.n >= chunk_size) break;
		}
	} else {
		while (kseq_read(fp->ks) >= 0) {
			bseq_read_fasta(fp, &seq, &mem);
			if (mem.n >= chunk_size) break;
		}
	}
	kv_pushm(uint8_t, mem, margin, 64);

	for (uint64_t i = 0; i < seq.n; i++) {
		seq.a[i].name += (ptrdiff_t)mem.a;
		seq.a[i].seq += (ptrdiff_t)mem.a;
		seq.a[i].qual += (ptrdiff_t)mem.a;
		seq.a[i].tag += (ptrdiff_t)mem.a;
	}
	if (seq.n == 0) free(mem.a), mem.a = 0, fp->is_eof = 1;
	fp->base_rid += seq.n;
	*n = seq.n;
	*base = (void*)mem.a;
	*size = mem.n;
	return seq.a;
}

static int bseq_eof(bseq_file_t *fp)
{
	return fp->is_eof;
}

/* end of bseq.c */

/* sketch.c */

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
		if (q+64 > t) { p->n = q-p->a; p->m = MAX2(128, p->m*2); p->a = realloc(p->a, p->m*sizeof(mm128_t)); q = p->a+p->n; t = p->a+p->m; }
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
			buf[l&0x0f] = (mm128_t){ .u64 = { h = hash64(k0 < k1? k0 : k1, mask), (k0 < k1? 0 : 0xffffffff) ^ i } };
			if (h <= min) min = h, min_pos = l & 0x0f;
			l++;
		}
		for (uint64_t j = kk; j < kk + w; ++j) if (buf[j&0x0f].u64[0] == min) *q++ = buf[j&0x0f];
		q--;
		while (1) {
			uint64_t c, h;
			if ((c = seq4[++i]) == 0) goto _loop_tail;
			_push_kmer(c); if (k0 == k1) continue;
			buf[l&0x0f] = (mm128_t){ .u64 = { h = hash64(k0 < k1? k0 : k1, mask), (k0 < k1? 0 : 0xffffffff) ^ i } };

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

#define MM_BLAST6		( 0x01ULL<<56 )
#define MM_BLASR1		( 0x02ULL<<56 )
#define MM_BLASR4		( 0x03ULL<<56 )
#define MM_PAF			( 0x04ULL<<56 )

typedef struct {
	uint32_t sidx, eidx, n_threads, min, k, w, b;
	uint64_t flag;
	double min_ratio;
	int32_t m, x, gi, ge, xdrop, llim, hlim, elim, blim;
	uint32_t n_frq;
	float frq[16];
	uint16_v tags;
	uint64_t batch_size, outbuf_size;
	char *rg_line, *rg_id;
} mm_mapopt_t;

static void mm_mapopt_destroy(mm_mapopt_t *opt)
{
	free(opt->rg_line); free(opt->rg_id); free(opt);
	return;
}

static mm_mapopt_t *mm_mapopt_init(void)
{
	set_info("[mm_mapopt_init] initialize mapopt object");
	mm_mapopt_t *opt = calloc(1, sizeof(mm_mapopt_t));
	*opt = (mm_mapopt_t){
		/* -f, -k, -w, -b, -T */ .k = 15, .w = 16, .b = 14, .flag = 0,
		/* -a, -b, -p, -q */ .m = 1, .x = 1, .gi = 1, .ge = 1, .xdrop = 50,
		/* -s, -m */ .min = 50, .min_ratio = 0.3,
		/* -f */ .n_frq = 3, .frq[0] = 0.05, .frq[1] = 0.01, .frq[2] = 0.001,
		/* -t */ .n_threads = 1,
		/* -R */ .rg_line = NULL, .rg_id = NULL,

		/* -S, -E */.sidx = 0, .eidx = 3,
		.hlim = 5000, .llim = 5000, .blim = 0, .elim = 200,
		.batch_size = 1024 * 1024,
		.outbuf_size = 512 * 1024
	};
	return opt;
}

static int mm_mapopt_check(mm_mapopt_t *opt, int (*_fprintf)(FILE*,const char*,...), FILE *_fp)
{
	int ret = 0;
	if (opt->k >= 32) _fprintf(_fp, "[M::%s] ERROR: k must be inside [1,32).\n", __func__), ret = 1;
	if (ret) return(ret);
	if (opt->w >= 16) _fprintf(_fp, "[M::%s] ERROR: w must be inside [1,16).\n", __func__), ret = 1;
	if (ret) return(ret);

	if (opt->m < 1 || opt->m > 5) _fprintf(_fp, "[M::%s] ERROR: Match award must be inside [1,5].\n", __func__), ret = 1;
	if (opt->x < 1 || opt->x > 5) _fprintf(_fp, "[M::%s] ERROR: Mismatch penalty must be inside [1,5].\n", __func__), ret = 1;
	if (opt->gi < 1 || opt->gi > 5) _fprintf(_fp, "[M::%s] ERROR: Gap open penalty must be inside [1,5].\n", __func__), ret = 1;
	if (opt->ge < 1 || opt->ge > 5) _fprintf(_fp, "[M::%s] ERROR: Gap extension penalty must be inside [1,5].\n", __func__), ret = 1;
	if (ret) return(ret);

	if (opt->x >= (opt->gi + opt->ge))
		_fprintf(_fp, "[M::%s] info: Large mismatch penalty with respect to the gap open/extend penalty may cause SEGV or broken CIGAR. [issue #2]\n", __func__);
	if (opt->m + 2*(opt->gi + opt->ge) > 10)
		_fprintf(_fp, "[M::%s] info: Large match award or large gap open/extend penalty may cause SEGV or broken CIGAR. [issue #7]\n", __func__);
	if (opt->xdrop < 10 || opt->xdrop > 100) _fprintf(_fp, "[M::%s] ERROR: Xdrop cutoff must be inside [10,100].\n", __func__), ret = 1;
	if (opt->min > INT32_MAX) _fprintf(_fp, "[M::%s] ERROR: Minimum alignment score must be > 0.\n", __func__), ret = 1;
	if (opt->min_ratio < 0.0 || opt->min_ratio > 1.0) _fprintf(_fp, "[M::%s] ERROR: Minimum alignment score ratio must be inside [0.0,1.0].\n", __func__), ret = 1;
	if (opt->n_frq >= 16) _fprintf(_fp, "[M::%s] ERROR: Frequency thresholds must be fewer than 16.\n", __func__), ret = 1;
	for (uint64_t i = 0; i < opt->n_frq; ++i)
		if (opt->frq[i] < 0.0 || opt->frq[i] > 1.0 || (i != 0 && opt->frq[i-1] < opt->frq[i]))
			_fprintf(_fp, "[M::%s] ERROR: Frequency thresholds must be inside [0.0,1.0] and descending.\n", __func__), ret = 1;
	if (ret) return(ret);

	if (opt->n_threads < 1) _fprintf(_fp, "[M::%s] ERROR: Thread counts must be > 0.\n", __func__), ret = 1;
	if (opt->batch_size < 64 * 1024) _fprintf(_fp, "[M::%s] ERROR: Batch size must be > 64k.\n", __func__), ret = 1;
	if (opt->outbuf_size < 64 * 1024) _fprintf(_fp, "[M::%s] ERROR: Output buffer size must be > 64k.\n", __func__), ret = 1;
	if (ret) return(ret);

	if (opt->sidx >= 16) _fprintf(_fp, "[M::%s] ERROR: sidx must be inside [0,16).\n", __func__), ret = 1;
	if (opt->eidx >= 16) _fprintf(_fp, "[M::%s] ERROR: eidx must be inside [0,16).\n", __func__), ret = 1;
	if (opt->hlim < 100 || opt->hlim >= 100000) _fprintf(_fp, "[M::%s] ERROR: hlim must be inside [100,100000).\n", __func__), ret = 1;
	if (opt->llim < 100 || opt->llim >= 100000) _fprintf(_fp, "[M::%s] ERROR: llim must be inside [100,100000).\n", __func__), ret = 1;
	if (opt->blim >= 10000) _fprintf(_fp, "[M::%s] ERROR: blim must be inside [0,10000).\n", __func__), ret = 1;
	if (opt->elim >= 10000) _fprintf(_fp, "[M::%s] ERROR: elim must be inside [0,10000).\n", __func__), ret = 1;
	if (ret) return(ret);

	#define _dup(x)	({ char *_p = malloc(strlen(x)+1); memcpy(_p, (x), strlen(x)); _p[strlen(x)] = '\0'; _p; })
	if (opt->flag&(0x01ULL<<MM_RG)) {
		if (opt->rg_line == NULL && opt->rg_id == NULL) {
			opt->rg_line = _dup("@RG\tID:1\tSM:default"); opt->rg_id = _dup("1");
		}
	}
	#undef _dup
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
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

typedef struct {
	uint64_t mask;
	mm_idx_bucket_t *B;
	mm_idx_seq_v s;
	uint8_t b, w, k, circular;
	uint32_t base_rid;

	// work
	mm128_v a;
	uint64_v size;
	ptr_v base;
} mm_idx_t;

static mm_idx_t *mm_idx_init(uint32_t w, uint32_t k, uint32_t b)
{
	mm_idx_t *mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w<1? 1 : w, mi->k = k; mi->b = MIN2(k*2, b); mi->mask = (1<<mi->b) - 1;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	return mi;
}

static void mm_idx_destroy(mm_idx_t *mi)
{
	if (mi == 0) return;
	for (uint64_t i = 0; i < 1ULL<<mi->b; ++i) {
		free(mi->B[i].p);
		free(mi->B[i].a.a);
		kh_destroy((kh_t*)mi->B[i].h);
	}
	free(mi->B);
	for (uint64_t i = 0; i < mi->base.n; ++i) free(mi->base.a[i]);
	free(mi->base.a); free(mi->size.a);
	free(mi->s.a); free(mi->a.a); free(mi);
}

static const v2u32_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, uint64_t *n)
{
	mm_idx_bucket_t *b = &mi->B[minier & mi->mask];
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

static uint32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	uint64_t n = 0;
	uint32_t thres;

	set_info("[mm_idx_cal_max_occ] calculate occurrence thresholds");
	if (f <= 0.) return UINT32_MAX;
	for (uint64_t i = (n = 0); i < 1ULL<<mi->b; ++i)
		if (mi->B[i].h) n += kh_cnt((kh_t*)mi->B[i].h);
	
	uint32_t *a = (uint32_t*)malloc(n * sizeof(uint32_t));
	for (uint64_t i = (n = 0); i < 1ULL<<mi->b; ++i) {
		kh_t *h = (kh_t*)mi->B[i].h;
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
	uint64_t batch_size;
	bseq_file_t *fp;
	mm_idx_t *mi;
} mm_idx_pipeline_t;

typedef struct {
    uint32_t n_seq;
	const bseq_t *seq;	// const!
	void *base;
	uint64_t size;
	mm128_v a;
} mm_idx_step_t;

static void *mm_idx_source(void *arg)
{
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t*)arg;
	mm_idx_step_t *s = (mm_idx_step_t*)calloc(1, sizeof(mm_idx_step_t));
	uint64_t size;
	void *base;

	set_info("[mm_idx_source] fetch sequence block");
	s->seq = bseq_read(q->fp, q->batch_size, &s->n_seq, &base, &size);
	if (s->seq == 0) { free(s), s = 0; return 0; }

	kv_push(void*, q->mi->base, base);
	kv_push(uint64_t, q->mi->size, size);
	return s;
}

static void *mm_idx_worker(void *arg, void *item)
{
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t*)arg;
	mm_idx_step_t *s = (mm_idx_step_t*)item;

	char buf[128], *p = buf;
	p += _pstr(p, "[mm_idx_worker] bin id "); p += _pnum(uint32_t, p, s->seq[0].rid); p += _pstr(p, ":"); p += _pnum(uint32_t, p, s->seq[s->n_seq-1].rid);
	set_info(buf);
	for (uint64_t i = 0; i < s->n_seq; ++i)
		mm_sketch(s->seq[i].seq, s->seq[i].l_seq, q->mi->w, q->mi->k, s->seq[i].rid, &s->a);
	return s;
}

static void mm_idx_drain(void *arg, void *item)
{
	set_info("[mm_idx_drain] dump minimizers to pool");
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t*)arg;
	mm_idx_step_t *s = (mm_idx_step_t*)item;
	q->mi->s.n = MAX2(q->mi->s.n, s->seq[s->n_seq-1].rid+1-q->mi->base_rid);
	kv_reserve(mm_idx_seq_t, q->mi->s, q->mi->s.n);
	for (uint64_t i = 0; i < s->n_seq; ++i) {
		const bseq_t *p = &s->seq[i];
		q->mi->s.a[s->seq[i].rid-q->mi->base_rid] = (mm_idx_seq_t){
			.l_seq = p->l_seq, .rid = p->rid, .l_name = p->l_name, .circular = q->mi->circular,
			.name = p->name, .seq = p->seq
		};
	}
	uint64_t mask = q->mi->mask;
	for (uint64_t i = 0; i < s->a.n; ++i) {
		mm128_v *p = &q->mi->B[s->a.a[i].u64[0]&mask].a;
		kv_push(mm128_t, *p, s->a.a[i]);
	}

	free((void*)s->seq); free(s->a.a); free(s);
	return;
}

typedef struct {
	mm_idx_t *mi;
	uint32_t from, to;
} mm_idx_post_t;

static void *mm_idx_post(void *arg, void *item)
{
	set_info("[mm_idx_post] indexing postprocess");
	mm_idx_post_t *q = (mm_idx_post_t*)arg;
	mm_idx_t *mi = q->mi;

	for (uint64_t i = q->from; i < q->to; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		if (b->a.n == 0) continue;

		// sort by minimizer
		radix_sort_128x(b->a.a, b->a.a + b->a.n);

		// count and preallocate
		uint64_t n_keys = 0; b->n = 0;
		for (uint64_t j = 1, n = 1; j <= b->a.n; ++j, ++n) {
			if (j == b->a.n || b->a.a[j].u64[0] != b->a.a[j-1].u64[0]) {
				b->n += (n > 1)? n : 0; ++n_keys; n = 0;
			}
		}
		kh_t *h = kh_init(n_keys);
		b->p = (uint64_t*)malloc(sizeof(uint64_t) * b->n);

		// create the hash table
		for (uint64_t j = 1, n = 1, sp = 0; j <= b->a.n; ++j, ++n) {
			if (j == b->a.n || b->a.a[j].u64[0] != b->a.a[j-1].u64[0]) {
				mm128_t *p = &b->a.a[j-n];
				uint64_t key = p->u64[0]>>mi->b, val = p->u64[1];
				if (n != 1) {
					b->p[sp] = val;	// k = 0
					for (uint64_t k = 1; k < n; ++k) b->p[sp+k] = b->a.a[j-n+k].u64[1];
					val = sp<<32 | n | 0x01ULL<<63;
				}
				kh_put(h, key, val); n = 0;
			}
		}
		b->h = h;

		// deallocate and clear b->a
		free(b->a.a);
		b->a.n = b->a.m = 0, b->a.a = 0;
	}
	return 0;
}

static mm_idx_t *mm_idx_gen(const mm_mapopt_t *opt, bseq_file_t *fp)
{
	mm_idx_pipeline_t pl = {0}, **p;

	set_info("[mm_idx_gen] initialize index object");
	pl.batch_size = opt->batch_size;
	pl.fp = fp;
	if (pl.fp == 0) return 0;
	pl.mi = mm_idx_init(opt->w, opt->k, opt->b); pl.mi->circular = (opt->flag&MM_CIRCULAR)!=0; pl.mi->base_rid = pl.fp->base_rid;

	p = (mm_idx_pipeline_t**)calloc(opt->n_threads, sizeof(mm_idx_pipeline_t*));
	for (uint64_t i = 0; i < opt->n_threads; ++i) p[i] = &pl;

	pt_t *pt = pt_init(opt->n_threads);
	pt_stream(pt, mm_idx_source, &pl, mm_idx_worker, (void**)p, mm_idx_drain, &pl);
	free(p);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	mm_idx_post_t *q = (mm_idx_post_t*)calloc(opt->n_threads, sizeof(mm_idx_post_t));
	mm_idx_post_t **qq = (mm_idx_post_t**)calloc(opt->n_threads, sizeof(mm_idx_post_t*));
	for (uint64_t i = 0; i < opt->n_threads; ++i) {
		q[i].mi = pl.mi;
		q[i].from = (1ULL<<pl.mi->b)*i/opt->n_threads;
		q[i].to = (1ULL<<pl.mi->b)*(i+1)/opt->n_threads;
		qq[i] = &q[i];
	}
	pt_parallel(pt, mm_idx_post, (void**)qq, NULL, NULL);
	pt_destoroy(pt);

	free(q); free(qq);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
	return pl.mi;
}

/*************
 * index I/O *
 *************/

#define MM_IDX_MAGIC "MAI\5"		/* minialign index version 5, differs from minimap index signature */

static void mm_idx_dump(gzFile fp, const mm_idx_t *mi)
{
	uint32_t x[3];
	uint64_t size = 0, y[2];
	set_info("[mm_idx_dump] dump index to file");
	for (uint64_t i = (size = 0); i < mi->size.n; ++i) size += mi->size.a[i];
	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b; y[0] = mi->s.n, y[1] = size;
	gzwrite(fp, MM_IDX_MAGIC, sizeof(char) * 4);
	gzwrite(fp, x, sizeof(uint32_t) * 3);
	gzwrite(fp, y, sizeof(uint64_t) * 2);
	for (uint64_t i = 0; i < 1ULL<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		gzwrite(fp, &b->n, sizeof(uint32_t) * 1);
		gzwrite(fp, b->p, sizeof(uint64_t) * b->n);
		kh_dump(fp, (kh_t*)b->h);
	}
	for (uint64_t i = 0; i < mi->base.n; ++i)
		gzwrite(fp, mi->base.a[i], sizeof(char) * mi->size.a[i]);
	for (uint64_t i = 0, j = 0, s = 0; i < mi->s.n; ++i) {
		if ((uintptr_t)mi->s.a[i].seq < (uintptr_t)mi->base.a[j]
		|| (uintptr_t)mi->s.a[i].seq >= (uintptr_t)mi->base.a[j] + (ptrdiff_t)mi->size.a[j]) {
			s += mi->size.a[j++];
		}
		mi->s.a[i].name -= (ptrdiff_t)mi->base.a[j], mi->s.a[i].seq -= (ptrdiff_t)mi->base.a[j];
		mi->s.a[i].name += (ptrdiff_t)s, mi->s.a[i].seq += (ptrdiff_t)s;
	}
	gzwrite(fp, mi->s.a, sizeof(mm_idx_seq_t) * mi->s.n);
	// restore pointers
	for (uint64_t i = 0, j = 0, s = 0; i < mi->s.n; ++i) {
		mi->s.a[i].name -= (ptrdiff_t)s, mi->s.a[i].seq -= (ptrdiff_t)s;
		mi->s.a[i].name += (ptrdiff_t)mi->base.a[j], mi->s.a[i].seq += (ptrdiff_t)mi->base.a[j];
		if ((uintptr_t)mi->s.a[i].name > s + mi->size.a[j]) s += mi->size.a[j++];
	}
	return;
}

static mm_idx_t *mm_idx_load(gzFile fp)
{
	char magic[4];
	uint32_t x[3];
	uint64_t bsize, y[2];

	mm_idx_t *mi;
	set_info("[mm_idx_load] load index from file");
	if (gzread(fp, magic, sizeof(char) * 4) != sizeof(char) * 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (gzread(fp, x, sizeof(uint32_t) * 3) != sizeof(uint32_t) * 3) return 0;
	if (gzread(fp, y, sizeof(uint64_t) * 2) != sizeof(uint64_t) * 2) return 0;
	mi = mm_idx_init(x[0], x[1], x[2]); mi->s.n = y[0]; bsize = y[1];
	for (uint64_t i = 0; i < 1ULL<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		if (gzread(fp, &b->n, sizeof(uint32_t) * 1) != sizeof(uint32_t) * 1) goto _mm_idx_load_fail;
		b->p = (uint64_t*)malloc(b->n * sizeof(uint64_t));
		if (gzread(fp, b->p, sizeof(uint64_t) * b->n) != sizeof(uint64_t) * (size_t)b->n) goto _mm_idx_load_fail;
		b->h = kh_load(fp);
	}
	mi->base.n = mi->size.n = 1;
	mi->base.a = malloc(sizeof(void*) * mi->base.n);
	mi->size.a = malloc(sizeof(uint64_t) * mi->size.n);
	mi->base.a[0] = malloc(sizeof(char) * bsize);
	mi->size.a[0] = bsize;
	if (gzread(fp, mi->base.a[0], sizeof(char) * mi->size.a[0]) != sizeof(char) * mi->size.a[0]) goto _mm_idx_load_fail;
	mi->s.a = malloc(sizeof(mm_idx_seq_t) * mi->s.n);
	if (gzread(fp, mi->s.a, sizeof(mm_idx_seq_t) * mi->s.n) != sizeof(mm_idx_seq_t) * mi->s.n) goto _mm_idx_load_fail;
	for (uint64_t i = 0; i < mi->s.n; ++i) mi->s.a[i].name += (ptrdiff_t)mi->base.a[0], mi->s.a[i].seq += (ptrdiff_t)mi->base.a[0];
	mi->base_rid = mi->s.a[0].rid;
	return mi;
_mm_idx_load_fail:
	mm_idx_destroy(mi);
	return 0;
}

/* end of index.c */

/* map.c */
/**************************
 * Multi-threaded mapping *
 **************************/
typedef struct mm_align_s mm_align_t;
typedef struct {
	void (*header)(mm_align_t*);
	void (*unmapped)(mm_align_t*, const bseq_t*);
	void (*mapped)(mm_align_t*, const bseq_t*, uint32_t, const mm128_t*);
} mm_printer_t;

struct mm_align_s {
	uint8_t *base, *tail, *p;
	const mm_idx_t *mi;
	const mm_mapopt_t *opt;
	uint64_t *occ;
	uint32_t n_occ;
	bseq_file_t *fp;
	mm_printer_t printer;
	gaba_t *gaba;
	void **t;	// mm_tbuf_t **
	pt_t *pt;
};

typedef struct {
	// query sequences
	uint32_t n_seq;
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
	// poshash_t *pos;
	kh_t *pos;
	gaba_dp_t *dp;	// alignment work
} mm_tbuf_t;

#define _s(x)		( (x)<0?-1:1)
#define _m(x)		( (((int32_t)(x))>>31)^(x) )
static void mm_expand(uint32_t n, const v2u32_t *r, uint32_t qid, int32_t qs, uint32_t thresh, mm128_v *coef)
{
	const int32_t ofs = 0x40000000;
	for (uint64_t i = 0; i < n; ++i) {	// iterate over all the collected minimizers
		if (ofs + r[i].x[1] - qid > thresh) continue;	// 0x3fffffff to skip self and dual, 0xffffffff to keep all
		int32_t rs = (int32_t)r[i].x[0], _qs = (rs>>31) ^ qs, _rs = (rs>>31) ^ rs;
		mm128_t *p;
		kv_pushp(mm128_t, *coef, &p);
		p->u32[0] = ofs + _rs - (_qs>>1); p->u32[1] = r[i].x[1]; p->u32[2] = _rs; p->u32[3] = _qs;
	}
	return;
}

static void mm_collect(const mm_idx_t *mi, uint32_t l_seq, const uint8_t *seq, uint32_t qid, uint32_t max_occ, uint32_t resc_occ, uint32_t thresh, mm128_v *mini, mm128_v *coef)
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
		mm_expand(n, r, qid, qs, thresh, coef);
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
	gaba_dp_t *dp, gaba_section_t *qf, gaba_section_t *qr, gaba_section_t *t, kh_t *pos/*poshash_t *pos*/, lmm_t *lmm)
{
	const uint32_t mask = 0x7fffffff;
	const uint8_t *lim = (const uint8_t*)0x800000000000;
	gaba_section_t *qu, *qd, *r, *q;
	gaba_section_t rf = { .id = ref->rid<<1, .len = ref->l_seq, .base = (const uint8_t*)ref->seq };
	gaba_section_t rr = { .id = (ref->rid<<1)+1, .len = ref->l_seq, .base = gaba_rev((const uint8_t*)ref->seq+ref->l_seq-1, lim) };
	gaba_pos_pair_t p = {0};
	gaba_dp_flush(dp, lim, lim);
	uint64_t key = (uint64_t)ref->rid<<32, *pval;
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
			if (f->status & GABA_STATUS_UPDATE_A) flag |= GABA_STATUS_UPDATE_A, r = t;
			if (f->status & GABA_STATUS_UPDATE_B) flag |= GABA_STATUS_UPDATE_B, q = t;
			if ((f = gaba_dp_fill(dp, f, r, q) ) == NULL) goto _abort;
			m = (f->max > m->max)? f : m;
		} while(!(flag & f->status));
		// find max
		p = gaba_dp_search_max(dp, m);
		// check duplicate
		key |= p.apos - (p.bpos>>1);
		if ((pval = kh_get_ptr(pos, key)) != NULL) { return 0; }	// already evaluated
		// downward extension from max
		gaba_dp_flush_stack(dp, stack);
		if ((m = f = gaba_dp_fill_root(dp, r = &rf, ref->l_seq-p.apos-1, q = qd, qd->len-p.bpos-1)) == NULL) goto _abort;
		flag = GABA_STATUS_TERM;
		do {
			if (f->status & GABA_STATUS_UPDATE_A) flag |= GABA_STATUS_UPDATE_A, r = t;
			if (f->status & GABA_STATUS_UPDATE_B) flag |= GABA_STATUS_UPDATE_B, q = t;
			if ((f = gaba_dp_fill(dp, f, r, q)) == NULL) goto _abort;
			m = (f->max > m->max)? f : m;
		} while(!(flag & f->status));
		if (m->max < min) { key &= 0xffffffff00000000; continue; }
		// convert alignment to cigar
		a = gaba_dp_trace(dp, NULL, m, GABA_TRACE_PARAMS( .lmm = lmm ));
		break;
	}
	// record head
	kh_put(pos, key, (uintptr_t)a);
	return a;
_abort:;
	oom_abort(__func__);
	return 0;
}

#define MAPQ_DEC	( 4 )
#define MAPQ_COEF	( 1<<MAPQ_DEC )
#define _clip(x)	MAX2(0, MIN2(((uint32_t)(x)), 60 * MAPQ_COEF))
#define _aln(x)		( (const gaba_alignment_t*)(x).u64[1] )
static uint32_t mm_post_map(const mm_mapopt_t *opt, uint32_t n_reg, mm128_t *reg)
{
	uint32_t i, score = _aln(reg[0])->score, min = score * opt->min_ratio;
	uint32_t ssc = (n_reg>1)? _aln(reg[1])->score : 0, bsc = (n_reg>1)? _aln(reg[n_reg-1])->score : 0, tsc = 0;
	double elen = (double)gaba_plen(_aln(reg[0])->sec) / 2.0, pid = 1.0 - (double)(elen * opt->m - score) / (double)(opt->m + opt->x) / elen;
	double ec = 2.0 / (pid * (double)(opt->m + opt->x) - (double)opt->x), ulen = ec * (score - ssc), pe = 1.0 / (ulen * ulen + (double)n_reg);

	for (i = 1; i < n_reg; ++i) tsc += _aln(reg[i])->score - bsc + 1;
	reg[0].u32[0] |= _clip(-10.0 * MAPQ_COEF * log10(pe));
	for (i = 1; i < n_reg && (score = _aln(reg[i])->score) >= min; ++i)
		reg[i].u32[0] |= (0x100<<16) | _clip(-10.0 * MAPQ_COEF * log10(1.0 - pe * (double)(score - bsc + 1) / (double)tsc));
	return i;
}

static uint32_t mm_post_ava(const mm_mapopt_t *opt, uint32_t n_reg, mm128_t *reg)
{
	uint32_t i, score, min = _aln(reg[0])->score * opt->min_ratio;
	for (i = 0; i < n_reg && (score = _aln(reg[i])->score) >= min; ++i) {
		double elen = (double)gaba_plen(_aln(reg[i])->sec) / 2.0, pid = 1.0 - (double)(elen * opt->m - score) / (double)(opt->m + opt->x) / elen;
		double ec = 2.0 / (pid * (double)(opt->m + opt->x) - (double)opt->x), ulen = ec * score, pe = 1.0 / (ulen + 1);
		reg[i].u32[0] |= _clip(-10.0 * MAPQ_COEF * log10(pe)) | ((i == 0? 0 : 0x800)<<16);
	}
	return i;
}
#undef _clip
#undef _aln

static const mm128_t *mm_align_seq(
	mm_tbuf_t *b, const mm_mapopt_t *opt, const mm_idx_t *mi, uint32_t l_seq, const uint8_t *seq, uint32_t qid,
	uint32_t n_occ, uint64_t *occ, lmm_t *lmm, uint64_t *n_reg)
{
	int32_t min = opt->min;
	const int32_t ofs = 0x40000000;
	const uint32_t mask = 0x7fffffff, thresh = opt->flag&MM_AVA? (uint32_t)ofs - 1 : 0xffffffff;
	mm128_v reg = {0};
	// prepare section info for alignment
	uint8_t tail[96];
	const uint8_t *lim = (const uint8_t*)0x800000000000;
	gaba_section_t qf, qr, mg;
	if (l_seq == 0) return 0;
	// qid = opt->flag&MM_AVA? qid : 0xffffffff;
	memset(tail, 0, 96);
	qf.id = 0; qf.len = l_seq; qf.base = (const uint8_t*)seq;
	qr.id = 1; qr.len = l_seq; qr.base = gaba_rev((const uint8_t*)seq+l_seq-1, lim);
	mg.id = 4; mg.len = 32; mg.base = tail+32;
	b->coef.n = b->resc.n = b->intv.n = 0; kh_clear(b->pos);
	for (uint64_t i = 0, j = 0; reg.n == 0 && i < n_occ; ++i) {
		if (i == 0) {
			mm_collect(mi, l_seq, seq, qid, occ[n_occ-1], occ[0], thresh, &b->resc, &b->coef);
		} else {
			if (i == 1) radix_sort_128x(b->resc.a, b->resc.a + b->resc.n);
			for (uint64_t k = 0; k < b->coef.n; ++k) b->coef.a[k].u32[2] &= mask;
			for (; j < b->resc.n && b->resc.a[j].u32[1] <= occ[i]; ++j)
				mm_expand(b->resc.a[j].u32[1], (const v2u32_t*)b->resc.a[j].u64[1], qid, b->resc.a[j].u32[0], thresh, &b->coef);
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
	if (reg.n == 0) return 0;
	radix_sort_128x(reg.a, reg.a + reg.n);
	*n_reg = ((opt->flag & MM_AVA)? mm_post_ava : mm_post_map)(opt, reg.n, reg.a);
	for (uint64_t i = *n_reg; i < reg.n; ++i) { lmm_free(lmm, (void*)reg.a[i].u64[1]); }
	return reg.a;
}

#undef _s
#undef _m

#define BUF_SIZE	( 512 * 1024 )
#define _put(_buf, _c) { \
	*(_buf)->p++ = (_c); \
	if ((uintptr_t)(_buf)->p >= (uintptr_t)(_buf)->tail) { \
		fwrite((_buf)->base, sizeof(uint8_t), (_buf)->p - (_buf)->base, stdout); \
		(_buf)->p = (_buf)->base; \
	} \
}
#define _putfi(type, _buf, _n, _c) ({ \
	uint8_t _b[16] = {0}; \
	type _m = (type)(_n); int64_t _i = 0; \
	while (_m) _b[_i++] = _m % 10, _m /= 10; \
	_i += (_i==0); \
	_i += (_c - _i + 1 > 0)? _c - _i + 1 : 0; \
	for (int64_t _j = _i; _j > (_c); _j--) { _put(_buf, _b[_j-1]+'0'); } \
	_put(_buf, '.'); \
	for (int64_t _j = (_c); _j > 0; _j--) { _put(_buf, _b[_j-1]+'0'); } \
	_i; \
})
#define _puti(type, _buf, _n) ({ \
	uint8_t _b[16] = {0}; \
	type _m = (type)(_n); int64_t _i = 0; \
	while (_m) _b[_i++] = _m % 10, _m /= 10; \
	_i += (_i==0); \
	for (int64_t _j = _i; _j > 0; _j--) { _put(_buf, _b[_j-1]+'0'); } \
	_i; \
})
#define _putn(_buf, _n) _puti(uint32_t, _buf, _n)
#define _puts(_buf, _s) { \
	for (const uint8_t *_q = (const uint8_t*)(_s); *_q; _q++) { _put(_buf, *_q); } \
}

#define MM_RG			( 0 )		// Z: read group
#define MM_CO			( 1 )		// Z: comment
#define MM_NH			( 2 )		// i: #hits
#define MM_IH 			( 3 )		// i: index of the record in the hits
#define MM_AS			( 4 )		// i: score
#define MM_XS			( 5 )		// i: suboptimal score
#define MM_NM 			( 6 )		// i: editdist to the reference
#define MM_SA			( 7 )		// Z: supplementary records

static uint64_t mm_print_num(mm_align_t *b, uint8_t type, const uint8_t *p)
{
	if (type == 'a') { _put(b, *p); return 1; }
	else if (type == 'c') { _puti(int8_t, b, (int8_t)*p); return 1; }
	else if (type == 'C') { _puti(uint8_t, b, (uint8_t)*p); return 1; }
	else if (type == 's') { _puti(int16_t, b, *((int16_t*)p)); return 2; }
	else if (type == 'S') { _puti(uint16_t, b, *((uint16_t*)p)); return 2; }
	else if (type == 'i') { _puti(int32_t, b, *((int32_t*)p)); return 4; }
	else if (type == 'I') { _puti(uint32_t, b, *((uint32_t*)p)); return 4; }
	else if (type == 'f') {
		char buf[32]; uint64_t l = sprintf(buf, "%f", *((float*)p));
		for (uint64_t i = 0; i < l; ++i) { _put(b, buf[i]); }
		return 4;
	}
	return 1;
}

#define _t(b)			_put(b, '\t')
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
	_puts(b, "@HD\tVN:1.0\tSO:unsorted\n");
	for (uint64_t i = 0; i < mi->s.n; ++i) {
		_puts(b, "@SQ\tSN:"); _puts(b, mi->s.a[i].name);
		_puts(b, "\tLN:"); _putn(b, mi->s.a[i].l_seq);
		_cr(b);
	}
	if (b->opt->flag & 0x01ULL<<MM_RG) { _puts(b, b->opt->rg_line); _cr(b); }
	_puts(b, "@PG\tID:minialign\tPN:minialign\n");
	return;
}

static void mm_print_unmapped_sam(mm_align_t *b, const bseq_t *t)
{
	_puts(b, t->name); _puts(b, "\t4\t*\t0\t0\t*\t*\t0\t0\t");
	for (uint64_t k = 0; k < t->l_seq; k++) _put(b, "NACMGRSVTWYHKDBN"[(uint8_t)t->seq[k]]);
	_t(b);
	if (b->opt->flag&MM_KEEP_QUAL && t->qual[0] != '\0') { _puts(b, t->qual); }
	else { _put(b, '*'); }
	mm_restore_tags(b, t); _cr(b);
	return;
}

static int mm_cigar_printer(void *_b, int64_t len, char c)
{
	mm_align_t *b = (mm_align_t*)_b;
	len = _putn(b, (uint32_t)len); _put(b, c);
	return len+1;
}

#define _aln(x)		( (const gaba_alignment_t*)(x).u64[1] )
#define _mapq(x)	( (x).u32[0] & 0xffff )
#define _flag(x)	( (x).u32[0]>>16 )
static void mm_print_mapped_sam(mm_align_t *b, const bseq_t *t, uint32_t n_reg, const mm128_t *reg)
{
	for (uint64_t j = 0; j < n_reg; ++j) {
		uint16_t mapq = _mapq(reg[j]), flag = _flag(reg[j]);
		const gaba_alignment_t *a = _aln(reg[j]);
		const mm_idx_t *mi = b->mi;
		const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];
		const gaba_path_section_t *s = &a->sec[0];
		uint32_t rs = r->l_seq-s->apos-s->alen, hl = t->l_seq-s->bpos-s->blen, tl = s->bpos;
		uint32_t qs = (flag&0x900)? hl : 0, qe = t->l_seq - ((flag&0x900)? tl : 0);

		_puts(b, t->name); _t(b); _putn(b, flag); _t(b); _puts(b, r->name); _t(b);
		_putn(b, rs+1); _t(b); _putn(b, mapq>>MAPQ_DEC); _t(b);
		if (hl) { _putn(b, hl); _put(b, (flag&0x900)? 'H' : 'S'); }
		gaba_dp_print_cigar_reverse(mm_cigar_printer, b, a->path->array, 0, a->path->len);
		if (tl) { _putn(b, tl); _put(b, (flag&0x900)? 'H' : 'S'); }
		_puts(b, "\t*\t0\t0\t");
		if (flag&0x10) { for (int64_t k = t->l_seq-qs; k > t->l_seq-qe; k--) _put(b, "NTGKCYSBAWRDMHVN"[(uint8_t)t->seq[k-1]]); }
		else { for (int64_t k = qs; k < qe; k++) _put(b, "NACMGRSVTWYHKDBN"[(uint8_t)t->seq[k]]); }
		_t(b);
		if (b->opt->flag&MM_KEEP_QUAL && t->qual[0] != '\0') {
			if (flag&0x10) { for (int64_t k = t->l_seq-qs; k > t->l_seq-qe; k--) _put(b, t->qual[k-1]); }
			else { for (int64_t k = qs; k < qe; k++) _put(b, t->qual[k]); }
		} else {
			_put(b, '*');
		}
		// mm_print_tags(b, t, a, flag, idx);
		// print tags
		uint64_t f = b->opt->flag;
		if (f & 0x01ULL<<MM_RG) { _puts(b, "\tRG:Z:"); _puts(b, b->opt->rg_id); }
		if (f & 0x01ULL<<MM_NH) { _puts(b, "\tNH:i:"); _putn(b, n_reg); }
		if (f & 0x01ULL<<MM_IH) { _puts(b, "\tIH:i:"); _putn(b, j); }
		if (f & 0x01ULL<<MM_AS) { _puts(b, "\tAS:i:"); _putn(b, a->score); }
		if ((f & 0x01ULL<<MM_XS) && j == 0) { _puts(b, "\tXS:i:"); _putn(b, n_reg>1? _aln(reg[1])->score : 0); }
		if (f & 0x01ULL<<MM_NM) { _puts(b, "\tNM:i:"); _putn(b, a->xcnt + a->gecnt); }
		if (f & 0x01ULL<<MM_SA) { }
		if (j == 0) mm_restore_tags(b, t);
		_cr(b);
	}
	return;
}

#define _putd(b, _id)		_put(b, (_id)&0x01? '+' : '-');
// qname rname idt len #x #gi qs qe rs re e-value bitscore
static void mm_print_mapped_blast6(mm_align_t *b, const bseq_t *t, uint32_t n_reg, const mm128_t *reg)
{
	for (uint64_t j = 0; j < n_reg; ++j) {
		const gaba_alignment_t *a = _aln(reg[j]);
		const mm_idx_t *mi = b->mi;
		const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];
		const gaba_path_section_t *s = &a->sec[0];
		int32_t dcnt = (a->path->len - a->gecnt)>>1, slen = dcnt + a->gecnt;
		int32_t mid = 100000.0 * (float)(dcnt - a->xcnt) / (float)slen;	// percent identity
		uint32_t rs = s->bid&0x01? r->l_seq-s->apos-s->alen+1 : r->l_seq-s->apos;
		uint32_t re = s->bid&0x01? r->l_seq-s->apos : r->l_seq-s->apos-s->alen+1;
		uint32_t qs = t->l_seq-s->bpos-s->blen+1, qe = t->l_seq-s->bpos;

		_puts(b, t->name); _t(b); _puts(b, r->name); _t(b);
		_putfi(int32_t, b, mid, 3); _t(b); _putn(b, slen); _t(b); _putn(b, a->xcnt); _t(b); _putn(b, a->gicnt); _t(b);
		_putn(b, qs); _t(b); _putn(b, qe); _t(b); _putn(b, rs); _t(b); _putn(b, re); _t(b);

		float bit = 1.85 * (float)a->score - 0.02;	// fixme, lambda and k should be calcd from the scoring params
		int32_t e = 1000.0 * (float)r->l_seq * (float)t->l_seq * pow(2.0, -bit);	// fixme, lengths are not corrected
		_putfi(int32_t, b, e, 3); _t(b); _putn(b, (int32_t)bit); _cr(b);
	}
	return;
}

// qname rname qd rd score idt rs re rl qs qe ql #cells
static void mm_print_mapped_blasr1(mm_align_t *b, const bseq_t *t, uint32_t n_reg, const mm128_t *reg)
{
	for (uint64_t j = 0; j < n_reg; ++j) {
		const gaba_alignment_t *a = _aln(reg[j]);
		const mm_idx_t *mi = b->mi;
		const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];
		const gaba_path_section_t *s = &a->sec[0];
		int32_t dcnt = (a->path->len - a->gecnt)>>1, slen = dcnt + a->gecnt;
		int32_t mid = 1000000.0 * (float)(dcnt - a->xcnt) / (float)slen;	// percent identity
		uint32_t rs = s->bid&0x01? r->l_seq-s->apos-s->alen : s->apos, re = rs+s->alen;
		uint32_t qs = s->bid&0x01? t->l_seq-s->bpos-s->blen : s->bpos, qe = qs+s->blen;

		_puts(b, t->name); _put(b, '/'); _put(b, '0'); _put(b, '_'); _putn(b, t->l_seq); _sp(b);
		_puts(b, r->name); _sp(b); _put(b, '0'); _sp(b); _put(b, s->bid&0x01? '0' : '1'); _sp(b);
		_put(b, '-'); _putn(b, a->score); _sp(b);	// score in negative
		_putfi(int32_t, b, mid, 4); _sp(b);
		_putn(b, rs); _sp(b); _putn(b, re); _sp(b); _putn(b, r->l_seq); _sp(b);
		_putn(b, qs); _sp(b); _putn(b, qe); _sp(b); _putn(b, t->l_seq); _sp(b);
		_put(b, '0'); _cr(b);	// #cells is always 0
	}
	return;
}

// qname rname score idt qd qs qe ql rd rs re rl mapq
static void mm_print_mapped_blasr4(mm_align_t *b, const bseq_t *t, uint32_t n_reg, const mm128_t *reg)
{
	for (uint64_t j = 0; j < n_reg; ++j) {
		uint16_t mapq = _mapq(reg[j]);
		const gaba_alignment_t *a = _aln(reg[j]);
		const mm_idx_t *mi = b->mi;
		const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];
		const gaba_path_section_t *s = &a->sec[0];
		int32_t dcnt = (a->path->len - a->gecnt)>>1, slen = dcnt + a->gecnt;
		int32_t mid = 1000000.0 * (float)(dcnt - a->xcnt) / (float)slen;	// percent identity
		uint32_t rs = s->bid&0x01? r->l_seq-s->apos-s->alen : s->apos, re = rs+s->alen;
		uint32_t qs = s->bid&0x01? t->l_seq-s->bpos-s->blen : s->bpos, qe = qs+s->blen;

		_puts(b, t->name); _put(b, '/'); _put(b, '0'); _put(b, '_'); _putn(b, t->l_seq); _sp(b);
		_puts(b, r->name); _sp(b);
		_put(b, '-'); _putn(b, a->score); _sp(b); _putfi(int32_t, b, mid, 4); _sp(b);
		_put(b, '0'); _sp(b); _putn(b, qs); _sp(b); _putn(b, qe); _sp(b); _putn(b, t->l_seq); _sp(b);
		_put(b, s->bid&0x01? '0' : '1'); _sp(b); _putn(b, rs); _sp(b); _putn(b, re); _sp(b); _putn(b, r->l_seq); _sp(b);
		_putn(b, (mapq>>(MAPQ_DEC - 2)) + ((mapq>>(MAPQ_DEC + 2)) & ~0x01)); _cr(b);
	}
	return;
}

// qname ql qs qe qd rname rl rs re #m block_len mapq
static void mm_print_mapped_paf(mm_align_t *b, const bseq_t *t, uint32_t n_reg, const mm128_t *reg)
{
	for (uint64_t j = 0; j < n_reg; ++j) {
		uint16_t mapq = _mapq(reg[j]);
		const gaba_alignment_t *a = _aln(reg[j]);
		const mm_idx_t *mi = b->mi;
		const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];
		const gaba_path_section_t *s = &a->sec[0];
		uint32_t dcnt = (a->path->len - a->gecnt)>>1;
		uint32_t rs = r->l_seq-s->apos-s->alen, re = r->l_seq-s->apos, qs = t->l_seq-s->bpos-s->blen, qe = t->l_seq-s->bpos;
		_puts(b, t->name); _t(b); _putn(b, t->l_seq); _t(b); _putn(b, qs); _t(b); _putn(b, qe); _t(b); _putd(b, s->bid); _t(b);
		_puts(b, r->name); _t(b); _putn(b, r->l_seq); _t(b); _putn(b, rs); _t(b); _putn(b, re); _t(b);
		_putn(b, dcnt - a->xcnt); _t(b); _putn(b, dcnt + a->gecnt); _t(b); _putn(b, mapq>>MAPQ_DEC); _cr(b);
	}
	return;
}
#undef _d
#undef _t
#undef _cr
#undef _put
#undef _putn
#undef _puts
#undef _aln

static void *mm_align_source(void *arg)
{
	set_info("[mm_align_source] fetch sequence block");
	mm_align_t *b = (mm_align_t*)arg;
	mm_align_step_t *s = (mm_align_step_t*)calloc(1, sizeof(mm_align_step_t));
	s->lmm = lmm_init(NULL, 512 * 1024);
	if (s->lmm == 0) return 0;
	s->seq = bseq_read(b->fp, b->opt->batch_size, &s->n_seq, &s->base, &s->size);
	if (s->seq == 0) lmm_clean(s->lmm), free(s), s = 0;
	return s;
}

static void *mm_align_worker(void *arg, void *item)
{
	mm_tbuf_t *t = (mm_tbuf_t*)arg;
	mm_align_step_t *s = (mm_align_step_t*)item;

	char buf[128], *p = buf;
	p += _pstr(p, "[mm_align_worker] bin id "); p += _pnum(uint32_t, p, s->seq[0].rid); p += _pstr(p, ":"); p += _pnum(uint32_t, p, s->seq[s->n_seq-1].rid);
	set_info(buf);

	for (uint64_t i = 0; i < s->n_seq; ++i) {
		mm128_t *r;
		kv_pushp(mm128_t, s->reg, &r);
		r->u64[0] = (uintptr_t)mm_align_seq(t, t->b->opt, t->b->mi, s->seq[i].l_seq, s->seq[i].seq, s->seq[i].rid, t->b->n_occ, t->b->occ, s->lmm, &r->u64[1]);
	}
	return s;
}

static void mm_align_drain(void *arg, void *item)
{
	mm_align_t *b = (mm_align_t*)arg;
	mm_align_step_t *s = (mm_align_step_t*)item;

	set_info("[mm_align_drain] dump alignments to sam records");
	for (uint64_t i = 0; i < s->n_seq; ++i) {
		mm128_t *r = (mm128_t*)s->reg.a[i].u64[0];
		uint32_t n_reg = s->reg.a[i].u64[1];
		if (r == 0) {
			if (b->printer.unmapped) b->printer.unmapped(b, &s->seq[i]);
			continue;
		}
		b->printer.mapped(b, &s->seq[i], n_reg, r);
		for (uint64_t j = 0; j < n_reg; ++j) { lmm_free(s->lmm, (void*)r[j].u64[1]); }
		free(r);
	}
	free(s->reg.a); free(s->base); free((void*)s->seq); lmm_clean(s->lmm); free(s);
	return;
}

static void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == 0) return;
	free(b->resc.a); free(b->coef.a); free(b->intv.a);
	kh_destroy(b->pos); gaba_dp_clean(b->dp);
	free(b);
}

static mm_tbuf_t *mm_tbuf_init(mm_align_t *b)
{
	mm_tbuf_t *t = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	if (t == 0) return 0;
	const uint8_t *lim = (const uint8_t *)0x800000000000;
	t->b = b;
	if ((t->pos = kh_init(0)) == 0) goto _fail;
	if ((t->dp = gaba_dp_init(b->gaba, lim, lim)) == 0) goto _fail;
	return t;
_fail:
	mm_tbuf_destroy(t);
	return 0;
}

static void mm_align_destroy(mm_align_t *b)
{
	assert(b != 0);
	fwrite(b->base, sizeof(uint8_t), b->p - b->base, stdout);
	for (uint64_t i = 0; i < b->opt->n_threads; ++i) mm_tbuf_destroy((mm_tbuf_t*)b->t[i]);
	free(b->base); free(b->occ); free(b->t); gaba_clean(b->gaba); pt_destoroy(b->pt); // ptask_clean(b->pt);
	free(b);
	return;
}

static mm_align_t *mm_align_init(const mm_mapopt_t *opt, const mm_idx_t *mi)
{
	assert(mi != 0);
	assert(opt != 0);

	static const mm_printer_t printer[] = {
		[0] = { .header = mm_print_header_sam, .unmapped = mm_print_unmapped_sam, .mapped = mm_print_mapped_sam },
		[MM_BLAST6>>56] = { .mapped = mm_print_mapped_blast6 },
		[MM_BLASR1>>56] = { .mapped = mm_print_mapped_blasr1 },
		[MM_BLASR4>>56] = { .mapped = mm_print_mapped_blasr4 },
		[MM_PAF>>56] = { .mapped = mm_print_mapped_paf }
	};

	// init output buf and printer
	mm_align_t *b = calloc(1, sizeof(mm_align_t));
	if (b == 0) return 0;
	b->tail = (b->base = b->p = malloc(sizeof(uint8_t) * opt->outbuf_size)) + opt->outbuf_size;
	if (b->base == 0) goto _fail;
	b->printer = printer[opt->flag>>56];
	if (!b->printer.mapped) goto _fail;

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
	if ((b->t = (void**)calloc(b->opt->n_threads, sizeof(mm_tbuf_t*))) == 0) goto _fail;
	for (uint64_t i = 0; i < b->opt->n_threads; ++i) { if ((b->t[i] = (void*)mm_tbuf_init(b)) == 0) goto _fail; }
	if ((b->pt = pt_init(b->opt->n_threads)) == 0) goto _fail;

	// print sam header
	if (b->printer.header) b->printer.header(b);
	return b;
_fail:
	mm_align_destroy(b);
	return 0;
}

static int mm_align_file(mm_align_t *b, bseq_file_t *fp)
{
	assert(b != 0);
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

static void mm_print_help(const mm_mapopt_t *opt)
{
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
					"    $ minialign -X -xava <reads.fa> [<reads.fa> ...] | samsplit <prefix>\n"
					"\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  Global:\n");
	fprintf(stderr, "    -x STR       load preset params {pacbio,ont,ava} [ont]\n");
	fprintf(stderr, "    -t INT       number of threads [%d]\n", opt->n_threads);
	fprintf(stderr, "    -X           switch to all-versus-all alignment mode\n");
	fprintf(stderr, "    -v           show version number [%s]\n", MM_VERSION);
	fprintf(stderr, "  Indexing:\n");
	// fprintf(stderr, "    -c STR,...   treat specified sequences as circular []\n");
	fprintf(stderr, "    -k INT       k-mer size [%d]\n", opt->k);
	fprintf(stderr, "    -w INT       minimizer window size [{-k}*2/3]\n");
	fprintf(stderr, "    -d FILE      dump index to FILE []\n");
	fprintf(stderr, "    -l FILE      load index from FILE [] (overriding -k and -w)\n");
	fprintf(stderr, "  Mapping:\n");
	// fprintf(stderr, "    -f FLOAT,... occurrence thresholds [0.05,0.01,0.001]\n");
	fprintf(stderr, "    -a INT       match award [%d]\n", opt->m);
	fprintf(stderr, "    -b INT       mismatch penalty [%d]\n", opt->x);
	fprintf(stderr, "    -p INT       gap open penalty [%d]\n", opt->gi);
	fprintf(stderr, "    -q INT       gap extension penalty [%d]\n", opt->ge);
	fprintf(stderr, "    -s INT       minimum alignment score [%d]\n", opt->min);
	fprintf(stderr, "    -m INT       minimum alignment score ratio [%1.2f]\n", opt->min_ratio);
	fprintf(stderr, "  Output:\n");
	fprintf(stderr, "    -O STR       output format {sam,blast6,blasr1,blasr4,paf} [sam]\n");
	fprintf(stderr, "    -Q           include quality string\n");
	fprintf(stderr, "    -R STR       read group header line, like \"@RG\\tID:1\" []\n");
	fprintf(stderr, "    -T STR,...   list of optional tags: {RG,AS,} []\n");
	fprintf(stderr, "                   (RG is also inferred from -R)\n");
	fprintf(stderr, "    -U STR,...   tags to be transferred from the input bam file []\n");
	fprintf(stderr, "\n");
	return;
}

static int mm_mapopt_load_preset(mm_mapopt_t *o, const char *arg)
{
	if (strcmp(arg, "pacbio") == 0) {
		o->k = 15; o->w = 10; o->m = 1; o->x = 2; o->gi = 2; o->ge = 1; o->xdrop = 50; o->min = 50; o->min_ratio = 0.3;
	} else if (strcmp(arg, "ont") == 0) {
		o->k = 15; o->w = 10; o->m = 1; o->x = 1; o->gi = 1; o->ge = 1; o->xdrop = 50; o->min = 50; o->min_ratio = 0.3;
	} else if (strcmp(arg, "ava") == 0) {
		o->k = 14; o->w = 10; o->m = 1; o->x = 1; o->gi = 1; o->ge = 1; o->xdrop = 50; o->min = 30; o->min_ratio = 0.3;
		o->flag |= MM_AVA | MM_PAF;
	} else {
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
	return 0;
}

static uint64_t mm_mapopt_parse_tags(const char *p, uint16_v *buf)
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

static uint64_t mm_mapopt_parse_format(const char *arg)
{
	if (strcmp(arg, "sam") == 0) return 0;
	else if (strcmp(arg, "blast6") == 0) return MM_BLAST6;
	else if (strcmp(arg, "blasr1") == 0) return MM_BLASR1;
	else if (strcmp(arg, "blasr4") == 0) return MM_BLASR4;
	else if (strcmp(arg, "paf") == 0) return MM_PAF;
	return 0;
}

static int mm_mapopt_parse(mm_mapopt_t *o, int argc, char *argv[], const char **fnr, const char **fnw, ptr_v *v)
{
	while (optind < argc) {
		int ch;
		if ((ch = getopt(argc, argv, "k:w:f:c:x:B:t:V:d:l:Xs:m:r:a:b:p:q:L:H:I:J:S:E:Y:O:QR:T:U:vh")) < 0) {
			kv_push(void*, *v, argv[optind]); optind++; continue;
		}

		if (ch == 'k') o->k = atoi(optarg);
		else if (ch == 'w') o->w = atoi(optarg);
		else if (ch == 'f') {
			if (mm_mapopt_parse_threshs(o, optarg) && mm_verbose >= 3)
				fprintf(stderr, "[M::%s] Warning: Invalid thresholds: `%s'.\n", __func__, optarg);
		}
		else if (ch == 'x') {
			if (mm_mapopt_load_preset(o, optarg) && mm_verbose >= 3)
				fprintf(stderr, "[M::%s] Warning: Unknown preset tag: `%s'.\n", __func__, optarg);
		}
		else if (ch == 'B') o->b = atoi(optarg);
		else if (ch == 't') o->n_threads = atoi(optarg);
		else if (ch == 'V') mm_verbose = atoi(optarg);
		else if (ch == 'd') *fnw = optarg;
		else if (ch == 'l') *fnr = optarg;
		else if (ch == 'X') o->flag |= MM_AVA;
		else if (ch == 's') o->min = atoi(optarg);
		else if (ch == 'm') o->min_ratio = atof(optarg);
		else if (ch == 'r') {
			if (mm_verbose >= 3)
				fprintf(stderr, "[M::%s] Warning: Minimum length threshold option is deprecated in version 0.4.0 and later, interpreted as score ratio.\n", __func__);
			o->min_ratio = atof(optarg);
		}
		else if (ch == 'a') o->m = atoi(optarg);
		else if (ch == 'b') o->x = atoi(optarg);
		else if (ch == 'p') o->gi = atoi(optarg);
		else if (ch == 'q') o->ge = atoi(optarg);
		else if (ch == 'L') o->llim = atoi(optarg);
		else if (ch == 'H') o->hlim = atoi(optarg);
		else if (ch == 'I') o->blim = atoi(optarg);
		else if (ch == 'J') o->elim = atoi(optarg);
		else if (ch == 'S') o->sidx = atoi(optarg);
		else if (ch == 'E') o->eidx = atoi(optarg);
		else if (ch == 'Y') o->xdrop = atoi(optarg);
		else if (ch == 'O') o->flag &= ~(0xffULL<<56); o->flag |= mm_mapopt_parse_format(optarg);
		else if (ch == 'Q') o->flag |= MM_KEEP_QUAL;
		else if (ch == 'R') mm_mapopt_parse_rg(o, optarg);
		else if (ch == 'T') o->flag |= mm_mapopt_parse_tags(optarg, NULL);
		else if (ch == 'U') mm_mapopt_parse_tags(optarg, &o->tags);
		else if (ch == 'v') { return 1; }
		else if (ch == 'h') { return 2; }
	}
	return 0;
}

int main(int argc, char *argv[])
{
	int ret = 1;
	uint32_t base_rid = 0, base_qid = 0;
	const char *fnr = 0, *fnw = 0;
	bseq_file_t *fp = 0;
	ptr_v v = {0};
	gzFile fpr = 0, fpw = 0;

	set_info("[main] parsing arguments");
	liftrlimit(); posixly_correct();
	mm_realtime0 = realtime();
	mm_mapopt_t *opt = mm_mapopt_init();
	switch (mm_mapopt_parse(opt, argc, argv, &fnr, &fnw, &v)) {
		case 1: puts(MM_VERSION); ret = 0; goto _final;
		case 2: mm_print_help(opt); ret = 0; goto _final;
	}
	if (!fnr && v.n == 0) { mm_print_help(opt); ret = 1; goto _final; }
	if (!fnw && ((fnr && v.n == 0) || (!fnr && v.n == 1 && !(opt->flag&MM_AVA)))) {
		fprintf(stderr, "[M::%s] query-side input redirected to stdin.\n", __func__);
		kv_push(void*, v, "-");
	}
	if (opt->w >= 16) opt->w = (int)(.6666667 * opt->k + .499);
	if (mm_mapopt_check(opt, fprintf, stderr)) {
		ret = 1; goto _final;
	}

	set_info("[main] open index file");
	if (fnr) fpr = gzopen(fnr, "rb");
	if (fnw) fpw = gzopen(fnw, "wb1");
	for (uint64_t i = 0; i < (fpr? 0x7fffffff : (((opt->flag&MM_AVA) || fpw)? v.n : 1)); ++i) {
		uint32_t qid = base_qid;
		mm_idx_t *mi = 0;
		mm_align_t *aln = 0;
		if (fpr) mi = mm_idx_load(fpr);
		else mi = mm_idx_gen(opt, (fp = bseq_open((const char *)v.a[i], base_rid, 0, 0, NULL))), base_rid = bseq_close(fp);
		if (mi == 0) {
			if (fpr && i > 0) break;
			fprintf(stderr, "[M::%s] ERROR: failed to %s `%s'. Please check %s.\n", __func__,
				fpr? "load index file" : "open sequence file",
				fpr? fnr : (const char*)v.a[i], fpr? "file path, format and its version" : "file path and format");
			ret = 1; goto _final;
		}
		if (mm_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built index for %lu target sequence(s)\n", __func__,
			realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->s.n);
		if (fpw) mm_idx_dump(fpw, mi);
		else aln = mm_align_init(opt, mi);
		for (uint64_t j = (!fpr && !(opt->flag&MM_AVA)); j < (fpw? 0 : v.n); ++j)
			mm_align_file(aln, (fp = bseq_open((const char*)v.a[j], qid, (opt->flag&MM_KEEP_QUAL)!=0, opt->tags.n, opt->tags.a))), qid = bseq_close(fp);
		if (!fpw) mm_align_destroy(aln);
		mm_idx_destroy(mi);
	}
	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (uint64_t i = 0; i < (uint64_t)argc; ++i) fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
	ret = 0;
_final:
	free(v.a); mm_mapopt_destroy(opt);
	if (fpr) gzclose(fpr);
	if (fpw) gzclose(fpw);
	return ret;
}

/* end of main.c */
