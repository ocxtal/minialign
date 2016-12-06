#include <unistd.h>
#include <zlib.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/types.h>

#include "kvec.h"
#include "ptask.h"
#include "gaba.h"
#include "sassert.h"

#define MM_VERSION "0.4.2"

#include "arch/arch.h"
#define _VECTOR_ALIAS_PREFIX		v16i8
#include "arch/vector_alias.h"

#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

// #define DEBUG
#include "log.h"

/* minimap.h */

typedef struct bseq_file_s {
	int is_eof;
	gzFile fp;
	kseq_t *ks;
} bseq_file_t;

typedef struct {
	uint32_t l_seq, rid;
	char *name;
	uint8_t *seq;
} bseq_t;
typedef struct { size_t n, m; bseq_t *a; } bseq_v;

typedef struct { uint32_t x[2]; } v2u32_t;
typedef struct { uint32_t x[4]; } v4u32_t;
typedef struct { uint64_t x[2]; } v2u64_t;

typedef union {
	uint64_t u64[2];
	uint32_t u32[4];
} mm128_t;

typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; v2u32_t *a; } v2u32_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; uint32_t *a; } uint32_v;
typedef struct { size_t n, m; uint8_t *a; } uint8_v;
typedef struct { size_t n, m; void **a; } ptr_v;

typedef struct {
	mm128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

typedef struct {
	uint32_t w, k, b;
	uint32_t n;  // number of reference sequences
	mm_idx_bucket_t *B;
	bseq_v s;

	// work
	mm128_v a;
	uint64_v size;
	ptr_v base;
} mm_idx_t;

typedef struct {
	uint64_t batch_size;
	uint32_t sidx, eidx, n_threads, min, k, w, b;
	double min_ratio;
	int32_t m, x, gi, ge, xdrop, llim, hlim, elim, blim;
	uint32_t n_frq;
	float frq[16];
} mm_mapopt_t;

/* end of minimap.h */

#include "ksort.h"
#define sort_key_64(a) ((a))
KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 4)
#define sort_key_64x(a) ((a).x[0])
KRADIX_SORT_INIT(64x, v2u32_t, sort_key_64x, 4)
#define sort_key_128x(a) ((a).u64[0])
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8) 
KSORT_INIT_GENERIC(uint32_t)

#include "khash.h"
#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;
#define pos_hash(a)	(a)
#define pos_eq(a, b)	((a)==(b))
KHASH_INIT(pos, uint64_t, uint64_t, 1, pos_hash, pos_eq)
typedef khash_t(pos) poshash_t;

/* misc.c */

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

/* end of misc.c */

/* bseq.c */

// ascii to 4bit conversion
static unsigned char seq_nt4_table_4bit[32] = {
	0, 1, 0, 2,  0, 0, 0, 4,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  8, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
};

static bseq_file_t *bseq_open(const char *fn)
{
	bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (f == 0) return 0;
	fp = (bseq_file_t*)calloc(1, sizeof(bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

static void bseq_close(bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

static bseq_t *bseq_read(bseq_file_t *fp, uint64_t chunk_size, uint32_t *n_, void **base, uint64_t *size)
{
	uint32_t i, n = *n_;
	kvec_t(uint8_t) mem = {0};
	kvec_t(bseq_t) seq = {0};
	static const uint8_t margin[64] = {0};

	kv_pushm(uint8_t, mem, margin, 64);
	kseq_t *ks = fp->ks;
	while (kseq_read(ks) >= 0) {
		bseq_t *s;
		kv_pushp(bseq_t, seq, &s);

		kv_reserve(uint8_t, mem, mem.n + ks->name.l + ks->seq.l + 2);
		s->l_seq = ks->seq.l;
		s->rid = seq.n + n - 1;
		s->name = (char *)mem.n;
		kv_pushm(uint8_t, mem, ks->name.s, ks->name.l);
		kv_push(uint8_t, mem, '\0');
		s->seq = (uint8_t *)mem.n;
		for (i = 0; i < ks->seq.l; ++i)
			kv_push(uint8_t, mem, seq_nt4_table_4bit[0x1f&ks->seq.s[i]]);
		kv_push(uint8_t, mem, '\0');
		if (mem.n >= chunk_size) break;
	}
	kv_pushm(uint8_t, mem, margin, 64);

	for (i = 0; i < seq.n; i++) {
		seq.a[i].name += (ptrdiff_t)mem.a;
		seq.a[i].seq += (ptrdiff_t)mem.a;
	}
	if (seq.n == 0) free(mem.a), mem.a = 0, fp->is_eof = 1;
	*n_ += seq.n;
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

/* end of sketch.c */

/* index.c */

static mm_idx_t *mm_idx_init(uint32_t w, uint32_t k, uint32_t b)
{
	mm_idx_t *mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w<1? 1 : w, mi->k = k; mi->b = MIN2(k*2, b);
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	return mi;
}

static void mm_idx_destroy(mm_idx_t *mi)
{
	uint32_t i;
	if (mi == 0) return;
	for (i = 0; i < 1<<mi->b; ++i) {
		free(mi->B[i].p);
		free(mi->B[i].a.a);
		kh_destroy(idx, (idxhash_t*)mi->B[i].h);
	}
	free(mi->B);
	for (i = 0; i < mi->base.n; ++i) free(mi->base.a[i]);
	free(mi->base.a); free(mi->size.a);
	free(mi->s.a); free(mi->a.a); free(mi);
}

static const v2u32_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, uint64_t *n)
{
	uint64_t mask = (1<<mi->b) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>mi->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) {
		*n = 1;
		return (const v2u32_t*)&kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return (const v2u32_t*)&b->p[kh_val(h, k)>>32];
	}
}

static uint32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	uint64_t i, n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.) return UINT32_MAX;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = n = 0; i < 1<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
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
	uint32_t n_processed;
	bseq_file_t *fp;
	mm_idx_t *mi;
} mm_idx_pipeline_t;

typedef struct {
    uint32_t n_seq;
	bseq_t *seq;
	void *base;
	uint64_t size;
	mm128_v a;
} mm_idx_step_t;

static void *mm_idx_source(void *arg)
{
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t*)arg;
	mm_idx_step_t *s = (mm_idx_step_t*)calloc(1, sizeof(mm_idx_step_t));
	uint64_t base_n_seq = q->mi->n, size;
	void *base;
	s->seq = bseq_read(q->fp, q->batch_size, &q->mi->n, &base, &size);
	s->n_seq = q->mi->n - base_n_seq;
	if (s->seq == 0) { free(s), s = 0; return 0; }

	kv_push(void*, q->mi->base, base);
	kv_push(uint64_t, q->mi->size, size);
	// kv_pushm(bseq_t, q->mi->s, seq, n_seq);
	return s;
}

static void *mm_idx_worker(void *arg, void *item)
{
	uint32_t i;
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t*)arg;
	mm_idx_step_t *s = (mm_idx_step_t*)item;
	for (i = 0; i < s->n_seq; ++i)
		mm_sketch(s->seq[i].seq, s->seq[i].l_seq, q->mi->w, q->mi->k, s->seq[i].rid, &s->a);
	return s;
}

static void mm_idx_drain(void *arg, void *item)
{
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t*)arg;
	mm_idx_step_t *s = (mm_idx_step_t*)item;
	uint64_t i, mask = (1<<q->mi->b) - 1;
	q->mi->s.n = MAX2(q->mi->s.n, s->seq[s->n_seq-1].rid+1);
	kv_reserve(bseq_t, q->mi->s, q->mi->s.n);
	for (i = 0; i < s->n_seq; ++i)
		q->mi->s.a[s->seq[i].rid] = s->seq[i];
	// kv_pushm(mm128_t, q->mi->a, s->a.a, s->a.n);

	for (i = 0; i < s->a.n; ++i) {
		mm128_v *p = &q->mi->B[s->a.a[i].u64[0]&mask].a;
		kv_push(mm128_t, *p, s->a.a[i]);
	}

	free(s->seq); free(s->a.a);
	return;
}

typedef struct {
	mm_idx_t *mi;
	uint32_t from, to;
} mm_idx_post_t;

static void *mm_idx_post(void *arg, void *item)
{
	uint32_t i, j, start_a, start_p, n, n_keys;
	mm_idx_post_t *q = (mm_idx_post_t*)arg;
	mm_idx_t *mi = q->mi;

	for (i = q->from; i < q->to; ++i) {
		idxhash_t *h;
		mm_idx_bucket_t *b = &mi->B[i];
		if (b->a.n == 0) continue;

		// sort by minimizer
		radix_sort_128x(b->a.a, b->a.a + b->a.n);

		// count and preallocate
		for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
			if (j == b->a.n || b->a.a[j].u64[0] != b->a.a[j-1].u64[0]) {
				++n_keys;
				if (n > 1) b->n += n;
				n = 1;
			} else ++n;
		}
		h = kh_init(idx);
		kh_resize(idx, h, n_keys);
		b->p = (uint64_t*)calloc(b->n, 8);

		// create the hash table
		for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
			if (j == b->a.n || b->a.a[j].u64[0] != b->a.a[j-1].u64[0]) {
				khint_t itr;
				int absent;
				mm128_t *p = &b->a.a[j-1];
				itr = kh_put(idx, h, p->u64[0]>>mi->b<<1, &absent);
				assert(absent && j - start_a == n);
				if (n == 1) {
					kh_key(h, itr) |= 1;
					kh_val(h, itr) = p->u64[1];
				} else {
					int k;
					for (k = 0; k < n; ++k)
						b->p[start_p + k] = b->a.a[start_a + k].u64[1];
					kh_val(h, itr) = (uint64_t)start_p<<32 | n;
					start_p += n;
				}
				start_a = j, n = 1;
			} else ++n;
		}
		b->h = h;
		assert(b->n == start_p);

		// deallocate and clear b->a
		free(b->a.a);
		b->a.n = b->a.m = 0, b->a.a = 0;
	}
	return 0;
}

static mm_idx_t *mm_idx_gen(bseq_file_t *fp, uint32_t w, uint32_t k, uint32_t b, uint64_t batch_size, uint32_t n_threads)
{
	uint64_t i;
	mm_idx_pipeline_t pl = {0}, **p;
	n_threads += (n_threads == 0);
	pl.batch_size = batch_size;
	pl.fp = fp;
	if (pl.fp == 0) return 0;
	pl.mi = mm_idx_init(w, k, b);

	p = (mm_idx_pipeline_t**)calloc(n_threads, sizeof(mm_idx_pipeline_t*));
	for (i = 0; i < n_threads; ++i) p[i] = &pl;
	ptask_t *pt = ptask_init(mm_idx_worker, (void**)p, n_threads, 256);
	ptask_stream(pt, mm_idx_source, &pl, mm_idx_drain, &pl, 256/n_threads);
	ptask_clean(pt);
	free(p);
	
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	mm_idx_post_t *q = (mm_idx_post_t*)calloc(n_threads, sizeof(mm_idx_post_t));
	mm_idx_post_t **qq = (mm_idx_post_t**)calloc(n_threads, sizeof(mm_idx_post_t*));
	for (i = 0; i < n_threads; ++i) {
		q[i].mi = pl.mi;
		q[i].from = (1ULL<<pl.mi->b)*i/n_threads;
		q[i].to = (1ULL<<pl.mi->b)*(i+1)/n_threads;
		qq[i] = &q[i];
	}
	pt = ptask_init(mm_idx_post, (void**)qq, n_threads, 256);
	ptask_parallel(pt, NULL, NULL);
	ptask_clean(pt);
	free(q); free(qq);

	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
	return pl.mi;
}

static mm_idx_t *mm_idx_build(const char *fn, uint32_t w, uint32_t k, uint32_t b, float *frq, uint32_t n_frq, uint32_t n_threads) // a simpler interface
{
	bseq_file_t *fp;
	mm_idx_t *mi;
	fp = bseq_open(fn);
	if (fp == 0) return 0;
	mi = mm_idx_gen(fp, w, k, b, 1<<18, n_threads);
	bseq_close(fp);
	return mi;
}

/*************
 * index I/O *
 *************/

#define MM_IDX_MAGIC "MAI\4"		/* minialign index version 4, differs from minimap index signature */

static void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint32_t x[3];
	uint64_t i, j, size = 0, y[2];
	for (i = 0; i < mi->size.n; ++i) size += mi->size.a[i];
	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b; y[0] = mi->s.n, y[1] = size;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 3, fp);
	fwrite(y, 8, 2, fp);
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&b->n, 4, 1, fp);
		fwrite(b->p, 8, b->n, fp);
		fwrite(&size, 4, 1, fp);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, fp);
		}
	}
	for (i = 0; i < mi->base.n; ++i)
		fwrite(mi->base.a[i], sizeof(char), mi->size.a[i], fp);
	for (i = j = size = 0; i < mi->s.n; ++i) {
		if ((uintptr_t)mi->s.a[i].seq < (uintptr_t)mi->base.a[j]
		|| (uintptr_t)mi->s.a[i].seq >= (uintptr_t)mi->base.a[j] + (ptrdiff_t)mi->size.a[j]) {
			size += mi->size.a[j++];
		}
		mi->s.a[i].name -= (ptrdiff_t)mi->base.a[j], mi->s.a[i].seq -= (ptrdiff_t)mi->base.a[j];
		mi->s.a[i].name += (ptrdiff_t)size, mi->s.a[i].seq += (ptrdiff_t)size;
	}
	fwrite(mi->s.a, sizeof(bseq_t), mi->s.n, fp);
	// restore pointers
	for (i = j = size = 0; i < mi->s.n; ++i) {
		mi->s.a[i].name -= (ptrdiff_t)size, mi->s.a[i].seq -= (ptrdiff_t)size;
		mi->s.a[i].name += (ptrdiff_t)mi->base.a[j], mi->s.a[i].seq += (ptrdiff_t)mi->base.a[j];
		if ((uintptr_t)mi->s.a[i].name > size + mi->size.a[j]) size += mi->size.a[j++];
	}
	return;
}

static mm_idx_t *mm_idx_load(FILE *fp)
{
	char magic[4];
	uint32_t x[3];
	uint64_t i, bsize, y[2];

	mm_idx_t *mi;
	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (fread(x, 4, 3, fp) != 3) return 0;
	if (fread(y, 8, 2, fp) != 2) return 0;
	mi = mm_idx_init(x[0], x[1], x[2]); mi->n = mi->s.n = y[0]; bsize = y[1];
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j, hsize;
		khint_t k;
		idxhash_t *h;
		if (fread(&b->n, 4, 1, fp) != 1) goto _mm_idx_load_fail;
		b->p = (uint64_t*)malloc(b->n * 8);
		if (fread(b->p, 8, b->n, fp) != b->n) goto _mm_idx_load_fail;
		if (fread(&hsize, 4, 1, fp) != 1) goto _mm_idx_load_fail;
		if (hsize == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, hsize);
		for (j = 0; j < hsize; ++j) {
			uint64_t x[2];
			int absent;
			if (fread(x, 8, 2, fp) != 2) goto _mm_idx_load_fail;
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}

	mi->base.n = mi->size.n = 1;
	mi->base.a = malloc(sizeof(void*) * mi->base.n);
	mi->size.a = malloc(sizeof(uint64_t) * mi->size.n);
	mi->base.a[0] = malloc(sizeof(char) * bsize);
	mi->size.a[0] = bsize;
	if (fread(mi->base.a[0], sizeof(char), mi->size.a[0], fp) != mi->size.a[0]) goto _mm_idx_load_fail;
	mi->s.a = malloc(sizeof(bseq_t) * mi->s.n);
	if ((i = fread(mi->s.a, sizeof(bseq_t), mi->s.n, fp)) != mi->s.n) goto _mm_idx_load_fail;
	for (i = 0; i < mi->s.n; ++i) mi->s.a[i].name += (ptrdiff_t)mi->base.a[0], mi->s.a[i].seq += (ptrdiff_t)mi->base.a[0];
	return mi;
_mm_idx_load_fail:
	mm_idx_destroy(mi);
	return 0;
}

/* end of index.c */

/* map.c */

static void mm_mapopt_init(mm_mapopt_t *opt)
{
	opt->k = 15;
	opt->w = 16;
	opt->b = 14;
	opt->sidx = 0;
	opt->eidx = 3;
	opt->hlim = 5000;
	opt->llim = 5000;
	opt->blim = 0;
	opt->elim = 200;
	opt->m = 1;
	opt->x = 2;
	opt->gi = 2;
	opt->ge = 1;
	opt->xdrop = 50;
	opt->min = 200;
	opt->min_ratio = 0.8;
	opt->frq[0] = 0.05, opt->frq[1] = 0.01, opt->frq[2] = 0.001;
	opt->n_frq = 3;
	opt->n_threads = 1;
	opt->batch_size = 10000000;
	return;
}

static const char *mm_mapopt_check(mm_mapopt_t *opt)
{
	uint64_t i;
	if (opt->w >= 16) return "w must be inside [1,16).";
	if (opt->k >= 32) return "k must be inside [1,32).";
	if (opt->sidx >= 16) return "sidx must be inside [0,16).";
	if (opt->eidx >= 16) return "eidx must be inside [0,16).";
	if (opt->hlim < 100 || opt->hlim >= 100000) return "hlim must be inside [100,100000).";
	if (opt->llim < 100 || opt->llim >= 100000) return "llim must be inside [100,100000).";
	if (opt->m < 1 || opt->m > 5) return "Match award must be inside [1,5].";
	if (opt->x < 1 || opt->x > 5) return "Mismatch penalty must be inside [1,5].";
	if (opt->gi < 1 || opt->gi > 5) return "Gap open penalty must be inside [1,5].";
	if (opt->ge < 1 || opt->ge > 5) return "Gap extension penalty must be inside [1,5].";
	if (opt->xdrop < 10 || opt->xdrop > 100) return "Xdrop cutoff must be inside [10,100].";
	if (opt->min > INT32_MAX) return  "Minimum alignment score must be > 0.";
	if (opt->min_ratio < 0.0 || opt->min_ratio > 1.0) return  "Minimum alignment score ratio must be inside [0.0,1.0].";
	if (opt->n_frq >= 16) return "Frequency thresholds must be fewer than 16.";
	for (i = 0; i < opt->n_frq; ++i) if (opt->frq[i] < 0.0 || opt->frq[i] > 1.0 || (i != 0 && opt->frq[i-1] < opt->frq[i])) return "Frequency thresholds must be inside [0.0,1.0] and descending.";
	if (opt->n_threads < 1) return "Thread counts must be > 0.";
	if (opt->batch_size > INT64_MAX) return "Batch size must be > 0.";
	return 0;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int32_t rs, qs, re, qe;
	int32_t rid;
	uint32_t score;
	char *cigar;	// 160907: alignment cigar string
} reg_t;
typedef struct { size_t n, m; reg_t *a; } reg_v;
_static_assert(offsetof(reg_t, re) == offsetof(mm128_t, u32[2]));
_static_assert(offsetof(reg_t, qe) == offsetof(mm128_t, u32[3]));

typedef struct {
	uint64_t buf_size;
	uint32_t n_processed;
	const mm_mapopt_t *opt;
	uint64_t *occ;
	uint32_t n_occ;
	bseq_file_t *fp;
	const mm_idx_t *mi;
	gaba_t *gaba;
	char *base, *tail, *p;
	v2u32_v sort;
} mm_align_pipeline_t;

typedef struct {
	// query sequences
    uint32_t n_seq;
	bseq_t *seq;
	void *base;
	uint64_t size;
	// results
	mm128_v reg;
} mm_align_step_t;

typedef struct { // per-thread buffer
	const mm_align_pipeline_t *p;
	gaba_section_t qf, qr, t;
	mm128_v mini; // query minimizers
	mm128_v resc;
	mm128_v coef; // Hough transform coefficient
	v2u32_v intv; // intervals on sorted coef
	poshash_t *pos;
	gaba_dp_t *dp;	// alignment work
} tbuf_t;

static tbuf_t *mm_tbuf_init(const mm_align_pipeline_t *pl)
{
	tbuf_t *b = (tbuf_t*)calloc(1, sizeof(tbuf_t));
	const uint8_t *lim = (const uint8_t *)0x800000000000;
	b->p = pl; b->pos = kh_init(pos); b->dp = gaba_dp_init(pl->gaba, lim, lim);
	return b;
}

static void mm_tbuf_destroy(tbuf_t *b)
{
	if (b == 0) return;
	free(b->mini.a); free(b->coef.a); free(b->intv.a);
	kh_destroy(pos, b->pos); gaba_dp_clean(b->dp);
	free(b);
}

#define _s(x)		( (x)<0?-1:1)
#define _m(x)		( (((int32_t)(x))>>31)^(x) )
#define _len(x)		( (x)->re+(x)->qe-(x)->rs-(x)->qs )
#define _clip(x)	MAX2(0, MIN2(((uint32_t)(x)), 60))
static void mm_expand(const v2u32_t *r, uint32_t n, int32_t qs, mm128_v *coef)
{
	uint64_t i;
	const int32_t ofs = 0x40000000;
	for (i = 0; i < n; ++i) {	// iterate over all the collected minimizers
		int32_t rs = (int32_t)r[i].x[0], _qs = (rs>>31) ^ qs, _rs = (rs>>31) ^ rs;
		mm128_t *p;
		kv_pushp(mm128_t, *coef, &p);
		p->u32[0] = ofs + _rs - (_qs>>1); p->u32[1] = r[i].x[1];
		p->u32[2] = _rs; p->u32[3] = _qs;
	}
	return;
}

static uint64_t mm_collect(const mm_idx_t *mi, uint32_t l_seq, const uint8_t *seq, uint32_t max_occ, uint32_t resc_occ, mm128_v *mini, mm128_v *coef, mm128_v *resc)
{
	uint64_t i, n, cnt = 0;
	mini->n = coef->n = 0;
	mm_sketch(seq, l_seq, mi->w, mi->k, 0, mini);

	for (i = 0; i < mini->n; ++i) {
		const v2u32_t *r;
		int32_t qs = (int32_t)mini->a[i].u32[2];
		r = mm_idx_get(mi, mini->a[i].u64[0], &n);	// get minimizer at current pos
		if (n > max_occ) continue;	// skip if exceeds repetitive threshold
		cnt += n;
		if (n > resc_occ) {
			mm128_t *q;
			kv_pushp(mm128_t, *resc, &q);
			q->u32[0] = qs; q->u32[1] = n; q->u64[1] = (uint64_t)r;
			continue;
		}
		mm_expand(r, n, qs, coef);
	}
	return(cnt);
}

#if 0
static void mm_chain(mm128_v *coef, uint32_t llim, uint32_t hlim, uint32_t min, v2u32_v *intv)
{
	uint64_t i, j, k, n;
	const int32_t ofs = 0x40000000;
	const uint32_t chained = 0x80000000, mask = 0x7fffffff;
	for (intv->n = i = 0; i < coef->n; i = MIN2(j, k)) {
		uint32_t rid = coef->a[i].u32[1];
		int32_t rs = coef->a[i].u32[2], qs = coef->a[i].u32[3], re = rs, qe = qs;
		int32_t l, h, lub = coef->a[i].u32[0] + llim, hub = rs - (qs<<1) + ofs, hlb = hub - hlim;
		uint32_t len = 0;
		for (j = i+1, k = UINT64_MAX, n = i; j < coef->n && (l = (int32_t)coef->a[j].u32[0]) < lub; j++) {
			if ((int32_t)coef->a[j].u32[2] < 0) continue;
			re = coef->a[j].u32[2]/* & mask*/, qe = coef->a[j].u32[3];	/* tagged seeds are skipped, no need to mask */
			if (rid != coef->a[j].u32[1] || (qs^qe)&chained || (h = ofs + re - (qe<<1)) < hlb || h > hub) { k = MIN2(j, k); continue; }	// out of range, skip
			lub = l + llim; hub = h; hlb = h - hlim;
			coef->a[j].u32[2] |= chained; n = j;
		}
		re = coef->a[n].u32[2] & mask; qe = coef->a[n].u32[3]; qs = _m(qs); qe = _m(qe);
		// len = _s(re-rs)*(re-rs)+_s(qe-qs)*(qe-qs);
		len = re-rs+_s(qe-qs)*(qe-qs);
		if (len < min) continue;
		v2u32_t *p;
		kv_pushp(v2u32_t, *intv, &p);
		p->x[0] = (uint32_t)ofs - len; p->x[1] = i;
	}
	return;
}
#endif

static void mm_chain(uint64_t l_coef, mm128_t *coef, uint32_t llim, uint32_t hlim, uint32_t min, v2u32_v *intv)
{
	uint64_t i, j, k, n;
	const int32_t ofs = 0x40000000;
	const uint32_t chained = 0x80000000, mask = 0x7fffffff;
	for (intv->n = i = 0; i < l_coef; i = MIN2(j, k)) {
		uint32_t rid = coef[i].u32[1];
		int32_t rs = coef[i].u32[2], qs = coef[i].u32[3], re = rs, qe = qs;
		int32_t l, h, lub = coef[i].u32[0] + llim, hub = rs - (qs<<1) + ofs, hlb = hub - hlim;
		uint32_t len = 0;
		for (j = i+1, k = UINT64_MAX, n = i; j < l_coef && (l = (int32_t)coef[j].u32[0]) < lub; j++) {
			if ((int32_t)coef[j].u32[2] < 0) continue;
			re = coef[j].u32[2]/* & mask*/, qe = coef[j].u32[3];	/* tagged seeds are skipped, no need to mask */
			if (rid != coef[j].u32[1] || (qs^qe)&chained || (h = ofs + re - (qe<<1)) < hlb || h > hub) { k = MIN2(j, k); continue; }	// out of range, skip
			lub = l + llim; hub = h; hlb = h - hlim;
			coef[j].u32[2] |= chained; n = j;
		}
		re = coef[n].u32[2] & mask; qe = coef[n].u32[3]; qs = _m(qs); qe = _m(qe);
		// len = _s(re-rs)*(re-rs)+_s(qe-qs)*(qe-qs);
		len = re-rs+_s(qe-qs)*(qe-qs);
		debug("len(%u), re(%d), rs(%d), qe(%d), qs(%d)", len, re, rs, qe, qs);
		if (len < min) continue;
		v2u32_t *p;
		kv_pushp(v2u32_t, *intv, &p);
		p->x[0] = (uint32_t)ofs - len; p->x[1] = i;
	}
	return;
}

static uint64_t mm_short_chain(uint64_t l_coef, mm128_t *coef, uint32_t llim, uint32_t hlim, uint32_t rid, uint32_t eidx)
{
	uint64_t j, k;
	const int32_t ofs = 0x40000000;
	const uint32_t chained = 0x80000000, mask = 0x7fffffff;
	// uint32_t rid = coef->u32[1];
	int32_t rs = coef->u32[2] & mask, qs = coef->u32[3], re = rs, qe = qs;
	int32_t l, h, lub = coef->u32[0] + llim, hub = rs - (qs<<1) + ofs, hlb = hub - hlim;
	uint32_t len = 0;
	debug("short chain, rid(%u), s(%d, %d), lub(%d), hub(%d), hlb(%d)", rid, rs, qs, lub, hub, hlb);
	if (rid != coef->u32[1]) return 1;
	rid = coef->u32[1];
	for (j = 1, k = UINT64_MAX; j < l_coef && (l = (int32_t)coef[j].u32[0]) < lub; j++) {
		re = coef[j].u32[2] & mask, qe = coef[j].u32[3];
		debug("test (%d, %d)", re, qe);
		if (rid != coef[j].u32[1] || (qs^qe)&chained || (h = ofs + re - (qe<<1)) < hlb || h > hub) { k = MIN2(j, k); continue; }
		// hlb = h - hlim; n = j;
		debug("chained, len(%u)", len+1);
		lub = l + llim; hub = h; hlb = h - hlim;	// n = j;
		if (++len > eidx) break;
	}
	return MIN2(j, k);
}

static uint64_t mm_rescue(const bseq_t *ref, uint64_t l_coef, mm128_t *coef, mm128_t *r, uint32_t llim, uint32_t hlim, uint32_t elim, uint32_t blim, uint32_t eidx)
{
	uint64_t i;
	const int32_t ofs = 0x40000000, mask = 0x7fffffff;
	int32_t re = r->u32[2], qe = r->u32[3], l, h;
	int32_t lt = ofs + re - (qe>>1), ht = ofs + re - (qe<<1), lb = lt - blim, le = lt + elim, hb = ht + blim, he = ht - elim, lub = lt + llim;
	debug("tail(%d, %d), lt(%d), ht(%d), lb(%d), le(%d), hb(%d), he(%d), lub(%d), llim(%d), hlim(%d), elim(%d), blim(%d)", re, qe, lt, ht, lb, le, hb, he, lub, llim, hlim, elim, blim);
	debug("coef(%d), lb(%d)", coef[0].u32[0], lb);
	for (i = 0; i < l_coef && coef[i].u32[0] < lb; ++i) { /*debug("coef(%u), lb(%u)", coef[i].u32[0], lb);*/ }
	for (i = 0; i < l_coef && (l = (int32_t)coef[i].u32[0]) < lub; ++i) {
		uint32_t rs = coef[i].u32[2] & mask, qs = coef[i].u32[3];
		debug("test (%d, %d), l(%u), h(%u)", rs, qs, l, ofs + rs - (qs<<1));
		if ((h = ofs + rs - (qs<<1)) > hb || (l < le && h > he && ((l > lt) ^ (h < ht)))) continue;
		debug("prev tail(%d, %d), chain head(%d, %d)", re, qe, rs, qs);
		return mm_short_chain(l_coef - i, coef + i, llim, hlim, ref->rid, eidx);
	}
	return 1;
}

static const gaba_alignment_t *mm_extend(const bseq_t *ref, uint32_t l_coef, mm128_t *coef, uint32_t k, uint32_t min, uint32_t sidx, uint32_t eidx, gaba_dp_t *dp, gaba_section_t *qf, gaba_section_t *qr, gaba_section_t *t, poshash_t *pos, khint_t *pitr)
{
	uint64_t i, j;
	khint_t itr;
	int absent;
	const uint32_t mask = 0x7fffffff;
	const uint8_t *lim = (const uint8_t*)0x800000000000;
	gaba_section_t rf, rr, *qu, *qd, *r, *q;
	gaba_pos_pair_t p = {0};
	rf.id = 2; rf.len = ref->l_seq; rf.base = (const uint8_t*)ref->seq;
	rr.id = 3; rr.len = ref->l_seq; rr.base = gaba_rev((const uint8_t*)ref->seq+ref->l_seq-1, lim);
	gaba_dp_flush(dp, lim, lim);
	gaba_alignment_t *a[3] = {0};
	for (i = sidx, j = 0; i < eidx && i < l_coef; ++i) {
		if (i != 0 && (int32_t)coef[i].u32[2] >= 0) continue;	// skip head
		const gaba_stack_t *stack = gaba_dp_save_stack(dp);
		int32_t rs = coef[i].u32[2] & mask, qs = coef[i].u32[3];
		uint64_t rev = qs<0; qu = rev? qf : qr; qd = rev? qr : qf;
		qs = _m(qs); qs = rev? qf->len-qs+k-1 : qs;
		// upward extension
		gaba_fill_t *f = gaba_dp_fill_root(dp, r = &rr, ref->l_seq-rs-1, q = qu, qu->len-qs-1), *m = f;
		uint32_t flag = GABA_STATUS_TERM;
		do {
			if (f->status & GABA_STATUS_UPDATE_A) flag |= GABA_STATUS_UPDATE_A, r = t;
			if (f->status & GABA_STATUS_UPDATE_B) flag |= GABA_STATUS_UPDATE_B, q = t;
			f = gaba_dp_fill(dp, f, r, q);
			m = (f->max > m->max)? f : m;
		} while(!(flag & f->status));
		// find max
		p = gaba_dp_search_max(dp, m);
		// check duplicate
		// if ((itr = kh_get(pos, pos, ((uint64_t)ref->rid<<32) | p.apos)) != kh_end(pos) && (uint32_t)kh_val(pos, itr) == p.bpos) continue;
		debug("test hash, (%u, %u, %u)", ref->rid, p.apos, p.bpos);
		if ((itr = kh_get(pos, pos, ((uint64_t)ref->rid<<32) | p.apos)) != kh_end(pos)) {
			debug("hit hash, (%lu, %lu)", kh_val(pos, itr)>>32, kh_val(pos, itr) & 0xffffffff);
			// if ((uint32_t)kh_val(pos, itr) == p.bpos) continue;
			if ((uint32_t)kh_val(pos, itr) == p.bpos) return 0;
		}
		// downward extension from max
		gaba_dp_flush_stack(dp, stack);
		m = f = gaba_dp_fill_root(dp, r = &rf, ref->l_seq-p.apos-1, q = qd, qd->len-p.bpos-1);
		flag = GABA_STATUS_TERM;
		do {
			if (f->status & GABA_STATUS_UPDATE_A) flag |= GABA_STATUS_UPDATE_A, r = t;
			if (f->status & GABA_STATUS_UPDATE_B) flag |= GABA_STATUS_UPDATE_B, q = t;
			f = gaba_dp_fill(dp, f, r, q);
			m = (f->max > m->max)? f : m;
		} while(!(flag & f->status));
		// convert alignment to cigar
		a[j++] = gaba_dp_trace(dp, NULL, m, NULL);
		debug("score(%ld), min(%u)", a[j-1]->score, min);
		if (a[j-1]->score < min) continue;	// search again
		break;
	}
	if (j == 0) return 0;
	// collect longest
	for (i = 1; i < j; ++i) if (a[0]->score < a[i]->score) a[0] = a[i];
	// record head
	if ((itr = kh_get(pos, pos, ((uint64_t)ref->rid<<32) | (a[0]->rapos-1))) == kh_end(pos)) {
		itr = kh_put(pos, pos, (((uint64_t)ref->rid<<32) | (a[0]->rapos-1)), &absent);
		debug("record hash, (%u, %u, %u)", ref->rid, a[0]->rapos-1, a[0]->rbpos-1);
		kh_val(pos, itr) = 0xffffffff00000000 | (a[0]->rbpos-1);
	}
	*pitr = itr;

	#if 0
	if ((itr = kh_get(pos, pos, ((uint64_t)ref->rid<<32) | (a[0]->rapos-1))) == kh_end(pos)) {
		itr = kh_put(pos, pos, (((uint64_t)ref->rid<<32) | (a[0]->rapos-1)), &absent);
		debug("record hash, (%u, %u, %u)", ref->rid, a[0]->rapos-1, a[0]0>rbpos-1);
		kh_val(pos, itr) = (reg->n<<32) | (a[0]->rbpos-1);
		kv_pushp(reg_t, *reg, &s);
	} else {
		s = &reg->a[kh_val(pos, itr)>>32];
		if (s->score >= a->score) return 0;
		kh_val(pos, itr) &= 0xffffffff00000000; kh_val(pos, itr) |= rqs;	// replace qpos
		debug("replace hash, (%lu, %u, %u)", kh_val(pos, itr)>>32, rrs, rqs);
		free(s->cigar);	// replace cigar
	}
	#endif
	return a[0];	// one alignment is returned regardless of the score
}

static const reg_t *mm_record(const bseq_t *ref, uint32_t l_seq, const gaba_alignment_t *a, poshash_t *pos, khint_t itr, reg_v *reg)
{
	debug("apos(%u, %u), bpos(%u, %u)", a->sec->apos, a->rapos, a->sec->bpos, a->rbpos);
	reg_t *s;
	if ((int64_t)kh_val(pos, itr) < 0) {
		kh_val(pos, itr) = (reg->n<<32) | (0xffffffff & kh_val(pos, itr));
		kv_pushp(reg_t, *reg, &s);
	} else {
		s = &reg->a[kh_val(pos, itr)>>32];
		free(s->cigar);
	}
	s->rid = ((a->sec->bid&0x01)? 0 : 0xffffffff) ^ ref->rid; s->score = a->score;
	s->rs = ref->l_seq-a->sec->apos-a->sec->alen; s->re = ref->l_seq-a->sec->apos;
	s->qs = l_seq-a->sec->bpos-a->sec->blen; s->qe = l_seq-a->sec->bpos;
	s->cigar = (char *)calloc(a->path->len < 512 ? 1024 : 2*a->path->len, 1);
	gaba_dp_dump_cigar_reverse(s->cigar, 2*a->path->len, a->path->array, 0, a->path->len);
	return s;
}

static reg_t *mm_align(tbuf_t *b, const mm_idx_t *mi, uint32_t l_seq, const uint8_t *seq, uint32_t n_occ, uint64_t *occ, const mm_mapopt_t *opt, uint64_t *n_reg)
{
	uint64_t i, j, k;
	int32_t min = opt->min;
	const int32_t ofs = 0x40000000;
	const uint32_t mask = 0x7fffffff;
	reg_v reg = {0};
	// prepare section info for alignment
	uint8_t tail[96];
	const uint8_t *lim = (const uint8_t*)0x800000000000;
	gaba_section_t qf, qr, mg;
	memset(tail, 0, 96);
	qf.id = 0; qf.len = l_seq; qf.base = (const uint8_t*)seq;
	qr.id = 1; qr.len = l_seq; qr.base = gaba_rev((const uint8_t*)seq+l_seq-1, lim);
	mg.id = 4; mg.len = 32; mg.base = tail+32;
	b->mini.n = b->coef.n = b->resc.n = b->intv.n = 0; kh_clear(pos, b->pos);
	for (i = j = 0; reg.n == 0 && i < n_occ; ++i) {
		if (i == 0) {
			mm_collect(mi, l_seq, seq, occ[n_occ-1], occ[0], &b->mini, &b->coef, &b->resc);
		} else {
			if (i == 1) radix_sort_128x(b->resc.a, b->resc.a + b->resc.n);
			for (k = 0; k < b->coef.n; ++k) b->coef.a[k].u32[2] &= mask;
			for (; j < b->resc.n && b->resc.a[j].u32[1] <= occ[i]; ++j)
				mm_expand((const v2u32_t*)b->resc.a[j].u64[1], b->resc.a[j].u32[1], b->resc.a[j].u32[0], &b->coef);
		}
		radix_sort_128x(b->coef.a, b->coef.a + b->coef.n);
		mm_chain(b->coef.n, b->coef.a, opt->llim, opt->hlim, min>>2, &b->intv);
		radix_sort_64x(b->intv.a, b->intv.a + b->intv.n);
		if (b->intv.a == 0) continue;
		debug("chain count, n(%lu), qlen(%u)", b->intv.n, l_seq);
		for (k = 0; k < b->intv.n && ofs - b->intv.a[k].x[0] > min>>2; ++k) {
			const gaba_alignment_t *a = 0;
			mm128_t *t = &b->coef.a[b->coef.n], r;
			uint64_t p = b->coef.n-b->intv.a[k].x[1], q = p, l = 0, m = ofs - b->intv.a[k].x[0];
			debug("ofs(%u), x(%u), m(%lu)", ofs, b->intv.a[k].x[0], m);
			bseq_t *ref = &mi->s.a[(t-p)->u32[1]];
			do {
				debug("k(%lu), l(%lu), m(%lu)", k, l, m);
				khint_t itr;
				a = mm_extend(ref, q, t-q, mi->k, min, opt->sidx, opt->eidx, b->dp, &qf, &qr, &mg, b->pos, &itr);
				debug("a(%p)", a);
				if (p == q && a == 0) break;	// skip chain if first extension did not result in meaningful alignment
				r.u32[2] = (t-p)->u32[2] & mask, r.u32[3] = (t-p)->u32[3];
				if (a != 0) {
					debug("l(%u), k(%lu), r(%u, %u), q(%u, %u)", l_seq, k, ref->l_seq-a->sec->apos-a->sec->alen, ref->l_seq-a->sec->apos, l_seq-a->sec->bpos-a->sec->blen, l_seq-a->sec->bpos);
					
					debug("record alignment, score(%ld), r(%u, %u), q(%u, %u), root(%u, %u), len(%u, %u)",
						a->score,
						a->sec->apos, a->sec->apos + a->sec->alen,
						a->sec->bpos, a->sec->bpos + a->sec->blen,
						a->rapos, a->rbpos,
						ofs - b->intv.a[k].x[0], a->sec->alen + a->sec->blen);
					
					r.u32[2] = ref->l_seq-a->sec->apos; r.u32[3] = (a->sec->bid&0x01)? l_seq-a->sec->bpos : -a->sec->bpos+mi->k+1;
					debug("r(%u, %d, %d), root(%u, %d, %d)", r.u32[1], r.u32[2], r.u32[3], ref->rid, ref->l_seq-a->rapos, l_seq-a->rbpos);
					min = MAX2(min, opt->min_ratio*a->score);
					l += a->sec->alen+a->sec->blen;
					if (a->score > min) mm_record(ref, l_seq, a, b->pos, itr, &reg);
					debug("l(%lu), m(%lu)", l, m);
				} else {
					debug("l(%u), k(%lu), r(%d, %d), q(%d, %d)", l_seq, k, r.u32[2], r.u32[2], r.u32[3], r.u32[3]);
				}
				q = p;
			} while (l < m && (p -= mm_rescue(ref, p, t-p, &r, opt->llim, opt->hlim, opt->elim, opt->blim, opt->eidx)) < q - 1);
		}
	}
	*n_reg = reg.n;
	return reg.a;
}

static void mm_est_mapq(uint32_t n_reg, const reg_t *reg, v2u32_t *map, uint32_t m, uint32_t x)
{
	uint64_t i;
	const reg_t *preg = &reg[map[0].x[1]], *sreg = (n_reg>1)? &reg[map[1].x[1]] : NULL, *breg = (n_reg>1)? &reg[map[n_reg-1].x[1]] : NULL;
	uint32_t ssc = (n_reg>1)? sreg->score : 0, bsc = (n_reg>1)? breg->score : 0, tsc = 0;
	double elen = (double)_len(preg) / 2.0, pid = 1.0 - (double)(elen * m - preg->score) / (double)(m + x) / elen;
	double ec = 2.0 / (pid * (double)(m + x) - (double)x), ulen = ec * (preg->score - ssc), pe = 1.0 / (ulen * ulen + (double)n_reg);

	for (i = 0; i < n_reg; ++i) {
		const reg_t *r = &reg[i];
		if (r != preg) tsc += r->score - bsc + 1;
	}
	map[0].x[0] = _clip(-10.0 * log10(pe));
	for (i = 1; i < n_reg; ++i) {
		const reg_t *r = &reg[map[i].x[1]];
		map[i].x[0] = _clip(-10.0 * log10(1.0 - pe * (double)(r->score - bsc + 1) / (double)tsc));
	}
	return;
}
#undef _s
#undef _m
#undef _len
#undef _clip

static void *mm_align_source(void *arg)
{
	mm_align_pipeline_t *p = (mm_align_pipeline_t*)arg;
	mm_align_step_t *s = (mm_align_step_t*)calloc(1, sizeof(mm_align_step_t));
	uint64_t base_n_seq = p->n_processed;
	s->seq = bseq_read(p->fp, p->opt->batch_size, &p->n_processed, &s->base, &s->size);
	s->n_seq = p->n_processed - base_n_seq;
	if (s->seq == 0) free(s), s = 0;
	return s;
}

static void *mm_align_worker(void *arg, void *item)
{
	uint32_t i;
	tbuf_t *t = (tbuf_t*)arg;
	mm_align_step_t *s = (mm_align_step_t*)item;
	for (i = 0; i < s->n_seq; ++i) {
		mm128_t *r;
		kv_pushp(mm128_t, s->reg, &r);
		r->u64[0] = (uintptr_t)mm_align(t, t->p->mi, s->seq[i].l_seq, s->seq[i].seq, t->p->n_occ, t->p->occ, t->p->opt, &r->u64[1]);
	}
	return s;
}

static void mm_align_drain(void *arg, void *item)
{
	uint32_t i, j;
	int32_t k;
	mm_align_pipeline_t *p = (mm_align_pipeline_t*)arg;
	mm_align_step_t *s = (mm_align_step_t*)item;

	#define _put(_p, _c) { \
		*(_p)->p++ = (_c); \
		if ((_p)->p == (_p)->tail) { \
			fwrite((_p)->base, sizeof(char), (_p)->buf_size, stdout); \
			(_p)->p = (_p)->base; \
		} \
	}
	#define _putn(_p, _n) { \
		char _b[16] = {0}; \
		int _m = (_n), _i = 0, _j; \
		while (_m) _b[_i++] = _m % 10, _m /= 10; \
		for (_j = _i - (_i!=0); _j >= 0; _j--) _put(_p, _b[_j]+'0'); \
	}

	for (i = 0; i < s->n_seq; ++i) {
		reg_t *r = (reg_t*)s->reg.a[i].u64[0];
		uint32_t n_reg = s->reg.a[i].u64[1];
		bseq_t *t = &s->seq[i];
		char *q;
		if (r == 0) {
			const char *c1 = "\t4\t*\t0\t0\t*\t*\t0\t0\t", *c2 = "\t*\n";
			for (q = t->name; *q; q++) _put(p, *q);
			for (k = 0; k < strlen(c1); k++) _put(p, c1[k]);
			for (k = 0; k < t->l_seq; k++) _put(p, "NACMGRSVTWYHKDBN"[(uint8_t)t->seq[k]]);
			for (k = 0; k < strlen(c2); k++) _put(p, c2[k]);
			continue;
		}
		p->sort.n = 0;
		for (j = 0; j < n_reg; ++j) {
			v2u32_t *s;
			kv_pushp(v2u32_t, p->sort, &s);
			s->x[0] = -r[j].score; s->x[1] = j;
		}
		radix_sort_64x(p->sort.a, p->sort.a + p->sort.n);
		mm_est_mapq(n_reg, r, p->sort.a, p->opt->m, p->opt->x);

		for (j = 0; j < n_reg && r[p->sort.a[j].x[1]].score > p->opt->min_ratio*r[p->sort.a[0].x[1]].score; ++j) {
			const char *c1 = "\t*\t0\t0\t", *c2 = "\t*\tRG:Z:1\n";
			int qs = j==0? 0 : r[p->sort.a[j].x[1]].qs, qe = j==0? t->l_seq : r[p->sort.a[j].x[1]].qe;
			for (q = t->name; *q; q++) _put(p, *q);
			_put(p, '\t'); _putn(p, (r[p->sort.a[j].x[1]].rid<0?16:0) | (j==0?0:256)); _put(p, '\t');
			for (q = p->mi->s.a[(r[p->sort.a[j].x[1]].rid>>31) ^ r[p->sort.a[j].x[1]].rid].name; *q; q++) _put(p, *q);
			_put(p, '\t'); _putn(p, r[p->sort.a[j].x[1]].rs+1);
			_put(p, '\t'); _putn(p, p->sort.a[j].x[0]); _put(p, '\t');
			if (r[p->sort.a[j].x[1]].qs) { _putn(p, r[p->sort.a[j].x[1]].qs); _put(p, j==0? 'S' : 'H'); }
			for (q = r[p->sort.a[j].x[1]].cigar; *q; q++) _put(p, *q);
			debug("head(%u), tail(%u)", r[p->sort.a[j].x[1]].qs, t->l_seq-r[p->sort.a[j].x[1]].qe);
			if (t->l_seq-r[p->sort.a[j].x[1]].qe) { _putn(p, t->l_seq-r[p->sort.a[j].x[1]].qe); _put(p, j==0? 'S' : 'H'); }
			for (k = 0; k < strlen(c1); k++) _put(p, c1[k]);
			if (r[p->sort.a[j].x[1]].rid < 0) { for (k = t->l_seq-qs; k > t->l_seq-qe; k--) _put(p, "NTGKCYSBAWRDMHVN"[(uint8_t)t->seq[k-1]]); }
			else { for (k = qs; k < qe; k++) _put(p, "NACMGRSVTWYHKDBN"[(uint8_t)t->seq[k]]); }
			for (k = 0; k < strlen(c2); k++) _put(p, c2[k]);
			free(r[p->sort.a[j].x[1]].cigar);
		}
		for (; j < n_reg; ++j) free(r[p->sort.a[j].x[1]].cigar);
		free(r);
	}
	free(s->reg.a); free(s->base); free(s->seq); free(s);
	return;
}

static void mm_print_header(const mm_idx_t *idx)
{
	uint32_t i;
	printf("@HD\tVN:1.0\tSO:unsorted\n");
	for (i = 0; i < idx->n; ++i) printf("@SQ\tSN:%s\tLN:%d\n", idx->s.a[i].name, idx->s.a[i].l_seq);
	printf("@RG\tID:1\n");
	printf("@PG\tID:6\tPN:minialign\n");
	return;
}

static int mm_align_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt)
{
	uint32_t i;
	mm_align_pipeline_t pl = {0};
	tbuf_t **t;
	pl.fp = bseq_open(fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.mi = idx;

	// calc occ
	pl.occ = calloc(pl.opt->n_frq, sizeof(uint64_t));
	pl.n_occ = pl.opt->n_frq;
	for (i = 0; i < pl.n_occ; ++i) pl.occ[i] = mm_idx_cal_max_occ(idx, pl.opt->frq[i]);

	// init output buf
	pl.buf_size = 512 * 1024;
	pl.base = pl.p = malloc(sizeof(char) * pl.buf_size);
	pl.tail = pl.base + pl.buf_size;

	// initialize alignment context
	int m = opt->m, x = -opt->x, gi = opt->gi, ge = opt->ge;
	struct gaba_score_s sc = {{{m,x,x,x},{x,m,x,x},{x,x,m,x},{x,x,x,m}},gi,ge,gi,ge};
	struct gaba_params_s p = {0,0,0,opt->xdrop,&sc};
	pl.gaba = gaba_init(&p);

	t = (tbuf_t**)calloc(pl.opt->n_threads, sizeof(tbuf_t*));
	for (i = 0; i < pl.opt->n_threads; ++i) t[i] = mm_tbuf_init(&pl);
	ptask_t *pt = ptask_init(mm_align_worker, (void**)t, pl.opt->n_threads, 1024);
	ptask_stream(pt, mm_align_source, &pl, mm_align_drain, &pl, 1024/pl.opt->n_threads);
	ptask_clean(pt);
	for (i = 0; i < pl.opt->n_threads; ++i) mm_tbuf_destroy(t[i]);
	free(t);

	fwrite(pl.base, sizeof(char), pl.p - pl.base, stdout);

	bseq_close(pl.fp);
	gaba_clean(pl.gaba);
	free(pl.base); free(pl.occ);
	return 0;
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

int main(int argc, char *argv[])
{
	int i, ch, is_idx = 0;
	mm_mapopt_t opt;
	const char *err;
	bseq_file_t *fp = 0;
	char *fnw = 0;
	FILE *fpr = 0, *fpw = 0;

	liftrlimit();
	mm_realtime0 = realtime();
	mm_mapopt_init(&opt);

	while ((ch = getopt(argc, argv, "k:w:f:B:t:V:d:ls:m:r:a:b:p:q:L:H:I:J:S:E:X:v")) >= 0) {
		if (ch == 'k') opt.k = atoi(optarg);
		else if (ch == 'w') opt.w = atoi(optarg);
		else if (ch == 'f') {
			const char *p = optarg; opt.n_frq = 0;
			while (*p) {
				opt.frq[opt.n_frq++] = atof(p);
				while (*p && *p != ',') p++;
				if (!*p || opt.n_frq >= 16) break;
				p++;
			}
		}
		else if (ch == 'B') opt.b = atoi(optarg);
		else if (ch == 't') opt.n_threads = atoi(optarg);
		else if (ch == 'V') mm_verbose = atoi(optarg);
		else if (ch == 'd') fnw = optarg;
		else if (ch == 'l') is_idx = 1;
		else if (ch == 's') opt.min = atoi(optarg);
		else if (ch == 'm') opt.min_ratio = atof(optarg);
		else if (ch == 'r') {
			if (mm_verbose >= 3) fprintf(stderr, "Minimum length threshold option is deprecated in version 0.4.0. It is interpreted as score ratio.).\n");
			opt.min_ratio = atof(optarg);
		}
		else if (ch == 'a') opt.m = atoi(optarg);
		else if (ch == 'b') opt.x = atoi(optarg);
		else if (ch == 'p') opt.gi = atoi(optarg);
		else if (ch == 'q') opt.ge = atoi(optarg);
		else if (ch == 'L') opt.llim = atoi(optarg);
		else if (ch == 'H') opt.hlim = atoi(optarg);
		else if (ch == 'I') opt.blim = atoi(optarg);
		else if (ch == 'J') opt.elim = atoi(optarg);
		else if (ch == 'S') opt.sidx = atoi(optarg);
		else if (ch == 'E') opt.eidx = atoi(optarg);
		else if (ch == 'X') opt.xdrop = atoi(optarg);
		else if (ch == 'v') { puts(MM_VERSION); return 0; }
	}
	if (opt.w >= 16) opt.w = (int)(.6666667 * opt.k + .499);

	if (argc == optind) {
		fprintf(stderr, "\n"
						"  minialign - fast aligner for PacBio and Nanopore long reads\n"
						"\n"
						"minialign is a fast long-read (nucleotide sequence) alignment tool built on\n"
						"the top of minimap long-read overlapper adopting libgaba SIMD-parallelized\n"
						"Smith-Waterman extension algorithm.\n"
						"\n");
		fprintf(stderr, "Usage:\n"
						"  without index:\n"
						"    $ minialign [options] <ref.fa> <reads.fa> > mapping.sam\n"
						"\n"
						"  with prebuilt index (saves ~1min for human genome per run):\n"
						"    $ minialign [options] -d <index.mai> <ref.fa>\n"
						"    $ minialign -l <index.mai> <reads.fa> > mapping.sam\n"
						"\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Indexing:\n");
		fprintf(stderr, "    -k INT       k-mer size [%d]\n", opt.k);
		fprintf(stderr, "    -w INT       minimizer window size [{-k}*2/3]\n");
		fprintf(stderr, "    -d FILE      dump index to FILE []\n");
		fprintf(stderr, "    -l FILE      load index from FILE [] (overriding -k, -w, and -f)\n");
		fprintf(stderr, "  Mapping:\n");
		fprintf(stderr, "    -f FLOAT,... occurrence thresholds [0.05,0.01,0.001]\n");
		fprintf(stderr, "    -a INT       match award [%d]\n", opt.m);
		fprintf(stderr, "    -b INT       mismatch penalty [%d]\n", opt.x);
		fprintf(stderr, "    -p INT       gap open penalty [%d]\n", opt.gi);
		fprintf(stderr, "    -q INT       gap extension penalty [%d]\n", opt.ge);
		fprintf(stderr, "    -s INT       minimum alignment score [%d]\n", opt.min);
		fprintf(stderr, "    -m INT       minimum alignment score ratio [%f]\n", opt.min_ratio);
		fprintf(stderr, "  Misc:\n");
		fprintf(stderr, "    -t INT       number of threads [%d]\n", opt.n_threads);
		fprintf(stderr, "    -v           show version number\n");
		fprintf(stderr, "\n");
		return 1;
	}

	if ((err = mm_mapopt_check(&opt))) {
		fprintf(stderr, "[M::%s::] ERROR: %s\n", __func__, err);
		return -1;
	}

	if (is_idx) fpr = fopen(argv[optind], "rb");
	else fp = bseq_open(argv[optind]);
	if (fnw) fpw = fopen(fnw, "wb");
	while (fpr || fp) {
		mm_idx_t *mi = 0;
		if (fpr) mi = mm_idx_load(fpr);
		else if (!bseq_eof(fp))
			mi = mm_idx_gen(fp, opt.w, opt.k, opt.b, opt.batch_size, opt.n_threads);
		if (mi == 0) break;
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built index for %d target sequence(s)\n",
					__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n);
		if (fpw) mm_idx_dump(fpw, mi);
		if (argc - optind > 1) mm_print_header(mi);
		for (i = optind + 1; i < argc; ++i)
			mm_align_file(mi, argv[i], &opt);
		mm_idx_destroy(mi);
	}
	if (fpw) fclose(fpw);
	if (fpr) fclose(fpr);
	if (fp)  bseq_close(fp);

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
	return 0;
}

/* end of main.c */
