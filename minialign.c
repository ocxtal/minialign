#include <unistd.h>
#include <zlib.h>
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

#define MM_VERSION "0.2-unstable"

#include "arch/arch.h"
#define _VECTOR_ALIAS_PREFIX		v16i8
#include "arch/vector_alias.h"

#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

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
	float f;
	uint32_t n;  // number of reference sequences
	mm_idx_bucket_t *B;
	uint64_t max_occ;
	bseq_v s;

	// work
	mm128_v a;
	uint64_v size;
	ptr_v base;
} mm_idx_t;

typedef struct {
	uint32_t min_len;
	int32_t m, x, gi, ge, xdrop, min_score, ofs_llim, ofs_hlim;
	double min_ratio;
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

/* misc.c */

int mm_verbose = 3;
double mm_realtime0;

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

/* end of misc.c */

/* bseq.c */

// ascii to 4bit conversion
static uint8_t seq_nt4_table_4bit[256] = {
	0, 1, 2, 0,  4, 0, 0, 0,  8, 0, 0, 0,  0, 0, 0, 0, 
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	0, 1, 0, 2,  0, 0, 0, 4,  0, 0, 0, 0,  0, 0, 0, 0, 
	0, 0, 0, 0,  8, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	0, 1, 0, 2,  0, 0, 0, 4,  0, 0, 0, 0,  0, 0, 0, 0, 
	0, 0, 0, 0,  8, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
};

bseq_file_t *bseq_open(const char *fn)
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

void bseq_close(bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

bseq_t *bseq_read(bseq_file_t *fp, uint64_t chunk_size, uint32_t *n_, void **base, uint64_t *size)
{
	uint32_t i, n = *n_;
	kseq_t *ks = fp->ks;
	kvec_t(uint8_t) buf = {0};
	kvec_t(bseq_t) seq = {0};
	static const uint8_t margin[64] = {0};

	kv_pushm(uint8_t, buf, margin, 64);
	while (kseq_read(ks) >= 0) {
		bseq_t *s;
		kv_pushp(bseq_t, seq, &s);

		kv_reserve(uint8_t, buf, buf.n + ks->name.l + ks->seq.l + 2);
		s->l_seq = ks->seq.l;
		s->rid = seq.n + n - 1;
		s->name = (char *)buf.n;
		kv_pushm(uint8_t, buf, ks->name.s, ks->name.l);
		kv_push(uint8_t, buf, '\0');
		s->seq = (uint8_t *)buf.n;
		for (i = 0; i < ks->seq.l; ++i)
			kv_push(uint8_t, buf, seq_nt4_table_4bit[(uint8_t)ks->seq.s[i]]);
		kv_push(uint8_t, buf, '\0');
		if (buf.n >= chunk_size) break;
	}
	kv_pushm(uint8_t, buf, margin, 64);

	for (i = 0; i < seq.n; i++) {
		seq.a[i].name += (ptrdiff_t)buf.a;
		seq.a[i].seq += (ptrdiff_t)buf.a;
	}
	if (seq.n == 0) free(buf.a), buf.a = 0, fp->is_eof = 1;
	*n_ += seq.n;
	*base = (void*)buf.a;
	*size = buf.n;
	return seq.a;
}

int bseq_eof(bseq_file_t *fp)
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
void mm_sketch(const uint8_t *seq4, uint32_t len, uint32_t w, uint32_t k, uint32_t rid, mm128_v *p)
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

mm_idx_t *mm_idx_init(uint32_t w, uint32_t k, uint32_t b, float f)
{
	mm_idx_t *mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w<1? 1 : w, mi->k = k; mi->b = MIN2(k*2, b);
	mi->f = f<0.0? 0.0 : f;
	mi->max_occ = UINT32_MAX;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	return mi;
}

void mm_idx_destroy(mm_idx_t *mi)
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

const v2u32_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, uint64_t *n)
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

uint32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	int i;
	size_t n = 0;
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

mm_idx_t *mm_idx_gen(bseq_file_t *fp, uint32_t w, uint32_t k, uint32_t b, float f, uint64_t batch_size, uint32_t n_threads)
{
	uint64_t i;
	mm_idx_pipeline_t pl = {0}, **p;
	n_threads += (n_threads == 0);
	pl.batch_size = batch_size;
	pl.fp = fp;
	if (pl.fp == 0) return 0;
	pl.mi = mm_idx_init(w, k, b, f);

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

	pl.mi->max_occ = mm_idx_cal_max_occ(pl.mi, pl.mi->f);
	return pl.mi;
}

mm_idx_t *mm_idx_build(const char *fn, uint32_t w, uint32_t k, uint32_t b, float f, uint32_t n_threads) // a simpler interface
{
	bseq_file_t *fp;
	mm_idx_t *mi;
	fp = bseq_open(fn);
	if (fp == 0) return 0;
	mi = mm_idx_gen(fp, w, k, b, f, 1<<18, n_threads);
	bseq_close(fp);
	return mi;
}

/*************
 * index I/O *
 *************/

#define MM_IDX_MAGIC "MAI\2"		/* minialign index version 2, differs from minimap index signature */

void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	union { uint32_t i; float f; } x[6];
	uint64_t i, j, size = 0, y[2];
	for (i = 0; i < mi->size.n; ++i) size += mi->size.a[i];
	x[0].i = mi->w, x[1].i = mi->k, x[2].i = mi->b, x[3].f = mi->f, x[4].i = mi->n, x[5].i = mi->max_occ;
	y[0] = mi->s.n, y[1] = size;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 6, fp);
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
}

mm_idx_t *mm_idx_load(FILE *fp)
{
	char magic[4];
	union { uint32_t i; float f; } x[6];
	uint64_t i, bsize, y[2];

	mm_idx_t *mi;
	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (fread(x, 4, 6, fp) != 6) return 0;
	if (fread(y, 8, 2, fp) != 2) return 0;
	mi = mm_idx_init(x[0].i, x[1].i, x[2].i, x[3].f);
	mi->n = x[4].i, mi->max_occ = x[5].i;
	mi->s.n = y[0]; bsize = y[1];

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
	if (fread(mi->s.a, sizeof(bseq_t), mi->s.n, fp) != mi->s.n) goto _mm_idx_load_fail;
	for (i = 0; i < mi->s.n; ++i)
		mi->s.a[i].name += (ptrdiff_t)mi->base.a[0], mi->s.a[i].seq += (ptrdiff_t)mi->base.a[0];
	return mi;
_mm_idx_load_fail:
	mm_idx_destroy(mi);
	return 0;
}

/* end of index.c */

/* map.c */

void mm_mapopt_init(mm_mapopt_t *opt)
{
	opt->ofs_hlim = 1500;
	opt->ofs_llim = 1500;
	opt->min_len = 2000;
	opt->m = 1;
	opt->x = 1;
	opt->gi = 1;
	opt->ge = 1;
	opt->xdrop = 50;
	opt->min_score = 100;
	opt->min_ratio = 0.8;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	uint32_t rid;
	uint32_t flag;
	int32_t qs, qe, rs, re;
	char *cigar;	// 160907: alignment cigar string
} reg_t;
typedef struct { size_t n, m; reg_t *a; } reg_v;

typedef struct {
	uint64_t batch_size, buf_size;
	uint32_t n_processed, n_threads;
	const mm_mapopt_t *opt;
	bseq_file_t *fp;
	const mm_idx_t *mi;
	gaba_t *gaba;
	char *base, *tail, *p;
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
	mm128_v mini; // query minimizers
	mm128_v coef; // Hough transform coefficient
	v2u32_v intv; // intervals on sorted coef
	gaba_dp_t *dp;	// alignment work
} tbuf_t;

tbuf_t *mm_tbuf_init(const mm_align_pipeline_t *pl)
{
	tbuf_t *b = (tbuf_t*)calloc(1, sizeof(tbuf_t));
	const uint8_t *lim = (const uint8_t *)0x800000000000;
	b->p = pl; b->dp = gaba_dp_init(pl->gaba, lim, lim);
	return b;
}

void mm_tbuf_destroy(tbuf_t *b)
{
	if (b == 0) return;
	free(b->mini.a); free(b->coef.a); free(b->intv.a);
	gaba_dp_clean(b->dp);
	free(b);
}

reg_t *mm_align(tbuf_t *b, const mm_idx_t *mi, uint32_t l_seq, const uint8_t *seq, const mm_mapopt_t *opt, uint64_t *n_reg)
{
	uint64_t i, j, k, l, n;
	uint64_t const chained = 0x8000000000000000;
	int32_t min_score;
	int32_t const ofs = 0x40000000;
	reg_v reg = {0};
	#define _s(x)	( (x)<0?-1:1)
	#define _m(x)	( (((int32_t)(x))>>31)^(x) )

	// collect minimizers
	b->mini.n = b->coef.n = 0;
	mm_sketch(seq, l_seq, mi->w, mi->k, 0, &b->mini);
	for (j = 0; j < b->mini.n; ++j) {
		const v2u32_t *r;
		int32_t qs = (int32_t)b->mini.a[j].u32[2];
		r = mm_idx_get(mi, b->mini.a[j].u64[0], &n);	// get minimizer at current pos
		if (n > mi->max_occ) continue;	// skip if exceeds repetitive threshold
		for (k = 0; k < n; ++k) {	// iterate over all the collected minimizers
			int32_t rs = (int32_t)r[k].x[0];
			mm128_t *p;
			kv_pushp(mm128_t, b->coef, &p);
			p->u32[0] = ofs + ((rs>>31)^(rs - (qs>>1))); p->u32[1] = r[k].x[1];
			p->u32[2] = rs; p->u32[3] = qs;
		}
	}
	if (b->coef.n == 0) return 0;
	radix_sort_128x(b->coef.a, b->coef.a + b->coef.n);

	// chain (modified lis algorithm)
	for (b->intv.n = i = 0; i < b->coef.n; i = MIN2(j, k)) {
		uint32_t rid = b->coef.a[i].u32[1];
		int32_t rs = b->coef.a[i].u32[2], qs = b->coef.a[i].u32[3], re, qe;
		int32_t lub = b->coef.a[i].u32[0] + opt->ofs_llim, hub = ((rs>>31)^(rs - (qs<<1))) + ofs, hlb = hub - opt->ofs_hlim;
		uint32_t cnt = 0, len = 0;
		for (j = i+1, k = UINT64_MAX, l = i; j < b->coef.n && (int32_t)b->coef.a[j].u32[0] < lub; j++) {
			if ((int64_t)b->coef.a[j].u64[0] < 0) continue;
			re = b->coef.a[j].u32[2], qe = b->coef.a[j].u32[3];
			int32_t h = ofs + ((re>>31)^(re - (qe<<1)));
			if (rid != b->coef.a[j].u32[1] || (h > hub && h < hlb)) { k = MIN2(j, k); continue; }	// out of range, skip
			lub = b->coef.a[j].u32[0] + opt->ofs_llim; hlb = h - opt->ofs_hlim;
			b->coef.a[j].u64[0] |= chained; l = j; cnt++;
		}
		re = b->coef.a[l].u32[2], qe = b->coef.a[l].u32[3];
		rs = _m(rs); re = _m(re); qs = _m(qs); qe = _m(qe);
		len = _s(re-rs)*(re-rs)+_s(qe-qs)*(qe-qs);
		if (len < opt->min_len) continue;
		v2u32_t *p;
		kv_pushp(v2u32_t, b->intv, &p);
		p->x[0] = (uint32_t)ofs - len; p->x[1] = i;
	}
	if (b->intv.n == 0) return 0;
	radix_sort_64x(b->intv.a, b->intv.a + b->intv.n);

	// prepare section info for alignment
	uint8_t tail[96];
	const uint8_t *lim = (const uint8_t*)0x800000000000;
	gaba_section_t qf, qr, rf, rr, t, *qu, *qd, *r, *q;
	memset(tail, 0, 96);
	qf.id = 0; qf.len = l_seq; qf.base = (const uint8_t*)seq;
	qr.id = 1; qr.len = l_seq; qr.base = gaba_rev((const uint8_t*)seq+l_seq-1, lim);
	t.id = 4; t.len = 32; t.base = tail+32;

	// extend
	min_score = opt->min_score;
	for (i = 0; i < b->intv.n && b->intv.a[i].x[0] < ofs - opt->min_len; ++i) {
		char *cig;
		int64_t rev = 0;
		uint32_t rid = b->coef.a[b->intv.a[i].x[1]].u32[1] & 0x7fffffff;
		bseq_t *ref = &mi->s.a[rid];
		rf.id = 2; rf.len = ref->l_seq; rf.base = (const uint8_t*)ref->seq;
		rr.id = 3; rr.len = ref->l_seq; rr.base = gaba_rev((const uint8_t*)ref->seq+ref->l_seq-1, lim);

		gaba_dp_flush(b->dp, lim, lim);
		gaba_alignment_t *a[3] = {0};
		for (j = b->intv.a[i].x[1], k = l = 0; k < 3 && j+k < b->coef.n; ++k) {
			if (k != 0 && (int64_t)b->coef.a[j+k].u64[0] >= 0) continue;	// skip head
			int32_t rs = b->coef.a[j+k].u32[2], qs = b->coef.a[j+k].u32[3];
			rev = (qs^rs)<0; qu = rev? &qf : &qr; qd = rev? &qr : &qf;
			rs = _m(rs); qs = _m(qs); qs = rev? l_seq-qs+mi->k : qs;
			// upward extension
			gaba_fill_t *f = gaba_dp_fill_root(b->dp, r = &rr, rf.len-rs, q = qu, qf.len-qs), *m = f;
			uint32_t mask = GABA_STATUS_TERM;
			do {
				if (f->status & GABA_STATUS_UPDATE_A) mask |= GABA_STATUS_UPDATE_A, r = &t;
				if (f->status & GABA_STATUS_UPDATE_B) mask |= GABA_STATUS_UPDATE_B, q = &t;
				f = gaba_dp_fill(b->dp, f, r, q);
				m = (f->max > m->max)? f : m;
			} while(!(mask & f->status));
			// find max
			gaba_pos_pair_t p = gaba_dp_search_max(b->dp, m);
			// downward extension from max
			m = f = gaba_dp_fill_root(b->dp, r = &rf, rf.len-p.apos, q = qd, qf.len-p.bpos);
			mask = GABA_STATUS_TERM;
			do {
				if (f->status & GABA_STATUS_UPDATE_A) mask |= GABA_STATUS_UPDATE_A, r = &t;
				if (f->status & GABA_STATUS_UPDATE_B) mask |= GABA_STATUS_UPDATE_B, q = &t;
				f = gaba_dp_fill(b->dp, f, r, q);
				m = (f->max > m->max)? f : m;
			} while(!(mask & f->status));
			// convert alignment to cigar
			a[l++] = gaba_dp_trace(b->dp, NULL, m, NULL);
			if (a[l-1]->score < min_score || gaba_plen(a[l-1]->sec) < ofs - b->intv.a[i].x[0]) continue;
			break;
		}
		if (l == 0) continue;
		// collect longest
		for (j = 1; j < l; ++j) if (a[0]->score < a[j]->score) a[0] = a[j];
		if (a[0]->score < min_score) continue;
		min_score = MAX2(min_score, opt->min_ratio*a[0]->score);

		reg_t *r;
		kv_pushp(reg_t, reg, &r);
		r->rid = rid; r->flag = (rev? 0x10 : 0) | (reg.n==1? 0 : 0x100);
		r->rs = rf.len-a[0]->sec->apos-a[0]->sec->alen; r->re = rf.len-a[0]->sec->apos;
		r->qs = qf.len-a[0]->sec->bpos-a[0]->sec->blen; r->qe = qf.len-a[0]->sec->bpos;
		if(r->qs < 0) printf("%c\n", *((char*)0));
		r->cigar = cig = (char *)calloc(a[0]->path->len < 512 ? 1024 : 2*a[0]->path->len, 1);
		if (r->qs) cig += sprintf(cig, "%d%c", r->qs, (reg.n==1)? 'S' : 'H');
		cig += gaba_dp_dump_cigar_reverse(cig, 2*a[0]->path->len, a[0]->path->array, 0, a[0]->path->len);
		if (qf.len-r->qe) cig += sprintf(cig, "%d%c", qf.len-r->qe, (reg.n==1)? 'S' : 'H');
	}
	*n_reg = reg.n;
	return reg.a;
}

static void *mm_align_source(void *arg)
{
	mm_align_pipeline_t *p = (mm_align_pipeline_t*)arg;
	mm_align_step_t *s = (mm_align_step_t*)calloc(1, sizeof(mm_align_step_t));
	uint64_t base_n_seq = p->n_processed;
	s->seq = bseq_read(p->fp, p->batch_size, &p->n_processed, &s->base, &s->size);
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
		r->u64[0] = (uintptr_t)mm_align(t, t->p->mi, s->seq[i].l_seq, s->seq[i].seq, t->p->opt, &r->u64[1]);
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
		for (j = 0; j < s->reg.a[i].u64[1]; ++j) {
			const char *c1 = "\t255\t", *c2 = "\t*\t0\t0\t", *c3 = "\t*\tRG:Z:1\n";
			int qs = j==0? 0 : r[j].qs, qe = j==0? t->l_seq : r[j].qe;
			for (q = t->name; *q; q++) _put(p, *q);
			_put(p, '\t'); _putn(p, r[j].flag); _put(p, '\t');
			for (q = p->mi->s.a[r[j].rid].name; *q; q++) _put(p, *q);
			_put(p, '\t'); _putn(p, r[j].rs+1); 
			for (k = 0; k < strlen(c1); k++) _put(p, c1[k]);
			for (q = r[j].cigar; *q; q++) _put(p, *q);
			for (k = 0; k < strlen(c2); k++) _put(p, c2[k]);
			if (r[j].flag&0x10) { for (k = qe-1; k >= qs; k--) _put(p, "NTGKCYSBAWRDMHVN"[(uint8_t)t->seq[k]]); }
			else { for (k = qs; k < qe; k++) _put(p, "NACMGRSVTWYHKDBN"[(uint8_t)t->seq[k]]); }
			for (k = 0; k < strlen(c3); k++) _put(p, c3[k]);
			free(r[j].cigar);
		}
		free(r);
	}
	free(s->reg.a); free(s->base); free(s->seq); free(s);
	return;
}

int mm_align_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, uint32_t n_threads, uint32_t tbatch_size)
{
	uint32_t i;
	mm_align_pipeline_t pl = {0};
	tbuf_t **t;
	n_threads += (n_threads == 0);
	pl.fp = bseq_open(fn);
	if (pl.fp == 0) return -1;

	// initialize alignment context
	pl.opt = opt, pl.mi = idx;
	pl.n_threads = n_threads, pl.batch_size = tbatch_size;
	int m = opt->m, x = -opt->x, gi = opt->gi, ge = opt->ge;
	struct gaba_score_s sc = {{{m,x,x,x},{x,m,x,x},{x,x,m,x},{x,x,x,m}},gi,ge,gi,ge};
	struct gaba_params_s p = {0,0,0,opt->xdrop,&sc};
	pl.gaba = gaba_init(&p);
	pl.buf_size = 512 * 1024;
	pl.base = pl.p = malloc(sizeof(char) * pl.buf_size);
	pl.tail = pl.base + pl.buf_size;

	// print sam header
	printf("@HD\tVN:1.0\tSO:unsorted\n");
	for (i = 0; i < idx->n; ++i) printf("@SQ\tSN:%s\tLN:%d\n", idx->s.a[i].name, idx->s.a[i].l_seq);
	printf("@RG\tID:1\n");
	printf("@PG\tID:6\tPN:minialign\n");

	t = (tbuf_t**)calloc(n_threads, sizeof(tbuf_t*));
	for (i = 0; i < n_threads; ++i) t[i] = mm_tbuf_init(&pl);
	ptask_t *pt = ptask_init(mm_align_worker, (void**)t, n_threads, 256);
	ptask_stream(pt, mm_align_source, &pl, mm_align_drain, &pl, 256/n_threads);
	ptask_clean(pt);
	for (i = 0; i < n_threads; ++i) mm_tbuf_destroy(t[i]);
	free(t);

	fwrite(pl.base, sizeof(char), pl.p - pl.base, stdout);

	bseq_close(pl.fp);
	gaba_clean(pl.gaba);
	free(pl.base);
	return 0;
}

/* end of map.c */

/* main.c */

void liftrlimit()
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
	mm_mapopt_t opt;
	int i, ch, k = 15, w = -1, b = 14, n_threads = 3, is_idx = 0;
	int tbatch_size = 10000000;
	float f = 0.001;
	bseq_file_t *fp = 0;
	char *fnw = 0;
	FILE *fpr = 0, *fpw = 0;

	liftrlimit();
	mm_realtime0 = realtime();
	mm_mapopt_init(&opt);

	while ((ch = getopt(argc, argv, "k:w:f:B:t:v:d:lM:m:a:b:p:q:L:H:X:Vx:")) >= 0) {
		if (ch == 'k') k = atoi(optarg);
		else if (ch == 'w') w = atoi(optarg);
		else if (ch == 'f') f = atof(optarg);
		else if (ch == 'B') b = atoi(optarg);
		else if (ch == 't') n_threads = atoi(optarg);
		else if (ch == 'v') mm_verbose = atoi(optarg);
		else if (ch == 'd') fnw = optarg;
		else if (ch == 'l') is_idx = 1;
		else if (ch == 'M') opt.min_len = atoi(optarg);
		else if (ch == 'm') opt.min_score = atoi(optarg);
		else if (ch == 'a') opt.m = atoi(optarg);
		else if (ch == 'b') opt.x = atoi(optarg);
		else if (ch == 'p') opt.gi = atoi(optarg);
		else if (ch == 'q') opt.ge = atoi(optarg);
		else if (ch == 'L') opt.ofs_llim = atoi(optarg);
		else if (ch == 'H') opt.ofs_hlim = atoi(optarg);
		else if (ch == 'X') opt.xdrop = atoi(optarg);
		else if (ch == 'V') {
			puts(MM_VERSION);
			return 0;
		} else if (ch == 'x') {
			if (strcmp(optarg, "ava10k") == 0) {
				opt.min_len = 100;
				opt.m = 1; opt.x = 1; opt.gi = 1; opt.ge = 1;
				w = 5;
			}
		}
	}
	if (w < 0) w = (int)(.6666667 * k + .499);

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
						"    $ minialign [options] <ref.fa> -d <index.mai>\n"
						"    $ minialign -l <index.mai> <reads.fa> > mapping.sam\n"
						"\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Indexing:\n");
		fprintf(stderr, "    -k INT     k-mer size [%d]\n", k);
		fprintf(stderr, "    -w INT     minimizer window size [{-k}*2/3]\n");
		fprintf(stderr, "    -f FLOAT   filter out top FLOAT fraction of repetitive minimizers [%.3f]\n", f);
		fprintf(stderr, "    -d FILE    dump index to FILE []\n");
		fprintf(stderr, "    -l FILE    load index from FILE [] (overriding -k, -w, and -f)\n");
		fprintf(stderr, "  Mapping:\n");
		fprintf(stderr, "    -a INT     match award [%d]\n", opt.m);
		fprintf(stderr, "    -b INT     mismatch penalty [%d]\n", opt.x);
		fprintf(stderr, "    -p INT     gap open penalty [%d]\n", opt.gi);
		fprintf(stderr, "    -q INT     gap extension penalty [%d]\n", opt.ge);
		fprintf(stderr, "    -m INT     minimum alignment score [%d]\n", opt.min_score);
		fprintf(stderr, "  Misc:\n");
		fprintf(stderr, "    -t INT     number of threads [%d]\n", n_threads);
		fprintf(stderr, "    -V         show version number\n");
		fprintf(stderr, "\nSee minialign.1 for detailed description of the command-line options.\n");
		return 1;
	}

	if (is_idx) fpr = fopen(argv[optind], "rb");
	else fp = bseq_open(argv[optind]);
	if (fnw) fpw = fopen(fnw, "wb");
	for (;;) {
		mm_idx_t *mi = 0;
		if (fpr) mi = mm_idx_load(fpr);
		else if (!bseq_eof(fp))
			mi = mm_idx_gen(fp, w, k, b, f, tbatch_size, n_threads);
		if (mi == 0) break;
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n);
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s] max occurrences of a minimizer to consider: %lu\n", __func__, mi->max_occ);
		if (fpw) mm_idx_dump(fpw, mi);
		for (i = optind + 1; i < argc; ++i)
			mm_align_file(mi, argv[i], &opt, n_threads, tbatch_size);
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
