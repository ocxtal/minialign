#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "bseq.h"
#include "kvec.h"
#include "minimap.h"

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	int batch_size, n_processed, n_threads;
	int radius, max_gap, min_cnt;
	uint32_t thres;
	float f;
	bseq_file_t *fp;
	const mm_idx_t *mi;
} pipeline_t;

typedef struct {
	uint32_t cnt:31, rev:1;
	int32_t rid;
	int32_t qs, qe, rs, re;
} mm_reg1_t;

typedef struct { // per-thread buffer
	mm128_v mini, coef, intv;
	uint64_v stack;
	kvec_t(mm_reg1_t) reg;
	uint32_t n, m;
	uint64_t *a;
	size_t *b, *p;
} tbuf_t;

typedef struct {
	const pipeline_t *p;
    int n_seq;
	bseq1_t *seq;
	int *n_reg;
	mm_reg1_t **reg;
	tbuf_t *buf;
} step_t;

#include "ksort.h"
#define sort_key_64(a) (a)
KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 8) 
#define lt_low32(a, b) ((uint32_t)(a) < (uint32_t)(b))
KSORT_INIT(low32lt, uint64_t, lt_low32)
#define gt_low32(a, b) ((uint32_t)(a) > (uint32_t)(b))
KSORT_INIT(low32gt, uint64_t, gt_low32)

static void proc_intv(tbuf_t *b, int which, int k, int min_cnt, int max_gap)
{
	int i, j, l_lis, rid = -1, rev = 0, start = b->intv.a[which].y, end = start + b->intv.a[which].x;
	if (end - start > b->m) {
		b->m = end - start;
		kv_roundup32(b->m);
		b->a = (uint64_t*)realloc(b->a, b->m * 8);
		b->b = (size_t*)realloc(b->b, b->m * sizeof(size_t));
		b->p = (size_t*)realloc(b->p, b->m * sizeof(size_t));
	}
	b->n = 0;
	for (i = start; i < end; ++i)
		if (b->coef.a[i].x != UINT64_MAX)
			b->a[b->n++] = b->coef.a[i].y, rid = b->coef.a[i].x << 1 >> 33, rev = b->coef.a[i].x >> 63;
	if (b->n < min_cnt) return;
	radix_sort_64(b->a, b->a + b->n);
	l_lis = rev? ks_lis_low32gt(b->n, b->a, b->b, b->p) : ks_lis_low32lt(b->n, b->a, b->b, b->p);
	if (l_lis < min_cnt) return;
	for (i = 1, j = 1; i < l_lis; ++i)
		if (b->a[b->b[i]]>>32 != b->a[b->b[i-1]]>>32)
			b->a[b->b[j++]] = b->a[b->b[i]];
	l_lis = j;
	if (l_lis < min_cnt) return;
	for (i = 1, start = 0; i <= l_lis; ++i) {
		if (i == l_lis || ((b->a[b->b[i]]>>32) - (b->a[b->b[i-1]]>>32) > max_gap && abs((int32_t)b->a[b->b[i]] - (int32_t)b->a[b->b[i-1]]) > max_gap)) {
			if (i - start >= min_cnt) {
				mm_reg1_t *r;
				kv_pushp(mm_reg1_t, b->reg, &r);
				r->rid = rid, r->rev = rev, r->cnt = i - start;
				r->qs = (b->a[b->b[start]]>>32) - (k - 1);
				r->qe = (b->a[b->b[i-1]]>>32) + 1;
				r->rs = rev? (uint32_t)b->a[b->b[i-1]] : (uint32_t)b->a[b->b[start]];
				r->re = rev? (uint32_t)b->a[b->b[start]] : (uint32_t)b->a[b->b[i-1]];
				r->rs -= k - 1;
				r->re += 1;
			}
			start = i;
		}
	}
}

static inline void push_intv(mm128_v *intv, int start, int end)
{
	mm128_t *p;
	if (intv->n > 0) {
		int last_start, last_end, min;
		p = &intv->a[intv->n-1];
		last_start = p->y, last_end = p->x + last_start;
		min = end - start < last_end - last_start? end - start : last_end - last_start;
		if (last_end > start && last_end - start > (min>>1) + (min>>2)) {
			p->x = end - last_start;
			return;
		}
	}
	kv_pushp(mm128_t, *intv, &p);
	p->x = end - start, p->y = start;
}

static void get_reg(tbuf_t *b, int radius, int k, int min_cnt, int max_gap)
{
	mm128_v *c = &b->coef;
	int i, start = 0;
	if (c->n < min_cnt) return;
	b->intv.n = 0;
	for (i = 1; i < c->n; ++i) { // identify all (possibly overlapping) clusters within _radius_
		if (c->a[i].x - c->a[start].x > radius) {
			if (i - start >= min_cnt) push_intv(&b->intv, start, i);
			for (++start; start < i && c->a[i].x - c->a[start].x > radius; ++start);
		}
	}
	if (i - start >= min_cnt) push_intv(&b->intv, start, i);
	radix_sort_128x(b->intv.a, b->intv.a + b->intv.n); // sort by the size of the cluster
	for (i = b->intv.n - 1; i >= 0; --i) proc_intv(b, i, k, min_cnt, max_gap);
}

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *step = (step_t*)_data;
	bseq1_t *t = &step->seq[i];
	tbuf_t *b = &step->buf[tid];
	const mm_idx_t *mi = step->p->mi;
	int j;

	b->mini.n = b->coef.n = 0;
	mm_sketch(t->seq, t->l_seq, mi->w, mi->k, t->rid, &b->mini);
	for (j = 0; j < b->mini.n; ++j) {
		int k, n;
		const uint64_t *r;
		int32_t qpos = (uint32_t)b->mini.a[j].y>>1, strand = b->mini.a[j].y&1;
		r = mm_idx_get(mi, b->mini.a[j].x, &n);
		if (n > step->p->thres) continue;
		for (k = 0; k < n; ++k) {
			int32_t rpos = (uint32_t)r[k] >> 1;
			mm128_t *p;
			kv_pushp(mm128_t, b->coef, &p);
			if ((r[k]&1) == strand) { // forward strand
				p->x = (uint64_t)r[k] >> 32 << 32 | (0x80000000U + rpos - qpos);
				p->y = (uint64_t)qpos << 32 | rpos;
			} else { // reverse strand
				p->x = (uint64_t)r[k] >> 32 << 32 | (rpos + qpos) | 1ULL<<63;
				p->y = (uint64_t)qpos << 32 | rpos;
			}
		}
	}
	radix_sort_128x(b->coef.a, b->coef.a + b->coef.n);
	if (mm_verbose >= 5) { // NB: the following block may be corrupted by multi-threading
		printf(">%s\n", t->name);
		for (j = 0; j < b->coef.n; ++j) {
			uint64_t x = b->coef.a[j].x;
			uint32_t off = x>>63 == 0? 0x80000000U : 0;
			printf("%d\t%d\t%c\t%d\n", j, (uint32_t)(x<<1>>33), "+-"[x>>63], (int32_t)((uint32_t)x - off));
		}
		printf("//\n");
	}
	b->reg.n = 0;
	get_reg(b, step->p->radius, mi->k, step->p->min_cnt, step->p->max_gap);
	step->n_reg[i] = b->reg.n;
	step->reg[i] = (mm_reg1_t*)malloc(b->reg.n * sizeof(mm_reg1_t));
	memcpy(step->reg[i], b->reg.a, b->reg.n * sizeof(mm_reg1_t));
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i, j;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = bseq_read(p->fp, p->batch_size, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].rid = p->n_processed++;
			s->buf = (tbuf_t*)calloc(p->n_threads, sizeof(tbuf_t));
			s->n_reg = (int*)calloc(s->n_seq, sizeof(int));
			s->reg = (mm_reg1_t**)calloc(s->n_seq, sizeof(mm_reg1_t**));
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: map
		kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_seq);
		return in;
    } else if (step == 2) { // step 2: output
        step_t *s = (step_t*)in;
		const mm_idx_t *mi = p->mi;
		for (i = 0; i < p->n_threads; ++i) { // free temporary data
			tbuf_t *b = &s->buf[i];
			free(b->mini.a); free(b->coef.a); free(b->intv.a);
			free(b->stack.a); free(b->reg.a);
			free(b->a); free(b->b); free(b->p);
		}
		free(s->buf);
		for (i = 0; i < s->n_seq; ++i) {
			bseq1_t *t = &s->seq[i];
			for (j = 0; j < s->n_reg[i]; ++j) {
				mm_reg1_t *r = &s->reg[i][j];
				printf("%s\t%d\t%d\t%d\t%c\t", t->name, t->l_seq, r->qs, r->qe, "+-"[r->rev]);
				if (mi->name) fputs(mi->name[r->rid], stdout);
				else printf("%d", r->rid + 1);
				printf("\t%d\t%d\t%d\t%d\n", mi->len[r->rid], r->rs, r->re, r->cnt);
			}
			free(s->reg[i]);
			free(s->seq[i].seq); free(s->seq[i].name);
		}
		free(s->reg); free(s->n_reg); free(s->seq);
		free(s);
	}
    return 0;
}

int mm_map(const mm_idx_t *idx, const char *fn, int radius, int max_gap, int min_cnt, float f, int n_threads, int batch_size)
{
	pipeline_t pl;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.fp = bseq_open(fn);
	if (pl.fp == 0) return -1;
	pl.mi = idx, pl.radius = radius, pl.min_cnt = min_cnt, pl.max_gap = max_gap;
	pl.thres = mm_idx_thres(idx, f);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s] max occurrences of a minimizer to consider: %d\n", __func__, pl.thres);
	pl.n_threads = n_threads, pl.batch_size = batch_size;
	kt_pipeline(n_threads == 1? 1 : 2, worker_pipeline, &pl, 3);
	bseq_close(pl.fp);
	return 0;
}
