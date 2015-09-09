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
	int d, m;
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
	mm128_v a, c[2];
	uint64_v stack;
	kvec_t(mm_reg1_t) reg;
} tbuf_t;

typedef struct {
	const pipeline_t *p;
    int n_seq;
	bseq1_t *seq;
	kvec_t(mm_reg1_t) *reg;
	tbuf_t *buf;
} step_t;

static void get_reg(tbuf_t *b, int radius, int min_cnt, int k, int rev)
{
	mm128_v *c = &b->c[rev];
	mm128_t *p;
	int i, j, start = 0;
	if (c->n < min_cnt) return;
	b->a.n = 0;
	for (i = 1; i < c->n; ++i) { // identify all (possibly overlapping) clusters within _radius_
		if (c->a[i].x - c->a[start].x > radius) {
			if (i - start >= min_cnt) {
				kv_pushp(mm128_t, b->a, &p);
				p->x = i - start, p->y = start;
			}
			for (++start; start < i && c->a[i].x - c->a[start].x > radius; ++start);
		}
	}
	if (i - start >= min_cnt) { // the last cluster
		kv_pushp(mm128_t, b->a, &p);
		p->x = i - start, p->y = start;
	}
	radix_sort_128x(b->a.a, b->a.a + b->a.n); // sort by the size of the cluster
	for (i = b->a.n - 1; i >= 0; --i) { // starting from the largest cluster
		int start = b->a.a[i].y;
		int end = start + b->a.a[i].x;
		int cnt = 0;
		mm_reg1_t *r;
		for (j = start; j < end; ++j) // exclude minimizer hits that have been used in larger clusters
			if (c->a[j].x != UINT64_MAX) ++cnt;
		if (cnt < min_cnt) continue;
		while (c->a[start].x == UINT64_MAX) ++start; // we do this as we need the first to be present
		kv_pushp(mm_reg1_t, b->reg, &r);
		r->cnt = cnt, r->rid = c->a[start].x>>32, r->rev = rev;
		r->qs = r->rs = INT32_MAX; r->qe = r->re = 0;
		for (j = start; j < end; ++j) {
			int qpos, rpos;
			if (c->a[j].x == UINT64_MAX) continue;
			qpos = c->a[j].y>>32;
			rpos = (int32_t)c->a[j].y;
			r->qs = r->qs < qpos? r->qs : qpos;
			r->rs = r->rs < rpos? r->rs : rpos;
			r->qe = r->qe > qpos? r->qe : qpos;
			r->re = r->re > rpos? r->re : rpos;
			c->a[j].x = UINT64_MAX; // mark the minimizer hit having been used
		}
		if (rev == 0) r->qs -= k - 1, ++r->qe; // we need to correct for k
		else r->qe += k; // and a hit to the reverse strand is different
		r->rs -= k - 1, ++r->re;
	}
}

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *step = (step_t*)_data;
	bseq1_t *t = &step->seq[i];
	tbuf_t *b = &step->buf[tid];
	const mm_idx_t *mi = step->p->mi;
	int j, N = 0;

	b->a.n = b->c[0].n = b->c[1].n = 0;
	mm_sketch(t->seq, t->l_seq, mi->w, mi->k, t->rid, &b->a);
	for (j = 0; j < b->a.n; ++j) {
		int k, n;
		const uint64_t *r;
		int32_t qpos = (uint32_t)b->a.a[j].y>>1, strand = b->a.a[j].y&1;
		r = mm_idx_get(mi, b->a.a[j].x, &n);
		N += n;
		if (n > step->p->thres) continue;
		for (k = 0; k < n; ++k) {
			int32_t rpos = (uint32_t)r[k] >> 1;
			mm128_t *p;
			if ((r[k]&1) == strand) {
				kv_pushp(mm128_t, b->c[0], &p);
				p->x = (uint64_t)r[k] >> 32 << 32 | (0x80000000U + rpos - qpos);
				p->y = (uint64_t)qpos << 32 | rpos;
			} else {
				kv_pushp(mm128_t, b->c[1], &p);
				p->x = (uint64_t)r[k] >> 32 << 32 | (rpos + qpos - mi->k + 1);
				p->y = (uint64_t)(qpos - mi->k + 1) << 32 | rpos;
			}
		}
	}
	printf(">%s\t%ld\n", t->name, b->reg.n);
	radix_sort_128x(b->c[0].a, b->c[0].a + b->c[0].n);
	radix_sort_128x(b->c[1].a, b->c[1].a + b->c[1].n);
	/*
	for (j = 0; j < 2; ++j) {
		int k;
		for (k = 0; k < b->c[j].n; ++k) {
			uint64_t x = b->c[j].a[k].x;
			uint32_t off = j == 0? 0x80000000U : 0;
			printf("%d\t%d\t%c\t%d\n", k, (uint32_t)(x>>32), "+-"[j], (int32_t)((uint32_t)x - off));
		}
	}
	*/
	b->reg.n = 0;
	get_reg(b, step->p->d, step->p->m, mi->k, 0);
	get_reg(b, step->p->d, step->p->m, mi->k, 1);
	for (j = 0; j < b->reg.n; ++j) {
		mm_reg1_t *r = &b->reg.a[j];
		printf("%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\n", t->name, r->qs, r->qe, "+-"[r->rev], r->rid, r->rs, r->re, r->cnt);
	}
	printf("//\n");
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
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
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: map
		kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_seq);
		return in;
    } else if (step == 2) { // step 2: output
        step_t *s = (step_t*)in;
		for (i = 0; i < p->n_threads; ++i) {
			free(s->buf[i].a.a);
			free(s->buf[i].c[0].a);
			free(s->buf[i].c[1].a);
		}
		free(s->buf);
		free(s);
	}
    return 0;
}

int mm_map(const mm_idx_t *idx, const char *fn, int d, int m, float f, int n_threads, int batch_size)
{
	pipeline_t pl;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.fp = bseq_open(fn);
	if (pl.fp == 0) return -1;
	pl.mi = idx, pl.d = d, pl.m = m;
	pl.thres = mm_idx_thres(idx, f);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s] max occurrences of a minimizer to consider: %d\n", __func__, pl.thres);
	pl.n_threads = n_threads, pl.batch_size = batch_size;
	kt_pipeline(n_threads == 1? 1 : 2, worker_pipeline, &pl, 3);
	bseq_close(pl.fp);
	return 0;
}
