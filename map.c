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

typedef struct { // per-thread buffer
	mm128_v a, c[2];
} tbuf_t;

typedef struct {
	const pipeline_t *p;
    int n_seq;
	bseq1_t *seq;
	tbuf_t *buf;
} step_t;

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
				p->x = (uint64_t)r[k] >> 32 << 32 | (rpos + qpos - mi->k);
				p->y = (uint64_t)(qpos - mi->k) << 32 | rpos;
			}
		}
	}
	radix_sort_128x(b->c[0].a, b->c[0].a + b->c[0].n);
	radix_sort_128x(b->c[1].a, b->c[1].a + b->c[1].n);
	printf(">%s\n", t->name);
	for (j = 0; j < 2; ++j) {
		int k;
		for (k = 0; k < b->c[j].n; ++k) {
			uint64_t x = b->c[j].a[k].x;
			uint32_t off = j == 0? 0x80000000U : 0;
			printf("%d\t%c\t%d\n", (uint32_t)(x>>32), "+-"[j], (int32_t)((uint32_t)x - off));
		}
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
