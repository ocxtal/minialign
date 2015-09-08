#include <stdlib.h>
#include <stdio.h>
#include "bseq.h"
#include "minimap.h"

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	int batch_size, n_processed, n_threads;
	int d, m, thres;
	float f;
	bseq_file_t *fp;
	const mm_idx_t *mi;
} pipeline_t;

typedef struct { // per-thread buffer
	mm128_v a;
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
	b->a.n = 0;
	mm_sketch(t->seq, t->l_seq, step->p->mi->w, step->p->mi->k, t->rid, &b->a);
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
		}
		free(s->buf);
		free(s);
	}
    return 0;
}

int mm_map(const mm_idx_t *idx, const char *fn, int d, int m, float f, int n_threads, int batch_size)
{
	pipeline_t pl;
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
