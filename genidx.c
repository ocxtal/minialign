#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "minimap.h"
#include "bseq.h"

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	int batch_size, n_processed;
	bseq_file_t *fp;
	mm128_v a;
	mm_idx_t *mi;
} pipeline_t;

typedef struct {
    int n_seq;
	bseq1_t *seq;
} step_t;

extern int mm_verbose;
extern double mm_realtime0;
double cputime(void);
double realtime(void);

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = bseq_read(p->fp, p->batch_size, &s->n_seq);
		if (s->seq) {
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].rid = p->n_processed++;
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: compute sketch
        step_t *s = (step_t*)in;
		p->a.n = 0;
		for (i = 0; i < s->n_seq; ++i) {
			bseq1_t *t = &s->seq[i];
			mm_sketch(t->seq, t->l_seq, p->mi->w, p->mi->k, t->rid, &p->a);
			free(t->seq); free(t->name);
		}
		free(s->seq);
		return s;
    } else if (step == 2) { // dispatch sketch to buckets
        step_t *s = (step_t*)in;
		mm_idx_add(p->mi, p->a.n, p->a.a);
		p->a.n = 0;
		free(s);
	}
    return 0;
}

int main_index(int argc, char *argv[])
{
	int c, k = 16, w = 16, b = 14, n_threads = 2, debug_print = 0;
	mm_idx_t *mi = 0;
	pipeline_t pl;

	memset(&pl, 0, sizeof(pipeline_t));
	pl.batch_size = 10000000;

	while ((c = getopt(argc, argv, "w:k:B:b:t:p")) >= 0) {
		if (c == 'w') w = atoi(optarg);
		else if (c == 'k') k = atoi(optarg);
		else if (c == 'b') b = atoi(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'B') pl.batch_size = atoi(optarg);
		else if (c == 'p') debug_print = 1;
	}

	if (argc == optind) {
		fprintf(stderr, "Usage: minimap index [options] <in.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", k);
		fprintf(stderr, "  -w INT     minizer window size [%d]\n", w);
		fprintf(stderr, "  -b INT     bucket bits [%d]\n", b);
		fprintf(stderr, "  -t INT     number of threads for post-processing [%d]\n", n_threads);
		fprintf(stderr, "  -B INT     batch size [%d]\n", pl.batch_size);
		return 1;
	}

	pl.fp = bseq_open(argv[optind]);
	assert(pl.fp);
	pl.mi = mm_idx_init(w, k, b);

	kt_pipeline(3, worker_pipeline, &pl, 3);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] finished collecting minimizers...\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	bseq_close(pl.fp);
	if (debug_print) {
		int i;
		for (i = 0; i < pl.a.n; ++i)
			printf("%ld\t%ld\t%lx\n", (long)(pl.a.a[i].y>>32), (long)pl.a.a[i].y, (long)pl.a.a[i].x);
	}
	free(pl.a.a);

	mi = pl.mi;
	mm_idx_post(mi, n_threads);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] finished sorting minimizers...\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	mm_idx_destroy(mi);
	return 0;
}
