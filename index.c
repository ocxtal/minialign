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
	int w, k;
	int batch_size, n_processed;
	bseq_file_t *fp;
	mm128_v a;
} pipeline_t;

typedef struct {
    int n_seq;
	bseq1_t *seq;
} step_t;

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
		int i;
        step_t *s = (step_t*)in;
		for (i = 0; i < s->n_seq; ++i) {
			bseq1_t *t = &s->seq[i];
			mm_sketch(t->seq, t->l_seq, p->w, p->k, t->rid, &p->a);
			free(t->seq); free(t->name);
		}
		free(s->seq); free(s);
        return 0;
    }
    return 0;
}

int main_index(int argc, char *argv[])
{
	int c, debug_print = 0;
	pipeline_t pl;

	memset(&pl, 0, sizeof(pipeline_t));
	pl.w = 20, pl.k = 16, pl.batch_size = 10000000;

	while ((c = getopt(argc, argv, "w:k:B:p")) >= 0) {
		if (c == 'w') pl.w = atoi(optarg);
		else if (c == 'k') pl.k = atoi(optarg);
		else if (c == 'B') pl.batch_size = atoi(optarg);
		else if (c == 'p') debug_print = 1;
	}

	if (argc == optind) {
		fprintf(stderr, "Usage: minimap index [options] <in.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", pl.k);
		fprintf(stderr, "  -w INT     minizer window size [%d]\n", pl.w);
		fprintf(stderr, "  -B INT     batch size [%d]\n", pl.batch_size);
		return 1;
	}

	pl.fp = bseq_open(argv[optind]);
	assert(pl.fp);
	kt_pipeline(2, worker_pipeline, &pl, 2);
	bseq_close(pl.fp);
	if (debug_print) {
		int i;
		for (i = 0; i < pl.a.n; ++i)
			printf("%ld\t%ld\t%lx\n", (long)(pl.a.a[i].y>>32), (long)pl.a.a[i].y, (long)pl.a.a[i].x);
	}
	free(pl.a.a);
	return 0;
}
