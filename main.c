#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "minimap.h"

#define MM_VERSION "r49"

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
	int i, c, k = 15, w = -1, b = 14, radius = 500, max_gap = 10000, min_cnt = 4, n_threads = 3, batch_size = 10000000, keep_name = 1;
	float f = 0.001;
	mm_idx_t *mi = 0;

	liftrlimit();
	mm_realtime0 = realtime();

	while ((c = getopt(argc, argv, "w:k:B:b:t:r:c:f:Vv:Ng:")) >= 0) {
		if (c == 'w') w = atoi(optarg);
		else if (c == 'k') k = atoi(optarg);
		else if (c == 'b') b = atoi(optarg);
		else if (c == 'r') radius = atoi(optarg);
		else if (c == 'c') min_cnt = atoi(optarg);
		else if (c == 'f') f = atof(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'B') batch_size = atoi(optarg);
		else if (c == 'v') mm_verbose = atoi(optarg);
		else if (c == 'g') max_gap = atoi(optarg);
		else if (c == 'N') keep_name = 0;
		else if (c == 'V') {
			puts(MM_VERSION);
			return 0;
		}
	}
	if (w < 0) w = k;

	if (argc == optind) {
		fprintf(stderr, "Usage: minimap [options] <target.fa> [query.fa]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", k);
		fprintf(stderr, "  -w INT     minizer window size [same as -k]\n");
		fprintf(stderr, "  -b INT     bucket bits [%d]\n", b);
		fprintf(stderr, "  -t INT     number of threads [%d]\n", n_threads);
		fprintf(stderr, "  -r INT     bandwidth [%d]\n", radius);
		fprintf(stderr, "  -c INT     min minimizer match count [%d]\n", min_cnt);
		fprintf(stderr, "  -f FLOAT   minimizer filteration threshold [%.3f]\n", f);
		fprintf(stderr, "  -g INT     break a chain if there is a gap longer than INT [%d]\n", max_gap);
		fprintf(stderr, "  -B INT     batch size [%d]\n", batch_size);
		fprintf(stderr, "  -v INT     verbose level [%d]\n", mm_verbose);
		fprintf(stderr, "  -N         use integer as target names\n");
		fprintf(stderr, "  -V         show version number\n");
		return 1;
	}

	mi = mm_idx_gen(argv[optind], w, k, b, batch_size, n_threads, keep_name);
	if (argc - optind >= 2)
		mm_map(mi, argv[optind+1], radius, max_gap, min_cnt, f, n_threads, batch_size);
	mm_idx_destroy(mi);

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
	return 0;
}
