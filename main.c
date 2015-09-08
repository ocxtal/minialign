#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "minimap.h"

#define MM_VERSION "r25"

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
	int i, c, k = 16, w = 16, b = 14, d = 100, m = 3, n_threads = 3, batch_size = 10000000;
	mm_idx_t *mi = 0;

	liftrlimit();
	mm_realtime0 = realtime();

	while ((c = getopt(argc, argv, "w:k:B:b:t:d:m:")) >= 0) {
		if (c == 'w') w = atoi(optarg);
		else if (c == 'k') k = atoi(optarg);
		else if (c == 'b') b = atoi(optarg);
		else if (c == 'd') d = atoi(optarg);
		else if (c == 'm') m = atoi(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'B') batch_size = atoi(optarg);
	}

	if (argc == optind) {
		fprintf(stderr, "Usage: minimap [options] <target.fa> [query.fa]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", k);
		fprintf(stderr, "  -w INT     minizer window size [%d]\n", w);
		fprintf(stderr, "  -b INT     bucket bits [%d]\n", b);
		fprintf(stderr, "  -t INT     number of threads [%d]\n", n_threads);
		fprintf(stderr, "  -d INT     bandwidth [%d]\n", d);
		fprintf(stderr, "  -m INT     min minimizer matches [%d]\n", m);
		fprintf(stderr, "  -B INT     batch size [%d]\n", batch_size);
		return 1;
	}

	mi = mm_idx_gen(argv[optind], w, k, b, batch_size, n_threads);
	if (argc - optind >= 2)
		mm_map(mi, argv[optind+1], d, m, n_threads, batch_size);
	mm_idx_destroy(mi);

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
	return 0;
}
