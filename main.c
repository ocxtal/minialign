#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "minimap.h"

#define MM_VERSION "r20"

void liftrlimit()
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}

int main_index(int argc, char *argv[])
{
	int c, k = 16, w = 16, b = 14, n_threads = 3, batch_size = 10000000;
	mm_idx_t *mi = 0;

	while ((c = getopt(argc, argv, "w:k:B:b:t:")) >= 0) {
		if (c == 'w') w = atoi(optarg);
		else if (c == 'k') k = atoi(optarg);
		else if (c == 'b') b = atoi(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'B') batch_size = atoi(optarg);
	}

	if (argc == optind) {
		fprintf(stderr, "Usage: minimap index [options] <in.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", k);
		fprintf(stderr, "  -w INT     minizer window size [%d]\n", w);
		fprintf(stderr, "  -b INT     bucket bits [%d]\n", b);
		fprintf(stderr, "  -t INT     number of threads [%d]\n", n_threads);
		fprintf(stderr, "  -B INT     batch size [%d]\n", batch_size);
		return 1;
	}

	mi = mm_idx_gen(argv[optind], w, k, b, batch_size, n_threads);
	mm_idx_destroy(mi);
	return 0;
}

int main(int argc, char *argv[])
{
	int ret = 0, i;
	double t_start;
	liftrlimit();
	if (argc == 1) {
		fprintf(stderr, "Usage: minimap <command> [arguments]\n");
		fprintf(stderr, "Commands:\n");
		fprintf(stderr, "  index    create minimap index\n");
		return 1;
	}
	mm_realtime0 = t_start = realtime();
	if (strcmp(argv[1], "index") == 0) ret = main_index(argc-1, argv+1);
	else {
		fprintf(stderr, "[E::%s] unknown command\n", __func__);
		return 1;
	}
	if (ret == 0) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_start, cputime());
	}
	return ret;
}
