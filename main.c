#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "minimap.h"

#define MM_VERSION "r69"

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
	int i, c, k = 15, w = -1, b = MM_IDX_DEF_B, radius = 500, max_gap = 10000, min_cnt = 4, n_threads = 3, keep_name = 1, is_idx = 0, flag = 0;
	int tbatch_size = 10000000;
	uint64_t ibatch_size = 10000000000ULL;
	float f = 0.001;
	bseq_file_t *fp = 0;
	char *fnw = 0;
	FILE *fpr = 0, *fpw = 0;

	liftrlimit();
	mm_realtime0 = realtime();

	while ((c = getopt(argc, argv, "w:k:B:b:t:r:c:f:Vv:Ng:I:d:lRS")) >= 0) {
		if (c == 'w') w = atoi(optarg);
		else if (c == 'k') k = atoi(optarg);
		else if (c == 'b') b = atoi(optarg);
		else if (c == 'r') radius = atoi(optarg);
		else if (c == 'c') min_cnt = atoi(optarg);
		else if (c == 'f') f = atof(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'v') mm_verbose = atoi(optarg);
		else if (c == 'g') max_gap = atoi(optarg);
		else if (c == 'N') keep_name = 0;
		else if (c == 'd') fnw = optarg;
		else if (c == 'l') is_idx = 1;
		else if (c == 'R') flag |= MM_F_WITH_REP;
		else if (c == 'S') flag |= MM_F_NO_SELF;
		else if (c == 'V') {
			puts(MM_VERSION);
			return 0;
		} else if (c == 'B' || c == 'I') {
			double x;
			char *p;
			x = strtod(optarg, &p);
			if (*p == 'G' || *p == 'g') x *= 1e9;
			else if (*p == 'M' || *p == 'm') x *= 1e6;
			else if (*p == 'K' || *p == 'k') x *= 1e3;
			if (c == 'B') tbatch_size = (uint64_t)(x + .499);
			else ibatch_size = (uint64_t)(x + .499);
		}
	}
	if (w < 0) w = k;

	if (argc == optind) {
		fprintf(stderr, "Usage: minimap [options] <target.fa> [query.fa] [...]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  Indexing:\n");
		fprintf(stderr, "    -k INT     k-mer size [%d]\n", k);
		fprintf(stderr, "    -w INT     minizer window size [same as -k]\n");
		fprintf(stderr, "    -d FILE    dump index to FILE []\n");
		fprintf(stderr, "    -l         the 1st argument is a index file\n");
//		fprintf(stderr, "    -b INT     bucket bits [%d]\n", b); // most users would care about this
		fprintf(stderr, "  Mapping:\n");
		fprintf(stderr, "    -f FLOAT   filter out top FLOAT fraction of repetitive minimizers [%.3f]\n", f);
		fprintf(stderr, "    -r INT     bandwidth [%d]\n", radius);
		fprintf(stderr, "    -c INT     retain a mapping if it consists of >=INT minimizers [%d]\n", min_cnt);
		fprintf(stderr, "    -g INT     split a mapping if there is a gap longer than INT [%d]\n", max_gap);
		fprintf(stderr, "    -R         skip post-mapping repeat filtering\n");
		fprintf(stderr, "  Input/Output:\n");
		fprintf(stderr, "    -t INT     number of threads [%d]\n", n_threads);
		fprintf(stderr, "    -B NUM     process ~NUM bp in each batch [10M]\n");
		fprintf(stderr, "    -I NUM     create an index for every ~NUM bp [10G]\n");
		fprintf(stderr, "    -v INT     verbose level [%d]\n", mm_verbose);
		fprintf(stderr, "    -S         skip self mapping\n");
		fprintf(stderr, "    -N         use integer as target names\n");
		fprintf(stderr, "    -V         show version number\n");
		return 1;
	}

	if (is_idx) fpr = fopen(argv[optind], "rb");
	else fp = bseq_open(argv[optind]);
	if (fnw) fpw = fopen(fnw, "wb");
	for (;;) {
		mm_idx_t *mi = 0;
		if (fpr) mi = mm_idx_load(fpr);
		else if (!bseq_eof(fp))
			mi = mm_idx_gen(fp, w, k, b, tbatch_size, n_threads, ibatch_size, keep_name);
		if (mi == 0) break;
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n);
		mm_idx_set_max_occ(mi, f);
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s] max occurrences of a minimizer to consider: %d\n", __func__, mi->max_occ);
		if (fpw) mm_idx_dump(fpw, mi);
		for (i = optind + 1; i < argc; ++i)
			mm_map_file(mi, argv[i], radius, max_gap, min_cnt, flag, n_threads, tbatch_size);
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
