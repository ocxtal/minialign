#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include "minimap.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int c, w = 20, k = 17;
	uint64_t shift, mask;
	mm128_v a = {0,0,0};
	while ((c = getopt(argc, argv, "w:k:")) >= 0) {
		if (c == 'w') w = atoi(optarg);
		else if (c == 'k') k = atoi(optarg);
	}
	shift = 2 * k, mask = (1ULL<<shift) - 1;
	fp = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		int i;
		mm_sketch128(seq->seq.s, seq->seq.l, w, k, &a);
		printf(">%s\n", seq->name.s);
		for (i = 0; i < a.n; ++i)
			printf("%llu\t%llx\n", a.a[i].y, a.a[i].x);
		printf("//\n");
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
