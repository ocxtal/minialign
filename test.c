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
	int c, w = 30, k = 14;
	uint64_t shift, mask;
	while ((c = getopt(argc, argv, "w:k:")) >= 0) {
		if (c == 'w') w = atoi(optarg);
		else if (c == 'k') k = atoi(optarg);
	}
	shift = 2 * k, mask = (1ULL<<shift) - 1;
	fp = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		int i, n;
		uint64_t *a;
		a = mm_sketch(seq->seq.s, seq->seq.l, w, k, &n);
		printf(">%s\n", seq->name.s);
		for (i = 0; i < n; ++i)
			printf("%llu\t%llx\n", a[i]>>shift, a[i]&mask);
		printf("//\n");
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
