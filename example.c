#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include "minimap.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	bseq_file_t *fp;
	mm_mapopt_t opt;
	mm_idx_t *mi;
	gzFile f;
	kseq_t *ks;
	mm_tbuf_t *tbuf;
	int n_threads = 4, w = 10, k = 15;

	if (argc < 3) {
		fprintf(stderr, "Usage: minimap-lite <target.fa> <query.fa>\n");
		return 1;
	}
	
	// initialize parameters
	mm_mapopt_init(&opt);

	// open query file for reading
	f = gzopen(argv[2], "r"); // FIXME: test file existence
	ks = kseq_init(f);

	// create index for target; we are creating one index for all target sequence
	fp = bseq_open(argv[1]); // FIXME: test file existence
	mi = mm_idx_gen(fp, w, k, MM_IDX_DEF_B, 1<<18, n_threads, UINT64_MAX, 1);
	mm_idx_set_max_occ(mi, 0.001);
	bseq_close(fp);

	// mapping
	tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
	while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
		const mm_reg1_t *reg;
		int j, n_reg;
		// get all hits for the query
		reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &opt, 0);
		// print them out
		for (j = 0; j < n_reg; ++j) {
			const mm_reg1_t *r = &reg[j];
			printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
			if (mi->name) fputs(mi->name[r->rid], stdout);
			else printf("%d", r->rid + 1);
			printf("\t%d\t%d\t%d\t%d\t%d\n", mi->len[r->rid], r->rs, r->re, r->len, r->cnt);
		}
	}
	mm_tbuf_destroy(tbuf);

	// deallocate index and close the query file
	mm_idx_destroy(mi);
	kseq_destroy(ks);
	gzclose(f);
	return 0;
}
