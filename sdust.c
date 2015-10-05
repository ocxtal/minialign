#include <string.h>
#include "kdq.h"
#include "kvec.h"

#define SD_WLEN 3
#define SD_WTOT (1<<(SD_WLEN<<1))
#define SD_WMSK (SD_WTOT - 1)

typedef struct {
	int start, finish, S;
} perf_intv_t;

typedef kvec_t(perf_intv_t) perf_intv_v;
typedef kvec_t(uint64_t) uint64_v;

KDQ_INIT(int)

#ifndef HAVE_NT4_TBL
unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
#else
extern unsigned char seq_nt4_table[256];
#endif

static inline void shift_window(int t, kdq_t(int) *w, double T, int W, int *L, int rw, int rv, int *cw, int *cv)
{
	int s;
	if (kdq_size(w) >= W - SD_WLEN + 1) { // TODO: is this right for SD_WLEN!=3?
		s = *kdq_shift(int, w);
		rw -= --cw[s];
		if (*L > kdq_size(w))
			--*L, rv -= --cv[s];
	}
	kdq_push(int, w, t);
	++*L;
	rw += cw[t]++;
	rv += cv[t]++;
	if (cv[t] > 2. * T) {
		do {
			s = kdq_at(w, kdq_size(w) - *L);
			rv -= --cv[s];
			--*L;
		} while (s == t);
	}
}

static inline void save_masked_regions(uint64_v *res, perf_intv_v *P, int start)
{
	int i, saved = 0;
	perf_intv_t *p;
	if (P->n == 0 || P->a[P->n - 1].start >= start) return;
	p = &P->a[P->n - 1];
	if (res->n) {
		int s = res->a[res->n - 1]>>32, f = (uint32_t)res->a[res->n - 1];
		if (p->start <= f + 1)
			saved = 1, res->a[res->n - 1] = (uint64_t)s<<32 | (f > p->finish? f : p->finish);
	}
	if (!saved) kv_push(uint64_t, *res, (uint64_t)p->start<<32|p->finish);
	for (i = P->n - 1; i >= 0 && P->a[i].start < start; --i); // remove perfect intervals that have falled out of the window
	P->n = i + 1;
}

static void find_perfect(perf_intv_v *P, const kdq_t(int) *w, double T, int start, int L, int rv, const int *cv)
{
	int c[SD_WTOT], r = rv, i;
	double max_score = 0.;
	memcpy(c, cv, SD_WTOT * sizeof(int));
	for (i = (long)kdq_size(w) - L - 1; i >= 0; --i) {
		int j, t = kdq_at(w, i);
		double new_score;
		r += c[t]++;
		new_score = (double)r / (kdq_size(w) - i - 1);
		if (new_score > T) {
			for (j = 0; j < P->n && P->a[j].start >= i + start; ++j) // find insertion position
				max_score = max_score > P->a[j].S? max_score : P->a[j].S;
			if (new_score >= max_score) { // then insert
				max_score = new_score;
				if (P->n == P->m) kv_resize(perf_intv_t, *P, P->n + 1);
				memmove(&P->a[j+1], &P->a[j], (P->n - j) * sizeof(perf_intv_t)); // make room
				++P->n;
				P->a[j].start = i + start, P->a[j].finish = kdq_size(w) + 1 + start, P->a[j].S = new_score;
			}
		}
	}
}

uint64_t *sdust(const uint8_t *seq, int l_seq, double T, int W, int *n)
{
	int rv = 0, rw = 0, L = 0, cv[SD_WTOT], cw[SD_WTOT];
	uint64_v res = {0,0,0};  // the result
	perf_intv_v P = {0,0,0}; // _P_ keeps the list of perfect intervals for the current window, sorted by descending start and then by ascending finish
	kdq_t(int) *w; // this caches previous words
	int i, start, l; // _start_: start of the current window; _l_: length of a contiguous A/C/G/T (sub)sequence
	unsigned t; // current word

	memset(cv, 0, SD_WTOT * sizeof(int));
	memset(cw, 0, SD_WTOT * sizeof(int));
	w = kdq_init(int);
	if (l_seq < 0) l_seq = strlen((const char*)seq);

	for (i = l = t = 0; i <= l_seq; ++i) {
		int b = i < l_seq? seq_nt4_table[seq[i]] : 4;
		if (b < 3) { // an A/C/G/T base
			++l, t = (t<<2 | b) & SD_WMSK;
			if (l >= SD_WLEN) { // we have seen a word
				start = (l - W + 1 > 0? l - W + 1 : 0) + (i - l); // set the start of the current window
				save_masked_regions(&res, &P, start); // save intervals falling out of the current window?
				shift_window(t, w, T, W, &L, rw, rv, cw, cv);
				if (rw > L * T)
					find_perfect(&P, w, T, start, L, rv, cv);
			}
		} else { // N or the end of sequence; N effectively breaks input into pieces of independent sequences
			start = (l - W + 1 > 0? l - W + 1 : 0) + (i - l);
			while (P.n) save_masked_regions(&res, &P, start++); // clear up unsaved perfect intervals
			l = t = 0;
		}
	}

	kdq_destroy(int, w); free(P.a);
	*n = res.n;
	return res.a;
}

#ifdef SDUST_MAIN
#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *ks;
	int W = 64, c;
	double T = 2.0;

	while ((c = getopt(argc, argv, "w:t:")) >= 0) {
		if (c == 'w') W = atoi(optarg);
		else if (c == 't') T = atof(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: sdust [-w %d] [-t %.1f] <in.fa>\n", W, T);
		return 1;
	}
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		uint64_t *r;
		int i, n;
		r = sdust((uint8_t*)ks->seq.s, -1, T, W, &n);
		for (i = 0; i < n; ++i)
			printf("%s\t%d\t%d\n", ks->name.s, (int)(r[i]>>32), (int)r[i]);
		free(r);
	}
	kseq_destroy(ks);
	gzclose(fp);
	return 0;
}
#endif
