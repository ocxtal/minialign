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

static inline void shift_window(int t, kdq_t(int) *w, double T, int W, int L, int rw, int rv, int *cw, int *cv)
{
	int s;
	if (kdq_size(w) >= W - 2) {
		s = kdq_shift(int, w);
		rw -= --cw[t];
		if (L > kdq_size(w))
			--L, rv -= --cv[t];
	}
	kdq_push(int, w, t);
	++L;
	rw += cw[t]++;
	rv += cv[t]++;
	if (cv[t] > 2. * T) {
		do {
			s = kdq_at(w, kdq_size(w) - L);
			rv -= --cv[t];
			--L;
		} while (s == t);
	}
}

static inline void save_masked_region(uint64_v *res, perf_intv_v *P, int start)
{
	int i, saved = 0;
	if (P->n == 0 || P->a[P->n - 1].start >= start) return;
	if (res->n) {
		int s = res->a[res->n - 1]>>32, f = (uint32_t)res->a[res->n - 1];
		perf_intv_t *p = &P->a[P->n - 1];
		if (p->start <= f + 1)
			saved = 1, res->a[res->n - 1] = (uint64_t)s<<32 | (f > p->finish? f : p->finish);
	}
	if (!saved) kv_push(uint64_t, *res, (uint64_t)p->start<<32|p->finish);
	for (i = P->n - 1; i >= 0 && P->a[i].start < start; --i);
	P->n = i + 1;
}

static void find_perfect(perf_intv_v *P, const kdq_t(int) *w, double T, int start, int L, int rv, const int *cv)
{
	int c[SD_WTOT], r = rv, i;
	double max_score = 0.;
	memcpy(c, cv, SD_WTOT * sizeof(int));
	for (i = kdq_size(w) - L - 1; i >= 0; --i) {
		int j, t = kdq_at(w, i);
		double new_score;
		r += c[t]++;
		new_score = (double)r / (kdq_size(w) - i - 1);
		if (new_score > T) {
			for (j = 0; j < P->n && P->a[j].start >= i + start; ++j) // find insert point
				max_score = max_score > P->a[j].S? max_score : P->a[j].S;
			if (new_score >= max_score) { // then insert
				if (P->n == P->m)
					kv_resize(perf_intv_t, *P, P->n + 1);
				memmove(&P->a[j+1], &P->a[j], (P->n - j) * sizeof(perf_intv_t));
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
	perf_intv_v P = {0,0,0}; // _P_ keeps the list of perfect intervals for the current window, sorted by ascending start and then by descending finish
	kdq_t(int) *w;
	int i, start, l;
	unsigned t;

	memset(cv, 0, SD_WTOT * sizeof(int));
	memset(cw, 0, SD_WTOT * sizeof(int));
	w = kdq_init(int);

	for (i = l = t = 0; i < l_seq; ++i) {
		int b = seq_nt4_table[seq[i]];
		if (b < 3) {
			++l, t = (t<<2 | b) & SD_WMSK;
			if (l >= SD_WLEN) {
				start = i - W + 1 > 0? i - W + 1 : 0;
				save_masked_regions(&res, &P, start);
				shift_window(t, w, T, start, L, rv, cv);
			}
		} else l = t = 0;
	}
	start = l_seq - W + 1 > 0? l_seq - W + 1 : 0;
	while (P->n)
		save_masked_regions(&res, &P, start++);

	kdq_destroy(int, w); free(P.a);
	*n = res.n;
	return res.a;
}
