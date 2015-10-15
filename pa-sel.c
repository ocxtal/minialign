#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

#include "khash.h"
KHASH_MAP_INIT_STR(str, int)
typedef khash_t(str) shash_t;

typedef struct {
	int min_match;
	int min_len;
	float min_match_frac;
} pas_opt_t;

typedef struct {
	uint64_t id_start;
	uint32_t end, match:31, rev:1;
} pas_reg_t;

typedef struct {
	char *name;
	int len, start, end, type;
} pas_seq_t;

typedef struct {
	int n_seq, m_seq;
	pas_seq_t *seq;
	shash_t *h;
	uint64_t n_reg, m_reg;
	pas_reg_t *reg;
} pas_t;

#include "ksort.h"
#define sort_reg(a) ((a).id_start)
KRADIX_SORT_INIT(reg, pas_reg_t, sort_reg, 8) 

void pas_opt_init(pas_opt_t *opt)
{
	opt->min_match = 40;
	opt->min_len = 1000;
	opt->min_match_frac = .1;
}

void pas_destroy(pas_t *p)
{
	uint64_t i;
	if (p == 0) return;
	for (i = 0; i < p->n_seq; ++i) free(p->seq[i].name);
	kh_destroy(str, p->h);
	free(p->seq); free(p->reg);
	free(p);
}

void pas_process1(const pas_opt_t *opt, pas_t *p, kstring_t *s)
{
	char *q, *r;
	int i, k, id = -1, start = -1, end = -1, rev = -1, score = -1;
	for (i = k = 0, q = s->s; i <= s->l; ++i) {
		if (i == s->l || s->s[i] == '\t') {
			s->s[i] = 0;
			if (k == 1) { // length
				khint_t itr;
				int len, absent;
				len = strtol(q, &r, 10);
				if (len < opt->min_len) break;
				itr = kh_put(str, p->h, s->s, &absent);
				if (absent) {
					if (p->n_seq == p->m_seq) {
						p->m_seq = p->m_seq? p->m_seq<<1 : 16;
						p->seq = (pas_seq_t*)realloc(p->seq, sizeof(pas_seq_t) * p->m_seq);
					}
					kh_val(p->h, itr) = id = p->n_seq++;
					kh_key(p->h, itr) = p->seq[id].name = strdup(s->s);
					p->seq[id].len = p->seq[id].end = strtol(q, &r, 10);
					p->seq[id].start = 0;
				} else id = kh_val(p->h, itr);
			} else if (k == 2) { // start
				start = strtol(q, &r, 10);
			} else if (k == 3) { // end
				end = strtol(q, &r, 10);
				if (end - start < opt->min_len) break;
			} else if (k == 4) { // strand
				rev = (*q == '-');
			} else if (k == 9) { // match length
				pas_reg_t *t;
				score = strtol(q, &r, 10);
				if (score < opt->min_match) break; // do nothing
				if ((float)score / (end - start) < opt->min_match_frac) break;
				if (p->n_reg == p->m_reg) {
					p->m_reg = p->m_reg? p->m_reg<<1 : 16;
					p->reg = (pas_reg_t*)realloc(p->reg, sizeof(pas_reg_t) * p->m_reg);
				}
				t = &p->reg[p->n_reg++];
				t->id_start = (uint64_t)id<<32 | start;
				t->end = end, t->match = score, t->rev = rev;
			}
			++k, q = &s->s[i+1];
		}
	}
}

pas_t *pas_read(const char *fn, const pas_opt_t *opt)
{
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	pas_t *p;
	int dret;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	p = (pas_t*)calloc(1, sizeof(pas_t));
	p->h = kh_init(str);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0)
		pas_process1(opt, p, &str);
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	radix_sort_reg(p->reg, p->reg + p->n_reg);
	return p;
}

int main(int argc, char *argv[])
{
	int c;
	pas_opt_t opt;
	pas_t *p;
	char *r;

	pas_opt_init(&opt);
	while ((c = getopt(argc, argv, "m:l:")) >= 0) {
		if (c == 'm') {
			opt.min_match = strtol(optarg, &r, 10);
			if (*r == ',') opt.min_match_frac = strtod(r, &r);
		} else if (c == 'l') opt.min_len = atoi(optarg);
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: pa-sel [options] <in.pmf>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -m INT[,FLOAT]   min match length and fraction [%d,%.2f]\n", opt.min_match, opt.min_match_frac);
		fprintf(stderr, "  -l INT           min length [%d]\n", opt.min_len);
		return 1;
	}
	p = pas_read(argv[optind], &opt);
	fprintf(stderr, "[M::%s] read %d sequences and %ld hits\n", __func__, p->n_seq, (long)p->n_reg);
	pas_destroy(p);
	return 0;
}
