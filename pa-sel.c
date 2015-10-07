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
} pas_opt_t;

typedef struct {
	uint64_t id_start;
	uint32_t end, match:31, rev:1;
} pas_reg_t;

typedef struct {
	int n_seq, m_seq;
	char **name;
	int *len;
	shash_t *h;
	uint64_t n_reg, m_reg;
	pas_reg_t *reg;
} pas_t;

void pas_opt_init(pas_opt_t *opt)
{
	opt->min_match = 40;
}

void pas_destroy(pas_t *p)
{
	uint64_t i;
	if (p == 0) return;
	for (i = 0; i < p->n_seq; ++i) free(p->name[i]);
	kh_destroy(str, p->h);
	free(p->name); free(p->len); free(p->reg);
	free(p);
}

void pas_process1(const pas_opt_t *opt, pas_t *p, kstring_t *s)
{
	char *q, *r;
	int i, k, id = -1, start = -1, end = -1, rev = -1, score = -1;
	for (i = k = 0, q = s->s; i <= s->l; ++i) {
		if (i == s->l || s->s[i] == '\t') {
			s->s[i] = 0;
			if (k == 0) { // name
				khint_t itr;
				int absent;
				itr = kh_put(str, p->h, q, &absent);
				if (absent) {
					if (p->n_seq == p->m_seq) {
						p->m_seq = p->m_seq? p->m_seq<<1 : 16;
						p->name = (char**)realloc(p->name, sizeof(char*) * p->m_seq);
						p->len = (int*)realloc(p->len, sizeof(int) * p->m_seq);
					}
					id = p->n_seq++;
					kh_key(p->h, itr) = p->name[id] = strdup(q);
				} else id = kh_val(p->h, itr);
			} else if (k == 1) { // length
				p->len[id] = strtol(q, &r, 10);
			} else if (k == 2) { // start
				start = strtol(q, &r, 10);
			} else if (k == 3) { // end
				end = strtol(q, &r, 10);
			} else if (k == 4) { // strand
				rev = (*q == '-');
			} else if (k == 9) { // match length
				pas_reg_t *t;
				score = strtol(q, &r, 10);
				if (score < opt->min_match) break; // do nothing
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
	ks_destroy(ks);
	gzclose(fp);
	return p;
}

int main(int argc, char *argv[])
{
	int c;
	pas_opt_t opt;
	pas_t *p;

	pas_opt_init(&opt);
	while ((c = getopt(argc, argv, "m:")) >= 0) {
		if (c == 'm') opt.min_match = atoi(optarg);
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: pa-sel [options] <in.pmf>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -m INT     min match length [%d]\n", opt.min_match);
		return 1;
	}
	p = pas_read(argv[optind], &opt);
	fprintf(stderr, "[M::%s] read %d sequences and %ld hits\n", __func__, p->n_seq, (long)p->n_reg);
	pas_destroy(p);
	return 0;
}
