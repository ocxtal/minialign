#ifndef MINIMAP_H
#define MINIMAP_H

#include <stdint.h>
#include <sys/types.h>

typedef struct {
 	uint64_t x, y;
} mm128_t;

typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;

typedef struct {
	mm128_v a;   // (minimizer, position) array
	int n;       // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

typedef struct {
	int b, w, k;
	uint32_t n;  // number of reference sequences
	mm_idx_bucket_t *B;
	int *len;    // length of each reference sequence
	char **name; // TODO: if this uses too much RAM, switch one concatenated string
} mm_idx_t;

typedef struct {
	uint32_t cnt:31, rev:1;
	int32_t rid;
	int32_t qs, qe, rs, re;
} mm_reg1_t;

extern int mm_verbose;
extern double mm_realtime0;

struct mm_tbuf_s;
typedef struct mm_tbuf_s mm_tbuf_t;

#ifdef __cplusplus
extern "C" {
#endif

void mm_sketch(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p);

mm_idx_t *mm_idx_init(int w, int k, int b);
void mm_idx_destroy(mm_idx_t *mi);
mm_idx_t *mm_idx_gen(const char *fn, int w, int k, int b, int batch_size, int n_threads, int keep_name);
uint32_t mm_idx_thres(const mm_idx_t *mi, float f);
const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n);

const mm_reg1_t *mm_map(const mm_idx_t *mi, int l_seq, const char *seq, int *n_regs, mm_tbuf_t *b, int max_occ, int radius, int min_cnt, int max_gap);
int mm_map_file(const mm_idx_t *idx, const char *fn, int radius, int max_gap, int min_cnt, float f, int n_threads, int batch_size);

double cputime(void);
double realtime(void);
void radix_sort_128x(mm128_t *beg, mm128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);

#ifdef __cplusplus
}
#endif

#endif
