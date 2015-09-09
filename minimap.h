#ifndef MINIMAP_H
#define MINIMAP_H

#include <stdint.h>

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
	mm_idx_bucket_t *B;
} mm_idx_t;

extern int mm_verbose;
extern double mm_realtime0;

#ifdef __cplusplus
extern "C" {
#endif

void mm_sketch(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p);

mm_idx_t *mm_idx_init(int w, int k, int b);
void mm_idx_destroy(mm_idx_t *mi);
mm_idx_t *mm_idx_gen(const char *fn, int w, int k, int b, int batch_size, int n_threads);
uint32_t mm_idx_thres(const mm_idx_t *mi, float f);
const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n);

int mm_map(const mm_idx_t *idx, const char *fn, int radius, int min_cnt, float f, int n_threads, int batch_size);

double cputime(void);
double realtime(void);
void radix_sort_128x(mm128_t *beg, mm128_t *end);

#ifdef __cplusplus
}
#endif

#endif
