#ifndef MINIMAP_H
#define MINIMAP_H

#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>

typedef struct {
 	uint64_t x, y;
} mm128_t;

typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; uint32_t *a; } uint32_v;
typedef struct { size_t n, m; uint8_t *a; } uint8_v;
typedef struct { size_t n, m; void **a; } ptr_v;

typedef struct {
	int w, k;
	float f;
	uint32_t n;  // number of reference sequences
	uint32_t max_occ;
	uint64_t *h, *m;
	bseq_v b;

	// work
	mm128_v a;
	uint64_v size;
	ptr_v base;
} mm_idx_t;

typedef struct {
	double lbc, ubc;	// lower/upper bound coefficiencies
	int max_gap, min_cnt;
	int m, x, gi, ge, xdrop, min_score, ofs_lim;
	double min_ratio;
} mm_mapopt_t;

extern int mm_verbose;
extern double mm_realtime0;

#if 0
struct mm_tbuf_s;
typedef struct mm_tbuf_s mm_tbuf_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

// compute minimizers
void mm_sketch(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p);

// minimizer indexing
mm_idx_t *mm_idx_init(int w, int k, float f);
void mm_idx_destroy(mm_idx_t *mi);
mm_idx_t *mm_idx_gen(bseq_file_t *fp, int w, int k, float f, int batch_size, int n_threads);
const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n);
mm_idx_t *mm_idx_build(const char *fn, int w, int k, float f, int n_threads);

// minimizer index I/O
void mm_idx_dump(FILE *fp, const mm_idx_t *mi);
mm_idx_t *mm_idx_load(FILE *fp);

// mapping
void mm_mapopt_init(mm_mapopt_t *opt);
#if 0
mm_tbuf_t *mm_tbuf_init(const mm_pipline_t *pl);
void mm_tbuf_destroy(mm_tbuf_t *b);
mm_reg_t *mm_align(mm_tbuf_t *b, const mm_idx_t *mi, int l_seq, const char *seq, const mm_mapopt_t *opt, uint64_t *n_reg);
#endif

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads, int tbatch_size);

// private functions (may be moved to a "mmpriv.h" in future)
double cputime(void);
double realtime(void);
void radix_sort_128x(mm128_t *beg, mm128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);

#ifdef __cplusplus
}
#endif

#endif
