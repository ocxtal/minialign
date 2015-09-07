#ifndef MINIMAP_H
#define MINIMAP_H

#include <stdint.h>

typedef struct {
 	uint64_t x, y;
} mm128_t;

typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;

typedef struct {
	int b, w, k;
	mm128_v *a;
} mm_idx_t;

#ifdef __cplusplus
extern "C" {
#endif

void mm_sketch(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p);

mm_idx_t *mm_idx_init(int w, int k, int b);
void mm_idx_destroy(mm_idx_t *mi);
void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a);
void mm_idx_post(mm_idx_t *mi, int n_threads);

#ifdef __cplusplus
}
#endif

#endif
