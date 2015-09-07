#ifndef MINIMAP_H
#define MINIMAP_H

#include <stdint.h>

typedef struct {
 	uint64_t x, y;
} mm128_t;

typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;

#ifdef __cplusplus
extern "C" {
#endif

void mm_sketch(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p);

#ifdef __cplusplus
}
#endif

#endif
