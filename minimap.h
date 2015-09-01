#ifndef MINIMAP_H
#define MINIMAP_H

#include <stdint.h>

typedef struct { size_t n, m; uint64_t *a; } uint64_v;

#ifdef __cplusplus
extern "C" {
#endif

uint64_t *mm_sketch(const char *str, int len, int w, int k, int *n);

#ifdef __cplusplus
}
#endif

#endif
