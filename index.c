#include <stdlib.h>
#include "minimap.h"
#include "kvec.h"

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

mm_idx_t *mm_idx_init(int w, int k, int b)
{
	mm_idx_t *mi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b;
	mi->a = (mm128_v*)calloc(1<<b, sizeof(mm128_v));
	return mi;
}

void mm_idx_destroy(mm_idx_t *mi)
{
	int i;
	if (mi == 0) return;
	for (i = 0; i < 1<<mi->b; ++i)
		free(mi->a[i].a);
	free(mi->a); free(mi);
}

void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
	int i, mask = (1<<mi->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &mi->a[a[i].x&mask];
		kv_push(mm128_t, *p, a[i]);
	}
}

#include "ksort.h"
#define sort_key(a) ((a).x)
KRADIX_SORT_INIT(128, mm128_t, sort_key, 8) 

static void worker_post(void *g, long i, int tid)
{
	mm_idx_t *mi = (mm_idx_t*)g;
	radix_sort_128(mi->a[i].a, mi->a[i].a + mi->a[i].n);
}
 
void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1<<mi->b);
}
