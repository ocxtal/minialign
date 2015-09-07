#include <stdlib.h>
#include "minimap.h"
#include "kvec.h"
#include "khash.h"

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

mm_idx_t *mm_idx_init(int w, int k, int b)
{
	mm_idx_t *mi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	return mi;
}

void mm_idx_destroy(mm_idx_t *mi)
{
	int i;
	if (mi == 0) return;
	for (i = 0; i < 1<<mi->b; ++i) {
		free(mi->B[i].a.a);
		kh_destroy(idx, mi->B[i].h);
	}
	free(mi->B); free(mi);
}

void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
	int i, mask = (1<<mi->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &mi->B[a[i].x&mask].a;
		kv_push(mm128_t, *p, a[i]);
	}
}

#include "ksort.h"
#define sort_key(a) ((a).x)
KRADIX_SORT_INIT(128, mm128_t, sort_key, 8) 

static inline void insert_key(idxhash_t *h, const mm128_t *p, int start, int n, int b)
{
	int absent;
	khint_t k;
	k = kh_put(idx, h, p->x>>b<<1, &absent);
	if (n == 1) {
		kh_key(h, k) |= 1;
		kh_val(h, k) = p->y;
	} else kh_val(h, k) = (uint64_t)start<<32 | n;
}

static void worker_post(void *g, long i, int tid)
{
	int j, start, n;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)g;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;
	radix_sort_128(b->a.a, b->a.a + b->a.n);
	for (j = 1, n = 0; j < b->a.n; ++j) // count the number of keys to preallocate the hash table
		if (b->a.a[j].x != b->a.a[j-1].x) ++n;
	h = kh_init(idx);
	kh_resize(idx, h, n + 1);
	for (j = 1, n = 1, start = 0; j < b->a.n; ++j) {
		if (b->a.a[j].x != b->a.a[j-1].x) {
			insert_key(h, &b->a.a[j-1], start, n, mi->b);
			start = j, n = 1;
		} else ++n;
	}
	insert_key(h, &b->a.a[j-1], start, n, mi->b);
	b->h = h;
}
 
void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1<<mi->b);
}
