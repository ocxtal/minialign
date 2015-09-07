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
		free(mi->B[i].p);
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

const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n)
{
	int mask = (1<<mi->b) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	k = kh_get(idx, h, minier>>mi->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) {
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

#include "ksort.h"
#define sort_key(a) ((a).x)
KRADIX_SORT_INIT(128, mm128_t, sort_key, 8) 

static void worker_post(void *g, long i, int tid)
{
	int j, start_a, start_p, n, n_keys;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)g;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer
	radix_sort_128(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t*)calloc(b->n, 8);

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>mi->b<<1, &absent);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += k;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;

	// deallocate and clear b->a
	free(b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}
 
void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1<<mi->b);
}
