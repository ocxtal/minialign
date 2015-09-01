#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "kvec.h"
#include "minimap.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

uint64_t *mm_sketch(const char *str, int len, int w, int k, int *n)
{
	uint64_t shift1 = 2 * (k - 1), shift = 2 * k, mask = (1ULL<<shift) - 1;
	int i, j, l, buf_pos, n_min;
	uint64_t *buf, *min, kmer[2] = {0,0};
	uint64_v a = {0,0,0};

	assert(len <= UINT64_MAX>>2*k);
	buf = (uint64_t*)alloca(w * 8);
	memset(buf, 0xff, w * 8);
	min = (uint64_t*)alloca(w * 8);
	memset(min, 0xff, w * 8);

	for (i = l = n_min = buf_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		uint64_t info = UINT64_MAX;
		if (c < 4) { // not an ambiguous base
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (++l >= k)
				info = (uint64_t)i << shift | hash64(kmer[(kmer[0] > kmer[1])], mask); // hash the smaller k-mer
		} else l = 0;
//		printf("X\t%d\t%d\t%d\t%ld\t%llx\t%llx\n", i, l, n_min, a.n, info&mask, min[0]&mask);
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (n_min == 0) { // directly write to the min buffer if it is empty
			if (info != UINT64_MAX) min[n_min++] = info;
		} else if ((info&mask) < (min[0]&mask)) { // a new minimum; then write the old min(s)
			if (l >= w + k)
				for (j = 0; j < n_min; ++j) kv_push(uint64_t, a, min[j]);
			min[0] = info, n_min = 1;
		} else if ((info&mask) == (min[0]&mask)) {
			min[n_min++] = info;
		} else if (buf_pos == (min[n_min-1]>>shift)%w) { // all min(s) have moved outside the window
			uint64_t m;
			if (l >= w + k) // write the old min(s)
				for (j = 0; j < n_min; ++j) kv_push(uint64_t, a, min[j]);
			n_min = 0;
			for (j = 0, m = UINT64_MAX; j < w; ++j) // find the new min
				if (m > (buf[j]&mask)) m = buf[j]&mask;
			if (m != UINT64_MAX)
				for (j = 0; j < w; ++j) 
					if (m == (buf[j]&mask)) min[n_min++] = buf[j];
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	for (j = 0; j < n_min; ++j) kv_push(uint64_t, a, min[j]);
	*n = a.n;
	return a.a;
}
