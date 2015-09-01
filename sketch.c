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
	int i, j, l, buf_pos;
	uint64_t *buf, min, kmer[2] = {0,0};
	uint64_v a = {0,0,0};

	assert(len <= UINT64_MAX>>2*k);
	buf = (uint64_t*)alloca(w * 8);
	memset(buf, 0xff, w * 8);

	for (i = l = buf_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		uint64_t info = UINT64_MAX;
		if (c < 4) { // not an ambiguous base
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (++l >= k)
				info = (uint64_t)i << shift | hash64(kmer[(kmer[0] > kmer[1])], mask); // hash the smaller k-mer
		} else l = 0;
		printf("X\t%d\t%d\t%ld\t%llx\t(%lld,%llx)\n", i, l, a.n, info&mask, min>>shift, min&mask);
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if ((info&mask) <= (min&mask)) { // a new minimum; then write the old min
			if (l >= w + k || (l == w + k - 1 && (info&mask) == (min&mask)))
				kv_push(uint64_t, a, min);
			min = info;
		} else if (buf_pos == (min>>shift)%w) { // all min(s) have moved outside the window
			if (l >= w + k - 1) kv_push(uint64_t, a, min);
			for (j = buf_pos + 1, min = UINT64_MAX; j < w; ++j)
				if ((min&mask) >= (buf[j]&mask)) min = buf[j];
			for (j = 0; j <= buf_pos; ++j)
				if ((min&mask) >= (buf[j]&mask)) min = buf[j];
			if (l >= w + k - 1) {
				for (j = buf_pos + 1; j < w; ++j)
					if ((min&mask) == (buf[j]&mask) && min != buf[j])
						kv_push(uint64_t, a, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if ((min&mask) == (buf[j]&mask) && min != buf[j])
						kv_push(uint64_t, a, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	kv_push(uint64_t, a, min);
	*n = a.n;
	return a.a;
}
