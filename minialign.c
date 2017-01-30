#include <getopt.h>
#include <unistd.h>
#include <zlib.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/types.h>

#include "kvec.h"
#include "ptask.h"
#include "gaba.h"
#include "lmm.h"
#include "sassert.h"

#define MM_VERSION "0.4.4"


#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

/* minimap.h */
typedef struct { uint32_t x[2]; } v2u32_t;
typedef struct { uint32_t x[4]; } v4u32_t;
typedef struct { uint64_t x[2]; } v2u64_t;
_static_assert(sizeof(v2u32_t) == 8);
_static_assert(sizeof(v4u32_t) == 16);
_static_assert(sizeof(v2u64_t) == 16);

typedef union {
	uint64_t u64[2];
	uint32_t u32[4];
} mm128_t;
_static_assert(sizeof(mm128_t) == 16);

typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; v2u32_t *a; } v2u32_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; uint32_t *a; } uint32_v;
typedef struct { size_t n, m; uint16_t *a; } uint16_v;
typedef struct { size_t n, m; uint8_t *a; } uint8_v;
typedef struct { size_t n, m; void **a; } ptr_v;

#include "ksort.h"
#define sort_key_64(a) ((a))
KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 4)
#define sort_key_64x(a) ((a).x[0])
KRADIX_SORT_INIT(64x, v2u32_t, sort_key_64x, 4)
#define sort_key_128x(a) ((a).u64[0])
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8) 
KSORT_INIT_GENERIC(uint32_t)

#include "khash.h"
#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;
#define pos_hash(a)	(a)
#define pos_eq(a, b)	((a)==(b))
KHASH_INIT(pos, uint64_t, uint64_t, 1, pos_hash, pos_eq)
typedef khash_t(pos) poshash_t;

/* misc.c */

int mm_verbose = 3;
double mm_realtime0;

static double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

static double realtime()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

/* end of misc.c */

/* bamlite.c in bwa */
typedef struct {
	int32_t n_targets;
	char **target_name;
	uint32_t *target_len;
	size_t l_text, n_text;
	char *text;
} bam_header_t;

typedef struct {
	int32_t tid;
	int32_t pos;
	uint8_t l_qname, qual;
	uint16_t bin, n_cigar, flag;
	int32_t l_qseq;
	int32_t mtid;
	int32_t mpos;
	int32_t isize;
} bam_core_t;
_static_assert(sizeof(bam_core_t) == 32);

static void bam_header_destroy(bam_header_t *header)
{
	if (header == 0) return;
	if (header->target_name) {
		for (uint64_t i = 0; i < (uint64_t)header->n_targets; ++i)
			if (header->target_name[i]) free(header->target_name[i]);
		if (header->target_len) free(header->target_len);
		free(header->target_name);
	}
	if (header->text) free(header->text);
	free(header);
}

static bam_header_t *bam_read_header(gzFile fp)
{
	bam_header_t *header;
	char buf[4];
	int magic_len;
	// read "BAM1"
	magic_len = gzread(fp, buf, 4);
	if (magic_len != 4 || strncmp(buf, "BAM\001", 4) != 0) {
		// fprintf(stderr, "[bam_read_header] invalid BAM binary header (this is not a BAM file).\n");
		return NULL;
	}
	header = calloc(1, sizeof(bam_header_t));
	// read plain text and the number of reference sequences
	if (gzread(fp, &header->l_text, 4) != 4) goto _bam_read_header_fail;
	header->text = (char*)calloc(header->l_text + 1, 1);
	if ((size_t)gzread(fp, header->text, header->l_text) != header->l_text) goto _bam_read_header_fail;
	if (gzread(fp, &header->n_targets, 4) != 4) goto _bam_read_header_fail;
	// read reference sequence names and lengths
	header->target_name = (char**)calloc(header->n_targets, sizeof(char*));
	header->target_len = (uint32_t*)calloc(header->n_targets, 4);
	for (uint64_t i = 0; i != (uint64_t)header->n_targets; ++i) {
		int32_t name_len;
		if (gzread(fp, &name_len, 4) != 4) goto _bam_read_header_fail;
		header->target_name[i] = (char*)calloc(name_len, 1);
		if (gzread(fp, header->target_name[i], name_len) != name_len) {
			goto _bam_read_header_fail;
		}
		if (gzread(fp, &header->target_len[i], 4) != 4) goto _bam_read_header_fail;
	}
	return header;
_bam_read_header_fail:
	bam_header_destroy(header);
	return NULL;
}

/* end of bamlite.c */

/* bseq.c */
typedef struct {
	uint32_t is_eof, base_rid, keep_qual, keep_tag;
	gzFile fp;
	kseq_t *ks;
	bam_header_t *bh;
	uint8_v buf;
} bseq_file_t;

// sam optional tags
typedef struct {
	uint32_t size;	// size of data array
	char tag[2], type;
	uint8_t flag;	// 1 if only contained in primary
	uint8_t data[];
} bseq_tag_t;

typedef struct {
	// lengths
	uint32_t l_seq, rid, l_name, l_tag;
	// pointers
	char *name;
	uint8_t *seq, *qual, *tag;
} bseq_t;
_static_assert(sizeof(bseq_t) == 48);
typedef struct { size_t n, m; bseq_t *a; } bseq_v;

// ascii to 4bit conversion
static unsigned char seq_nt4_table_4bit[32] = {
	0, 1, 0, 2,  0, 0, 0, 4,  0, 0, 0, 0,  0, 0, 0, 0,
	0, 0, 0, 0,  8, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
};

static bseq_file_t *bseq_open(const char *fn, uint32_t base_rid, uint32_t keep_qual, uint32_t keep_tag)
{
	int c;
	bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (f == 0) return 0;
	fp = (bseq_file_t*)calloc(1, sizeof(bseq_file_t));
	fp->base_rid = base_rid; fp->keep_qual = keep_qual; fp->keep_tag = keep_tag;
	fp->fp = f;
	if ((c = gzgetc(fp->fp)) == 'B')	// test bam signature
		gzungetc(c, fp->fp), fp->bh = bam_read_header(fp->fp);
	else if (c == '>' || c == '@')
		gzungetc(c, fp->fp), fp->ks = kseq_init(fp->fp);
	else free(fp), fp = 0;
	return fp;
}

static uint32_t bseq_close(bseq_file_t *fp)
{
	if (fp == 0) return 0;
	uint32_t base_rid = fp->base_rid;
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
	return base_rid;
}

static void bseq_read_bam(bseq_file_t *fp, uint64_t size, bseq_v *seq, uint8_v *mem)
{
	// parse bam
	bseq_t *s;
	uint8_t *sname, *sseq, *squal, *stag;
	uint32_t l_tag;
	bam_core_t *c;

	if (fp->buf.n < size) kv_reserve(uint8_t, fp->buf, size);
	gzread(fp->fp, fp->buf.a, size);
	if ((c = (bam_core_t*)fp->buf.a)->flag&0x900) return;	// skip supp/secondary
	sname = fp->buf.a + sizeof(bam_core_t);
	sseq = sname + c->l_qname + sizeof(uint32_t)*c->n_cigar;
	squal = sseq + (c->l_qseq+1) / 2;
	stag = squal + c->l_qseq; l_tag = ((uint8_t*)fp->buf.a + size) - stag;

	kv_pushp(bseq_t, *seq, &s);
	kv_reserve(uint8_t, *mem, mem->n + c->l_qname + c->l_qseq + ((fp->keep_qual && *squal != 0xff)? c->l_qseq : 0) + (fp->keep_tag? l_tag : 0) + 3);
	s->l_seq = c->l_qseq;
	s->rid = seq->n + fp->base_rid - 1;
	s->l_name = c->l_qname - 1;		// remove tail '\0'
	s->l_tag = fp->keep_tag? l_tag : 0;

	s->name = (char *)mem->n;
	memcpy(mem->a + mem->n, sname, c->l_qname); mem->n += c->l_qname;

	s->seq = (uint8_t *)mem->n;
	if (c->flag&0x10) {
		for (int64_t i = c->l_qseq-1; i >= 0; --i) mem->a[mem->n++] = "\x0\x8\x4\x0\x2\x0\x0\x0\x1"[0x0f & (sseq[i>>1]>>((~i&0x01)<<2))];
	} else {
		for (int64_t i = 0; i < c->l_qseq; ++i) mem->a[mem->n++] = 0x0f & (sseq[i>>1]>>((~i&0x01)<<2));
	}
	mem->a[mem->n++] = '\0';

	s->qual = (uint8_t *)mem->n;
	if (fp->keep_qual && *squal != 0xff) {
		if (c->flag&0x10) {
			for (int64_t i = c->l_qseq-1; i >= 0; --i) mem->a[mem->n++] = squal[i] + 33;
		} else {
			for (int64_t i = 0; i < c->l_qseq; ++i) mem->a[mem->n++] = squal[i] + 33;
		}
	}
	mem->a[mem->n++] = '\0';

	s->tag = (uint8_t *)mem->n;
	if (fp->keep_tag && l_tag) { memcpy(mem->a + mem->n, stag, l_tag); mem->n += l_tag; }
	mem->a[mem->n++] = '\0';
	return;
}

static void bseq_read_fasta(bseq_file_t *fp, bseq_v *seq, uint8_v *mem)
{
	bseq_t *s;
	kseq_t *ks = fp->ks;

	kv_pushp(bseq_t, *seq, &s);
	kv_reserve(uint8_t, *mem, mem->n + ks->name.l + ks->seq.l + ((fp->keep_qual && ks->qual.l)? ks->seq.l : 0) + 4);
	s->l_seq = ks->seq.l;
	s->rid = seq->n + fp->base_rid - 1;
	s->l_name = ks->name.l;
	s->l_tag = 0;	// no tags available in fasta/q

	s->name = (char *)mem->n;
	memcpy(mem->a + mem->n, ks->name.s, ks->name.l); mem->n += ks->name.l;
	mem->a[mem->n++] = '\0';

	s->seq = (uint8_t *)mem->n;
	for (uint64_t i = 0; i < ks->seq.l; ++i) mem->a[mem->n++] = seq_nt4_table_4bit[0x1f&ks->seq.s[i]];
	mem->a[mem->n++] = '\0';

	s->qual = (uint8_t *)mem->n;
	if (fp->keep_qual && ks->qual.l) { memcpy(mem->a + mem->n, ks->qual.s, ks->qual.l); mem->n += ks->qual.l; }
	mem->a[mem->n++] = '\0';

	s->tag = (uint8_t *)mem->n;
	mem->a[mem->n++] = '\0';
	return;
}

static bseq_t *bseq_read(bseq_file_t *fp, uint64_t chunk_size, uint32_t *n, void **base, uint64_t *size)
{
	uint8_v mem = {0};
	bseq_v seq = {0};
	static const uint8_t margin[64] = {0};

	kv_reserve(uint8_t, mem, chunk_size + 128);
	kv_pushm(uint8_t, mem, margin, 64);
	if (fp->bh) {
		uint32_t size;
		while (gzread(fp->fp, &size, 4) == 4) {
			bseq_read_bam(fp, size, &seq, &mem);
			if (mem.n >= chunk_size) break;
		}
	} else {
		while (kseq_read(fp->ks) >= 0) {
			bseq_read_fasta(fp, &seq, &mem);
			if (mem.n >= chunk_size) break;
		}
	}
	kv_pushm(uint8_t, mem, margin, 64);

	for (uint64_t i = 0; i < seq.n; i++) {
		seq.a[i].name += (ptrdiff_t)mem.a;
		seq.a[i].seq += (ptrdiff_t)mem.a;
		seq.a[i].qual += (ptrdiff_t)mem.a;
		seq.a[i].tag += (ptrdiff_t)mem.a;
	}
	if (seq.n == 0) free(mem.a), mem.a = 0, fp->is_eof = 1;
	fp->base_rid += seq.n;
	*n = seq.n;
	*base = (void*)mem.a;
	*size = mem.n;
	return seq.a;
}

static int bseq_eof(bseq_file_t *fp)
{
	return fp->is_eof;
}

/* end of bseq.c */

/* sketch.c */

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

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param p      minimizers; p->a[i].x is the 2k-bit hash value;
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
static void mm_sketch(const uint8_t *seq4, uint32_t len, uint32_t w, uint32_t k, uint32_t rid, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, k0 = 0, k1 = 0;
	uint32_t i, j, l, buf_pos, min_pos;
	mm128_t buf[w], min = { .u64 = { UINT64_MAX, UINT64_MAX } };

	assert(len > 0 && w > 0 && k > 0);
	memset(buf, 0xff, w * sizeof(mm128_t));

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		uint64_t _c = seq4[i];
		mm128_t info = { .u64 = { UINT64_MAX, UINT64_MAX } };
		if (_c != 0) { // not an ambiguous base
			uint64_t c = 0x03 & ((_c>>1) - (_c>>3));
			k0 = (k0 << 2 | c) & mask;           // forward k-mer
			k1 = (k1 >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (k0 == k1) continue; // skip "symmetric k-mers" as we don't know it strand
			if (++l >= k) {
				info.u64[0] = hash64(k0 < k1? k0 : k1, mask);
				info.u32[2] = k0 < k1? i : ~i; info.u32[3] = rid;
			}
		} else l = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.u64[0] == buf[j].u64[0] && buf[j].u64[1] != min.u64[1]) kv_push(mm128_t, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.u64[0] == buf[j].u64[0] && buf[j].u64[1] != min.u64[1]) kv_push(mm128_t, *p, buf[j]);
		}
		if (info.u64[0] <= min.u64[0]) { // a new minimum; then write the old min
			if (l >= w + k) kv_push(mm128_t, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1) kv_push(mm128_t, *p, min);
			for (j = buf_pos + 1, min.u64[0] = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.u64[0] >= buf[j].u64[0]) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.u64[0] >= buf[j].u64[0]) min = buf[j], min_pos = j;
			if (l >= w + k - 1) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.u64[0] == buf[j].u64[0] && min.u64[1] != buf[j].u64[1]) kv_push(mm128_t, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.u64[0] == buf[j].u64[0] && min.u64[1] != buf[j].u64[1]) kv_push(mm128_t, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.u64[0] != UINT64_MAX)
		kv_push(mm128_t, *p, min);
}

/* end of sketch.c */

/* map.c, options */
#define MM_RG			( 0 )
#define MM_AS			( 1 )
#define MM_XS			( 2 )
#define MM_SA			( 3 )
#define MM_MD			( 4 )

#define MM_AVA			( 0x01ULL<<56 )
#define MM_KEEP_QUAL	( 0x02ULL<<56 )
#define MM_CIRCULAR		( 0x04ULL<<56 )

typedef struct {
	uint32_t sidx, eidx, n_threads, min, k, w, b;
	uint64_t flag;
	double min_ratio;
	int32_t m, x, gi, ge, xdrop, llim, hlim, elim, blim;
	uint32_t n_frq;
	float frq[16];
	uint16_v tags;
	uint64_t batch_size, outbuf_size;
	char *rg_line, *rg_id;
} mm_mapopt_t;

static void mm_mapopt_destroy(mm_mapopt_t *opt)
{
	free(opt->rg_line); free(opt->rg_id); free(opt);
	return;
}

static mm_mapopt_t *mm_mapopt_init(void)
{
	mm_mapopt_t *opt = calloc(1, sizeof(mm_mapopt_t));
	opt->k = 15;
	opt->w = 16;
	opt->b = 14;
	opt->flag = 0;
	opt->sidx = 0;
	opt->eidx = 3;
	opt->hlim = 5000;
	opt->llim = 5000;
	opt->blim = 0;
	opt->elim = 200;
	opt->m = 1;
	opt->x = 1;
	opt->gi = 1;
	opt->ge = 1;
	opt->xdrop = 50;
	opt->min = 50;
	opt->min_ratio = 0.3;
	opt->frq[0] = 0.05, opt->frq[1] = 0.01, opt->frq[2] = 0.001;
	opt->n_frq = 3;
	opt->n_threads = 1;
	opt->batch_size = 1024 * 1024;
	opt->outbuf_size = 512 * 1024;
	opt->rg_line = opt->rg_id = NULL;
	return opt;
}

static int mm_mapopt_check(mm_mapopt_t *opt, int (*_fprintf)(FILE*,const char*,...), FILE *_fp)
{
	int ret = 0;
	if (opt->k >= 32) _fprintf(_fp, "[M::%s] ERROR: k must be inside [1,32).\n", __func__), ret = 1;
	if (ret) return(ret);
	if (opt->w >= 16) _fprintf(_fp, "[M::%s] ERROR: w must be inside [1,16).\n", __func__), ret = 1;
	if (ret) return(ret);

	if (opt->m < 1 || opt->m > 5) _fprintf(_fp, "[M::%s] ERROR: Match award must be inside [1,5].\n", __func__), ret = 1;
	if (opt->x < 1 || opt->x > 5) _fprintf(_fp, "[M::%s] ERROR: Mismatch penalty must be inside [1,5].\n", __func__), ret = 1;
	if (opt->gi < 1 || opt->gi > 5) _fprintf(_fp, "[M::%s] ERROR: Gap open penalty must be inside [1,5].\n", __func__), ret = 1;
	if (opt->ge < 1 || opt->ge > 5) _fprintf(_fp, "[M::%s] ERROR: Gap extension penalty must be inside [1,5].\n", __func__), ret = 1;
	if (ret) return(ret);

	if (opt->x >= (opt->gi + opt->ge)) _fprintf(_fp, "[M::%s] info: Large mismatch penalty with respect to the gap open/extend penalty may cause SEGV or broken CIGAR. [issue #2]\n", __func__);
	if (opt->m + 2*(opt->gi + opt->ge) > 10) _fprintf(_fp, "[M::%s] info: Large match award or large gap open/extend penalty may cause SEGV or broken CIGAR. [issue #7]\n", __func__);
	if (opt->xdrop < 10 || opt->xdrop > 100) _fprintf(_fp, "[M::%s] ERROR: Xdrop cutoff must be inside [10,100].\n", __func__), ret = 1;
	if (opt->min > INT32_MAX) _fprintf(_fp, "[M::%s] ERROR: Minimum alignment score must be > 0.\n", __func__), ret = 1;
	if (opt->min_ratio < 0.0 || opt->min_ratio > 1.0) _fprintf(_fp, "[M::%s] ERROR: Minimum alignment score ratio must be inside [0.0,1.0].\n", __func__), ret = 1;
	if (opt->n_frq >= 16) _fprintf(_fp, "[M::%s] ERROR: Frequency thresholds must be fewer than 16.\n", __func__), ret = 1;
	for (uint64_t i = 0; i < opt->n_frq; ++i) if (opt->frq[i] < 0.0 || opt->frq[i] > 1.0 || (i != 0 && opt->frq[i-1] < opt->frq[i])) _fprintf(_fp, "[M::%s] ERROR: Frequency thresholds must be inside [0.0,1.0] and descending.\n", __func__), ret = 1;
	if (ret) return(ret);

	if (opt->n_threads < 1) _fprintf(_fp, "[M::%s] ERROR: Thread counts must be > 0.\n", __func__), ret = 1;
	if (opt->batch_size < 64 * 1024) _fprintf(_fp, "[M::%s] ERROR: Batch size must be > 64k.\n", __func__), ret = 1;
	if (opt->outbuf_size < 64 * 1024) _fprintf(_fp, "[M::%s] ERROR: Output buffer size must be > 64k.\n", __func__), ret = 1;
	if (ret) return(ret);

	if (opt->sidx >= 16) _fprintf(_fp, "[M::%s] ERROR: sidx must be inside [0,16).\n", __func__), ret = 1;
	if (opt->eidx >= 16) _fprintf(_fp, "[M::%s] ERROR: eidx must be inside [0,16).\n", __func__), ret = 1;
	if (opt->hlim < 100 || opt->hlim >= 100000) _fprintf(_fp, "[M::%s] ERROR: hlim must be inside [100,100000).\n", __func__), ret = 1;
	if (opt->llim < 100 || opt->llim >= 100000) _fprintf(_fp, "[M::%s] ERROR: llim must be inside [100,100000).\n", __func__), ret = 1;
	if (opt->blim >= 10000) _fprintf(_fp, "[M::%s] ERROR: blim must be inside [0,10000).\n", __func__), ret = 1;
	if (opt->elim >= 10000) _fprintf(_fp, "[M::%s] ERROR: elim must be inside [0,10000).\n", __func__), ret = 1;
	if (ret) return(ret);

	#define _dup(x)	({ char *_p = malloc(strlen(x)+1); memcpy(_p, (x), strlen(x)); _p[strlen(x)] = '\0'; _p; })
	if (opt->flag&(0x01ULL<<MM_RG)) {
		if (opt->rg_line == NULL && opt->rg_id == NULL) {
			opt->rg_line = _dup("@RG\tID:1\tSM:default"); opt->rg_id = _dup("1");
		}
	}
	#undef _dup
	return ret;
}

/* end of map.c, options */

/* index.c */
typedef struct {
	// lengths and flag
	uint32_t l_seq, rid, l_name, circular;
	// pointers
	char *name;
	uint8_t *seq;
} mm_idx_seq_t;
_static_assert(sizeof(mm_idx_seq_t) == 32);
typedef struct { size_t n, m; mm_idx_seq_t *a; } mm_idx_seq_v;

typedef struct {
	mm128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

typedef struct {
	uint8_t w, k, b, circular;
	uint32_t base_rid;
	mm_idx_bucket_t *B;
	mm_idx_seq_v s;

	// work
	mm128_v a;
	uint64_v size;
	ptr_v base;
} mm_idx_t;

static mm_idx_t *mm_idx_init(uint32_t w, uint32_t k, uint32_t b)
{
	mm_idx_t *mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w<1? 1 : w, mi->k = k; mi->b = MIN2(k*2, b);
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	return mi;
}

static void mm_idx_destroy(mm_idx_t *mi)
{
	if (mi == 0) return;
	for (uint64_t i = 0; i < 1ULL<<mi->b; ++i) {
		free(mi->B[i].p);
		free(mi->B[i].a.a);
		kh_destroy(idx, (idxhash_t*)mi->B[i].h);
	}
	free(mi->B);
	for (uint64_t i = 0; i < mi->base.n; ++i) free(mi->base.a[i]);
	free(mi->base.a); free(mi->size.a);
	free(mi->s.a); free(mi->a.a); free(mi);
}

static const v2u32_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, uint64_t *n)
{
	uint64_t mask = (1<<mi->b) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>mi->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) {
		*n = 1;
		return (const v2u32_t*)&kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return (const v2u32_t*)&b->p[kh_val(h, k)>>32];
	}
}

static uint32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	uint64_t n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.) return UINT32_MAX;
	for (uint64_t i = (n = 0); i < 1ULL<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (uint64_t i = (n = 0); i < 1ULL<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

/******************
 * Generate index *
 ******************/

typedef struct {
	uint64_t batch_size;
	bseq_file_t *fp;
	mm_idx_t *mi;
} mm_idx_pipeline_t;

typedef struct {
    uint32_t n_seq;
	const bseq_t *seq;	// const!
	void *base;
	uint64_t size;
	mm128_v a;
} mm_idx_step_t;

static void *mm_idx_source(void *arg)
{
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t*)arg;
	mm_idx_step_t *s = (mm_idx_step_t*)calloc(1, sizeof(mm_idx_step_t));
	uint64_t size;
	void *base;
	s->seq = bseq_read(q->fp, q->batch_size, &s->n_seq, &base, &size);
	if (s->seq == 0) { free(s), s = 0; return 0; }

	kv_push(void*, q->mi->base, base);
	kv_push(uint64_t, q->mi->size, size);
	return s;
}

static void *mm_idx_worker(void *arg, void *item)
{
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t*)arg;
	mm_idx_step_t *s = (mm_idx_step_t*)item;

	for (uint64_t i = 0; i < s->n_seq; ++i)
		mm_sketch(s->seq[i].seq, s->seq[i].l_seq, q->mi->w, q->mi->k, s->seq[i].rid, &s->a);
	return s;
}

static void mm_idx_drain(void *arg, void *item)
{
	mm_idx_pipeline_t *q = (mm_idx_pipeline_t*)arg;
	mm_idx_step_t *s = (mm_idx_step_t*)item;
	uint64_t mask = (1<<q->mi->b) - 1;
	q->mi->s.n = MAX2(q->mi->s.n, s->seq[s->n_seq-1].rid+1-q->mi->base_rid);
	kv_reserve(bseq_t, q->mi->s, q->mi->s.n);
	for (uint64_t i = 0; i < s->n_seq; ++i) {
		const bseq_t *p = &s->seq[i];
		q->mi->s.a[s->seq[i].rid-q->mi->base_rid] = (mm_idx_seq_t){
			.l_seq = p->l_seq, .rid = p->rid, .l_name = p->l_name, .circular = q->mi->circular,
			.name = p->name, .seq = p->seq
		};
	}
	for (uint64_t i = 0; i < s->a.n; ++i) {
		mm128_v *p = &q->mi->B[s->a.a[i].u64[0]&mask].a;
		kv_push(mm128_t, *p, s->a.a[i]);
	}

	free((void*)s->seq); free(s->a.a); free(s);
	return;
}

typedef struct {
	mm_idx_t *mi;
	uint32_t from, to;
} mm_idx_post_t;

static void *mm_idx_post(void *arg, void *item)
{
	mm_idx_post_t *q = (mm_idx_post_t*)arg;
	mm_idx_t *mi = q->mi;

	for (uint64_t i = q->from; i < q->to; ++i) {
		uint64_t n_keys;
		idxhash_t *h;
		mm_idx_bucket_t *b = &mi->B[i];
		if (b->a.n == 0) continue;

		// sort by minimizer
		radix_sort_128x(b->a.a, b->a.a + b->a.n);

		// count and preallocate
		for (uint64_t j = 1, n = (n_keys = 0, b->n = 0, 1); j <= b->a.n; ++j) {
			if (j == b->a.n || b->a.a[j].u64[0] != b->a.a[j-1].u64[0]) {
				++n_keys;
				if (n > 1) b->n += n;
				n = 1;
			} else ++n;
		}
		h = kh_init(idx);
		kh_resize(idx, h, n_keys);
		b->p = (uint64_t*)calloc(b->n, 8);

		// create the hash table
		for (uint64_t j = 1, n = 1, start_a = 0, start_p = 0; j <= b->a.n; ++j) {
			if (j == b->a.n || b->a.a[j].u64[0] != b->a.a[j-1].u64[0]) {
				khint_t itr;
				int absent;
				mm128_t *p = &b->a.a[j-1];
				itr = kh_put(idx, h, p->u64[0]>>mi->b<<1, &absent);
				assert(absent && j - start_a == n);
				if (n == 1) {
					kh_key(h, itr) |= 1;
					kh_val(h, itr) = p->u64[1];
				} else {
					for (uint64_t k = 0; k < n; ++k)
						b->p[start_p + k] = b->a.a[start_a + k].u64[1];
					kh_val(h, itr) = (uint64_t)start_p<<32 | n;
					start_p += n;
				}
				start_a = j, n = 1;
			} else ++n;
		}
		b->h = h;
		// assert(b->n == start_p);

		// deallocate and clear b->a
		free(b->a.a);
		b->a.n = b->a.m = 0, b->a.a = 0;
	}
	return 0;
}

static mm_idx_t *mm_idx_gen(const mm_mapopt_t *opt, bseq_file_t *fp)
{
	mm_idx_pipeline_t pl = {0}, **p;
	pl.batch_size = opt->batch_size;
	pl.fp = fp;
	if (pl.fp == 0) return 0;
	pl.mi = mm_idx_init(opt->w, opt->k, opt->b); pl.mi->circular = (opt->flag&MM_CIRCULAR)!=0; pl.mi->base_rid = pl.fp->base_rid;

	p = (mm_idx_pipeline_t**)calloc(opt->n_threads, sizeof(mm_idx_pipeline_t*));
	for (uint64_t i = 0; i < opt->n_threads; ++i) p[i] = &pl;
	ptask_t *pt = ptask_init(mm_idx_worker, (void**)p, opt->n_threads, 64);
	ptask_stream(pt, mm_idx_source, &pl, mm_idx_drain, &pl, 8*opt->n_threads);
	ptask_clean(pt);
	free(p);
	
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	mm_idx_post_t *q = (mm_idx_post_t*)calloc(opt->n_threads, sizeof(mm_idx_post_t));
	mm_idx_post_t **qq = (mm_idx_post_t**)calloc(opt->n_threads, sizeof(mm_idx_post_t*));
	for (uint64_t i = 0; i < opt->n_threads; ++i) {
		q[i].mi = pl.mi;
		q[i].from = (1ULL<<pl.mi->b)*i/opt->n_threads;
		q[i].to = (1ULL<<pl.mi->b)*(i+1)/opt->n_threads;
		qq[i] = &q[i];
	}
	pt = ptask_init(mm_idx_post, (void**)qq, opt->n_threads, 64);
	ptask_parallel(pt, NULL, NULL);
	ptask_clean(pt);
	free(q); free(qq);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
	return pl.mi;
}

/*************
 * index I/O *
 *************/

#define MM_IDX_MAGIC "MAI\5"		/* minialign index version 5, differs from minimap index signature */

static void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint32_t x[3];
	uint64_t size = 0, y[2];
	for (uint64_t i = (size = 0); i < mi->size.n; ++i) size += mi->size.a[i];
	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b; y[0] = mi->s.n, y[1] = size;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 3, fp);
	fwrite(y, 8, 2, fp);
	for (uint64_t i = 0; i < 1ULL<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&b->n, 4, 1, fp);
		fwrite(b->p, 8, b->n, fp);
		fwrite(&size, 4, 1, fp);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, fp);
		}
	}
	for (uint64_t i = 0; i < mi->base.n; ++i)
		fwrite(mi->base.a[i], sizeof(char), mi->size.a[i], fp);
	for (uint64_t i = 0, j = 0, s = 0; i < mi->s.n; ++i) {
		if ((uintptr_t)mi->s.a[i].seq < (uintptr_t)mi->base.a[j]
		|| (uintptr_t)mi->s.a[i].seq >= (uintptr_t)mi->base.a[j] + (ptrdiff_t)mi->size.a[j]) {
			s += mi->size.a[j++];
		}
		mi->s.a[i].name -= (ptrdiff_t)mi->base.a[j], mi->s.a[i].seq -= (ptrdiff_t)mi->base.a[j];
		mi->s.a[i].name += (ptrdiff_t)s, mi->s.a[i].seq += (ptrdiff_t)s;
	}
	fwrite(mi->s.a, sizeof(mm_idx_seq_t), mi->s.n, fp);
	// restore pointers
	for (uint64_t i = 0, j = 0, s = 0; i < mi->s.n; ++i) {
		mi->s.a[i].name -= (ptrdiff_t)s, mi->s.a[i].seq -= (ptrdiff_t)s;
		mi->s.a[i].name += (ptrdiff_t)mi->base.a[j], mi->s.a[i].seq += (ptrdiff_t)mi->base.a[j];
		if ((uintptr_t)mi->s.a[i].name > s + mi->size.a[j]) s += mi->size.a[j++];
	}
	return;
}

static mm_idx_t *mm_idx_load(FILE *fp)
{
	char magic[4];
	uint32_t x[3];
	uint64_t bsize, y[2];

	mm_idx_t *mi;
	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (fread(x, 4, 3, fp) != 3) return 0;
	if (fread(y, 8, 2, fp) != 2) return 0;
	mi = mm_idx_init(x[0], x[1], x[2]); mi->s.n = y[0]; bsize = y[1];
	for (uint64_t i = 0; i < 1ULL<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t hsize;
		khint_t k;
		idxhash_t *h;
		if (fread(&b->n, 4, 1, fp) != 1) goto _mm_idx_load_fail;
		b->p = (uint64_t*)malloc(b->n * 8);
		if (fread(b->p, 8, b->n, fp) != (size_t)b->n) goto _mm_idx_load_fail;
		if (fread(&hsize, 4, 1, fp) != 1) goto _mm_idx_load_fail;
		if (hsize == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, hsize);
		for (uint64_t j = 0; j < hsize; ++j) {
			uint64_t x[2];
			int absent;
			if (fread(x, 8, 2, fp) != 2) goto _mm_idx_load_fail;
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}

	mi->base.n = mi->size.n = 1;
	mi->base.a = malloc(sizeof(void*) * mi->base.n);
	mi->size.a = malloc(sizeof(uint64_t) * mi->size.n);
	mi->base.a[0] = malloc(sizeof(char) * bsize);
	mi->size.a[0] = bsize;
	if (fread(mi->base.a[0], sizeof(char), mi->size.a[0], fp) != mi->size.a[0]) goto _mm_idx_load_fail;
	mi->s.a = malloc(sizeof(mm_idx_seq_t) * mi->s.n);
	if (fread(mi->s.a, sizeof(mm_idx_seq_t), mi->s.n, fp) != mi->s.n) goto _mm_idx_load_fail;
	for (uint64_t i = 0; i < mi->s.n; ++i) mi->s.a[i].name += (ptrdiff_t)mi->base.a[0], mi->s.a[i].seq += (ptrdiff_t)mi->base.a[0];
	mi->base_rid = mi->s.a[0].rid;
	return mi;
_mm_idx_load_fail:
	mm_idx_destroy(mi);
	return 0;
}

/* end of index.c */

/* map.c */
/**************************
 * Multi-threaded mapping *
 **************************/
typedef struct {
	int32_t rs, re, qs, qe;
	int32_t rid, score;
	uint32_t l_cigar;
	uint8_t cigar[];
} reg_t;

typedef struct {
	uint8_t *base, *tail, *p;
	const mm_idx_t *mi;
	const mm_mapopt_t *opt;
	uint64_t *occ;
	uint32_t n_occ;
	gaba_t *gaba;
	void **t;	// mm_tbuf_t **
	ptask_t *pt;
	bseq_file_t *fp;
} mm_align_t;

typedef struct {
	// query sequences
	uint32_t n_seq;
	const bseq_t *seq;	// const!
	void *base;
	uint64_t size;
	// results
	lmm_t *lmm;		// alignment result container
	mm128_v reg;
} mm_align_step_t;

typedef struct mm_tbuf_s { // per-thread buffer
	const mm_align_t *b;
	gaba_section_t qf, qr, t;
	mm128_v resc; // also used as query minimizer vector
	mm128_v coef; // Hough transform coefficient
	v2u32_v intv; // intervals on sorted coef
	poshash_t *pos;
	gaba_dp_t *dp;	// alignment work
} mm_tbuf_t;

#define _s(x)		( (x)<0?-1:1)
#define _m(x)		( (((int32_t)(x))>>31)^(x) )
#define _len(x)		( (x)->re+(x)->qe-(x)->rs-(x)->qs )
#define _clip(x)	MAX2(0, MIN2(((uint32_t)(x)), 60))
static void mm_expand(uint32_t n, const v2u32_t *r, uint32_t qid, int32_t qs, uint32_t thresh, mm128_v *coef)
{
	const int32_t ofs = 0x40000000;
	for (uint64_t i = 0; i < n; ++i) {	// iterate over all the collected minimizers
		if (ofs + r[i].x[1] - qid > thresh) continue;	// 0x3fffffff to skip self and dual, 0xffffffff to keep all
		int32_t rs = (int32_t)r[i].x[0], _qs = (rs>>31) ^ qs, _rs = (rs>>31) ^ rs;
		mm128_t *p;
		kv_pushp(mm128_t, *coef, &p);
		p->u32[0] = ofs + _rs - (_qs>>1); p->u32[1] = r[i].x[1]; p->u32[2] = _rs; p->u32[3] = _qs;
	}
	return;
}

static void mm_collect(const mm_idx_t *mi, uint32_t l_seq, const uint8_t *seq, uint32_t qid, uint32_t max_occ, uint32_t resc_occ, uint32_t thresh, mm128_v *mini, mm128_v *coef)
{
	uint64_t n_resc = 0;
	mini->n = coef->n = 0;
	mm_sketch(seq, l_seq, mi->w, mi->k, 0, mini);

	for (uint64_t i = 0; i < mini->n; ++i) {
		uint64_t n;
		const v2u32_t *r;
		int32_t qs = (int32_t)mini->a[i].u32[2];
		r = mm_idx_get(mi, mini->a[i].u64[0], &n);	// get minimizer at current pos
		if (n > max_occ) continue;	// skip if exceeds repetitive threshold
		if (n > resc_occ) {
			mm128_t *q = &mini->a[n_resc++];
			q->u32[0] = qs; q->u32[1] = n; q->u64[1] = (uintptr_t)r;
			continue;
		}
		mm_expand(n, r, qid, qs, thresh, coef);
	}
	mini->n = n_resc;	// write back rescued array
	return;
}

static void mm_chain(uint64_t l_coef, mm128_t *coef, uint32_t llim, uint32_t hlim, uint32_t min, v2u32_v *intv)
{
	const int32_t ofs = 0x40000000;
	const uint32_t chained = 0x80000000, mask = 0x7fffffff;
	for (uint64_t i = intv->n = 0, j, k, n; i < l_coef; i = MIN2(j, k)) {
		uint32_t rid = coef[i].u32[1];
		int32_t rs = coef[i].u32[2], qs = coef[i].u32[3], re = rs, qe = qs;
		int32_t l, h, lub = coef[i].u32[0] + llim, hub = rs - (qs<<1) + ofs, hlb = hub - hlim;
		uint32_t len = 0;
		for (j = i+1, k = UINT64_MAX, n = i; j < l_coef && (l = (int32_t)coef[j].u32[0]) < lub; j++) {
			if ((int32_t)coef[j].u32[2] < 0) continue;
			re = coef[j].u32[2], qe = coef[j].u32[3];	/* tagged seeds are skipped, no need to mask */
			if (rid != coef[j].u32[1] || (qs^qe)&chained || (h = ofs + re - (qe<<1)) < hlb || h > hub) { k = MIN2(j, k); continue; }	// out of range, skip
			lub = l + llim; hub = h; hlb = h - hlim;
			coef[j].u32[2] |= chained; n = j;
		}
		re = coef[n].u32[2] & mask; qe = coef[n].u32[3]; qs = _m(qs); qe = _m(qe);
		len = re-rs+_s(qe-qs)*(qe-qs);
		if (len < min) continue;
		v2u32_t *p;
		kv_pushp(v2u32_t, *intv, &p);
		p->x[0] = (uint32_t)ofs - len; p->x[1] = i;
	}
	return;
}

static uint64_t mm_short_chain(uint64_t l_coef, mm128_t *coef, uint32_t llim, uint32_t hlim, uint32_t rid, uint32_t eidx)
{
	const int32_t ofs = 0x40000000;
	const uint32_t chained = 0x80000000, mask = 0x7fffffff;
	int32_t rs = coef->u32[2] & mask, qs = coef->u32[3], re = rs, qe = qs;
	int32_t l, h, lub = coef->u32[0] + llim, hub = rs - (qs<<1) + ofs, hlb = hub - hlim;
	uint32_t len = 0;
	if (rid != coef->u32[1]) return 1;
	rid = coef->u32[1];
	uint64_t j, k;
	for (j = 1, k = UINT64_MAX; j < l_coef && (l = (int32_t)coef[j].u32[0]) < lub; j++) {
		re = coef[j].u32[2] & mask, qe = coef[j].u32[3];
		if (rid != coef[j].u32[1] || (qs^qe)&chained || (h = ofs + re - (qe<<1)) < hlb || h > hub) { k = MIN2(j, k); continue; }
		lub = l + llim; hub = h; hlb = h - hlim;	// n = j;
		if (++len > eidx) break;
	}
	return MIN2(j, k);
}

static uint64_t mm_rescue(const mm_idx_seq_t *ref, uint64_t l_coef, mm128_t *coef, mm128_t *r, uint32_t llim, uint32_t hlim, uint32_t elim, uint32_t blim, uint32_t eidx)
{
	const int32_t ofs = 0x40000000, mask = 0x7fffffff;
	int32_t re = r->u32[2], qe = r->u32[3], l, h;
	int32_t lt = ofs + re - (qe>>1), ht = ofs + re - (qe<<1), lb = lt - blim, le = lt + elim, hb = ht + blim, he = ht - elim, lub = lt + llim;
	for (uint64_t i = 0; i < l_coef && (int32_t)coef[i].u32[0] < lb; ++i) {}
	for (uint64_t i = 0; i < l_coef && (l = (int32_t)coef[i].u32[0]) < lub; ++i) {
		uint32_t rs = coef[i].u32[2] & mask, qs = coef[i].u32[3];
		if ((h = ofs + rs - (qs<<1)) > hb || (l < le && h > he && ((l > lt) ^ (h < ht)))) continue;
		return mm_short_chain(l_coef - i, coef + i, llim, hlim, ref->rid, eidx);
	}
	return 1;
}

static const gaba_alignment_t *mm_extend(
	const mm_idx_seq_t *ref, uint32_t l_coef, mm128_t *coef, uint32_t k, uint32_t min, uint32_t sidx, uint32_t eidx,
	gaba_dp_t *dp, gaba_section_t *qf, gaba_section_t *qr, gaba_section_t *t, poshash_t *pos, lmm_t *lmm)
{
	khint_t itr;
	int absent;
	const uint32_t mask = 0x7fffffff;
	const uint8_t *lim = (const uint8_t*)0x800000000000;
	gaba_section_t *qu, *qd, *r, *q;
	gaba_section_t rf = { .id = ref->rid<<1, .len = ref->l_seq, .base = (const uint8_t*)ref->seq };
	gaba_section_t rr = { .id = (ref->rid<<1)+1, .len = ref->l_seq, .base = gaba_rev((const uint8_t*)ref->seq+ref->l_seq-1, lim) };
	gaba_pos_pair_t p = {0};
	gaba_dp_flush(dp, lim, lim);
	uint64_t key = (uint64_t)ref->rid<<32;
	gaba_alignment_t *a = NULL;
	for (uint64_t i = sidx; i < eidx && i < l_coef; ++i) {
		if (i != 0 && (int32_t)coef[i].u32[2] >= 0) continue;	// skip head
		const gaba_stack_t *stack = gaba_dp_save_stack(dp);
		int32_t rs = coef[i].u32[2] & mask, qs = coef[i].u32[3];
		uint64_t rev = qs<0; qu = rev? qf : qr; qd = rev? qr : qf;
		qs = _m(qs); qs = rev? (int32_t)(qf->len-qs+k-1) : qs;
		// upward extension
		gaba_fill_t *f = gaba_dp_fill_root(dp, r = &rr, ref->l_seq-rs-1, q = qu, qu->len-qs-1), *m = f;
		uint32_t flag = GABA_STATUS_TERM;
		do {
			if (f->status & GABA_STATUS_UPDATE_A) flag |= GABA_STATUS_UPDATE_A, r = t;
			if (f->status & GABA_STATUS_UPDATE_B) flag |= GABA_STATUS_UPDATE_B, q = t;
			f = gaba_dp_fill(dp, f, r, q);
			m = (f->max > m->max)? f : m;
		} while(!(flag & f->status));
		// find max
		p = gaba_dp_search_max(dp, m);
		// check duplicate
		key |= p.apos - (p.bpos>>1);
		if ((itr = kh_get(pos, pos, key)) != kh_end(pos)) { return 0; }	// already evaluated
		// downward extension from max
		gaba_dp_flush_stack(dp, stack);
		m = f = gaba_dp_fill_root(dp, r = &rf, ref->l_seq-p.apos-1, q = qd, qd->len-p.bpos-1);
		flag = GABA_STATUS_TERM;
		do {
			if (f->status & GABA_STATUS_UPDATE_A) flag |= GABA_STATUS_UPDATE_A, r = t;
			if (f->status & GABA_STATUS_UPDATE_B) flag |= GABA_STATUS_UPDATE_B, q = t;
			f = gaba_dp_fill(dp, f, r, q);
			m = (f->max > m->max)? f : m;
		} while(!(flag & f->status));
		if (m->max < min) { key &= 0xffffffff00000000; continue; }
		// convert alignment to cigar
		a = gaba_dp_trace(dp, NULL, m, GABA_TRACE_PARAMS( .lmm = lmm ));
		break;
	}
	// record head
	itr = kh_put(pos, pos, key, &absent);
	kh_val(pos, itr) = (uintptr_t)a;
	return a;
}

#if 0
static void mm_record(const mm_idx_seq_t *ref, uint32_t l_seq, const gaba_alignment_t *a, poshash_t *pos, khint_t itr, mm128_v *reg)
{
	uint64_t size;
	mm128_t *s;
	reg_t *r;
	const int32_t ofs = 0x40000000;
	if ((int64_t)kh_val(pos, itr) < 0) {
		kh_val(pos, itr) = (reg->n<<32) | (0xffffffff & kh_val(pos, itr));
		kv_pushp(mm128_t, *reg, &s);
	} else {
		s = &reg->a[kh_val(pos, itr)>>32];
		// free((void*)s->u64[1]);
	}
	/*
	r = (reg_t*)malloc((size = a->path->len < 512 ? 1024 : 2*a->path->len) + a->slen*sizeof(reg_t));
	char *cigar = 
	for (const gaba_section_t *s = a; s < a->sec+a->slen; s++, r++) {
		r->rs = ref->l_seq-s->apos-s->alen; r->re = ref->l_seq-s->apos;
		r->qs = l_seq-s->bpos-s->blen; r->qe = l_seq-s->bpos;
		r->rid = ref->rid; r->score = a->score;
		r->l_cigar = gaba_dp_dump_cigar_reverse((char*)r->cigar, size, a->path->array, s->ppos, gaba_plen(s));
		// r->l_cigar = gaba_dp_dump_cigar_reverse((char*)r->cigar, size - sizeof(reg_t), a->path->array, 0, a->path->len);
	}
	*/
	s->u32[0] = (a->sec->bid&0x01? 0 : 0x10)<<16;
	s->u32[1] = ofs - a->score;
	s->u64[1] = (uintptr_t)a;
	// s->u32[1] = ofs - r->score;
	// s->u64[1] = (uintptr_t)r;
	return;
}
#endif 

#define _reg(x)		( (reg_t*)(x).u64[1] )
static uint64_t mm_post_map(const mm_mapopt_t *opt, uint32_t n_reg, mm128_t *reg)
{
	uint64_t i;
	const reg_t *preg = _reg(reg[0]);
	uint32_t ssc = (n_reg>1)? _reg(reg[1])->score : 0, bsc = (n_reg>1)? _reg(reg[n_reg-1])->score : 0, tsc = 0;
	double elen = (double)_len(preg) / 2.0, pid = 1.0 - (double)(elen * opt->m - preg->score) / (double)(opt->m + opt->x) / elen;
	double ec = 2.0 / (pid * (double)(opt->m + opt->x) - (double)opt->x), ulen = ec * (preg->score - ssc), pe = 1.0 / (ulen * ulen + (double)n_reg);

	for (i = 1; i < n_reg; ++i) tsc += _reg(reg[i])->score - bsc + 1;
	reg[0].u32[0] |= _clip(-10.0 * log10(pe));
	for (i = 1; i < n_reg; ++i)
		reg[i].u32[0] |= (0x100<<16) | _clip(-10.0 * log10(1.0 - pe * (double)(_reg(reg[i])->score - bsc + 1) / (double)tsc));
	return n_reg;
}

static uint64_t mm_post_ava(const mm_mapopt_t *opt, uint32_t n_reg, mm128_t *reg)
{
	for (uint64_t i = 0; i < n_reg; ++i) {
		int32_t score = _reg(reg[i])->score;
		double elen = (double)_len((reg_t*)reg[i].u64[1]) / 2.0, pid = 1.0 - (double)(elen * opt->m - score) / (double)(opt->m + opt->x) / elen;
		double ec = 2.0 / (pid * (double)(opt->m + opt->x) - (double)opt->x), ulen = ec * score, pe = 1.0 / (ulen + 1);
		reg[i].u32[0] |= _clip(-10.0 * log10(pe)) | ((i == 0? 0 : 0x800)<<16);
	}
	return n_reg;
}
#undef _reg

static const mm128_t *mm_align_seq(mm_tbuf_t *b, const mm_mapopt_t *opt, const mm_idx_t *mi, uint32_t l_seq, const uint8_t *seq, uint32_t qid, uint32_t n_occ, uint64_t *occ, lmm_t *lmm, uint64_t *n_reg)
{
	int32_t min = opt->min;
	const int32_t ofs = 0x40000000;
	const uint32_t mask = 0x7fffffff, thresh = opt->flag&MM_AVA? (uint32_t)ofs - 1 : 0xffffffff;
	mm128_v reg = {0};
	// prepare section info for alignment
	uint8_t tail[96];
	const uint8_t *lim = (const uint8_t*)0x800000000000;
	gaba_section_t qf, qr, mg;
	if (l_seq == 0) return 0;
	// qid = opt->flag&MM_AVA? qid : 0xffffffff;
	memset(tail, 0, 96);
	qf.id = 0; qf.len = l_seq; qf.base = (const uint8_t*)seq;
	qr.id = 1; qr.len = l_seq; qr.base = gaba_rev((const uint8_t*)seq+l_seq-1, lim);
	mg.id = 4; mg.len = 32; mg.base = tail+32;
	b->coef.n = b->resc.n = b->intv.n = 0; kh_clear(pos, b->pos);
	for (uint64_t i = 0, j = 0; reg.n == 0 && i < n_occ; ++i) {
		if (i == 0) {
			mm_collect(mi, l_seq, seq, qid, occ[n_occ-1], occ[0], thresh, &b->resc, &b->coef);
		} else {
			if (i == 1) radix_sort_128x(b->resc.a, b->resc.a + b->resc.n);
			for (uint64_t k = 0; k < b->coef.n; ++k) b->coef.a[k].u32[2] &= mask;
			for (; j < b->resc.n && b->resc.a[j].u32[1] <= occ[i]; ++j)
				mm_expand(b->resc.a[j].u32[1], (const v2u32_t*)b->resc.a[j].u64[1], qid, b->resc.a[j].u32[0], thresh, &b->coef);
		}
		radix_sort_128x(b->coef.a, b->coef.a + b->coef.n);
		mm_chain(b->coef.n, b->coef.a, opt->llim, opt->hlim, min>>2, &b->intv);
		radix_sort_64x(b->intv.a, b->intv.a + b->intv.n);
		if (b->intv.a == 0) continue;
		for (uint64_t k = 0; k < b->intv.n && ofs - (int32_t)b->intv.a[k].x[0] > min>>2; ++k) {
			const gaba_alignment_t *a = 0;
			mm128_t *t = &b->coef.a[b->coef.n], r;
			uint64_t p = b->coef.n-b->intv.a[k].x[1], q = p, l = 0, m = ofs - b->intv.a[k].x[0];
			mm_idx_seq_t *ref = &mi->s.a[(t-p)->u32[1]-mi->base_rid];
			do {
				a = mm_extend(ref, q, t-q, mi->k, min, opt->sidx, opt->eidx, b->dp, &qf, &qr, &mg, b->pos, lmm);
				if (p == q && a == 0) break;	// skip chain if first extension did not result in meaningful alignment
				r.u32[2] = (t-p)->u32[2] & mask, r.u32[3] = (t-p)->u32[3];
				if (a != 0) {
					mm128_t *s;
					kv_pushp(mm128_t, reg, &s);
					s->u32[0] = (a->sec->bid&0x01? 0 : 0x10)<<16;
					s->u32[1] = ofs - a->score;
					s->u64[1] = (uintptr_t)a;
					r.u32[2] = ref->l_seq-a->sec->apos; r.u32[3] = (a->sec->bid&0x01)? l_seq-a->sec->bpos : -a->sec->bpos+mi->k+1;
					min = MAX2(min, opt->min_ratio*a->score);
					l += a->sec->alen+a->sec->blen;
				}
				q = p;
			} while (l < m && (p -= mm_rescue(ref, p, t-p, &r, opt->llim, opt->hlim, opt->elim, opt->blim, opt->eidx)) < q - 1);
		}
	}
	*n_reg = reg.n;
	if (reg.n == 0) return 0;
	radix_sort_128x(reg.a, reg.a + reg.n);
	((opt->flag & MM_AVA)? mm_post_ava : mm_post_map)(opt, reg.n, reg.a);
	return reg.a;
}

#undef _s
#undef _m
#undef _len
#undef _clip

#define BUF_SIZE	( 512 * 1024 )

#define _put(_buf, _c) { \
	*(_buf)->p++ = (_c); \
	if ((uintptr_t)(_buf)->p >= (uintptr_t)(_buf)->tail) { \
		fwrite((_buf)->base, sizeof(uint8_t), (_buf)->p - (_buf)->base, stdout); \
		(_buf)->p = (_buf)->base; \
	} \
}
#define _puti(type, _buf, _n) ({ \
	uint8_t _b[16] = {0}; \
	type _m = (type)(_n); int64_t _i = 0; \
	while (_m) _b[_i++] = _m % 10, _m /= 10; \
	for (int64_t _j = _i - (_i!=0); _j >= 0; _j--) { _put(_buf, _b[_j]+'0'); } \
	_i - (_i!=0) + 1; \
})
#define _putn(_buf, _n) _puti(uint32_t, _buf, _n)
#define _puts(_buf, _s) { \
	for (const uint8_t *_q = (const uint8_t*)(_s); *_q; _q++) { _put(_buf, *_q); } \
}

static void mm_print_header(mm_align_t *b)
{
	const mm_idx_t *mi = b->mi;
	_puts(b, "@HD\tVN:1.0\tSO:unsorted\n");
	for (uint64_t i = 0; i < mi->s.n; ++i) {
		_puts(b, "@SQ\tSN:"); _puts(b, mi->s.a[i].name);
		_puts(b, "\tLN:"); _putn(b, mi->s.a[i].l_seq);
		_put(b, '\n');
	}
	if (b->opt->flag & 0x01ULL<<MM_RG) { _puts(b, b->opt->rg_line); _put(b, '\n'); }
	_puts(b, "@PG\tID:minialign\tPN:minialign\n");
	return;
}

static void mm_print_tags(mm_align_t *b, const bseq_t *t, const gaba_alignment_t *a, uint16_t flag)
{
	uint64_t f = b->opt->flag;
	if (f & 0x01ULL<<MM_RG) { _puts(b, "\tRG:Z:"); _puts(b, b->opt->rg_id); }
	if (a == 0) return;

	if (f & 0x01ULL<<MM_AS) { _puts(b, "\tAS:i:"); _putn(b, a->score); }
	return;
}

static int32_t mm_search_tag(mm_align_t *b, char tag1, char tag2)
{
	uint16_t tag = (uint16_t)tag1 | ((uint16_t)tag2<<8);
	for (uint64_t i = 0; i < b->opt->tags.n; ++i) {
		if (b->opt->tags.a[i] == tag) return i;
	}
	return -1;
}

static uint64_t mm_print_num(mm_align_t *b, uint8_t type, const uint8_t *p)
{
	if (type == 'a') { _put(b, *p); return 1; }
	else if (type == 'c') { _puti(int8_t, b, (int8_t)*p); return 1; }
	else if (type == 'C') { _puti(uint8_t, b, (uint8_t)*p); return 1; }
	else if (type == 's') { _puti(int16_t, b, *((int16_t*)p)); return 2; }
	else if (type == 'S') { _puti(uint16_t, b, *((uint16_t*)p)); return 2; }
	else if (type == 'i') { _puti(int32_t, b, *((int32_t*)p)); return 4; }
	else if (type == 'I') { _puti(uint32_t, b, *((uint32_t*)p)); return 4; }
	else if (type == 'f') {
		char buf[32]; uint64_t l = sprintf(buf, "%f", *((float*)p));
		for (uint64_t i = 0; i < l; ++i) { _put(b, buf[i]); }
		return 4;
	}
	return 1;
}

static void mm_restore_tags(mm_align_t *b, const bseq_t *t)
{
	static const uint8_t tag_size[32] = {
		0, 1, 0xfe, 1,  0, 0, 4, 0,  0xfe, 4, 0, 0,  0, 0, 0, 0,
		0, 0, 0, 2,     0, 0, 0, 0,  0, 0, 0xff, 0,  0, 0, 0, 0,
	};
	uint8_t size;
	const uint8_t *p = t->tag, *tail = p + t->l_tag;
	while (p < tail) {
		if (mm_search_tag(b, p[0], p[1]) < 0) {	// skip
			if ((size = tag_size[p[2]&0x1f]) < 0x10) p += 3 + size;
			else if(size == 0xfe) p += 4 + tag_size[p[3]&0x1f] * *((uint32_t*)&p[4]);
			else while (*p++) {}
			continue;
		}
		_put(b, '\t'); _put(b, p[0]); _put(b, p[1]); _put(b, ':'); _put(b, p[2]); _put(b, ':');
		if (p[2] == 'Z') {
			p += 3; while(*p) { _put(b, *p); p++; } p++;
		} else if (tag_size[p[2]&0x1f] == 0xfe) {
			uint8_t type = p[3];
			uint32_t len = *((uint32_t*)&p[4]); p += 8;
			for (uint64_t i = 0; i < len; ++i) { p += mm_print_num(b, type, p); _put(b, ','); }
		} else {
			p += mm_print_num(b, p[2], p + 3) + 3;
		}
	}
	return;
}

static void mm_print_unmapped(mm_align_t *b, const bseq_t *t)
{
	_puts(b, t->name); _puts(b, "\t4\t*\t0\t0\t*\t*\t0\t0\t");
	for (uint64_t k = 0; k < t->l_seq; k++) _put(b, "NACMGRSVTWYHKDBN"[(uint8_t)t->seq[k]]);
	_put(b, '\t');
	if (b->opt->flag&MM_KEEP_QUAL && t->qual[0] != '\0') { _puts(b, t->qual); }
	else { _put(b, '*'); }
	mm_restore_tags(b, t); _put(b, '\n');
	return;
}

static int mm_cigar_printer(void *_b, int64_t len, char c)
{
	mm_align_t *b = (mm_align_t*)_b;
	len = _putn(b, (uint32_t)len); _put(b, c);
	return len+1;
}

static void mm_print_mapped(mm_align_t *b, const bseq_t *t, const gaba_alignment_t *a, uint16_t mapq, uint16_t flag)
{
	const mm_idx_t *mi = b->mi;
	const mm_idx_seq_t *r = &mi->s.a[(a->sec->aid>>1) - mi->base_rid];
	const gaba_path_section_t *s = &a->sec[0];
	uint32_t rs = r->l_seq - s->apos - s->alen;
	uint32_t qs = (flag&0x900)? t->l_seq-s->bpos-s->blen : 0, qe = (flag&0x900)? t->l_seq-s->bpos : t->l_seq;

	_puts(b, t->name); _put(b, '\t'); _putn(b, flag); _put(b, '\t'); _puts(b, r->name); _put(b, '\t');
	_putn(b, rs+1); _put(b, '\t'); _putn(b, mapq); _put(b, '\t');
	if (qs) { _putn(b, qs); _put(b, (flag&0x900)? 'H' : 'S'); }
	gaba_dp_print_cigar_reverse(mm_cigar_printer, b, a->path->array, 0, a->path->len);//s->ppos, gaba_plen(s));
	if (t->l_seq-qe) { _putn(b, t->l_seq-qe); _put(b, (flag&0x900)? 'H' : 'S'); }
	_puts(b, "\t*\t0\t0\t");
	if (flag&0x10) { for (int64_t k = t->l_seq-qs; k > t->l_seq-qe; k--) _put(b, "NTGKCYSBAWRDMHVN"[(uint8_t)t->seq[k-1]]); }
	else { for (int64_t k = qs; k < qe; k++) _put(b, "NACMGRSVTWYHKDBN"[(uint8_t)t->seq[k]]); }
	_put(b, '\t');
	if (b->opt->flag&MM_KEEP_QUAL && t->qual[0] != '\0') {
		if (flag&0x10) { for (int64_t k = t->l_seq-qs; k > t->l_seq-qe; k--) _put(b, t->qual[k-1]); }
		else { for (int64_t k = qs; k < qe; k++) _put(b, t->qual[k]); }
	} else {
		_put(b, '*');
	}
	mm_print_tags(b, t, a, flag); mm_restore_tags(b, t); _put(b, '\n');
	return;
}

#undef _put
#undef _putn
#undef _puts

static void *mm_align_source(void *arg)
{
	mm_align_t *b = (mm_align_t*)arg;
	mm_align_step_t *s = (mm_align_step_t*)calloc(1, sizeof(mm_align_step_t));
	s->lmm = lmm_init(NULL, 1024 * 1024);
	if (s->lmm == 0) return 0;
	s->seq = bseq_read(b->fp, b->opt->batch_size, &s->n_seq, &s->base, &s->size);
	if (s->seq == 0) lmm_clean(s->lmm), free(s), s = 0;
	return s;
}

static void *mm_align_worker(void *arg, void *item)
{
	mm_tbuf_t *t = (mm_tbuf_t*)arg;
	mm_align_step_t *s = (mm_align_step_t*)item;
	for (uint64_t i = 0; i < s->n_seq; ++i) {
		mm128_t *r;
		kv_pushp(mm128_t, s->reg, &r);
		r->u64[0] = (uintptr_t)mm_align_seq(t, t->b->opt, t->b->mi, s->seq[i].l_seq, s->seq[i].seq, s->seq[i].rid, t->b->n_occ, t->b->occ, s->lmm, &r->u64[1]);
	}
	return s;
}

static void mm_align_drain(void *arg, void *item)
{
	mm_align_t *b = (mm_align_t*)arg;
	mm_align_step_t *s = (mm_align_step_t*)item;

	for (uint64_t i = 0; i < s->n_seq; ++i) {
		mm128_t *r = (mm128_t*)s->reg.a[i].u64[0];
		uint32_t n_reg = s->reg.a[i].u64[1];
		if (r == 0) {
			mm_print_unmapped(b, &s->seq[i]);
			continue;
		}
		gaba_alignment_t *aln;
		int32_t min = ((gaba_alignment_t*)r[0].u64[1])->score * b->opt->min_ratio;
		for (uint64_t j = 0; j < n_reg && (aln = (gaba_alignment_t*)r[j].u64[1])->score >= min; ++j) {
			mm_print_mapped(b, &s->seq[i], aln, r[j].u32[0] & 0xffff, r[j].u32[0]>>16);
			lmm_free(s->lmm, aln);
		}
		free(r);
	}
	free(s->reg.a); free(s->base); free((void*)s->seq); lmm_clean(s->lmm); free(s);
	return;
}

static void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == 0) return;
	free(b->resc.a); free(b->coef.a); free(b->intv.a);
	kh_destroy(pos, b->pos); gaba_dp_clean(b->dp);
	free(b);
}

static mm_tbuf_t *mm_tbuf_init(mm_align_t *b)
{
	mm_tbuf_t *t = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	if (t == 0) return 0;
	const uint8_t *lim = (const uint8_t *)0x800000000000;
	t->b = b;
	if ((t->pos = kh_init(pos)) == 0) goto _fail;
	if ((t->dp = gaba_dp_init(b->gaba, lim, lim)) == 0) goto _fail;
	return t;
_fail:
	mm_tbuf_destroy(t);
	return 0;
}

static void mm_align_destroy(mm_align_t *b)
{
	assert(b != 0);
	fwrite(b->base, sizeof(uint8_t), b->p - b->base, stdout);
	for (uint64_t i = 0; i < b->opt->n_threads; ++i) mm_tbuf_destroy((mm_tbuf_t*)b->t[i]);
	free(b->base); free(b->occ); free(b->t); gaba_clean(b->gaba); ptask_clean(b->pt);
	free(b);
	return;
}

static mm_align_t *mm_align_init(const mm_mapopt_t *opt, const mm_idx_t *mi)
{
	assert(mi != 0);
	assert(opt != 0);

	// init output buf
	mm_align_t *b = calloc(1, sizeof(mm_align_t));
	if (b == 0) return 0;
	b->tail = (b->base = b->p = malloc(sizeof(uint8_t) * opt->outbuf_size)) + opt->outbuf_size;
	if (b->base == 0) goto _fail;

	// set consts
	b->mi = mi, b->opt = opt;

	// calc occ
	if ((b->occ = calloc(b->opt->n_frq, sizeof(uint64_t))) == 0) goto _fail;
	b->n_occ = b->opt->n_frq;
	for (uint64_t i = 0; i < b->n_occ; ++i) b->occ[i] = mm_idx_cal_max_occ(mi, b->opt->frq[i]);

	// initialize alignment context
	int m = opt->m, x = -opt->x, gi = opt->gi, ge = opt->ge;
	struct gaba_score_s sc = {{{m,x,x,x},{x,m,x,x},{x,x,m,x},{x,x,x,m}},gi,ge,gi,ge};
	struct gaba_params_s p = {0,0,7,opt->xdrop,&sc};
	if ((b->gaba = gaba_init(&p)) == 0) goto _fail;

	// initialize threads
	if ((b->t = (void**)calloc(b->opt->n_threads, sizeof(mm_tbuf_t*))) == 0) goto _fail;
	for (uint64_t i = 0; i < b->opt->n_threads; ++i) { if ((b->t[i] = (void*)mm_tbuf_init(b)) == 0) goto _fail; }
	if ((b->pt = ptask_init(mm_align_worker, (void**)b->t, b->opt->n_threads, 64)) == 0) goto _fail;

	// print sam header
	mm_print_header(b);
	return b;
_fail:
	mm_align_destroy(b);
	return 0;
}

static int mm_align_file(mm_align_t *b, bseq_file_t *fp)
{
	assert(b != 0);
	if (fp == 0) return -1;
	b->fp = fp;
	return ptask_stream(b->pt, mm_align_source, b, mm_align_drain, b, 8*b->opt->n_threads);
}

/* end of map.c */

/* main.c */

static void liftrlimit()
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}

static void posixly_correct()
{
#ifdef __linux__
	setenv("POSIXLY_CORRECT", "1", 0);
#endif
}

static void mm_print_help(const mm_mapopt_t *opt)
{
	fprintf(stderr, "\n"
					"  minialign - fast aligner for PacBio and Nanopore long reads\n"
					"\n"
					"minialign is a fast long-read (nucleotide sequence) alignment tool built on\n"
					"the top of minimap long-read overlapper adopting libgaba SIMD-parallelized\n"
					"Smith-Waterman extension algorithm.\n"
					"\n");
	fprintf(stderr, "Usage:\n"
					"  first trial:\n"
					"    $ minialign -t4 -xont <ref.fa> <ont2d.{fa,fq,bam}> > mapping.sam\n"
					"\n"
					"  mapping on a prebuilt index (saves ~1min for human genome per run):\n"
					"    $ minialign [indexing options] -d <index.mai> <ref.fa>\n"
					"    $ minialign -l <index.mai> <reads.{fa,fq,bam}> > mapping.sam\n"
					"\n"
					"  all-versus-all alignment in a read set:\n"
					"    $ minialign -X -xava <reads.fa> [<reads.fa> ...] | samsplit <prefix>\n"
					"\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  Global:\n");
	fprintf(stderr, "    -x STR       load preset params {pacbio, ont, ava} [ont]\n");
	fprintf(stderr, "    -t INT       number of threads [%d]\n", opt->n_threads);
	fprintf(stderr, "    -X           switch to all-versus-all alignment mode\n");
	fprintf(stderr, "    -v           show version number [%s]\n", MM_VERSION);
	fprintf(stderr, "  Indexing:\n");
	fprintf(stderr, "    -c STR,...   treat specified sequences as circular []\n");
	fprintf(stderr, "    -k INT       k-mer size [%d]\n", opt->k);
	fprintf(stderr, "    -w INT       minimizer window size [{-k}*2/3]\n");
	fprintf(stderr, "    -d FILE      dump index to FILE []\n");
	fprintf(stderr, "    -l FILE      load index from FILE [] (overriding -k and -w)\n");
	fprintf(stderr, "  Mapping:\n");
	// fprintf(stderr, "    -f FLOAT,... occurrence thresholds [0.05,0.01,0.001]\n");
	fprintf(stderr, "    -a INT       match award [%d]\n", opt->m);
	fprintf(stderr, "    -b INT       mismatch penalty [%d]\n", opt->x);
	fprintf(stderr, "    -p INT       gap open penalty [%d]\n", opt->gi);
	fprintf(stderr, "    -q INT       gap extension penalty [%d]\n", opt->ge);
	fprintf(stderr, "    -s INT       minimum alignment score [%d]\n", opt->min);
	fprintf(stderr, "    -m INT       minimum alignment score ratio [%1.2f]\n", opt->min_ratio);
	fprintf(stderr, "  Output:\n");
	fprintf(stderr, "    -Q           include quality string\n");
	fprintf(stderr, "    -R STR       read group header line, like \"@RG\\tID:1\" []\n");
	fprintf(stderr, "    -T STR,...   list of optional tags: {RG, AS} []\n");
	fprintf(stderr, "                   (RG is also inferred from -R)\n");
	fprintf(stderr, "    -U STR,...   tags to be transferred from the input bam file []\n");
	fprintf(stderr, "\n");
	return;
}

static uint64_t mm_mapopt_parse_tags(const char *p, uint16_v *buf)
{
	uint32_t flag = 0;
	while (*p != '\0') {
		const char *q = p;
		while (*q != ',' && *q != '\0') q++;
		if (q - p == 2) {
			if (buf) kv_push(uint16_t, *buf, (uint16_t)p[0] | ((uint16_t)p[1]<<8));
			if (strncmp(p, "RG", 2) == 0) flag |= 0x01ULL<<MM_RG;
			else if (strncmp(p, "AS", 2) == 0) flag |= 0x01ULL<<MM_AS;
			else if (strncmp(p, "XS", 2) == 0) flag |= 0x01ULL<<MM_XS;
			else if (strncmp(p, "SA", 2) == 0) flag |= 0x01ULL<<MM_SA;
			else if (strncmp(p, "MD", 2) == 0) flag |= 0x01ULL<<MM_MD;
		}
		if (*q == '\0') break;
		p = q + 1;
	}
	return flag;
}

static void mm_mapopt_parse_rg(mm_mapopt_t *o, const char *arg)
{
	free(o->rg_line); free(o->rg_id);
	o->rg_line = NULL; o->rg_id = NULL;
	o->flag &= ~(0x01ULL<<MM_RG);

	const char *id, *p = arg;
	char b[256], *q = b;
	while (q < b+256 && *p != '\0') {
		if (p[0] == '\\' && p[1] == 't') { if (p[1] == '\0') { break; } *q++ = '\t'; p += 2; }
		else { *q++ = *p++; }
	}
	*q = '\0';
	if (strstr(b, "@RG") != b) return;
	if ((id = strstr(b, "\tID:")) == 0) return;

	#define _dup(x)	({ char *_p = malloc(strlen(x)+1); memcpy(_p, (x), strlen(x)); _p[strlen(x)] = '\0'; _p; })
	o->rg_line = _dup(b); o->rg_id = _dup(id+4);
	char *t = o->rg_id;
	while (*t != '\t' && *t != '\n' && *t != '\0') { t++; }
	*t = '\0';
	o->flag |= 0x01ULL<<MM_RG;
	#undef _dup
	return;
}



static int mm_mapopt_parse(mm_mapopt_t *o, int argc, char *argv[], const char **fnr, const char **fnw, ptr_v *v)
{
	while (optind < argc) {
		int ch;
		if ((ch = getopt(argc, argv, "k:w:f:c:x:B:t:V:d:l:Xs:m:r:a:b:p:q:L:H:I:J:S:E:Y:QR:T:U:vh")) < 0) {
			kv_push(void*, *v, argv[optind]); optind++; continue;
		}

		if (ch == 'k') o->k = atoi(optarg);
		else if (ch == 'w') o->w = atoi(optarg);
		else if (ch == 'f') {
			const char *p = optarg; o->n_frq = 0;
			while (*p) {
				o->frq[o->n_frq++] = atof(p);
				while (*p && *p != ',') p++;
				if (!*p || o->n_frq >= 16) break;
				p++;
			}
		}
		else if (ch == 'x') {
			if (strcmp(optarg, "pacbio") == 0) {
				o->k = 15; o->w = 10; o->m = 1; o->x = 2; o->gi = 2; o->ge = 1; o->xdrop = 50; o->min = 50; o->min_ratio = 0.3;
			} else if (strcmp(optarg, "ont") == 0) {
				o->k = 15; o->w = 10; o->m = 1; o->x = 1; o->gi = 1; o->ge = 1; o->xdrop = 50; o->min = 50; o->min_ratio = 0.3;
			} else if (strcmp(optarg, "ava") == 0) {
				o->k = 14; o->w = 10; o->m = 1; o->x = 1; o->gi = 1; o->ge = 1; o->xdrop = 50; o->min = 30; o->min_ratio = 0.3; o->flag |= MM_AVA;
			} else {
				if (mm_verbose >= 3) fprintf(stderr, "[M::%s] Warning: Unknown preset tag: `%s'.\n", __func__, optarg);
			}
		}
		else if (ch == 'B') o->b = atoi(optarg);
		else if (ch == 't') o->n_threads = atoi(optarg);
		else if (ch == 'V') mm_verbose = atoi(optarg);
		else if (ch == 'd') *fnw = optarg;
		else if (ch == 'l') *fnr = optarg;
		else if (ch == 'X') o->flag |= MM_AVA;
		else if (ch == 's') o->min = atoi(optarg);
		else if (ch == 'm') o->min_ratio = atof(optarg);
		else if (ch == 'r') {
			if (mm_verbose >= 3) fprintf(stderr, "[M::%s] Warning: Minimum length threshold option is deprecated in version 0.4.0 and later, interpreted as score ratio.\n", __func__);
			o->min_ratio = atof(optarg);
		}
		else if (ch == 'a') o->m = atoi(optarg);
		else if (ch == 'b') o->x = atoi(optarg);
		else if (ch == 'p') o->gi = atoi(optarg);
		else if (ch == 'q') o->ge = atoi(optarg);
		else if (ch == 'L') o->llim = atoi(optarg);
		else if (ch == 'H') o->hlim = atoi(optarg);
		else if (ch == 'I') o->blim = atoi(optarg);
		else if (ch == 'J') o->elim = atoi(optarg);
		else if (ch == 'S') o->sidx = atoi(optarg);
		else if (ch == 'E') o->eidx = atoi(optarg);
		else if (ch == 'Y') o->xdrop = atoi(optarg);
		else if (ch == 'Q') o->flag |= MM_KEEP_QUAL;
		else if (ch == 'R') mm_mapopt_parse_rg(o, optarg);
		else if (ch == 'T') o->flag |= mm_mapopt_parse_tags(optarg, NULL);
		else if (ch == 'U') mm_mapopt_parse_tags(optarg, &o->tags);
		else if (ch == 'v') { return 1; }
		else if (ch == 'h') { return 2; }
	}
	return 0;
}

int main(int argc, char *argv[])
{
	int ret = 1;
	uint32_t base_rid = 0, base_qid = 0;
	const char *fnr = 0, *fnw = 0;
	bseq_file_t *fp = 0;
	ptr_v v = {0};
	FILE *fpr = 0, *fpw = 0;

	liftrlimit(); posixly_correct();
	mm_realtime0 = realtime();
	mm_mapopt_t *opt = mm_mapopt_init();
	switch (mm_mapopt_parse(opt, argc, argv, &fnr, &fnw, &v)) {
		case 1: puts(MM_VERSION); ret = 0; goto _final;
		case 2: mm_print_help(opt); ret = 0; goto _final;
	}
	if (!fnr && v.n == 0) { mm_print_help(opt); ret = 1; goto _final; }
	if ((fnr && v.n == 0) || (!fnr && v.n == 1 && !(opt->flag&MM_AVA))) {
		fprintf(stderr, "[M::%s] query-side input redirected to stdin.\n", __func__);
		kv_push(void*, v, "-");
	}
	if (opt->w >= 16) opt->w = (int)(.6666667 * opt->k + .499);
	if (mm_mapopt_check(opt, fprintf, stderr)) {
		ret = 1; goto _final;
	}

	if (fnr) fpr = fopen(fnr, "rb");
	if (fnw) fpw = fopen(fnw, "wb");
	for (uint64_t i = 0; i < (fpr? 0x7fffffff : (((opt->flag&MM_AVA) || fpw)? v.n : 1)); ++i) {
		uint32_t qid = base_qid;
		mm_idx_t *mi = 0;
		mm_align_t *aln = 0;
		if (fpr) mi = mm_idx_load(fpr);
		else mi = mm_idx_gen(opt, (fp = bseq_open((const char *)v.a[i], base_rid, 0, 0))), base_rid = bseq_close(fp);
		if (mi == 0) {
			if (fpr && i > 0) break;
			fprintf(stderr, "[M::%s] ERROR: failed to %s `%s'. Please check %s.\n", __func__, fpr? "load index file" : "open sequence file", fpr? fnr : (const char*)v.a[i], fpr? "file path, format and its version" : "file path and format");
			ret = 1; goto _final;
		}
		if (mm_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built index for %lu target sequence(s)\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->s.n);
		if (fpw) mm_idx_dump(fpw, mi);
		else aln = mm_align_init(opt, mi);
		for (uint64_t j = (!fpr && !(opt->flag&MM_AVA)); j < (fpw? 0 : v.n); ++j) mm_align_file(aln, (fp = bseq_open((const char*)v.a[j], qid, (opt->flag&MM_KEEP_QUAL)!=0, opt->tags.n!=0))), qid = bseq_close(fp);
		if (!fpw) mm_align_destroy(aln);
		mm_idx_destroy(mi);
	}
	fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (uint64_t i = 0; i < (uint64_t)argc; ++i) fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
	ret = 0;
_final:
	free(v.a); mm_mapopt_destroy(opt);
	if (fpr) fclose(fpr);
	if (fpw) fclose(fpw);
	return ret;
}

/* end of main.c */
