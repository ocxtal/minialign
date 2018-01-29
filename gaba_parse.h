#define TEMP		0
/**
 * @file gaba_parse.h
 *
 * @brief libgaba utility function implementations
 *
 * @author Hajime Suzuki
 * @date 2018/1/20
 * @license Apache v2
 *
 * NOTE: Including this header is not recommended since this header internally include
 *   arch/arch.h, which floods a lot of SIMD-related defines. The same function will
 *   be found in the libgaba.a. The functions can be overridden by this header by
 *   compiling the file with a proper SIMD-enabling flag, such as `-mavx2 -mbmi -mpopcnt'
 *   (see Makefile for the details), which result in a slightly better performace
 *   with the functions inlined.
 *
 * CIGAR printers:
 *   gaba_dump_cigar_forward:
 *   gaba_dump_cigar_reverse:
 *   gaba_print_cigar_forward:
 *   gaba_print_cigar_reverse:
 *
 * MAF-style gapped alignment printers:
 *   gaba_dump_gapped_forward
 *   gaba_dump_gapped_reverse
 *
 * match and gap counter (also calculates score of subalignment):
 *   gaba_calc_score
 */

#ifndef _GABA_PARSE_H_INCLUDED
#define _GABA_PARSE_H_INCLUDED
#include <stdint.h>				/* uint32_t, uint64_t, ... */
#include <stddef.h>				/* ptrdiff_t */
#include "gaba.h"
#include "arch/arch.h"


/**
 * @type gaba_printer_t
 * @brief callback to CIGAR element printer
 * called with a pair of cigar operation (c) and its length (len).
 * void *fp is an opaque pointer to the context of the printer.
 */
typedef int (*gaba_printer_t)(void *, uint64_t, char);	/* moved to gaba.h */


/* macros */
#ifndef _GABA_PARSE_EXPORT_LEVEL
#  define _GABA_PARSE_EXPORT_LEVEL		static inline	/* hidden */
#endif

#define _gaba_parse_min2(_x, _y)		( (_x) < (_y) ? (_x) : (_y) )
#define _gaba_parse_ptr(_p)				( (uint64_t const *)((uint64_t)(_p) & ~(sizeof(uint64_t) - 1)) )
#define _gaba_parse_ofs(_p)				( ((uint64_t)(_p) & sizeof(uint32_t)) ? 32 : 0 )


#if BIT == 2
/* 2bit encoding */
static uint8_t const gaba_parse_ascii_fw[16] __attribute__(( aligned(16) )) = {
	'A', 'C', 'G', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};
static uint8_t const gaba_parse_ascii_rv[16] __attribute__(( aligned(16) )) = {
	'T', 'G', 'C', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};
#else
/* 4bit encoding */
static uint8_t const gaba_parse_ascii_fw[16] __attribute__(( aligned(16) )) = {
	'N', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
};
static uint8_t const gaba_parse_ascii_rv[16] __attribute__(( aligned(16) )) = {
	'N', 'T', 'G', 'K', 'C', 'Y', 'S', 'B', 'A', 'W', 'R', 'D', 'M', 'H', 'V', 'N'
};
#endif


/* dp context (opaque) */
#ifndef _GABA_H_INCLUDED
typedef struct gaba_dp_context_s gaba_dp_t;
#endif
struct gaba_parse_dp_context_s {
	uint64_t _pad[12];					/** (96) margin */

	/* score constants */
	// struct gaba_score_vec_s scv;		/** (80) substitution matrix and gaps */

	/* scores */
	float imx, xmx;						/** (8) 1 / (M - X), X / (M - X) (precalculated constants) */
	int8_t tx;							/** (1) xdrop threshold */
	int8_t tf;							/** (1) filter threshold */
	int8_t gi, ge, gfa, gfb;			/** (4) negative integers */
	uint8_t aflen, bflen;				/** (2) short-gap length thresholds */
	/** 192; 64byte aligned */
};

/**
 * @fn gaba_parse_u64
 */
static inline
uint64_t gaba_parse_u64(
	uint64_t const *ptr,
	int64_t pos)
{
	int64_t rem = pos & 63;
	uint64_t a = (ptr[pos>>6]>>rem) | ((ptr[(pos>>6) + 1]<<(63 - rem))<<1);
	return(a);
}

/**
 * @macro _gaba_diag
 * @brief convert contiguous diagonal transitions to zeros
 */
#define _gaba_diag(_x)			( (_x) ^ 0x5555555555555555 )

/**
 * @macro _parser_init_fw, _parser_init_rv
 * @brief initialize variables for the parser template
 */
#define _parser_init_fw(_path, _offset, _len) \
	uint64_t const *p = _gaba_parse_ptr(_path); \
	uint64_t lim = (_offset) + _gaba_parse_ofs(_path) + (_len), ridx = (_len);
#define _parser_init_rv(_path, _offset, _len) \
	uint64_t const *p = _gaba_parse_ptr(_path); \
	uint64_t ofs = (int64_t)(_offset) + _gaba_parse_ofs(_path) - 64, idx = (_len);

/**
 * @macro _parser_loop_fw, _parser_loop_rv
 * @brief parser templates
 */
#define _parser_loop_fw(_del_dump, _ins_dump, _diag_step, _diag_end) { \
	while((int64_t)ridx > 0) { \
		ZCNT_RESULT uint64_t m; uint64_t c; \
		/* deletions */ \
		m = tzcnt(gaba_parse_u64(p, lim - ridx)); \
		c = _gaba_parse_min2(ridx, m); \
		ridx -= c; _del_dump(c);		/* at most 64bp */ \
		/* insertions */ \
		m = tzcnt(~gaba_parse_u64(p, lim - ridx)); \
		c = _gaba_parse_min2(ridx, m - 1); \
		ridx -= c; _ins_dump(c);		/* at most 64bp */ \
		/* diagonals */ \
		uint64_t sridx = ridx; \
		do { \
			m = tzcnt(_gaba_diag(gaba_parse_u64(p, lim - ridx))); \
			c = _gaba_parse_min2(ridx, m) & ~0x01; \
			ridx -= c; _diag_step(c>>1); \
		} while(c == 64); \
		_diag_end((sridx - ridx)>>1); \
	} \
}
#define _parser_loop_rv(_del_dump, _ins_dump, _diag_dump, _diag_end) { \
	while((int64_t)idx > 0) { \
		ZCNT_RESULT uint64_t m; uint64_t c; \
		/* insertions */ \
		m = lzcnt(~gaba_parse_u64(p, ofs + idx)); \
		c = _gaba_parse_min2(idx, m); \
		idx -= c; _ins_dump(c);		/* at most 64bp */ \
		/* deletions */ \
		m = lzcnt(gaba_parse_u64(p, ofs + idx)); \
		c = _gaba_parse_min2(idx, m - 1); \
		idx -= c; _del_dump(c);		/* at most 64bp */ \
		/* diagonals */ \
		uint64_t sidx = idx; \
		do { \
			m = lzcnt(_gaba_diag(gaba_parse_u64(p, ofs + idx))); \
			c = _gaba_parse_min2(idx, m) & ~0x01; \
			idx -= c; _diag_dump(c>>1); \
		} while(c == 64); \
		_diag_end((sidx - idx)>>1); \
	} \
}

/**
 * @fn gaba_parse_dump_num
 */
static inline
uint64_t gaba_parse_dump_num(
	char *buf,
	uint64_t len,
	char ch)
{
	if(len < 64) {
		static uint8_t const conv[64] = {
			0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09,
			0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19,
			0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29,
			0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39,
			0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49,
			0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59,
			0x60, 0x61, 0x62, 0x63
		};
		char *p = buf;
		*p = (conv[len]>>4) + '0'; p += (conv[len] & 0xf0) != 0;
		*p++ = (conv[len] & 0x0f) + '0';
		*p++ = ch;
		return(p - buf);
	} else {
		uint64_t adv;
		uint8_t b[16] = { (uint8_t)ch, '0' }, *p = &b[1];
		while(len != 0) { *p++ = (len % 10) + '0'; len /= 10; }
		for(p -= (p != &b[1]), adv = (uint64_t)((ptrdiff_t)(p - b)) + 1; p >= b; p--) { *buf++ = *p; }
		return(adv);
	}
}

/**
 * @fn gaba_print_cigar_forward
 * @brief parse path string and print cigar to file
 */
_GABA_PARSE_EXPORT_LEVEL
uint64_t gaba_print_cigar_forward(
	gaba_printer_t printer,
	void *fp,
	uint32_t const *path,
	uint64_t offset,
	uint64_t len)
{
	#define _del(_c)	{ if(_c) { clen += printer(fp, _c, 'D'); } }
	#define _ins(_c)	{ if(_c) { clen += printer(fp, _c, 'I'); } }
	#define _match(_c)	{ if(_c) { clen += printer(fp, _c, 'M'); } }
	#define _nop(_c)	{}

	uint64_t clen = 0;
	_parser_init_fw(path, offset, len);
	_parser_loop_fw(_del, _ins, _nop, _match);
	return(clen);

	#undef _del
	#undef _ins
	#undef _match
	#undef _nop
}

/**
 * @fn gaba_dump_cigar_forward
 * @brief parse path string and store cigar to buffer
 */
_GABA_PARSE_EXPORT_LEVEL
uint64_t gaba_dump_cigar_forward(
	char *buf,
	uint64_t buf_size,
	uint32_t const *path,
	uint64_t offset,
	uint64_t len)
{
	#define _del(_c)	{ if(_c) { b += gaba_parse_dump_num(b, _c, 'D'); } }
	#define _ins(_c)	{ if(_c) { b += gaba_parse_dump_num(b, _c, 'I'); } }
	#define _match(_c)	{ if(_c) { b += gaba_parse_dump_num(b, _c, 'M'); } }
	#define _nop(_c)	{}

	char *b = buf;
	_parser_init_fw(path, offset, len);
	_parser_loop_fw(_del, _ins, _nop, _match);
	*b = '\0';
	return(b - buf);

	#undef _del
	#undef _ins
	#undef _match
	#undef _nop
}

/**
 * @fn gaba_print_cigar_reverse
 * @brief parse path string and print cigar to file
 */
_GABA_PARSE_EXPORT_LEVEL
uint64_t gaba_print_cigar_reverse(
	gaba_printer_t printer,
	void *fp,
	uint32_t const *path,
	uint64_t offset,
	uint64_t len)
{
	#define _del(_c)	{ if(_c) { clen += printer(fp, _c, 'D'); } }
	#define _ins(_c)	{ if(_c) { clen += printer(fp, _c, 'I'); } }
	#define _match(_c)	{ if(_c) { clen += printer(fp, _c, 'M'); } }
	#define _nop(_c)	{ (void)(_c); }

	uint64_t clen = 0;
	_parser_init_rv(path, offset, len);
	_parser_loop_rv(_del, _ins, _nop, _match);
	return(clen);

	#undef _del
	#undef _ins
	#undef _match
	#undef _nop
}

/**
 * @fn gaba_dump_cigar_reverse
 * @brief parse path string and store cigar to buffer
 */
_GABA_PARSE_EXPORT_LEVEL
uint64_t gaba_dump_cigar_reverse(
	char *buf,
	uint64_t buf_size,
	uint32_t const *path,
	uint64_t offset,
	uint64_t len)
{
	#define _del(_c)	{ if(_c) { b += gaba_parse_dump_num(b, _c, 'D'); } }
	#define _ins(_c)	{ if(_c) { b += gaba_parse_dump_num(b, _c, 'I'); } }
	#define _match(_c)	{ if(_c) { b += gaba_parse_dump_num(b, _c, 'M'); } }
	#define _nop(_c)	{ (void)(_c); }

	char *b = buf;
	_parser_init_rv(path, offset, len);
	_parser_loop_rv(_del, _ins, _nop, _match);
	*b = '\0';
	return(b - buf);

	#undef _del
	#undef _ins
	#undef _match
	#undef _nop
}

/**
 * @macro GABA_SEQ_FW, GABA_SEQ_RV, GABA_SEQ_A, GABA_SEQ_B
 */
#define GABA_SEQ_FW				( 0x00 )
#define GABA_SEQ_RV				( 0x01 )
#define GABA_SEQ_A				( 0x00 )
#define GABA_SEQ_B				( 0x02 )

/**
 * @fn gaba_dump_seq_forward
 * @brief traverse the path in the forward direction and dump the sequence with appropriate gaps
 */
_GABA_PARSE_EXPORT_LEVEL
uint64_t gaba_dump_seq_forward(
	char *buf,
	uint64_t buf_size,
	uint32_t conf,				/* { SEQ_A, SEQ_B } x { SEQ_FW, SEQ_RV } */
	uint32_t const *path,
	uint64_t offset,
	uint64_t len,
	uint8_t const *seq,			/* a->seq[s->alen] when SEQ_RV */
	char gap)					/* gap char, '-' */
{
	#define _gap(c) { \
		r += c; \
		for(uint64_t i = c; i > 0; i -= _gaba_parse_min2(i, 16)) { \
			_storeu_v16i8(r - i, gv); \
		} \
	}
	#define _fw(c) { \
		q += c, r += c; \
		for(uint64_t i = c; i > 0; i -= _gaba_parse_min2(i, 16)) { \
			_storeu_v16i8(r - i, _shuf_v16i8(cv, _loadu_v16i8(q - i))); \
		} \
	}
	#define _rv(c) { \
		q -= c; r += c; \
		for(uint64_t i = c; i > 0; i -= _gaba_parse_min2(i, 16)) { \
			uint64_t l = _gaba_parse_min2(i, 16); \
			_storeu_v16i8(r - i, _shuf_v16i8(cv, \
				_swapn_v16i8(_loadu_v16i8(q + i - l), l) \
			)); \
		} \
	}
	#define _nop(c) { (void)(c); }

	v16i8_t const cv = _load_v16i8((conf & GABA_SEQ_RV) ? gaba_parse_ascii_rv : gaba_parse_ascii_fw);
	v16i8_t const gv = _set_v16i8(gap);
	uint8_t const *q = seq;
	char *r = buf;

	_parser_init_fw(path, offset, len);
	switch(conf) {
		case GABA_SEQ_A | GABA_SEQ_FW: _parser_loop_fw(_fw, _gap, _fw, _nop); break;
		case GABA_SEQ_A | GABA_SEQ_RV: _parser_loop_fw(_rv, _gap, _rv, _nop); break;
		case GABA_SEQ_B | GABA_SEQ_FW: _parser_loop_fw(_gap, _fw, _fw, _nop); break;
		case GABA_SEQ_B | GABA_SEQ_RV: _parser_loop_fw(_gap, _rv, _rv, _nop); break;
	}
	*r = '\0';
	return(r - buf);

	#undef _gap
	#undef _fw
	#undef _rv
	#undef _nop
}

/**
 * @fn gaba_dump_seq_reverse
 * @brief traverse the path in the reverse direction and dump the sequence with appropriate gaps
 */
_GABA_PARSE_EXPORT_LEVEL
uint64_t gaba_dump_seq_reverse(
	char *buf,
	uint64_t buf_size,
	uint32_t conf,				/* { SEQ_A, SEQ_B } x { SEQ_FW, SEQ_RV } */
	uint32_t const *path,
	uint32_t offset,
	uint32_t len,
	uint8_t const *seq,			/* a->seq[s->alen] when SEQ_RV */
	char gap)					/* gap char, '-' */
{
	#define _gap(c) { \
		r += c; \
		for(uint64_t i = c; i > 0; i -= _gaba_parse_min2(i, 16)) { \
			_storeu_v16i8(r - i, gv); \
		} \
	}
	#define _fw(c) { \
		q += c, r += c; \
		for(uint64_t i = c; i > 0; i -= _gaba_parse_min2(i, 16)) { \
			_storeu_v16i8(r - i, _shuf_v16i8(cv, _loadu_v16i8(q - i))); \
		} \
	}
	#define _rv(c) { \
		q -= c; r += c; \
		for(uint64_t i = c; i > 0; i -= _gaba_parse_min2(i, 16)) { \
			uint64_t l = _gaba_parse_min2(i, 16); \
			_storeu_v16i8(r - i, _shuf_v16i8(cv, \
				_swapn_v16i8(_loadu_v16i8(q + i - l), l) \
			)); \
		} \
	}
	#define _nop(c) { (void)(c); }

	v16i8_t const cv = _load_v16i8((conf & GABA_SEQ_RV) ? gaba_parse_ascii_rv : gaba_parse_ascii_fw);
	v16i8_t const gv = _set_v16i8(gap);
	uint8_t const *q = seq;
	char *r = buf;

	_parser_init_rv(path, offset, len);
	switch(conf) {
		case GABA_SEQ_A | GABA_SEQ_FW: _parser_loop_rv(_fw, _gap, _fw, _nop); break;
		case GABA_SEQ_A | GABA_SEQ_RV: _parser_loop_rv(_rv, _gap, _rv, _nop); break;
		case GABA_SEQ_B | GABA_SEQ_FW: _parser_loop_rv(_gap, _fw, _fw, _nop); break;
		case GABA_SEQ_B | GABA_SEQ_RV: _parser_loop_rv(_gap, _rv, _rv, _nop); break;
		default: break;
	}
	*r = '\0';
	return(r - buf);

	#undef _gap
	#undef _fw
	#undef _rv
	#undef _nop
}

/**
 * @fn gaba_dump_seq_ref, gaba_dump_seq_query
 */
_GABA_PARSE_EXPORT_LEVEL
uint64_t gaba_dump_seq_ref(
	char *buf,
	uint64_t buf_size,
	uint32_t const *path,
	gaba_path_section_t const *s,
	gaba_section_t const *a)
{
	return(gaba_dump_seq_forward(
		buf, buf_size,
		GABA_SEQ_A | (a->base < GABA_EOU ? GABA_SEQ_FW : GABA_SEQ_RV),
		path, s->ppos, gaba_plen(s),
		a->base < GABA_EOU ? &a->base[s->apos] : gaba_mirror(&a->base[s->apos], 0),
		'-'
	));
}
_GABA_PARSE_EXPORT_LEVEL
uint64_t gaba_dump_seq_query(
	char *buf,
	uint64_t buf_size,
	uint32_t const *path,
	gaba_path_section_t const *s,
	gaba_section_t const *b)
{
	return(gaba_dump_seq_forward(
		buf, buf_size,
		GABA_SEQ_B | (b->base < GABA_EOU ? GABA_SEQ_FW : GABA_SEQ_RV),
		path, s->ppos, gaba_plen(s),
		b->base < GABA_EOU ? &b->base[s->apos] : gaba_mirror(&b->base[s->apos], 0),
		'-'
	));
}

#endif /* _GABA_PARSE_H_INCLUDED */
/**
 * end of gaba_parse.h
 */
