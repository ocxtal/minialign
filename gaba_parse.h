
/**
 * @file gaba_parse.h
 *
 * @brief libgaba utility function implementations
 *
 * @author Hajime Suzuki
 * @date 2018/1/20
 * @license Apache v2
 */

#ifndef _GABA_PARSE_H_INCLUDED
#define _GABA_PARSE_H_INCLUDED
#include <stdint.h>				/* uint32_t, uint64_t, ... */
#include <stddef.h>				/* ptrdiff_t */
#include "gaba.h"
#include "arch/arch.h"


#define _gaba_parse_min2(_x, _y)		( (_x) < (_y) ? (_x) : (_y) )

/**
 * @type gaba_dp_printer_t
 * @brief callback to CIGAR element printer
 * called with a pair of cigar operation (c) and its length (len).
 * void *fp is an opaque pointer to the context of the printer.
 */
typedef int (*gaba_dp_printer_t)(void *, uint64_t, char);

/**
 * @fn gaba_parse_load_uint64
 */
static inline
uint64_t gaba_parse_load_uint64(
	uint64_t const *ptr,
	int64_t pos)
{
	int64_t rem = pos & 63;
	uint64_t a = (ptr[pos>>6]>>rem) | ((ptr[(pos>>6) + 1]<<(63 - rem))<<1);
	return(a);
}

/**
 * @fn gaba_parse_dump_match_string
 */
static inline
int64_t gaba_parse_dump_match_string(
	char *buf,
	int64_t len)
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
		*p++ = 'M';
		return(p - buf);
	} else {
		int64_t adv;
		uint8_t b[16] = { 'M', '0' }, *p = &b[1];
		while(len != 0) { *p++ = (len % 10) + '0'; len /= 10; }
		for(p -= (p != &b[1]), adv = (int64_t)((ptrdiff_t)(p - b)) + 1; p >= b; p--) { *buf++ = *p; }
		return(adv);
	}
}

/**
 * @fn gaba_parse_dump_gap_string
 */
static inline
int64_t gaba_parse_dump_gap_string(
	char *buf,
	int64_t len,
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
		int64_t adv;
		uint8_t b[16] = { ch, '0' }, *p = &b[1];
		while(len != 0) { *p++ = (len % 10) + '0'; len /= 10; }
		for(p -= (p != &b[1]), adv = (int64_t)((ptrdiff_t)(p - b)) + 1; p >= b; p--) { *buf++ = *p; }
		return(adv);
	}
}

/**
 * @macro _gaba_parse_count_match_forward, _gaba_parse_count_gap_forward
 */
#define _gaba_parse_count_match_forward(_arr) ({ \
	tzcnt((_arr) ^ 0x5555555555555555); \
})
#define _gaba_parse_count_gap_forward(_arr) ({ \
	uint64_t _a = (_arr); \
	uint64_t mask = 0ULL - (_a & 0x01); \
	uint64_t gc = tzcnt(_a ^ mask) + (uint64_t)mask; \
	gc; \
})

/**
 * @fn gaba_dp_print_cigar_forward
 * @brief parse path string and print cigar to file
 */
static inline
uint64_t gaba_dp_print_cigar_forward(
	gaba_dp_printer_t printer,
	void *fp,
	uint32_t const *path,
	uint32_t offset,
	uint32_t len)
{
	uint64_t clen = 0;

	/* convert path to uint64_t pointer */
	uint64_t const *p = (uint64_t const *)((uint64_t)path & ~(sizeof(uint64_t) - 1));
	uint64_t lim = offset + (((uint64_t)path & sizeof(uint32_t)) ? 32 : 0) + len;
	uint64_t ridx = len;

	while(1) {
		uint64_t rsidx = ridx;
		while(1) {
			uint64_t m = _gaba_parse_count_match_forward(gaba_parse_load_uint64(p, lim - ridx));
			uint64_t a = _gaba_parse_min2(m, ridx) & ~0x01;
			ridx -= a;
			ZCNT_RESULT uint64_t c = a;
			if(c < 64) { break; }
		}
		uint64_t m = (rsidx - ridx)>>1;
		if(m > 0) {
			clen += printer(fp, m, 'M');
		}
		if(ridx == 0) { break; }

		uint64_t arr;
		uint64_t g = _gaba_parse_min2(
			_gaba_parse_count_gap_forward(arr = gaba_parse_load_uint64(p, lim - ridx)),
			ridx);
		if(g > 0) {
			clen += printer(fp, g, 'D' + ((char)(0ULL - (arr & 0x01)) & ('I' - 'D')));
		}
		if((ridx -= g) <= 1) { break; }
	}
	return(clen);
}

/**
 * @fn gaba_dp_dump_cigar_forward
 * @brief parse path string and store cigar to buffer
 */
static inline
uint64_t gaba_dp_dump_cigar_forward(
	char *buf,
	uint64_t buf_size,
	uint32_t const *path,
	uint32_t offset,
	uint32_t len)
{
	uint64_t const filled_len_margin = 5;
	char *b = buf, *blim = buf + buf_size - filled_len_margin;

	/* convert path to uint64_t pointer */
	uint64_t const *p = (uint64_t const *)((uint64_t)path & ~(sizeof(uint64_t) - 1));
	uint64_t lim = offset + (((uint64_t)path & sizeof(uint32_t)) ? 32 : 0) + len;
	uint64_t ridx = len;

	while(1) {
		uint64_t rsidx = ridx;
		while(1) {
			uint64_t m = _gaba_parse_count_match_forward(gaba_parse_load_uint64(p, lim - ridx));
			uint64_t a = _gaba_parse_min2(m, ridx) & ~0x01;
			ridx -= a;
			ZCNT_RESULT uint64_t c = a;
			if(c < 64) { break; }
		}
		uint64_t m = (rsidx - ridx)>>1;
		if(m > 0) {
			b += gaba_parse_dump_match_string(b, m);
		}
		if(ridx == 0 || b > blim) { break; }

		uint64_t arr;
		uint64_t g = _gaba_parse_min2(
			_gaba_parse_count_gap_forward(arr = gaba_parse_load_uint64(p, lim - ridx)),
			ridx);
		if(g > 0) {
			b += gaba_parse_dump_gap_string(b, g, 'D' + ((char)(0ULL - (arr & 0x01)) & ('I' - 'D')));
		}
		if((ridx -= g) <= 1 || b > blim) { break; }
	}
	*b = '\0';
	return(b - buf);
}

/**
 * @macro _gaba_parse_count_match_reverse, _gaba_parse_count_gap_reverse
 */
#define _gaba_parse_count_match_reverse(_arr) ({ \
	lzcnt((_arr) ^ 0x5555555555555555); \
})
#define _gaba_parse_count_gap_reverse(_arr) ({ \
	uint64_t _a = (_arr); \
	uint64_t mask = (uint64_t)(((int64_t)_a)>>63); \
	uint64_t gc = lzcnt(_a ^ mask) - ((int64_t)mask + 1); \
	gc; \
})

/**
 * @fn gaba_dp_print_cigar_reverse
 * @brief parse path string and print cigar to file
 */
static inline
uint64_t gaba_dp_print_cigar_reverse(
	gaba_dp_printer_t printer,
	void *fp,
	uint32_t const *path,
	uint32_t offset,
	uint32_t len)
{
	int64_t clen = 0;

	/* convert path to uint64_t pointer */
	uint64_t const *p = (uint64_t const *)((uint64_t)path & ~(sizeof(uint64_t) - 1));
	uint64_t ofs = (int64_t)offset + (((uint64_t)path & sizeof(uint32_t)) ? -32 : -64);
	uint64_t idx = len;

	while(1) {
		uint64_t sidx = idx;
		while(1) {
			uint64_t m = _gaba_parse_count_match_reverse(gaba_parse_load_uint64(p, idx + ofs));
			uint64_t a = _gaba_parse_min2(m, idx) & ~0x01;
			idx -= a;
			if(a < 64) { break; }
		}
		uint64_t m = (sidx - idx)>>1;
		if(m > 0) {
			clen += printer(fp, m, 'M');
		}
		if(idx == 0) { break; }

		uint64_t arr;
		uint64_t g = _gaba_parse_min2(
			_gaba_parse_count_gap_reverse(arr = gaba_parse_load_uint64(p, idx + ofs)),
			idx);
		if(g > 0) {
			clen += printer(fp, g, 'D' + ((char)(((int64_t)arr)>>63) & ('I' - 'D')));
		}
		if((idx -= g) <= 1) { break; }
	}
	return(clen);
}

/**
 * @fn gaba_dp_dump_cigar_reverse
 * @brief parse path string and store cigar to buffer
 */
static inline
uint64_t gaba_dp_dump_cigar_reverse(
	char *buf,
	uint64_t buf_size,
	uint32_t const *path,
	uint32_t offset,
	uint32_t len)
{
	uint64_t const filled_len_margin = 5;
	char *b = buf, *blim = buf + buf_size - filled_len_margin;

	/* convert path to uint64_t pointer */
	uint64_t const *p = (uint64_t const *)((uint64_t)path & ~(sizeof(uint64_t) - 1));
	uint64_t ofs = (int64_t)offset + (((uint64_t)path & sizeof(uint32_t)) ? -32 : -64);
	uint64_t idx = len;

	while(1) {
		uint64_t sidx = idx;
		while(1) {
			uint64_t m = _gaba_parse_count_match_reverse(gaba_parse_load_uint64(p, idx + ofs));
			uint64_t a = _gaba_parse_min2(m, idx) & ~0x01;
			idx -= a;
			if(a < 64) { break; }
		}
		uint64_t m = (sidx - idx)>>1;
		if(m > 0) {
			b += gaba_parse_dump_match_string(b, m);
		}
		if(idx == 0 || b > blim) { break; }

		uint64_t arr;
		uint64_t g = _gaba_parse_min2(
			_gaba_parse_count_gap_reverse(arr = gaba_parse_load_uint64(p, idx + ofs)),
			idx);
		if(g > 0) {
			b += gaba_parse_dump_gap_string(b, g, 'D' + ((char)(((int64_t)arr)>>63) & ('I' - 'D')));
		}
		if((idx -= g) <= 1 || b > blim) { break; }
	}
	*b = '\0';
	return(b - buf);
}

#endif /* _GABA_PARSE_H_INCLUDED */
/**
 * end of gaba_parse.h
 */
