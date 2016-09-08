
/**
 * @file v32i8.h
 *
 * @brief struct and _Generic based vector class implementation
 */
#ifndef _V32I8_H_INCLUDED
#define _V32I8_H_INCLUDED

/* include header for intel / amd sse2 instruction sets */
#include <x86intrin.h>

/* 8bit 32cell */
typedef struct v32i8_s {
	__m256i v1;
} v32i8_t;

/* expanders (without argument) */
#define _e_x_v32i8_1(u)

/* expanders (without immediate) */
#define _e_v_v32i8_1(a)				(a).v1
#define _e_vv_v32i8_1(a, b)			(a).v1, (b).v1
#define _e_vvv_v32i8_1(a, b, c)		(a).v1, (b).v1, (c).v1

/* expanders with immediate */
#define _e_i_v32i8_1(imm)			(imm)
#define _e_vi_v32i8_1(a, imm)		(a).v1, (imm)
#define _e_vvi_v32i8_1(a, b, imm)	(a).v1, (b).v1, (imm)

/* address calculation macros */
#define _addr_v32i8_1(imm)			( (__m256i *)(imm) )
#define _pv_v32i8(ptr)				( _addr_v32i8_1(ptr) )
/* expanders with pointers */
#define _e_p_v32i8_1(ptr)			_addr_v32i8_1(ptr)
#define _e_pv_v32i8_1(ptr, a)		_addr_v32i8_1(ptr), (a).v1

/* expand intrinsic name */
#define _i_v32i8(intrin) 			_mm256_##intrin##_epi8
#define _i_v32i8x(intrin)			_mm256_##intrin##_si256

/* apply */
#define _a_v32i8(intrin, expander, ...) ( \
	(v32i8_t) { \
		_i_v32i8(intrin)(expander##_v32i8_1(__VA_ARGS__)) \
	} \
)
#define _a_v32i8x(intrin, expander, ...) ( \
	(v32i8_t) { \
		_i_v32i8x(intrin)(expander##_v32i8_1(__VA_ARGS__)) \
	} \
)
#define _a_v32i8xv(intrin, expander, ...) { \
	_i_v32i8x(intrin)(expander##_v32i8_1(__VA_ARGS__)); \
}

/* load and store */
#define _load_v32i8(...)	_a_v32i8x(load, _e_p, __VA_ARGS__)
#define _loadu_v32i8(...)	_a_v32i8x(loadu, _e_p, __VA_ARGS__)
#define _store_v32i8(...)	_a_v32i8xv(store, _e_pv, __VA_ARGS__)
#define _storeu_v32i8(...)	_a_v32i8xv(storeu, _e_pv, __VA_ARGS__)

/* broadcast */
#define _set_v32i8(...)		_a_v32i8(set1, _e_i, __VA_ARGS__)
#define _zero_v32i8()		_a_v32i8x(setzero, _e_x, _unused)

/* swap (reverse) */
#define _swap_idx_v32i8() ( \
	_mm256_broadcastsi128_si256(_mm_set_epi8( \
		0, 1, 2, 3, 4, 5, 6, 7, \
		8, 9, 10, 11, 12, 13, 14, 15)) \
)
#define _swap_v32i8(a) ( \
	(v32i8_t) { \
		_mm256_permute2x128_si256( \
			_mm256_shuffle_epi8((a).v1, _swap_idx_v32i8()), \
			_mm256_shuffle_epi8((a).v1, _swap_idx_v32i8()), \
			0x01) \
	} \
)

/* logics */
#define _not_v32i8(...)		_a_v32i8x(not, _e_v, __VA_ARGS__)
#define _and_v32i8(...)		_a_v32i8x(and, _e_vv, __VA_ARGS__)
#define _or_v32i8(...)		_a_v32i8x(or, _e_vv, __VA_ARGS__)
#define _xor_v32i8(...)		_a_v32i8x(xor, _e_vv, __VA_ARGS__)
#define _andn_v32i8(...)	_a_v32i8x(andnot, _e_vv, __VA_ARGS__)

/* arithmetics */
#define _add_v32i8(...)		_a_v32i8(add, _e_vv, __VA_ARGS__)
#define _sub_v32i8(...)		_a_v32i8(sub, _e_vv, __VA_ARGS__)
#define _adds_v32i8(...)	_a_v32i8(adds, _e_vv, __VA_ARGS__)
#define _subs_v32i8(...)	_a_v32i8(subs, _e_vv, __VA_ARGS__)
#define _max_v32i8(...)		_a_v32i8(max, _e_vv, __VA_ARGS__)
#define _min_v32i8(...)		_a_v32i8(min, _e_vv, __VA_ARGS__)

/* shuffle */
#define _shuf_v32i8(...)	_a_v32i8(shuffle, _e_vv, __VA_ARGS__)

/* blend */
// #define _sel_v32i8(...)		_a_v32i8(blendv, _e_vvv, __VA_ARGS__)

/* compare */
#define _eq_v32i8(...)		_a_v32i8(cmpeq, _e_vv, __VA_ARGS__)
#define _lt_v32i8(...)		_a_v32i8(cmplt, _e_vv, __VA_ARGS__)
#define _gt_v32i8(...)		_a_v32i8(cmpgt, _e_vv, __VA_ARGS__)

/* insert and extract */
#define _ins_v32i8(a, val, imm) { \
	(a).v1 = _i_v32i8(insert)((a).v1, (val), (imm)); \
}
#define _ext_v32i8(a, imm) ( \
	(int8_t)_i_v32i8(extract)((a).v1, (imm)) \
)

/* shift */
#define _bsl_v32i8(a, imm) ( \
	(v32i8_t) { \
		_mm256_alignr_epi8( \
			(a).v1, \
			_mm256_permute2x128_si256((a).v1, (a).v1, 0x08), \
			15) \
	} \
)
#define _bsr_v32i8(a, imm) ( \
	(v32i8_t) { \
		_mm256_alignr_epi8( \
			_mm256_castsi128_si256( \
				_mm256_extracti128_si256((a).v1, 1)), \
			(a).v1, \
			1) \
	} \
)
#define _shl_v32i8(a, imm) ( \
	(v32i8_t) { \
		_mm256_slli_epi32((a).v1, (imm)) \
	} \
)
#define _shr_v32i8(a, imm) ( \
	(v32i8_t) { \
		_mm256_srli_epi32((a).v1, (imm)) \
	} \
)
#define _sal_v32i8(a, imm) ( \
	(v32i8_t) { \
		_mm256_slai_epi32((a).v1, (imm)) \
	} \
)
#define _sar_v32i8(a, imm) ( \
	(v32i8_t) { \
		_mm256_srai_epi32((a).v1, (imm)) \
	} \
)

/* mask */
#define _mask_v32i8(a) ( \
	(v32_mask_t) { \
		.m1 = _i_v32i8(movemask)((a).v1) \
	} \
)

/* debug print */
#ifdef _LOG_H_INCLUDED
#define _print_v32i8(a) { \
	debug("(v32i8_t) %s(%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, " \
				 "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d)", \
		#a, \
		_ext_v32i8(a, 31), \
		_ext_v32i8(a, 30), \
		_ext_v32i8(a, 29), \
		_ext_v32i8(a, 28), \
		_ext_v32i8(a, 27), \
		_ext_v32i8(a, 26), \
		_ext_v32i8(a, 25), \
		_ext_v32i8(a, 24), \
		_ext_v32i8(a, 23), \
		_ext_v32i8(a, 22), \
		_ext_v32i8(a, 21), \
		_ext_v32i8(a, 20), \
		_ext_v32i8(a, 19), \
		_ext_v32i8(a, 18), \
		_ext_v32i8(a, 17), \
		_ext_v32i8(a, 16), \
		_ext_v32i8(a, 15), \
		_ext_v32i8(a, 14), \
		_ext_v32i8(a, 13), \
		_ext_v32i8(a, 12), \
		_ext_v32i8(a, 11), \
		_ext_v32i8(a, 10), \
		_ext_v32i8(a, 9), \
		_ext_v32i8(a, 8), \
		_ext_v32i8(a, 7), \
		_ext_v32i8(a, 6), \
		_ext_v32i8(a, 5), \
		_ext_v32i8(a, 4), \
		_ext_v32i8(a, 3), \
		_ext_v32i8(a, 2), \
		_ext_v32i8(a, 1), \
		_ext_v32i8(a, 0)); \
}
#else
#define _print_v32i8(x)		;
#endif

#endif /* _V32I8_H_INCLUDED */
/**
 * end of v32i8.h
 */
