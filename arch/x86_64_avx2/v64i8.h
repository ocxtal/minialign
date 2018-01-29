
/**
 * @file v64i8.h
 *
 * @brief struct and _Generic based vector class implementation
 */
#ifndef _V64I8_H_INCLUDED
#define _V64I8_H_INCLUDED

/* include header for intel / amd sse2 instruction sets */
#include <x86intrin.h>

/* 8bit 64cell */
typedef struct v64i8_s {
	__m256i v1;
	__m256i v2;
} v64i8_t;

/* expanders (without argument) */
#define _e_x_v64i8_1(u)
#define _e_x_v64i8_2(u)

/* expanders (without immediate) */
#define _e_v_v64i8_1(a)				(a).v1
#define _e_v_v64i8_2(a)				(a).v2
#define _e_vv_v64i8_1(a, b)			(a).v1, (b).v1
#define _e_vv_v64i8_2(a, b)			(a).v2, (b).v2
#define _e_vvv_v64i8_1(a, b, c)		(a).v1, (b).v1, (c).v1
#define _e_vvv_v64i8_2(a, b, c)		(a).v2, (b).v2, (c).v2

/* expanders with immediate */
#define _e_i_v64i8_1(imm)			(imm)
#define _e_i_v64i8_2(imm)			(imm)
#define _e_vi_v64i8_1(a, imm)		(a).v1, (imm)
#define _e_vi_v64i8_2(a, imm)		(a).v2, (imm)
#define _e_vvi_v64i8_1(a, b, imm)	(a).v1, (b).v1, (imm)
#define _e_vvi_v64i8_2(a, b, imm)	(a).v2, (b).v2, (imm)

/* address calculation macros */
#define _addr_v64i8_1(imm)			( (__m256i *)(imm) )
#define _addr_v64i8_2(imm)			( (__m256i *)(imm) + 1 )
#define _pv_v64i8(ptr)				( _addr_v64i8_1(ptr) )

/* expanders with pointers */
#define _e_p_v64i8_1(ptr)			_addr_v64i8_1(ptr)
#define _e_p_v64i8_2(ptr)			_addr_v64i8_2(ptr)
#define _e_pv_v64i8_1(ptr, a)		_addr_v64i8_1(ptr), (a).v1
#define _e_pv_v64i8_2(ptr, a)		_addr_v64i8_2(ptr), (a).v2

/* expand intrinsic name */
#define _i_v64i8(intrin) 			_mm256_##intrin##_epi8
#define _i_v64i8x(intrin)			_mm256_##intrin##_si256

/* apply */
#define _a_v64i8(intrin, expander, ...) ( \
	(v64i8_t) { \
		_i_v64i8(intrin)(expander##_v64i8_1(__VA_ARGS__)), \
		_i_v64i8(intrin)(expander##_v64i8_2(__VA_ARGS__)) \
	} \
)
#define _a_v64i8x(intrin, expander, ...) ( \
	(v64i8_t) { \
		_i_v64i8x(intrin)(expander##_v64i8_1(__VA_ARGS__)), \
		_i_v64i8x(intrin)(expander##_v64i8_2(__VA_ARGS__)) \
	} \
)
#define _a_v64i8xv(intrin, expander, ...) { \
	_i_v64i8x(intrin)(expander##_v64i8_1(__VA_ARGS__)); \
	_i_v64i8x(intrin)(expander##_v64i8_2(__VA_ARGS__)); \
}

/* load and store */
#define _load_v64i8(...)	_a_v64i8x(load, _e_p, __VA_ARGS__)
#define _loadu_v64i8(...)	_a_v64i8x(loadu, _e_p, __VA_ARGS__)
#define _store_v64i8(...)	_a_v64i8xv(store, _e_pv, __VA_ARGS__)
#define _storeu_v64i8(...)	_a_v64i8xv(storeu, _e_pv, __VA_ARGS__)

/* broadcast */
#define _set_v64i8(...)		_a_v64i8(set1, _e_i, __VA_ARGS__)
#define _zero_v64i8()		_a_v64i8x(setzero, _e_x, _unused)

/* swap (reverse) */
#define _swap_idx_v64i8() ( \
	_mm256_broadcastsi128_si256(_mm_set_epi8( \
		0, 1, 2, 3, 4, 5, 6, 7, \
		8, 9, 10, 11, 12, 13, 14, 15)) \
)
#define _swap_v64i8(a) ( \
	(v64i8_t) { \
		_mm256_permute2x128_si256( \
			_mm256_shuffle_epi8((a).v2, _swap_idx_v64i8()), \
			_mm256_shuffle_epi8((a).v2, _swap_idx_v64i8()), \
			0x01), \
		_mm256_permute2x128_si256( \
			_mm256_shuffle_epi8((a).v1, _swap_idx_v64i8()), \
			_mm256_shuffle_epi8((a).v1, _swap_idx_v64i8()), \
			0x01) \
	} \
)

/* logics */
#define _not_v64i8(...)		_a_v64i8x(not, _e_v, __VA_ARGS__)
#define _and_v64i8(...)		_a_v64i8x(and, _e_vv, __VA_ARGS__)
#define _or_v64i8(...)		_a_v64i8x(or, _e_vv, __VA_ARGS__)
#define _xor_v64i8(...)		_a_v64i8x(xor, _e_vv, __VA_ARGS__)
#define _andn_v64i8(...)	_a_v64i8x(andnot, _e_vv, __VA_ARGS__)

/* arithmetics */
#define _add_v64i8(...)		_a_v64i8(add, _e_vv, __VA_ARGS__)
#define _sub_v64i8(...)		_a_v64i8(sub, _e_vv, __VA_ARGS__)
#define _adds_v64i8(...)	_a_v64i8(adds, _e_vv, __VA_ARGS__)
#define _subs_v64i8(...)	_a_v64i8(subs, _e_vv, __VA_ARGS__)
#define _max_v64i8(...)		_a_v64i8(max, _e_vv, __VA_ARGS__)
#define _min_v64i8(...)		_a_v64i8(min, _e_vv, __VA_ARGS__)

/* shuffle */
#define _shuf_v64i8(...)	_a_v64i8(shuffle, _e_vv, __VA_ARGS__)

/* blend */
#define _sel_v64i8(...)		_a_v64i8(blendv, _e_vvv, __VA_ARGS__)

/* compare */
#define _eq_v64i8(...)		_a_v64i8(cmpeq, _e_vv, __VA_ARGS__)
#define _gt_v64i8(...)		_a_v64i8(cmpgt, _e_vv, __VA_ARGS__)

/* insert and extract */
#define _ins_v64i8(a, val, imm) { \
	if((imm) < sizeof(__m256i)) { \
		(a).v1 = _i_v64i8(insert)((a).v1, (val), (imm)); \
	} else { \
		(a).v2 = _i_v64i8(insert)((a).v2, (val), (imm) - sizeof(__m256i)); \
	} \
}
#define _ext_v64i8(a, imm) ( \
	(int8_t)(((imm) < sizeof(__m256i)) ? ( \
		_i_v64i8(extract)((a).v1, (imm)) \
	) : ( \
		_i_v64i8(extract)((a).v2, (imm) - sizeof(__m256i)) \
	)) \
)

/* shift */
#define _bsl_v64i8(a, imm) ( \
	(v64i8_t) { \
		_mm256_alignr_epi8( \
			(a).v1, \
			_mm256_permute2x128_si256((a).v1, (a).v2, 0x04), \
			15), \
		_mm256_alignr_epi8( \
			(a).v2, \
			_mm256_permute2x128_si256((a).v1, (a).v2, 0x21), \
			15) \
	} \
)
#define _bsr_v64i8(a, imm) ( \
	(v64i8_t) { \
		_mm256_alignr_epi8( \
			_mm256_permute2x128_si256((a).v1, (a).v2, 0x21), \
			(a).v1, \
			1), \
		_mm256_alignr_epi8( \
			_mm256_castsi128_si256(_mm256_extracti128_si256((a).v2, 1)), \
			(a).v2, \
			1) \
	} \
)
#define _shl_v64i8(a, imm) ( \
	(v64i8_t) { \
		_mm256_slli_epi32((a).v1, (imm)), \
		_mm256_slli_epi32((a).v2, (imm)) \
	} \
)
#define _shr_v64i8(a, imm) ( \
	(v64i8_t) { \
		_mm256_srli_epi32((a).v1, (imm)), \
		_mm256_srli_epi32((a).v2, (imm)) \
	} \
)
#define _sal_v64i8(a, imm) ( \
	(v64i8_t) { \
		_mm256_slai_epi32((a).v1, (imm)), \
		_mm256_slai_epi32((a).v2, (imm)) \
	} \
)
#define _sar_v64i8(a, imm) ( \
	(v64i8_t) { \
		_mm256_srai_epi32((a).v1, (imm)), \
		_mm256_srai_epi32((a).v2, (imm)), \
	} \
)

/* mask */
#define _mask_v64i8(a) ( \
	(v64_mask_t) { \
		.m1 = _i_v64i8(movemask)((a).v1), \
		.m2 = _i_v64i8(movemask)((a).v2) \
	} \
)

/* convert */
#define _cvt_v64i16_v64i8(a) ( \
	(v64i8_t) { \
		_mm256_packs_epi16( \
			_mm256_permute2x128_si256((a).v1, (a).v2, 0x20), \
			_mm256_permute2x128_si256((a).v1, (a).v2, 0x31)), \
		_mm256_packs_epi16( \
			_mm256_permute2x128_si256((a).v3, (a).v4, 0x20), \
			_mm256_permute2x128_si256((a).v3, (a).v4, 0x31)) \
	} \
)

/* debug print */
// #ifdef _LOG_H_INCLUDED
#define _print_v64i8(a) { \
	debug("(v64i8_t) %s(%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, " \
				 "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, " \
				 "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, " \
				 "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d)", \
		#a, \
		_ext_v64i8(a, 32 + 31), \
		_ext_v64i8(a, 32 + 30), \
		_ext_v64i8(a, 32 + 29), \
		_ext_v64i8(a, 32 + 28), \
		_ext_v64i8(a, 32 + 27), \
		_ext_v64i8(a, 32 + 26), \
		_ext_v64i8(a, 32 + 25), \
		_ext_v64i8(a, 32 + 24), \
		_ext_v64i8(a, 32 + 23), \
		_ext_v64i8(a, 32 + 22), \
		_ext_v64i8(a, 32 + 21), \
		_ext_v64i8(a, 32 + 20), \
		_ext_v64i8(a, 32 + 19), \
		_ext_v64i8(a, 32 + 18), \
		_ext_v64i8(a, 32 + 17), \
		_ext_v64i8(a, 32 + 16), \
		_ext_v64i8(a, 32 + 15), \
		_ext_v64i8(a, 32 + 14), \
		_ext_v64i8(a, 32 + 13), \
		_ext_v64i8(a, 32 + 12), \
		_ext_v64i8(a, 32 + 11), \
		_ext_v64i8(a, 32 + 10), \
		_ext_v64i8(a, 32 + 9), \
		_ext_v64i8(a, 32 + 8), \
		_ext_v64i8(a, 32 + 7), \
		_ext_v64i8(a, 32 + 6), \
		_ext_v64i8(a, 32 + 5), \
		_ext_v64i8(a, 32 + 4), \
		_ext_v64i8(a, 32 + 3), \
		_ext_v64i8(a, 32 + 2), \
		_ext_v64i8(a, 32 + 1), \
		_ext_v64i8(a, 32 + 0), \
		_ext_v64i8(a, 31), \
		_ext_v64i8(a, 30), \
		_ext_v64i8(a, 29), \
		_ext_v64i8(a, 28), \
		_ext_v64i8(a, 27), \
		_ext_v64i8(a, 26), \
		_ext_v64i8(a, 25), \
		_ext_v64i8(a, 24), \
		_ext_v64i8(a, 23), \
		_ext_v64i8(a, 22), \
		_ext_v64i8(a, 21), \
		_ext_v64i8(a, 20), \
		_ext_v64i8(a, 19), \
		_ext_v64i8(a, 18), \
		_ext_v64i8(a, 17), \
		_ext_v64i8(a, 16), \
		_ext_v64i8(a, 15), \
		_ext_v64i8(a, 14), \
		_ext_v64i8(a, 13), \
		_ext_v64i8(a, 12), \
		_ext_v64i8(a, 11), \
		_ext_v64i8(a, 10), \
		_ext_v64i8(a, 9), \
		_ext_v64i8(a, 8), \
		_ext_v64i8(a, 7), \
		_ext_v64i8(a, 6), \
		_ext_v64i8(a, 5), \
		_ext_v64i8(a, 4), \
		_ext_v64i8(a, 3), \
		_ext_v64i8(a, 2), \
		_ext_v64i8(a, 1), \
		_ext_v64i8(a, 0)); \
}
// #else
// #define _print_v64i8(x)		;
// #endif

#endif /* _V64I8_H_INCLUDED */
/**
 * end of v64i8.h
 */
