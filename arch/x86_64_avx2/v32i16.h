
/**
 * @file v32i16.h
 *
 * @brief struct and _Generic based vector class implementation
 */
#ifndef _V32I16_H_INCLUDED
#define _V32I16_H_INCLUDED

/* include header for intel / amd sse2 instruction sets */
#include <x86intrin.h>

/* 16bit 32cell */
typedef struct v32i16_s {
	__m256i v1;
	__m256i v2;
} v32i16_t;

/* expanders (without argument) */
#define _e_x_v32i16_1(u)
#define _e_x_v32i16_2(u)

/* expanders (without immediate) */
#define _e_v_v32i16_1(a)				(a).v1
#define _e_v_v32i16_2(a)				(a).v2
#define _e_vv_v32i16_1(a, b)			(a).v1, (b).v1
#define _e_vv_v32i16_2(a, b)			(a).v2, (b).v2
#define _e_vvv_v32i16_1(a, b, c)		(a).v1, (b).v1, (c).v1
#define _e_vvv_v32i16_2(a, b, c)		(a).v2, (b).v2, (c).v2

/* expanders with immediate */
#define _e_i_v32i16_1(imm)			(imm)
#define _e_i_v32i16_2(imm)			(imm)
#define _e_vi_v32i16_1(a, imm)		(a).v1, (imm)
#define _e_vi_v32i16_2(a, imm)		(a).v2, (imm)
#define _e_vvi_v32i16_1(a, b, imm)	(a).v1, (b).v1, (imm)
#define _e_vvi_v32i16_2(a, b, imm)	(a).v2, (b).v2, (imm)

/* address calculation macros */
#define _addr_v32i16_1(imm)			( (__m256i *)(imm) )
#define _addr_v32i16_2(imm)			( (__m256i *)(imm) + 1 )
#define _pv_v32i16(ptr)				( _addr_v32i16_1(ptr) )
/* expanders with pointers */
#define _e_p_v32i16_1(ptr)			_addr_v32i16_1(ptr)
#define _e_p_v32i16_2(ptr)			_addr_v32i16_2(ptr)
#define _e_pv_v32i16_1(ptr, a)		_addr_v32i16_1(ptr), (a).v1
#define _e_pv_v32i16_2(ptr, a)		_addr_v32i16_2(ptr), (a).v2

/* expand intrinsic name */
#define _i_v32i16(intrin) 			_mm256_##intrin##_epi16
#define _i_v32i16x(intrin)			_mm256_##intrin##_si256

/* apply */
#define _a_v32i16(intrin, expander, ...) ( \
	(v32i16_t) { \
		_i_v32i16(intrin)(expander##_v32i16_1(__VA_ARGS__)), \
		_i_v32i16(intrin)(expander##_v32i16_2(__VA_ARGS__)) \
	} \
)
#define _a_v32i16x(intrin, expander, ...) ( \
	(v32i16_t) { \
		_i_v32i16x(intrin)(expander##_v32i16_1(__VA_ARGS__)), \
		_i_v32i16x(intrin)(expander##_v32i16_2(__VA_ARGS__)) \
	} \
)
#define _a_v32i16xv(intrin, expander, ...) { \
	_i_v32i16x(intrin)(expander##_v32i16_1(__VA_ARGS__)); \
	_i_v32i16x(intrin)(expander##_v32i16_2(__VA_ARGS__)); \
}

/* load and store */
#define _load_v32i16(...)	_a_v32i16x(load, _e_p, __VA_ARGS__)
#define _loadu_v32i16(...)	_a_v32i16x(loadu, _e_p, __VA_ARGS__)
#define _store_v32i16(...)	_a_v32i16xv(store, _e_pv, __VA_ARGS__)
#define _storeu_v32i16(...)	_a_v32i16xv(storeu, _e_pv, __VA_ARGS__)

/* broadcast */
#define _set_v32i16(...)	_a_v32i16(set1, _e_i, __VA_ARGS__)
#define _zero_v32i16()		_a_v32i16x(setzero, _e_x, _unused)

/* logics */
#define _not_v32i16(...)	_a_v32i16x(not, _e_v, __VA_ARGS__)
#define _and_v32i16(...)	_a_v32i16x(and, _e_vv, __VA_ARGS__)
#define _or_v32i16(...)		_a_v32i16x(or, _e_vv, __VA_ARGS__)
#define _xor_v32i16(...)	_a_v32i16x(xor, _e_vv, __VA_ARGS__)
#define _andn_v32i16(...)	_a_v32i16x(andnot, _e_vv, __VA_ARGS__)

/* arithmetics */
#define _add_v32i16(...)	_a_v32i16(add, _e_vv, __VA_ARGS__)
#define _sub_v32i16(...)	_a_v32i16(sub, _e_vv, __VA_ARGS__)
#define _max_v32i16(...)	_a_v32i16(max, _e_vv, __VA_ARGS__)
#define _min_v32i16(...)	_a_v32i16(min, _e_vv, __VA_ARGS__)

/* compare */
#define _eq_v32i16(...)		_a_v32i16(cmpeq, _e_vv, __VA_ARGS__)
#define _gt_v32i16(...)		_a_v32i16(cmpgt, _e_vv, __VA_ARGS__)

/* insert and extract */
#define _ins_v32i16(a, val, imm) { \
	if((imm) < sizeof(__m256i)/sizeof(int16_t)) { \
		(a).v1 = _i_v32i16(insert)((a).v1, (val), (imm)); \
	} else if((imm) < 2*sizeof(__m256i)/sizeof(int16_t)) { \
		(a).v2 = _i_v32i16(insert)((a).v2, (val), (imm) - sizeof(__m256i)/sizeof(int16_t)); \
	} \
}
#define _ext_v32i16(a, imm) ( \
	(int16_t)(((imm) < sizeof(__m256i)/sizeof(int16_t)) \
		? _i_v32i16(extract)((a).v1, (imm)) \
		: _i_v32i16(extract)((a).v2, (imm) - sizeof(__m256i)/sizeof(int16_t))) \
)

/* mask */
#define _mask_v32i16(a) ( \
	(v32_mask_t) { \
		.m1 = _mm256_movemask_epi8( \
			_mm256_permute4x64_epi64( \
				_mm256_packs_epi16((a).v1, (a).v2), 0xd8)) \
	} \
)

/* horizontal max (reduction max) */
#define _hmax_v32i16(a) ({ \
	__m256i _s = _mm256_max_epi16((a).v1, (a).v2); \
	__m128i _t = _mm_max_epi16( \
		_mm256_castsi256_si128(_s), \
		_mm256_extracti128_si256(_s, 1) \
	); \
	_t = _mm_max_epi16(_t, _mm_srli_si128(_t, 8)); \
	_t = _mm_max_epi16(_t, _mm_srli_si128(_t, 4)); \
	_t = _mm_max_epi16(_t, _mm_srli_si128(_t, 2)); \
	(int16_t)_mm_extract_epi16(_t, 0); \
})

#define _cvt_v32i8_v32i16(a) ( \
	(v32i16_t) { \
		_mm256_cvtepi8_epi16(_mm256_castsi256_si128((a).v1)), \
		_mm256_cvtepi8_epi16(_mm256_extracti128_si256((a).v1, 1)) \
	} \
)

/* debug print */
// #ifdef _LOG_H_INCLUDED
#define _print_v32i16(a) { \
	debug("(v32i16_t) %s(%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, " \
				  "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d)", \
		#a, \
		_ext_v32i16(a, 31), \
		_ext_v32i16(a, 30), \
		_ext_v32i16(a, 29), \
		_ext_v32i16(a, 28), \
		_ext_v32i16(a, 27), \
		_ext_v32i16(a, 26), \
		_ext_v32i16(a, 25), \
		_ext_v32i16(a, 24), \
		_ext_v32i16(a, 23), \
		_ext_v32i16(a, 22), \
		_ext_v32i16(a, 21), \
		_ext_v32i16(a, 20), \
		_ext_v32i16(a, 19), \
		_ext_v32i16(a, 18), \
		_ext_v32i16(a, 17), \
		_ext_v32i16(a, 16), \
		_ext_v32i16(a, 15), \
		_ext_v32i16(a, 14), \
		_ext_v32i16(a, 13), \
		_ext_v32i16(a, 12), \
		_ext_v32i16(a, 11), \
		_ext_v32i16(a, 10), \
		_ext_v32i16(a, 9), \
		_ext_v32i16(a, 8), \
		_ext_v32i16(a, 7), \
		_ext_v32i16(a, 6), \
		_ext_v32i16(a, 5), \
		_ext_v32i16(a, 4), \
		_ext_v32i16(a, 3), \
		_ext_v32i16(a, 2), \
		_ext_v32i16(a, 1), \
		_ext_v32i16(a, 0)); \
}
// #else
// #define _print_v32i16(x)	;
// #endif

#endif /* _V32I16_H_INCLUDED */
/**
 * end of v32i16.h
 */
