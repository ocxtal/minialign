
/**
 * @file v4i32.h
 *
 * @brief struct and _Generic based vector class implementation
 */
#ifndef _V4I32_H_INCLUDED
#define _V4I32_H_INCLUDED

/* include header for intel / amd sse2 instruction sets */
#include <x86intrin.h>

/* 8bit 32cell */
typedef struct v4i32_s {
	__m128i v1;
} v4i32_t;

/* expanders (without argument) */
#define _e_x_v4i32_1(u)
#define _e_x_v4i32_2(u)

/* expanders (without immediate) */
#define _e_v_v4i32_1(a)				(a).v1
#define _e_v_v4i32_2(a)				(a).v1
#define _e_vv_v4i32_1(a, b)			(a).v1, (b).v1
#define _e_vv_v4i32_2(a, b)			(a).v1, (b).v1
#define _e_vvv_v4i32_1(a, b, c)		(a).v1, (b).v1, (c).v1
#define _e_vvv_v4i32_2(a, b, c)		(a).v1, (b).v1, (c).v1

/* expanders with immediate */
#define _e_i_v4i32_1(imm)			(imm)
#define _e_i_v4i32_2(imm)			(imm)
#define _e_vi_v4i32_1(a, imm)		(a).v1, (imm)
#define _e_vi_v4i32_2(a, imm)		(a).v1, (imm)
#define _e_vvi_v4i32_1(a, b, imm)	(a).v1, (b).v1, (imm)
#define _e_vvi_v4i32_2(a, b, imm)	(a).v1, (b).v1, (imm)

/* address calculation macros */
#define _addr_v4i32_1(imm)			( (__m128i *)(imm) )
#define _addr_v4i32_2(imm)			( (__m128i *)(imm) )
#define _pv_v4i32(ptr)				( _addr_v4i32_1(ptr) )
/* expanders with pointers */
#define _e_p_v4i32_1(ptr)			_addr_v4i32_1(ptr)
#define _e_p_v4i32_2(ptr)			_addr_v4i32_2(ptr)
#define _e_pv_v4i32_1(ptr, a)		_addr_v4i32_1(ptr), (a).v1
#define _e_pv_v4i32_2(ptr, a)		_addr_v4i32_2(ptr), (a).v1

/* expand intrinsic name */
#define _i_v4i32(intrin) 			_mm_##intrin##_epi32
#define _i_v4u32(intrin) 			_mm_##intrin##_epu32
#define _i_v4i32x(intrin)			_mm_##intrin##_si128

/* apply */
#define _a_v4i32(intrin, expander, ...) ( \
	(v4i32_t) { \
		_i_v4i32(intrin)(expander##_v4i32_1(__VA_ARGS__)) \
	} \
)
#define _a_v4u32(intrin, expander, ...) ( \
	(v4i32_t) { \
		_i_v4u32(intrin)(expander##_v4i32_1(__VA_ARGS__)) \
	} \
)
#define _a_v4i32x(intrin, expander, ...) ( \
	(v4i32_t) { \
		_i_v4i32x(intrin)(expander##_v4i32_1(__VA_ARGS__)) \
	} \
)
#define _a_v4i32xv(intrin, expander, ...) { \
	_i_v4i32x(intrin)(expander##_v4i32_1(__VA_ARGS__)); \
}

/* load and store */
#define _load_v4i32(...)	_a_v4i32x(load, _e_p, __VA_ARGS__)
#define _loadu_v4i32(...)	_a_v4i32x(load, _e_p, __VA_ARGS__)
#define _store_v4i32(...)	_a_v4i32xv(store, _e_pv, __VA_ARGS__)
#define _storeu_v4i32(...)	_a_v4i32xv(store, _e_pv, __VA_ARGS__)

/* broadcast */
#define _set_v4i32(...)		_a_v4i32(set1, _e_i, __VA_ARGS__)
#define _zero_v4i32()		_a_v4i32x(setzero, _e_x, _unused)
#define _seta_v4i32(i4, i3, i2, i1) ( \
	(v4i32_t) { \
		_mm_set_epi32(i4, i3, i2, i1) \
	} \
)
#define _swap_v4i32(x) ( \
	(v4i32_t) { \
		_mm_shuffle_epi32((x).v1, 0xe1) \
	} \
)

/* logics */
#define _not_v4i32(...)		_a_v4i32x(not, _e_v, __VA_ARGS__)
#define _and_v4i32(...)		_a_v4i32x(and, _e_vv, __VA_ARGS__)
#define _or_v4i32(...)		_a_v4i32x(or, _e_vv, __VA_ARGS__)
#define _xor_v4i32(...)		_a_v4i32x(xor, _e_vv, __VA_ARGS__)
#define _andn_v4i32(...)	_a_v4i32x(andnot, _e_vv, __VA_ARGS__)

/* arithmetics */
#define _add_v4i32(...)		_a_v4i32(add, _e_vv, __VA_ARGS__)
#define _sub_v4i32(...)		_a_v4i32(sub, _e_vv, __VA_ARGS__)
#define _max_v4i32(...)		_a_v4i32(max, _e_vv, __VA_ARGS__)
#define _min_v4i32(...)		_a_v4i32(min, _e_vv, __VA_ARGS__)
#define _maxu_v4i32(...)	_a_v4u32(max, _e_vv, __VA_ARGS__)
#define _minu_v4i32(...)	_a_v4u32(min, _e_vv, __VA_ARGS__)

/* blend: mask == 1 ? a : b */
#define _sel_v4i32(mask, a, b) ( \
	(v4i32_t) { \
		_mm_blendv_epi8((b).v1, (a).v1, (mask).v1) \
	} \
)

/* compare */
#define _eq_v4i32(...)		_a_v4i32(cmpeq, _e_vv, __VA_ARGS__)
#define _gt_v4i32(...)		_a_v4i32(cmpgt, _e_vv, __VA_ARGS__)

/* test: take mask and test if all zero */
#define _test_v4i32(x, y)	_mm_test_all_zeros((x).v1, (y).v1)

/* insert and extract */
#define _ins_v4i32(a, val, imm) { \
	(a).v1 = _i_v4i32((a).v1, (val), (imm)); \
}
#define _ext_v4i32(a, imm) ( \
	(int32_t)_i_v4i32(extract)((a).v1, (imm)) \
)

/* element shift */
#define _bsl_v4i32(a, imm) ( \
	(v4i32_t) { \
		_mm_slli_si128((a).v1, sizeof(int32_t) * (imm)) \
	} \
)
#define _bsr_v4i32(a, imm) ( \
	(v4i32_t) { \
		_mm_srli_si128((a).v1, sizeof(int32_t) * (imm)) \
	} \
)

/* double shift (palignr) */
#define _bsld_v4i32(a, b, imm) ( \
	(v4i32_t) { \
		_mm_alignr_epi8((a).v1, (b).v1, sizeof(__m128i) - sizeof(int32_t) * (imm)) \
	} \
)
#define _bsrd_v4i32(a, b, imm) ( \
	(v4i32_t) { \
		_mm_alignr_epi8((a).v1, (b).v1, sizeof(int32_t) * (imm)) \
	} \
)

/* bit shift */
#define _sal_v4i32(a, imm) ( \
	(v4i32_t) {_i_v4i32(slli)((a).v1, (imm))} \
)
#define _sar_v4i32(a, imm) ( \
	(v4i32_t) {_i_v4i32(srai)((a).v1, (imm))} \
)
#define _shl_v4i32(a, imm) ( \
	(v4i32_t) {_i_v4i32(slli)((a).v1, (imm))} \
)
#define _shr_v4i32(a, imm) ( \
	(v4i32_t) {_i_v4i32(srli)((a).v1, (imm))} \
)

/* mask */
#define _mask_v4i32(a) ( \
	(uint32_t)(_mm_movemask_epi8((a).v1)) \
)

/* transpose */
#define _lo_v4i32(a, b) ( \
	(v4i32_t) { \
		_mm_unpacklo_epi32((a).v1, (b).v1) \
	} \
)
#define _hi_v4i32(a, b) ( \
	(v4i32_t) { \
		_mm_shuffle_epi32(_mm_unpacklo_epi32((a).v1, (b).v1), 0x0e) \
	} \
)
#define _shuf_v4i32(a, i) ( \
	(v4i32_t) { \
		_mm_shuffle_epi32((a).v1, i) \
	} \
)
#define _km_v4i32(i4, i3, i2, i1)	( (i1) | ((i2)<<2) | ((i3)<<4) | ((i4)<<6) )

/* debug print */
#ifdef _LOG_H_INCLUDED
#define _print_v4i32(a) { \
	debug("(v4i32_t) %s(%d, %d, %d, %d)", #a, _ext_v4i32(a, 3), _ext_v4i32(a, 2), _ext_v4i32(a, 1), _ext_v4i32(a, 0)); \
}
#else
#define _print_v4i32(x)		;
#endif

#endif /* _V4I32_H_INCLUDED */
/**
 * end of v4i32.h
 */
