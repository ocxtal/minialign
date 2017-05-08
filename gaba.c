
/**
 * @file gaba.c
 *
 * @brief libgaba (libsea3) DP routine implementation
 *
 * @author Hajime Suzuki
 * @date 2016/1/11
 * @license Apache v2
 */

/* gap penalty model (linear or affine) */
#define LINEAR 						1
#define AFFINE						2

#ifdef MODEL
#  if !(MODEL == LINEAR || MODEL == AFFINE)
#    error "MODEL must be LINEAR (1) or AFFINE (2)."
#  endif
#else
#  define MODEL 					AFFINE
// #  define MODEL 					LINEAR
#endif


/* import unittest */
#ifndef UNITTEST_UNIQUE_ID
#  if MODEL == LINEAR
#    define UNITTEST_UNIQUE_ID		34
#  else
#    define UNITTEST_UNIQUE_ID		35
#  endif
#endif
#include  "unittest.h"


#include <stdio.h>				/* sprintf in dump_path */
#include <stdint.h>				/* uint32_t, uint64_t, ... */
#include <stddef.h>				/* offsetof */
#include <string.h>				/* memset, memcpy */
#include "gaba.h"
#include "log.h"
#include "lmm.h"
#include "sassert.h"

/* architecture dependent */
#include "arch/arch.h"

/* aliasing vector macros */
#define _VECTOR_ALIAS_PREFIX		v32i8
#include "arch/vector_alias.h"


/* add suffix */
#ifdef SUFFIX
#  if MODEL == LINEAR
#    define suffix(_base)			_base##_linear
#  else
#    define suffix(_base)			_base##_affine
#  endif
#else
#  define suffix(_base)				_base
#endif


/* constants */
#define BW_BASE						( 5 )
#define BW 							( 0x01<<BW_BASE )
#define BLK_BASE					( 5 )
#define BLK 						( 0x01<<BLK_BASE )
#define mask_t						uint32_t


#define MIN_BULK_BLOCKS				( 32 )
#define MEM_ALIGN_SIZE				( 32 )		/* 32byte aligned for AVX2 environments */
#define MEM_INIT_SIZE				( (uint64_t)256 * 1024 * 1024 )
#define MEM_MARGIN_SIZE				( 2048 )
#define PSUM_BASE					( 1 )

_static_assert(V2I32_MASK_01 == GABA_STATUS_UPDATE_A);
_static_assert(V2I32_MASK_10 == GABA_STATUS_UPDATE_B);


/**
 * @macro _likely, _unlikely
 * @brief branch prediction hint for gcc-compatible compilers
 */
#define _likely(x)		__builtin_expect(!!(x), 1)
#define _unlikely(x)	__builtin_expect(!!(x), 0)

/**
 * @macro _force_inline
 * @brief inline directive for gcc-compatible compilers
 */
#define _force_inline	inline
// #define _force_inline


/** assume 64bit little-endian system */
_static_assert(sizeof(void *) == 8);

/** check size of structs declared in gaba.h */
_static_assert(sizeof(struct gaba_params_s) == 8);
_static_assert(sizeof(struct gaba_section_s) == 16);
_static_assert(sizeof(struct gaba_fill_s) == 64);
_static_assert(sizeof(struct gaba_path_section_s) == 32);
_static_assert(sizeof(struct gaba_path_s) == 8);
_static_assert(sizeof(struct gaba_alignment_s) == 80);
_static_assert(sizeof(vec_masku_t) == 4);

/**
 * @macro _plen
 * @brief extract plen from path_section_s
 */
#define _plen(sec)		( (sec)->alen + (sec)->blen )

/* forward declarations */
static int32_t gaba_dp_add_stack(struct gaba_dp_context_s *this, uint64_t size);
static void *gaba_dp_malloc(struct gaba_dp_context_s *this, uint64_t size);
struct gaba_dp_context_s;


/**
 * @union gaba_dir_u
 */
union gaba_dir_u {
	struct gaba_dir_dynamic {
		int8_t acc;				/** (1) accumulator (v[0] - v[BW-1]) */
		int8_t _pad[3];			/** (3) */
		uint32_t array;			/** (4) dynamic band */
	} dynamic;
	struct gaba_dir_guided {
		uint8_t *ptr;			/** (8) guided band */
	} guided;
};
_static_assert(sizeof(union gaba_dir_u) == 8);

/**
 * @struct gaba_small_delta_s
 */
struct gaba_small_delta_s {
	int8_t delta[BW];			/** (32) small delta */
	int8_t max[BW];				/** (32) max */
};
_static_assert(sizeof(struct gaba_small_delta_s) == 64);

/**
 * @struct gaba_middle_delta_s
 */
struct gaba_middle_delta_s {
	int16_t delta[BW];		/** (64) middle delta */
};
_static_assert(sizeof(struct gaba_middle_delta_s) == 64);

/**
 * @struct gaba_mask_pair_u
 */
#if MODEL == LINEAR
union gaba_mask_pair_u {
	struct gaba_mask_pair_s {
		vec_masku_t h;	/** (4) horizontal mask vector */
		vec_masku_t v;	/** (4) vertical mask vector */
	} pair;
	uint64_t all;
};
_static_assert(sizeof(union gaba_mask_pair_u) == 8);
#else
union gaba_mask_pair_u {
	struct gaba_mask_pair_s {
		vec_masku_t h;	/** (4) horizontal mask vector */
		vec_masku_t v;	/** (4) vertical mask vector */
		vec_masku_t e;	/** (4) e mask vector */
		vec_masku_t f;	/** (4) f mask vector */
	} pair;
	uint64_t all;
};
_static_assert(sizeof(union gaba_mask_pair_u) == 16);
#endif

/**
 * @struct gaba_diff_vec_s
 */
struct gaba_diff_vec_s {
	uint8_t dh[BW];				/** (32) dh in the lower 5bits, de in the higher 3bits */
	uint8_t dv[BW];				/** (32) dv in the lower 5bits, df in the higher 3bits */
};
_static_assert(sizeof(struct gaba_diff_vec_s) == 64);

/**
 * @struct gaba_char_vec_s
 */
struct gaba_char_vec_s {
	uint8_t w[BW];				/** (32) a in the lower 4bit, b in the higher 4bit */
};
_static_assert(sizeof(struct gaba_char_vec_s) == 32);

/**
 * @struct gaba_block_s
 */
struct gaba_block_s {
	union gaba_mask_pair_u mask[BLK];/** (256 / 512) mask vectors */
	struct gaba_diff_vec_s diff; 		/** (64) */
	struct gaba_small_delta_s sd;		/** (64) */
	union gaba_dir_u dir;				/** (8) */
	int64_t offset;						/** (8) */
	int32_t aridx, bridx;				/** (8) reverse index in the current section */
	struct gaba_middle_delta_s const *md;/** (8) pointer to the middle delta vectors */
	struct gaba_char_vec_s ch;			/** (32) char vector */
};
struct gaba_phantom_block_s {
	struct gaba_diff_vec_s diff; 		/** (64) */
	struct gaba_small_delta_s sd;		/** (64) */
	union gaba_dir_u dir;				/** (8) */
	int64_t offset;						/** (8) */
	int32_t aridx, bridx;				/** (8) reverse index in the current section */
	struct gaba_middle_delta_s const *md;/** (8) pointer to the middle delta vectors */
	struct gaba_char_vec_s ch;			/** (32) char vector */
};
#if MODEL == LINEAR
_static_assert(sizeof(struct gaba_block_s) == 448);
#else
_static_assert(sizeof(struct gaba_block_s) == 704);
#endif
_static_assert(sizeof(struct gaba_phantom_block_s) == 192);
#define _last_block(x)				( (struct gaba_block_s *)(x) - 1 )

/**
 * @struct gaba_joint_tail_s
 *
 * @brief (internal) init vector container.
 * sizeof(struct gaba_joint_tail_s) == 64
 */
struct gaba_joint_tail_s {
	/* coordinates */
	int64_t psum;				/** (8) global p-coordinate of the tail */
	int32_t p;					/** (4) local p-coordinate of the tail */
	uint32_t ssum;				/** (4) */

	/* max scores */
	int64_t max;				/** (8) max */

	/* status */
	uint32_t stat;				/** (4) */
	uint32_t rem_len;			/** (4) */

	/* section info */
	struct gaba_joint_tail_s const *tail;/** (8) */
	uint32_t apos, bpos;		/** (8) pos */
	uint32_t alen, blen;		/** (8) len */
	uint32_t aid, bid;			/** (8) id */
};
_static_assert(sizeof(struct gaba_joint_tail_s) == 64);
_static_assert(offsetof(struct gaba_joint_tail_s, psum) == offsetof(struct gaba_fill_s, psum));
_static_assert(offsetof(struct gaba_joint_tail_s, p) == offsetof(struct gaba_fill_s, p));
_static_assert(offsetof(struct gaba_joint_tail_s, max) == offsetof(struct gaba_fill_s, max));
#define _tail(x)				( (struct gaba_joint_tail_s *)(x) )
#define _fill(x)				( (struct gaba_fill_s *)(x) )

/**
 * @struct gaba_merge_tail_s
 */
struct gaba_merge_tail_s {
	/* coordinates */
	int64_t psum;				/** (8) global p-coordinate of the tail */
	int32_t p;					/** (4) local p-coordinate of the tail */
	uint32_t ssum;				/** (4) */

	/* max scores */
	int64_t max;				/** (8) max */

	/* status */
	uint32_t stat;				/** (4) */
	uint32_t rem_len;			/** (4) */

	/* section info */
	struct gaba_joint_tail_s const *tail;/** (8) */
	uint32_t apos, bpos;		/** (8) pos */
	uint32_t alen, blen;		/** (8) len */
	uint32_t aid, bid;			/** (8) id */

	/* tail array */
	uint8_t tail_idx[2][BW];	/** (64) array of index of joint_tail */
};
_static_assert(sizeof(struct gaba_merge_tail_s) == 128);

/**
 * @struct gaba_path_intl_s
 */
struct gaba_path_intl_s {
	/* path arrays */
	uint32_t *phead, *ptail;
	uint32_t phofs, ptofs;

	/* section arrays */
	struct gaba_path_section_s *shead, *stail;

	/* gap counts */
	uint32_t gic, gec;
};
_static_assert(sizeof(struct gaba_path_intl_s) == 48);

/**
 * @struct gaba_reader_work_s
 *
 * @brief (internal) abstract sequence reader
 * sizeof(struct gaba_reader_s) == 16
 * sizeof(struct gaba_reader_work_s) == 384
 */
struct gaba_reader_work_s {
	/** 64byte alidned */
	uint8_t const *alim, *blim;			/** (16) max index of seq array */
	uint8_t const *atail, *btail;		/** (16) tail of the current section */
	int32_t alen, blen;					/** (8) lengths of the current section */
	uint32_t aid, bid;					/** (8) ids */
	uint64_t plim;						/** (8) p limit coordinate */
	uint64_t _pad;						/** (8) */
	/** 64, 64 */

	/** 64byte aligned */
	uint8_t bufa[BW + BLK];				/** (64) */
	uint8_t bufb[BW + BLK];				/** (64) */
	/** 128, 192 */
};
_static_assert(sizeof(struct gaba_reader_work_s) == 192);

/**
 * @struct gaba_writer_work_s
 */
struct gaba_writer_work_s {
	/* reserved to avoid collision with alim and blim */
	uint8_t const *alim, *blim;			/** (16) unused in writer */

	/** local context */
	struct gaba_path_intl_s path;		/** (48) */
	/** 64, 64 */

	/** 64byte aligned */
	/** block pointers */
	struct gaba_block_s const *blk;		/** (8) current block */

	/** work */
	int32_t aidx, bidx;					/** (8) indices of the current trace */
	int32_t p;							/** (4) */
	int32_t q;							/** (4) */

	/** save */
	int32_t alen, blen;					/** (8) section lengths */
	uint32_t aid, bid;					/** (8) */
	int32_t asum, bsum;					/** (8) sum length from current tail to base of each section */
	int32_t asidx, bsidx;				/** (8) base indices of the current trace */
	uint64_t _pad1;						/** (8) */
	/** 64, 128 */

	/** 64byte aligned */
	/** save (contd) */
	int64_t psum;						/** (8) */
	int64_t pspos;						/** (8) */

	/** section info */
	struct gaba_joint_tail_s const *tail;/** (8) current tail */
	struct gaba_joint_tail_s const *atail;/** (8) */
	struct gaba_joint_tail_s const *btail;/** (8) */
	uint64_t _pad4[3];					/** (24) */
	/** 64, 192 */
};
_static_assert(sizeof(struct gaba_writer_work_s) == 192);

/**
 * @struct gaba_score_vec_s
 */
struct gaba_score_vec_s {
	int8_t v1[16];
	int8_t v2[16];
	int8_t v3[16];
	int8_t v4[16];
	int8_t v5[16];
};
_static_assert(sizeof(struct gaba_score_vec_s) == 80);

/**
 * @struct gaba_mem_block_s
 */
struct gaba_mem_block_s {
	struct gaba_mem_block_s *next;
	struct gaba_mem_block_s *prev;
	uint64_t size;
	uint64_t pad;
};
_static_assert(sizeof(struct gaba_mem_block_s) == 32);

/**
 * @struct gaba_stack_s
 * @brief save stack pointer
 */
struct gaba_stack_s {
	struct gaba_mem_block_s *mem;
	uint8_t *stack_top;
	uint8_t *stack_end;
	uint64_t pad;
};
_static_assert(sizeof(struct gaba_stack_s) == 32);

/**
 * @struct gaba_dp_context_s
 *
 * @brief (internal) container for dp implementations
 * sizeof(struct gaba_dp_context_s) == 640
 */
struct gaba_dp_context_s {
	/** API function pointers */
	void *api[6];						/** (48) */
	uint8_t *stack_top;					/** (8) dynamic programming matrix */
	uint8_t *stack_end;					/** (8) the end of dp matrix */
	/** 64, 64 */

	/** individually stored on init */

	/** 64byte aligned */
	union gaba_work_s {
		struct gaba_writer_work_s l;	/** (192) */
		struct gaba_reader_work_s r;	/** (192) */
	} w;
	/** 192, 256 */

	/** 64byte aligned */
	/** loaded on init */
	struct gaba_score_vec_s scv;		/** (80) substitution matrix and gaps */
	/** 80, 336 */

	/** 64byte aligned */
	int8_t m;							/** (1) match award */
	int8_t x;							/** (1) mismatch penalty (neg.int) */
	int8_t gi;							/** (1) gap open penalty */
	int8_t ge;							/** (1) gap extension penalty */
	int8_t tx;							/** (1) xdrop threshold */
	uint8_t tf;							/** (1) ungapped alignment filter threshold */

	/** output options */
	uint8_t head_margin;				/** (1) margin at the head of gaba_res_t */
	uint8_t tail_margin;				/** (1) margin at the tail of gaba_res_t */

	/* memory management */
	struct gaba_mem_block_s *curr_mem;	/** (8) */
	struct gaba_mem_block_s mem;		/** (32) */
	/** 48, 384 */

	/** phantom vectors */
	/** 64byte aligned */
	struct gaba_phantom_block_s blk;	/** (192) */
	/** 192, 576 */

	/** 64byte aligned */
	struct gaba_joint_tail_s tail;		/** (64) */
	/** 64, 640 */
};
_static_assert(sizeof(struct gaba_dp_context_s) == 640);
#define GABA_DP_CONTEXT_LOAD_OFFSET	( offsetof(struct gaba_dp_context_s, scv) )
#define GABA_DP_CONTEXT_LOAD_SIZE	( sizeof(struct gaba_dp_context_s) - GABA_DP_CONTEXT_LOAD_OFFSET )
_static_assert(GABA_DP_CONTEXT_LOAD_OFFSET == 256);
_static_assert(GABA_DP_CONTEXT_LOAD_SIZE == 384);

/**
 * @struct gaba_context_s
 *
 * @brief (API) an algorithmic context.
 *
 * @sa gaba_init, gaba_close
 * sizeof(struct gaba_context_s) = 704
 */
struct gaba_context_s {
	/** templates */
	/** 64byte aligned */
	struct gaba_dp_context_s k;		/** (640) */
	/** 640, 640 */

	/** 64byte aligned */
	struct gaba_middle_delta_s md;	/** (64) */
	/** 64, 704 */
};
_static_assert(sizeof(struct gaba_context_s) == 704);

/**
 * @enum _STATE
 */
enum _STATE {
	CONT 	= 0,
	UPDATE  = 0x0100,
	TERM 	= 0x0200
};
_static_assert((int32_t)CONT == (int32_t)GABA_STATUS_CONT);
_static_assert((int32_t)UPDATE == (int32_t)GABA_STATUS_UPDATE);
_static_assert((int32_t)TERM == (int32_t)GABA_STATUS_TERM);


/**
 * coordinate conversion macros
 */
#define _rev(pos, len)				( (len) + (uint64_t)(len) - (uint64_t)(pos) - 1 )
#define _roundup(x, base)			( ((x) + (base) - 1) & ~((base) - 1) )

/**
 * max and min
 */
#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MAX3(x,y,z) 	( MAX2(x, MAX2(y, z)) )
#define MAX4(w,x,y,z) 	( MAX2(MAX2(w, x), MAX2(y, z)) )

#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )
#define MIN3(x,y,z) 	( MIN2(x, MIN2(y, z)) )
#define MIN4(w,x,y,z) 	( MIN2(MIN2(w, x), MIN2(y, z)) )


/**
 * aligned malloc
 */
static inline
void *gaba_aligned_malloc(
	size_t size,
	size_t align)
{
	void *ptr = NULL;
	if(posix_memalign(&ptr, align, size) != 0) {
		debug("posix_memalign failed");
		return(NULL);
	}
	debug("posix_memalign(%p)", ptr);
	return(ptr);
}
static inline
void gaba_aligned_free(
	void *ptr)
{
	free(ptr);
	return;
}


/* matrix fill functions */
/* direction macros */
#define DYNAMIC 		( 1 )
#if DYNAMIC
#define direction_prefix 		dynamic_

/* direction determiners for the dynamic band algorithms */
/**
 * @macro _dir_fetch
 */
#define _dir_fetch(_dir) { \
	(_dir).dynamic.array <<= 1; \
	(_dir).dynamic.array |= (uint32_t)((_dir).dynamic.acc < 0); \
	debug("fetched dir(%x), %s", (_dir).dynamic.array, _dir_is_down(dir) ? "go down" : "go right"); \
}
/**
 * @macro _dir_update
 * @brief update direction determiner for the next band
 */
#define _dir_update(_dir, _vector, _sign) { \
	(_dir).dynamic.acc += (_sign) * (_ext(_vector, 0) - _ext(_vector, BW-1)); \
	/*debug("acc(%d), (%d, %d)", (_dir).dynamic.acc, _ext(_vector, 0), _ext(_vector, BW-1));*/ \
}
/**
 * @macro _dir_adjust_remainder
 * @brief adjust direction array when termination is detected in the middle of the block
 */
#define _dir_adjust_remainder(_dir, _filled_count) { \
	debug("adjust remainder, array(%x), shifted array(%x)", \
		(_dir).dynamic.array, \
		(_dir).dynamic.array<<(BLK - (_filled_count))); \
	(_dir).dynamic.array <<= (BLK - (_filled_count)); \
}
/**
 * @macro _dir_test_bound, _dir_test_bound_cap
 * @brief test if the bound of the direction array is invaded
 */
#define _dir_test_bound(_dir, k, p) 		( /* nothing to do */ )
#define _dir_test_bound_cap(_dir, k, p)		( /* nothing to do */ )
/**
 * @macro _dir_is_down, _dir_is_right
 * @brief direction indicator (_dir_is_down returns true if dir == down)
 */
#define _dir_is_down(_dir)					( ((uint64_t)(_dir).dynamic.array) & 0x01 )
#define _dir_is_right(_dir)					( ~((uint64_t)(_dir).dynamic.array) & 0x01 )
/**
 * @macro _dir_load
 */
#define _dir_load(_blk, _local_idx) ({ \
	union gaba_dir_u _dir = (_blk)->dir; \
	debug("load dir idx(%d), mask(%x), shifted mask(%x)", (int32_t)_local_idx, _dir.dynamic.array, _dir.dynamic.array>>(BLK - (_local_idx) - 1)); \
	_dir.dynamic.array >>= (BLK - (_local_idx) - 1); \
	_dir; \
})
/**
 * @macro _dir_bcnt
 */
#define _dir_bcnt(_dir) ( \
	popcnt((_dir).dynamic.array) \
)
/**
 * @macro _dir_windback
 */
#define _dir_windback(_dir) { \
	(_dir).dynamic.array >>= 1; \
}

#else

#define direction_prefix 		guided_

/* direction determiners for the guided band algorithms */
#define _dir_fetch(_dir)					{ /* nothing to do */ }
#define _dir_update(_dir, _vector) { \
	(_dir).guided.ptr++; \
}
#define _dir_adjust_remainder(_dir, i)		{ /* nothing to do */ }
#define _dir_test_bound(_dir, k, p) ( \
	((k)->tdr - (_dir).guided.ptr + 2) - (p) - BLK \
)
#define _dir_test_bound_cap(_dir, k, p) ( \
	((k)->tdr - (_dir).guided.ptr + 2) - (p) \
)
#define _dir_is_down(_dir)					( *(_dir).guided.ptr != 0 )
#define _dir_is_right(_dir)					( *(_dir).guided.ptr == 0 )
#define _dir_load(_blk, _local_idx) ({ \
	union gaba_dir_u _dir = (_blk)->dir; \
	_dir.guided.ptr -= (BLK - (_local_idx)); \
	_dir; \
})
#define _dir_windback(_dir) { \
	(_dir).guided.ptr--; \
}

#endif


/**
 * @macro _match
 * @brief alias to sequence matcher macro
 */
#ifndef _match
#  define _match				_and		/* 4bit encoded */
#  define _match_v16i8			_and_v16i8
#  define _match_v32i8			_and_v32i8
#endif /* _match */

/**
 * seqreader macros
 */
#define _rd_bufa_base(k)		( (k)->w.r.bufa + BLK + BW )
#define _rd_bufb_base(k)		( (k)->w.r.bufb )
#define _rd_bufa(k, pos, len)	( _rd_bufa_base(k) - (pos) - (len) )
#define _rd_bufb(k, pos, len)	( _rd_bufb_base(k) + (pos) )
#define _lo64(v)		_ext_v2i64(v, 0)
#define _hi64(v)		_ext_v2i64(v, 1)
#define _lo32(v)		_ext_v2i32(v, 0)
#define _hi32(v)		_ext_v2i32(v, 1)

/**
 * @fn transpose_section_pair
 */
static _force_inline
struct gaba_trans_section_s {
	v2i32_t id;
	v2i32_t len;
	v2i64_t base;
} transpose_section_pair(
	v2i64_t a,
	v2i64_t b)
{
	_print_v2i64(a);
	_print_v2i64(b);

	/* transpose */
	v2i32_t id_len_a = _cast_v2i64_v2i32(a);
	v2i32_t id_len_b = _cast_v2i64_v2i32(b);

	_print_v2i32(id_len_a);
	_print_v2i32(id_len_b);

	v2i32_t id = _lo_v2i32(id_len_a, id_len_b);
	v2i32_t len = _hi_v2i32(id_len_a, id_len_b);
	v2i64_t base = _hi_v2i64(a, b);

	_print_v2i32(id);
	_print_v2i32(len);
	_print_v2i64(base);

	return((struct gaba_trans_section_s){
		.id = id,
		.len = len,
		.base = base
	});
}

/*
 * @fn fill_load_section
 */
static _force_inline
void fill_load_section(
	struct gaba_dp_context_s *this,
	struct gaba_section_s const *a,
	struct gaba_section_s const *b,
	uint64_t plim)
{
	/* load current section lengths */
	struct gaba_trans_section_s c = transpose_section_pair(
		_loadu_v2i64(a), _loadu_v2i64(b));

	/* convert current section to (pos of the end, len) */
	v2i64_t c_tail = _add_v2i64(c.base, _cvt_v2i32_v2i64(c.len));

	_print_v2i32(c.id);
	_print_v2i32(c.len);
	_print_v2i64(c_tail);

	/* store sections */
	_store_v2i64(&this->w.r.atail, c_tail);
	_store_v2i32(&this->w.r.alen, c.len);
	_store_v2i32(&this->w.r.aid, c.id);
	this->w.r.plim = plim;

	return;
}

/**
 * @struct gaba_joint_block_s
 * @brief result container for block fill functions
 */
struct gaba_joint_block_s {
	struct gaba_block_s *blk;
	int64_t p;
	int32_t stat;
};

/**
 * @fn fill_load_seq_a
 */
static _force_inline
void fill_load_seq_a(
	struct gaba_dp_context_s *this,
	uint8_t const *pos,
	uint64_t len)
{
	if(pos < this->w.r.alim) {
		debug("reverse fetch a: pos(%p), len(%llu)", pos, len);
		/* reverse fetch: 2 * alen - (2 * alen - pos) + (len - 32) */
		vec_t a = _loadu(pos + (len - BW));
		_storeu(_rd_bufa(this, BW, len), _swap(a));
	} else {
		debug("forward fetch a: pos(%p), len(%llu)", pos, len);
		/* take complement */
		static uint8_t const comp[16] __attribute__(( aligned(16) )) = {
			0x00, 0x08, 0x04, 0x0c, 0x02, 0x0a, 0x06, 0x0e,
			0x01, 0x09, 0x05, 0x0d, 0x03, 0x0b, 0x07, 0x0f
		};
		vec_t const cv = _from_v16i8(_load_v16i8(comp));

		/* forward fetch: 2 * alen - pos */
		vec_t a = _loadu(_rev(pos, this->w.r.alim) - (len - 1));
		_storeu(_rd_bufa(this, BW, len), _shuf(cv, a));
	}
	return;
}

/**
 * @fn fill_load_seq_b
 */
static _force_inline
void fill_load_seq_b(
	struct gaba_dp_context_s *this,
	uint8_t const *pos,
	uint64_t len)
{
	if(pos < this->w.r.blim) {
		debug("forward fetch b: pos(%p), len(%llu)", pos, len);
		/* forward fetch: pos */
		vec_t b = _loadu(pos);
		_storeu(_rd_bufb(this, BW, len), b);
	} else {
		debug("reverse fetch b: pos(%p), len(%llu)", pos, len);
		/* take complement */
		static uint8_t const comp[16] __attribute__(( aligned(16) )) = {
			0x00, 0x08, 0x04, 0x0c, 0x02, 0x0a, 0x06, 0x0e,
			0x01, 0x09, 0x05, 0x0d, 0x03, 0x0b, 0x07, 0x0f
		};
		vec_t const cv = _from_v16i8(_load_v16i8(comp));

		/* reverse fetch: 2 * blen - pos + (len - 32) */
		vec_t b = _loadu(_rev(pos, this->w.r.blim) - (BW - 1));
		_storeu(_rd_bufb(this, BW, len), _shuf(cv, _swap(b)));
	}
	return;
}

/**
 * @fn fill_bulk_fetch
 *
 * @brief fetch 32bases from current section
 */
static _force_inline
void fill_bulk_fetch(
	struct gaba_dp_context_s *this,
	struct gaba_block_s *blk)
{
	/* load sequence from previous block */
	vec_t const mask = _set(0x0f);
	vec_t w = _load(&(blk - 1)->ch.w);
	vec_t a = _and(mask, w);
	vec_t b = _and(mask, _shr(w, 4));

	_print(w);
	_print(a);
	_print(b);

	debug("atail(%p), aridx(%d)", this->w.r.atail, (blk-1)->aridx);
	debug("btail(%p), bridx(%d)", this->w.r.btail, (blk-1)->bridx);

	/* fetch seq a */
	fill_load_seq_a(this, this->w.r.atail - (blk - 1)->aridx, BLK);
	_store(_rd_bufa(this, 0, BW), a);

	/* fetch seq b */
	_store(_rd_bufb(this, 0, BW), b);
	fill_load_seq_b(this, this->w.r.btail - (blk - 1)->bridx, BLK);
	return;
}

/**
 * @fn fill_cap_fetch
 *
 * @return the length of the sequence fetched.
 */
static _force_inline
void fill_cap_fetch(
	struct gaba_dp_context_s *this,
	struct gaba_block_s *blk)
{
	/* const */
	v2i32_t const z = _zero_v2i32();
	v2i32_t const tot = _set_v2i32(BLK);

	/* load lengths */
	v2i32_t ridx = _load_v2i32(&(blk - 1)->aridx);

	/* clip len with max load length */
	v2i32_t len = _max_v2i32(_min_v2i32(ridx, tot), z);

	/* load sequence from previous block */
	vec_t const mask = _set(0x0f);
	vec_t w = _load(&(blk - 1)->ch.w);
	vec_t a = _and(mask, w);
	vec_t b = _and(mask, _shr(w, 4));

	_print(w);
	_print(a);
	_print(b);

	/* fetch seq a */
	fill_load_seq_a(this, this->w.r.atail - _lo32(ridx), _lo32(len));
	_store(_rd_bufa(this, 0, BW), a);

	/* fetch seq b */
	_store(_rd_bufb(this, 0, BW), b);
	fill_load_seq_b(this, this->w.r.btail - _hi32(ridx), _hi32(len));
	return;
}

/**
 * @fn fill_init_fetch
 * @brief similar to cap fetch, updating ridx and rem
 */
static _force_inline
struct gaba_joint_block_s fill_init_fetch(
	struct gaba_dp_context_s *this,
	struct gaba_joint_tail_s const *prev_tail,
	struct gaba_phantom_block_s *blk,
	v2i32_t ridx)
{
	/* restore (brem, arem) */
	int64_t prem = -prev_tail->psum;
	int64_t arem = (prem + 1) / 2;
	int64_t brem = prem / 2;
	v2i32_t rem = _seta_v2i32(brem, arem);

	/* calc the next rem */
	v2i32_t const z = _zero_v2i32();
	v2i32_t const ofs = _seta_v2i32(1, 0);
	v2i32_t nrem = _max_v2i32(_sub_v2i32(rem, ridx), z);
	v2i32_t rrem = _sub_v2i32(_swap_v2i32(nrem), ofs);
	nrem = _max_v2i32(nrem, rrem);

	/* cap fetch */
	v2i32_t len = _sub_v2i32(rem, nrem);

	/* load sequence from the previous block */ {
		struct gaba_block_s *prev_blk = _last_block(prev_tail);
		vec_t const mask = _set(0x0f);
		vec_t w = _load(&prev_blk->ch.w);
		vec_t a = _and(mask, w);
		vec_t b = _and(mask, _shr(w, 4));

		debug("prev_blk(%p), prev_tail(%p)", prev_blk, prev_tail);
		_print(w);
		_print(a);
		_print(b);

		/* fetch seq a */
		fill_load_seq_a(this, this->w.r.atail - _lo32(ridx), _lo32(len));
		_store(_rd_bufa(this, 0, BW), a);

		/* fetch seq b */
		_store(_rd_bufb(this, 0, BW), b);
		fill_load_seq_b(this, this->w.r.btail - _hi32(ridx), _hi32(len));
	}

	/* store char vector to the current block */ {
		vec_t a = _loadu(_rd_bufa(this, _lo32(len), BW));
		vec_t b = _loadu(_rd_bufb(this, _hi32(len), BW));
		_store(&blk->ch.w, _or(a, _shl(b, 4)));

		_print(a);
		_print(b);
	}

	/* adjust and store ridx */
	_print_v2i32(ridx);
	_print_v2i32(len);
	ridx = _sub_v2i32(ridx, len);
	_store_v2i32(&blk->aridx, ridx);
	_print_v2i32(ridx);
	debug("blk(%p), p(%d), stat(%x)", blk + 1, _lo32(len) + _hi32(len), _mask_v2i32(_eq_v2i32(ridx, z)));

	return((struct gaba_joint_block_s){
		.blk = (struct gaba_block_s *)(blk + 1),
		.p = _lo32(len) + _hi32(len),
		.stat = (_mask_v2i32(_eq_v2i32(ridx, z)) == V2I32_MASK_00) ? CONT : UPDATE
	});
}

/**
 * @fn fill_restore_fetch
 * @brief fetch sequence from existing block
 */
static _force_inline
void fill_restore_fetch(
	struct gaba_dp_context_s *this,
	struct gaba_block_s const *blk)
{
	vec_t const mask = _set(0x0f);

	/* calc cnt */
	v2i32_t curr_len = _load_v2i32(&blk->aridx);
	v2i32_t prev_len = _load_v2i32(&(blk - 1)->aridx);
	v2i32_t cnt = _sub_v2i32(prev_len, curr_len);

	/* from the current block */
	vec_t cw = _load(&blk->ch.w);
	vec_t ca = _and(mask, cw);
	vec_t cb = _and(mask, _shr(cw, 4));
	_storeu(_rd_bufa(this, _lo32(cnt), BW), ca);
	_storeu(_rd_bufb(this, _hi32(cnt), BW), cb);

	/* from the previous block */
	vec_t pw = _load(&(blk - 1)->ch.w);
	vec_t pa = _and(mask, pw);
	vec_t pb = _and(mask, _shr(pw, 4));
	_store(_rd_bufa(this, 0, BW), pa);
	_store(_rd_bufb(this, 0, BW), pb);

	_print(pa);
	_print(pb);

	return;
}

/**
 * @fn fill_update_section
 */
static _force_inline
v2i32_t fill_update_section(
	struct gaba_dp_context_s *this,
	struct gaba_block_s *blk,
	v2i32_t cnt)
{
	/* update pos and ridx */
	v2i32_t ridx = _load_v2i32(&(blk - 1)->aridx);
	ridx = _sub_v2i32(ridx, cnt);
	_print_v2i32(ridx);

	/* store ridx to block */
	_store_v2i32(&blk->aridx, ridx);

	vec_t a = _loadu(_rd_bufa(this, _lo32(cnt), BW));
	vec_t b = _loadu(_rd_bufb(this, _hi32(cnt), BW));
	_store(&blk->ch.w, _or(a, _shl(b, 4)));

	_print(a);
	_print(b);
	return(ridx);
}

/**
 * @fn fill_gapless_filter
 */
static _force_inline
int32_t fill_gapless_filter(
	struct gaba_dp_context_s *this,
	struct gaba_block_s *blk,
	int32_t stat)
{
	v16i8_t const load_mask = _set_v16i8(0x0f);
	v16i8_t const match_mask = _bsl_v16i8(_set_v16i8(0xff), 1);

	/* load char vectors */
	v16i8_t a = _load_v16i8(&blk->ch.w[0]);
	v16i8_t b = _load_v16i8(&blk->ch.w[16]);
	v16i8_t a0 = _swap_v16i8(_and_v16i8(load_mask, a));
	v16i8_t b0 = _bsr_v16i8(_and_v16i8(load_mask, _shr_v16i8(b, 4)), 1);

	/* make shifted vectors */
	v16i8_t a1 = _bsr_v16i8(a0, 1); //, a2 = _bsr_v16i8(a0, 2);
	v16i8_t b1 = _bsr_v16i8(b0, 1); //, b2 = _bsr_v16i8(b0, 2);

	_print_v16i8(a0);
	_print_v16i8(a1);
	// _print_v16i8(a2);
	_print_v16i8(b0);
	_print_v16i8(b1);
	// _print_v16i8(b2);

	/* ungapped alignment */
	// v16i8_t m0 = _shuf_v16i8(match_mask, _match_v16i8(a0, b2));
	v16i8_t m1 = _shuf_v16i8(match_mask, _match_v16i8(a0, b1));
	v16i8_t m2 = _shuf_v16i8(match_mask, _match_v16i8(a0, b0));
	v16i8_t m3 = _shuf_v16i8(match_mask, _match_v16i8(a1, b0));
	// v16i8_t m4 = _shuf_v16i8(match_mask, _match_v16i8(a2, b0));

	// _print_v16i8(m0);
	_print_v16i8(m1);
	_print_v16i8(m2);
	_print_v16i8(m3);
	// _print_v16i8(m4);

	// v16i8_t cnt_mask = _or_v16i8(_or_v16i8(_or_v16i8(m0, m1), m2), _or_v16i8(m3, m4));
	v16i8_t cnt_mask = _or_v16i8(_or_v16i8(m1, m2), m3);
	int64_t cnt = popcnt(((v16i8_masku_t){ .mask = _mask_v16i8(cnt_mask) }).all);

	_print_v16i8(cnt_mask);
	debug("cnt(%lld), tf(%u)", cnt, this->tf);

	return((cnt > this->tf) ? CONT : TERM);
}

/**
 * @fn fill_create_phantom_block
 * @brief create joint_head on the stack to start block extension
 */
static _force_inline
struct gaba_joint_block_s fill_create_phantom_block(
	struct gaba_dp_context_s *this,
	struct gaba_joint_tail_s const *prev_tail)
{
	debug("create head prev_tail(%p), p(%d), psum(%lld), ssum(%d)",
		prev_tail, prev_tail->p, prev_tail->psum, prev_tail->ssum);

	/* init working stack */
	struct gaba_phantom_block_s *blk = (struct gaba_phantom_block_s *)this->stack_top;
	struct gaba_block_s *pblk = _last_block(prev_tail);
	debug("start stack_top(%p), stack_end(%p)", this->stack_top, this->stack_end);

	/* copy phantom vectors from the previous fragment */
	// _memcpy_blk_aa(&blk->diff, &pblk->diff, offsetof(struct gaba_phantom_block_s, sd.max));
	_memcpy_blk_aa(&blk->diff, &pblk->diff, offsetof(struct gaba_phantom_block_s, dir));

	/* fill max vector with zero */
	// _store(&blk->sd.max, _zero());
	// _store(&blk->sd.max, _set(-128));

	/* copy remaining */
	blk->dir = pblk->dir;
	blk->offset = pblk->offset;

	/* calc ridx */
	v2i32_t ridx = _sub_v2i32(
		_load_v2i32(&this->w.r.alen),
		_load_v2i32(&prev_tail->apos));
	_print_v2i32(ridx);

	/* check if init fetch is needed */
	if(prev_tail->psum >= 0) {
		/* store index on the current section */
		_store_v2i32(&blk->aridx, ridx);
		
		/* copy char vectors from prev_tail */
		_store(&blk->ch, _load(&pblk->ch));

		return((struct gaba_joint_block_s){
			.blk = (struct gaba_block_s *)(blk + 1),
			.p = 0,
			.stat = CONT
		});
	} else {
		/* init fetch */
		// return(fill_init_fetch(this, prev_tail, blk, ridx));
		struct gaba_joint_block_s stat = fill_init_fetch(this, prev_tail, blk, ridx);

		/* check if initial vector is filled */
		if(prev_tail->psum + stat.p >= 0) {
			stat.stat = fill_gapless_filter(this, stat.blk - 1, stat.stat);
		}
		return(stat);
	}
}

/**
 * @fn fill_create_tail
 * @brief create joint_tail at the end of the blocks
 */
static _force_inline
struct gaba_joint_tail_s *fill_create_tail(
	struct gaba_dp_context_s *this,
	struct gaba_joint_tail_s const *prev_tail,
	struct gaba_block_s *blk,
	int64_t p,
	int32_t stat)
{
	debug("create tail, p(%lld)", p);

	/* create joint_tail */
	struct gaba_joint_tail_s *tail = (struct gaba_joint_tail_s *)blk;
	this->stack_top = (void *)(tail + 1);	/* write back stack_top */
	(blk - 1)->md = _last_block(prev_tail)->md;	/* copy middle delta pointer */
	debug("end stack_top(%p), stack_end(%p), blocks(%lld)", this->stack_top, this->stack_end, (p + BLK - 1) / BLK);

	/* update p */
	int64_t prev_psum = prev_tail->psum;
	int32_t np = (prev_psum < 0) ? MAX2(p + prev_psum, 0) : p;

	/* save misc to joint_tail */
	tail->psum = p + prev_tail->psum;
	tail->p = np;
	debug("p(%lld), prev_tail->psum(%lld), prev_tail->p(%d), tail->psum(%lld), tail->p(%d)",
		p, prev_tail->psum, prev_tail->p, tail->psum, tail->p);
	tail->ssum = prev_tail->ssum + 1;
	tail->tail = prev_tail;					/* to treat tail chain as linked list */
	tail->rem_len = 0;

	/* search max section */
	v32i16_t sd = _cvt_v32i8_v32i16(_load(&(blk - 1)->sd.max));
	v32i16_t md = _load_v32i16(_last_block(prev_tail)->md);
	_print_v32i16(sd);
	_print_v32i16(md);

	/* extract max */
	md = _add_v32i16(md, sd);
	int16_t max = _hmax_v32i16(md);
	_print_v32i16(md);

	/* store */
	// tail->mask_max.mask = mask_max;
	tail->max = max + (blk - 1)->offset;

	debug("offset(%lld)", (blk - 1)->offset);
	debug("max(%d)", _hmax_v32i16(md));
	// debug("mask_max(%u)", tail->mask_max.all);

	/* store section lengths */
	v2i32_t const z = _zero_v2i32();
	v2i32_t ridx = _load_v2i32(&(blk - 1)->aridx);
	v2i32_t len = _load_v2i32(&this->w.r.alen);
	_print_v2i32(ridx);
	_print_v2i32(len);
	_print_v2i32(_sel_v2i32(
		_eq_v2i32(ridx, z),
		z, _sub_v2i32(len, ridx)));
	_store_v2i32(&tail->apos, _sel_v2i32(_eq_v2i32(ridx, z),
		z, _sub_v2i32(len, ridx)));
	_store_v2i32(&tail->alen, len);

	/* store section ids */
	v2i32_t id = _load_v2i32(&this->w.r.aid);
	_store_v2i32(&tail->aid, id);

	/* store status */
	tail->stat = stat | _mask_v2i32(_eq_v2i32(ridx, z));
	return(tail);
}

/**
 * @macro _fill_load_context
 * @brief load vectors onto registers
 */
#if MODEL == LINEAR
#define _fill_load_context(_blk) \
	debug("blk(%p)", (_blk)); \
	/* load sequence buffer offset */ \
	uint8_t const *aptr = _rd_bufa(this, 0, BW); \
	uint8_t const *bptr = _rd_bufb(this, 0, BW); \
	/* load mask pointer */ \
	union gaba_mask_pair_u *ptr = (_blk)->mask; \
	/* load vector registers */ \
	register vec_t dh = _load(((_blk) - 1)->diff.dh); \
	register vec_t dv = _load(((_blk) - 1)->diff.dv); \
	_print(_add(dh, _load_ofsh(this->scv))); \
	_print(_add(dv, _load_ofsv(this->scv))); \
	/* load delta vectors */ \
	register vec_t delta = _load(((_blk) - 1)->sd.delta); \
	register vec_t max = _load(((_blk) - 1)->sd.max); \
	_print(max); \
	_print_v32i16(_add_v32i16(_cvt_v32i8_v32i16(delta), _load_v32i16(_last_block(&this->tail)->md))); \
	/* load direction determiner */ \
	union gaba_dir_u dir = ((_blk) - 1)->dir; \
	/* load large offset */ \
	int64_t offset = ((_blk) - 1)->offset; \
	debug("offset(%lld)", offset);
#else
#define _fill_load_context(_blk) \
	debug("blk(%p)", (_blk)); \
	/* load sequence buffer offset */ \
	uint8_t const *aptr = _rd_bufa(this, 0, BW); \
	uint8_t const *bptr = _rd_bufb(this, 0, BW); \
	/* load mask pointer */ \
	union gaba_mask_pair_u *ptr = (_blk)->mask; \
	/* load vector registers */ \
	vec_t const mask = _set(0x07); \
	register vec_t dh = _load(((_blk) - 1)->diff.dh); \
	register vec_t dv = _load(((_blk) - 1)->diff.dv); \
	register vec_t de = _and(mask, dh); \
	register vec_t df = _and(mask, dv); \
	dh = _shr(_andn(mask, dh), 3); \
	dv = _shr(_andn(mask, dv), 3); \
	de = _add(dv, de); \
	df = _add(dh, df); \
	dh = _sub(_zero(), dh); \
	_print(_add(dh, _load_ofsh(this->scv))); \
	_print(_add(dv, _load_ofsv(this->scv))); \
	_print(_sub(_sub(de, dv), _load_adjh(this->scv))); \
	_print(_sub(_add(df, dh), _load_adjv(this->scv))); \
	/* load delta vectors */ \
	register vec_t delta = _load(((_blk) - 1)->sd.delta); \
	register vec_t max = _load(((_blk) - 1)->sd.max); \
	_print(max); \
	_print_v32i16(_add_v32i16(_cvt_v32i8_v32i16(delta), _load_v32i16(_last_block(&this->tail)->md))); \
	/* load direction determiner */ \
	union gaba_dir_u dir = ((_blk) - 1)->dir; \
	/* load large offset */ \
	int64_t offset = ((_blk) - 1)->offset; \
	debug("offset(%lld)", offset);
#endif

/**
 * @macro _fill_body
 * @brief update vectors
 */
#if MODEL == LINEAR
#define _fill_body() { \
	register vec_t t = _match(_loadu(aptr), _loadu(bptr)); \
	/*t = _shuf(_load_sc(this, sb), t);*/ \
	t = _shuf(_load_sb(this->scv), t); \
	_print(t); \
	t = _max(dh, t); \
	t = _max(dv, t); \
	ptr->pair.h.mask = _mask(_eq(t, dv)); \
	ptr->pair.v.mask = _mask(_eq(t, dh)); \
	debug("mask(%llx)", ptr->all); \
	ptr++; \
	vec_t _dv = _sub(t, dh); \
	dh = _sub(t, dv); \
	dv = _dv; \
	_print(_add(dh, _load_ofsh(this->scv))); \
	_print(_add(dv, _load_ofsv(this->scv))); \
}
#else /* MODEL == AFFINE */
#define _fill_body() { \
	register vec_t t = _match(_loadu(aptr), _loadu(bptr)); \
	_print(_loadu(aptr)); \
	_print(_loadu(bptr)); \
	t = _shuf(_load_sb(this->scv), t); \
	_print(t); \
	t = _max(de, t); \
	t = _max(df, t); \
	ptr->pair.h.mask = _mask(_eq(t, de)); \
	ptr->pair.v.mask = _mask(_eq(t, df)); \
	debug("mask(%llx)", ptr->all); \
	/* update de and dh */ \
	de = _add(de, _load_adjh(this->scv)); \
	vec_t te = _max(de, t); \
	ptr->pair.e.mask = _mask(_eq(te, de)); \
	de = _add(te, dh); \
	dh = _add(dh, t); \
	/* update df and dv */ \
	df = _add(df, _load_adjv(this->scv)); \
	vec_t tf = _max(df, t); \
	ptr->pair.f.mask = _mask(_eq(tf, df)); \
	df = _sub(tf, dv); \
	t = _sub(dv, t); \
	ptr++; \
	dv = dh; dh = t; \
	_print(_add(dh, _load_ofsh(this->scv))); \
	_print(_add(dv, _load_ofsv(this->scv))); \
	_print(_sub(_sub(de, dv), _load_adjh(this->scv))); \
	_print(_sub(_add(df, dh), _load_adjv(this->scv))); \
}
#endif /* MODEL */

/**
 * @macro _fill_update_delta
 * @brief update small delta vector and max vector
 */
#define _fill_update_delta(_op, _vector, _offset, _sign) { \
	delta = _op(delta, _add(_vector, _offset)); \
	max = _max(max, delta); \
	_dir_update(dir, _vector, _sign); \
	_print_v32i16(_add_v32i16(_set_v32i16(offset), _add_v32i16(_cvt_v32i8_v32i16(delta), _load_v32i16(_last_block(&this->tail)->md)))); \
	_print_v32i16(_add_v32i16(_set_v32i16(offset), _add_v32i16(_cvt_v32i8_v32i16(max), _load_v32i16(_last_block(&this->tail)->md)))); \
}

/**
 * @macro _fill_right, _fill_down
 * @brief wrapper of _fill_body and _fill_update_delta
 */
#define _fill_right_update_ptr() { \
	aptr--;				/* increment sequence buffer pointer */ \
}
#define _fill_right_windback_ptr() { \
	aptr++; \
}
#if MODEL == LINEAR
#define _fill_right() { \
	dh = _bsl(dh, 1);	/* shift left dh */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(_add, dh, _load_ofsh(this->scv), 1); \
}
#else
#define _fill_right() { \
	dh = _bsl(dh, 1);	/* shift left dh */ \
	df = _bsl(df, 1);	/* shift left df */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(_sub, dh, _load_ofsh(this->scv), -1); \
}
#endif /* MODEL */
#define _fill_down_update_ptr() { \
	bptr++;				/* increment sequence buffer pointer */ \
}
#define _fill_down_windback_ptr() { \
	bptr--; \
}
#if MODEL == LINEAR
#define _fill_down() { \
	dv = _bsr(dv, 1);	/* shift right dv */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(_add, dv, _load_ofsv(this->scv), 1); \
}
#else
#define _fill_down() { \
	dv = _bsr(dv, 1);	/* shift right dv */ \
	de = _bsr(de, 1);	/* shift right de */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(_add, dv, _load_ofsv(this->scv), 1); \
}
#endif /* MODEL */

/**
 * @macro _fill_update_offset
 * @brief update offset and max vector, reset the small delta
 */
#define _fill_update_offset() { \
	int8_t _cd = _ext(delta, BW/2); \
	offset += _cd; \
	delta = _sub(delta, _set(_cd)); \
	max = _sub(max, _set(_cd)); \
}

/**
 * @macro _fill_store_vectors
 * @brief store vectors at the end of the block
 */
#if MODEL == LINEAR
#define _fill_store_vectors(_blk) ({ \
	/* store diff vectors */ \
	_store((_blk)->diff.dh, dh); \
	_store((_blk)->diff.dv, dv); \
	_print(dh); \
	_print(dv); \
	/* store delta vectors */ \
	_store((_blk)->sd.delta, delta); \
	_store((_blk)->sd.max, max); \
	/* store direction array */ \
	(_blk)->dir = dir; \
	/* store large offset */ \
	(_blk)->offset = offset; \
	/* calc cnt */ \
	uint64_t acnt = _rd_bufa(this, 0, BW) - aptr; \
	uint64_t bcnt = bptr - _rd_bufb(this, 0, BW); \
	_seta_v2i32(bcnt, acnt); \
})
#else
#define _fill_store_vectors(_blk) ({ \
	/* store diff vectors */ \
	de = _sub(de, dv); \
	df = _add(df, dh); \
	dh = _sub(_zero(), dh); \
	_print(dh); \
	_print(dv); \
	_print(de); \
	_print(df); \
	dh = _shl(dh, 3); \
	dv = _shl(dv, 3); \
	_print(dh); \
	_print(dv); \
	_store((_blk)->diff.dh, _add(dh, de)); \
	_store((_blk)->diff.dv, _add(dv, df)); \
	_print(_add(dh, de)); \
	_print(_add(dv, df)); \
	/* store delta vectors */ \
	_store((_blk)->sd.delta, delta); \
	_store((_blk)->sd.max, max); \
	/* store direction array */ \
	(_blk)->dir = dir; \
	/* store large offset */ \
	(_blk)->offset = offset; \
	/* calc cnt */ \
	uint64_t acnt = _rd_bufa(this, 0, BW) - aptr; \
	uint64_t bcnt = bptr - _rd_bufb(this, 0, BW); \
	_seta_v2i32(bcnt, acnt); \
})
#endif

/**
 * @fn fill_test_xdrop
 * @brief returns negative if terminate-condition detected
 */
static _force_inline
int64_t fill_test_xdrop(
	struct gaba_dp_context_s const *this,
	struct gaba_block_s const *blk)
{
	return(this->tx - blk->sd.max[BW/2]);
}

/**
 * @fn fill_bulk_test_seq_bound
 * @brief returns negative if ij-bound (for the bulk fill) is invaded
 */
static _force_inline
int64_t fill_bulk_test_seq_bound(
	struct gaba_dp_context_s const *this,
	struct gaba_block_s const *blk)
{
	debug("test(%lld, %lld), len(%d, %d)",
		(int64_t)blk->aridx - BW,
		(int64_t)blk->bridx - BW,
		blk->aridx, blk->bridx);
	return(((int64_t)blk->aridx - BW)
		 | ((int64_t)blk->bridx - BW));
}

/**
 * @fn fill_cap_test_seq_bound
 * @brief returns negative if ij-bound (for the cap fill) is invaded
 */
#define _fill_cap_test_seq_bound_init(_blk) \
	uint8_t *alim = _rd_bufa(this, ((_blk) - 1)->aridx, BW); \
	uint8_t *blim = _rd_bufb(this, ((_blk) - 1)->bridx, BW);
#define _fill_cap_test_seq_bound() ( \
	((int64_t)aptr - (int64_t)alim) | ((int64_t)blim - (int64_t)bptr) \
)

/**
 * @fn fill_bulk_block
 * @brief fill a block
 */
static _force_inline
void fill_bulk_block(
	struct gaba_dp_context_s *this,
	struct gaba_block_s *blk)
{
	/* fetch sequence */
	fill_bulk_fetch(this, blk);

	/* load vectors onto registers */
	debug("blk(%p)", blk);
	_fill_load_context(blk);
	/**
	 * @macro _fill_block
	 * @brief an element of unrolled fill-in loop
	 */
	#define _fill_block(_direction, _label, _jump_to) { \
		_dir_fetch(dir); \
		if(_unlikely(!_dir_is_##_direction(dir))) { \
			goto _fill_##_jump_to; \
		} \
		_fill_##_label: \
		_fill_##_direction##_update_ptr(); \
		_fill_##_direction(); \
		if(--i == 0) { break; } \
	}

	/* update diff vectors */
	int64_t i = BLK;
	while(1) {					/* 4x unrolled loop */
		_fill_block(down, d1, r1);
		_fill_block(right, r1, d2);
		_fill_block(down, d2, r2);
		_fill_block(right, r2, d1);
	}

	/* update seq offset */
	_fill_update_offset();

	/* store vectors */
	v2i32_t cnt = _fill_store_vectors(blk);

	/* update section */
	fill_update_section(this, blk, cnt);

	return;
}

/**
 * @fn fill_bulk_predetd_blocks
 * @brief fill <blk_cnt> contiguous blocks without ij-bound test
 */
static _force_inline
struct gaba_joint_block_s fill_bulk_predetd_blocks(
	struct gaba_dp_context_s *this,
	struct gaba_block_s *blk,
	uint64_t blk_cnt)
{
	int32_t stat = CONT;
	uint64_t bc = 0;
	for(bc = 0; bc < blk_cnt; bc++) {
		/* check xdrop termination */
		if(fill_test_xdrop(this, blk - 1) < 0) {
			stat = TERM; break;
		}

		/* bulk fill */
		debug("blk(%p)", blk);
		fill_bulk_block(this, blk++);
	}
	return((struct gaba_joint_block_s){
		.blk = blk,
		.p = (int64_t)bc * BLK,
		.stat = stat
	});
}

/**
 * @fn fill_bulk_seq_bounded
 * @brief fill blocks with ij-bound test
 */
static _force_inline
struct gaba_joint_block_s fill_bulk_seq_bounded(
	struct gaba_dp_context_s *this,
	struct gaba_block_s *blk)
{
	int32_t stat = CONT;

	/* init local coordinate */
	int64_t p = 0;

	/* bulk fill loop */
	while(1) {
		/* check termination */
		if((fill_test_xdrop(this, blk - 1)
		  | fill_bulk_test_seq_bound(this, blk - 1)) < 0) {
			break;
		}
		
		/* bulk fill */
		debug("blk(%p)", blk);
		fill_bulk_block(this, blk++);
		
		/* update p-coordinate */
		p += BLK;
	}
	if(fill_test_xdrop(this, blk - 1) < 0) { stat = TERM; }
	return((struct gaba_joint_block_s){
		.blk = blk,
		.p = p,
		.stat = stat 		/* < 0 if stat == TERM */
	});
}

/**
 * @fn fill_cap_seq_bounded
 * @brief fill blocks with cap test
 */
static _force_inline
struct gaba_joint_block_s fill_cap_seq_bounded(
	struct gaba_dp_context_s *this,
	struct gaba_block_s *blk)
{
	int32_t stat = CONT;
	int64_t p = 0;

	while(1) {
		/* check xdrop termination */
		if(fill_test_xdrop(this, blk - 1) < 0) {
			stat = TERM; goto _fill_cap_seq_bounded_finish;
		}
		/* fetch sequence */
		fill_cap_fetch(this, blk);

		/* vectors on registers inside this block */ {
			_fill_cap_test_seq_bound_init(blk);
			_fill_load_context(blk);

			/* update diff vectors */
			uint64_t i = 0;
			for(i = 0; i < BLK; i++) {
				/* determine direction */
				_dir_fetch(dir);

				if(_dir_is_right(dir)) {
					/* update sequence coordinate and then check term */
					_fill_right_update_ptr();
					
					if(_fill_cap_test_seq_bound() < 0) {
						_fill_right_windback_ptr();
						_dir_windback(dir);
						break;
					}

					/* update band */
					_fill_right();
				} else {
					/* update sequence coordinate and then check term */
					_fill_down_update_ptr();
					
					if(_fill_cap_test_seq_bound() < 0) {
						_fill_down_windback_ptr();
						_dir_windback(dir);
						break;
					}

					/* update band */
					_fill_down();
				}
			}
			/* adjust dir remainder */
			_dir_adjust_remainder(dir, i);

			/* update seq offset */
			_fill_update_offset();
			
			/* store mask and vectors */
			v2i32_t cnt = _fill_store_vectors(blk);

			/* update section */
			fill_update_section(this, blk, cnt);

			/* update block pointer and p-coordinate */
			blk += (i != 0); p += i;

			/* break if not filled full length */
			if(i != BLK) { stat = UPDATE; break; }
		}
	}

	debug("blk(%p), p(%lld), stat(%x)", blk, p, stat);

_fill_cap_seq_bounded_finish:;
	return((struct gaba_joint_block_s){
		.blk = blk,
		.p = p,
		.stat = stat
	});
}

/**
 * @fn calc_max_bulk_blocks_mem
 * @brief calculate maximum number of blocks (limited by stack size)
 */
static _force_inline
uint64_t calc_max_bulk_blocks_mem(
	struct gaba_dp_context_s const *this)
{
	uint64_t mem_size = this->stack_end - this->stack_top;
	uint64_t blk_cnt = mem_size / sizeof(struct gaba_block_s);
	debug("calc_max_block_mem, stack_top(%p), stack_end(%p), mem_size(%llu), cnt(%llu)",
		this->stack_top, this->stack_end, mem_size, (blk_cnt > 3) ? (blk_cnt - 3) : 0);
	return(((blk_cnt > 3) ? blk_cnt : 3) - 3);
}

/**
 * @fn calc_min_expected_blocks_blk
 * @brief calc min #expected blocks from block,
 * used to determine #blocks which can be filled without seq boundary check.
 */
static _force_inline
uint64_t calc_min_expected_blocks_blk(
	struct gaba_dp_context_s const *this,
	struct gaba_block_s const *blk)
{
	uint64_t p = MIN2(
		(blk - 1)->aridx,
		(blk - 1)->bridx);
	return((p + p/2) / BLK);
}

/**
 * @fn calc_max_expected_blocks_tail
 * @brief calc max #expected blocks from block,
 * used to determine #blocks which can be filled without mem boundary check.
 */
static _force_inline
uint64_t calc_max_expected_blocks_tail(
	struct gaba_dp_context_s const *this,
	struct gaba_joint_tail_s const *tail)
{
	uint64_t p = MIN2(
		this->w.r.alen - tail->apos,
		this->w.r.blen - tail->bpos);
	return((2*p + p/2) / BLK);
}

/**
 * @fn calc_min_expected_blocks_tail
 * @brief calc min #expected blocks from block,
 * used to determine #blocks which can be filled without seq boundary check.
 */
static _force_inline
uint64_t calc_min_expected_blocks_tail(
	struct gaba_dp_context_s const *this,
	struct gaba_joint_tail_s const *tail)
{
	uint64_t p = MIN2(
		this->w.r.alen - tail->apos,
		this->w.r.blen - tail->bpos);
	return((p + p/2) / BLK);
}

/**
 * @fn fill_mem_bounded
 * @brief fill <blk_cnt> contiguous blocks without seq bound tests, adding head and tail
 */
static _force_inline
struct gaba_joint_tail_s *fill_mem_bounded(
	struct gaba_dp_context_s *this,
	struct gaba_joint_tail_s const *prev_tail,
	uint64_t blk_cnt)
{
	struct gaba_joint_block_s h = fill_create_phantom_block(this, prev_tail);
	if(h.stat != CONT) {
		return(fill_create_tail(this, prev_tail, h.blk, h.p, h.stat));
	}

	struct gaba_joint_block_s b = fill_bulk_predetd_blocks(this, h.blk, blk_cnt);
	return(fill_create_tail(this, prev_tail, b.blk, h.p + b.p, b.stat));
}

/**
 * @fn fill_seq_bounded
 * @brief fill blocks with seq bound tests, adding head and tail
 */
static _force_inline
struct gaba_joint_tail_s *fill_seq_bounded(
	struct gaba_dp_context_s *this,
	struct gaba_joint_tail_s const *prev_tail)
{
	struct gaba_joint_block_s stat = fill_create_phantom_block(this, prev_tail);
	int64_t psum = stat.p;

	/* check if term detected in init fetch */
	if(stat.stat != CONT) {
		goto _fill_seq_bounded_finish;
	}

	/* calculate block size */
	uint64_t seq_bulk_blocks = calc_min_expected_blocks_blk(this, stat.blk);
	while(seq_bulk_blocks > MIN_BULK_BLOCKS) {
		/* bulk fill without ij-bound test */
		psum += (stat = fill_bulk_predetd_blocks(this, stat.blk, seq_bulk_blocks)).p;
		if(stat.stat != CONT) {
			goto _fill_seq_bounded_finish;	/* skip cap */
		}
		seq_bulk_blocks = calc_min_expected_blocks_blk(this, stat.blk);
	}

	/* bulk fill with ij-bound test */
	psum += (stat = fill_bulk_seq_bounded(this, stat.blk)).p;

	if(stat.stat != CONT) {
		goto _fill_seq_bounded_finish;	/* skip cap */
	}

	/* cap fill (without p-bound test) */
	psum += (stat = fill_cap_seq_bounded(this, stat.blk)).p;

_fill_seq_bounded_finish:;
	return(fill_create_tail(this, prev_tail, stat.blk, psum, stat.stat));
}

/**
 * @fn fill_section_seq_bounded
 * @brief fill dp matrix inside section pairs
 */
static _force_inline
struct gaba_joint_tail_s *fill_section_seq_bounded(
	struct gaba_dp_context_s *this,
	struct gaba_joint_tail_s const *prev_tail,
	struct gaba_section_s const *a,
	struct gaba_section_s const *b)
{
	/* init section and restore sequence reader buffer */
	fill_load_section(this, a, b, INT64_MAX);

	/* init tail pointer */
	struct gaba_joint_tail_s *tail = _tail(prev_tail);

	/* calculate block sizes */
	uint64_t mem_bulk_blocks = calc_max_bulk_blocks_mem(this);
	uint64_t seq_bulk_blocks = calc_max_expected_blocks_tail(this, tail);

	/* extra large bulk fill (with stack allocation) */
	debug("mem_bulk_blocks(%llu), seq_bulk_blocks(%llu)", mem_bulk_blocks, seq_bulk_blocks);
	while(_unlikely(mem_bulk_blocks < seq_bulk_blocks)) {
		uint64_t bulk_cnt = MIN2(
			mem_bulk_blocks,
			calc_min_expected_blocks_tail(this, tail));

		if(bulk_cnt > MIN_BULK_BLOCKS) {
			debug("mem bounded fill");
			if((tail = fill_mem_bounded(this, tail, bulk_cnt))->stat != CONT) {
				return(tail);
			}

			/* fill-in area has changed */
			seq_bulk_blocks = calc_max_expected_blocks_tail(this, tail);
		}

		/* malloc the next stack and set pointers */
		debug("add stack");
		if(gaba_dp_add_stack(this, 0) != GABA_SUCCESS) {
			return(NULL);
		}

		/* stack size has changed */
		mem_bulk_blocks = calc_max_bulk_blocks_mem(this);
		debug("mem_bulk_blocks(%llu), seq_bulk_blocks(%llu)", mem_bulk_blocks, seq_bulk_blocks);
	}

	debug("v(%p), psum(%lld), p(%d)", _last_block(tail)->md, tail->psum, tail->p);

	/* bulk fill with seq bound check */
	return(fill_seq_bounded(this, tail));
}

/**
 * @fn gaba_dp_fill_root
 *
 * @brief build_root API
 */
struct gaba_fill_s *suffix(gaba_dp_fill_root)(
	struct gaba_dp_context_s *this,
	struct gaba_section_s const *a,
	uint32_t apos,
	struct gaba_section_s const *b,
	uint32_t bpos)
{
	/* store section info */
	this->tail.apos = apos;
	this->tail.bpos = bpos;
	return(_fill(fill_section_seq_bounded(this, &this->tail, a, b)));
}

/**
 * @fn gaba_dp_fill
 *
 * @brief fill API
 */
struct gaba_fill_s *suffix(gaba_dp_fill)(
	struct gaba_dp_context_s *this,
	struct gaba_fill_s const *prev_sec,
	struct gaba_section_s const *a,
	struct gaba_section_s const *b)
{
	struct gaba_joint_tail_s const *tail = _tail(prev_sec);
	return(_fill(fill_section_seq_bounded(this, tail, a, b)));
}


/* trace leaf search functions */
/**
 * @struct gaba_leaf_s
 */
struct gaba_leaf_s {
	struct gaba_joint_tail_s const *tail;
	struct gaba_block_s const *blk;
	uint32_t aridx, bridx;
	int32_t p, q;
};

/**
 * @fn leaf_load_max_mask
 */
struct leaf_max_mask_s {
	vec_t max;
	int64_t offset;
	uint32_t mask_max;
};
static _force_inline
struct leaf_max_mask_s leaf_load_max_mask(
	struct gaba_dp_context_s *this,
	struct gaba_joint_tail_s const *tail)
{
	struct gaba_block_s *blk = _last_block(tail);

	debug("p(%d), psum(%lld), ssum(%d), offset(%lld)",
		tail->p, tail->psum, tail->ssum, blk->offset);

	/* load max vector, create mask */
	vec_t max = _load(&blk->sd.max);
	int64_t offset = blk->offset;
	uint32_t mask_max = ((vec_masku_t){
		.mask = _mask_v32i16(_eq_v32i16(
			_set_v32i16(tail->max - offset),
			_add_v32i16(_load_v32i16(_last_block(tail)->md), _cvt_v32i8_v32i16(max))))
	}).all;
	debug("mask_max(%x)", mask_max);
	_print_v32i16(_set_v32i16(tail->max - offset));
	_print_v32i16(_add_v32i16(_load_v32i16(_last_block(tail)->md), _cvt_v32i8_v32i16(max)));

	return((struct leaf_max_mask_s){
		.max = max,
		.offset = offset,
		.mask_max = mask_max
	});
}

/**
 * @fn leaf_detect_max_block
 */
struct leaf_max_block_s {
	vec_t max;
	struct gaba_block_s *blk;
	int32_t p;
	uint32_t mask_max;
};
static _force_inline
struct leaf_max_block_s leaf_detect_max_block(
	struct gaba_dp_context_s *this,
	struct gaba_joint_tail_s const *tail,
	int64_t offset,
	uint32_t mask_max,
	vec_t max)
{
	/* scan blocks backward */
	struct gaba_block_s *blk = _last_block(tail);
	int32_t p = -1;

	/* b must be sined integer, in order to detect negative index. */
	for(int32_t b = (tail->p - 1)>>BLK_BASE; b >= 0; b--, blk--) {

		/* load the previous max vector and offset */
		vec_t prev_max = _load(&(blk - 1)->sd.max);
		int64_t prev_offset = (blk - 1)->offset;

		/* adjust offset */
		max = _add(max, _set(offset - prev_offset));

		/* take mask */
		uint32_t prev_mask_max = mask_max & ((vec_masku_t){
			.mask = _mask(_eq(prev_max, max))
		}).all;

		debug("scan block: b(%d), offset(%lld), mask_max(%u), prev_mask_max(%u)",
			b, offset, mask_max, prev_mask_max);

		if(prev_mask_max == 0) {
			debug("block found: blk(%p), p(%d), mask_max(%u)", blk, b * BLK, mask_max);
			p = b * BLK;
			break;
		}

		/* windback a block */
		max = prev_max;
		offset = prev_offset;
		mask_max = prev_mask_max;
	}

	debug("loop break: blk(%p), p(%d), mask_max(%u)", blk, p, mask_max);
	return((struct leaf_max_block_s){
		.max = max,
		.blk = blk,
		.p = p,
		.mask_max = mask_max
	});
}

/**
 * @fn leaf_refill_block
 */
static _force_inline
void leaf_refill_block(
	struct gaba_dp_context_s *this,
	vec_masku_t *mask_max_ptr,
	int64_t len,
	struct gaba_block_s *blk,
	vec_t compd_max)
{
	/* fetch from existing blocks */
	fill_restore_fetch(this, blk);

	/* vectors on registers inside this block */ {
		#define _fill_block_leaf(_mask_ptr) { \
			_dir_fetch(dir); \
			if(_dir_is_right(dir)) { \
				_fill_right_update_ptr(); \
				_fill_right(); \
			} else { \
				_fill_down_update_ptr(); \
				_fill_down(); \
			} \
			(_mask_ptr)++->mask = _mask(_eq(max, delta)); \
			debug("mask(%x)", ((vec_masku_t){ .mask = _mask(_eq(max, delta)) }).all); \
		}

		/* load contexts and overwrite max vector */
		_fill_load_context(blk);
		max = compd_max;		/* overwrite with compensated max vector */
		(void)offset;			/* to avoid warning */

		for(int64_t i = 0; i < len; i++) {
			_fill_block_leaf(mask_max_ptr);
		}
	}
	return;
}

/**
 * @fn leaf_detect_max_pos
 */
struct leaf_max_pos_s {
	int32_t p;
	int32_t q;
};
static _force_inline
struct leaf_max_pos_s leaf_detect_max_pos(
	struct gaba_dp_context_s *this,
	vec_masku_t *mask_max_ptr,
	int64_t len,
	uint32_t mask_max)
{
	for(int64_t i = 0; i < len; i++) {
		uint32_t mask_update = mask_max_ptr[i].all & mask_max;
		if(mask_update != 0) {
			return((struct leaf_max_pos_s){
				.p = i,
				.q = tzcnt(mask_update)
			});
		}
	}

	/* not found in the block (never reaches here) */
	debug("max pos NOT found.");
	return((struct leaf_max_pos_s){
		.p = 0,
		.q = 0
	});
}

/**
 * @fn leaf_save_coordinates
 */
static _force_inline
void leaf_save_coordinates(
	struct gaba_dp_context_s *this,
	struct gaba_joint_tail_s const *tail,
	struct gaba_leaf_s *leaf,
	struct gaba_block_s const *blk,
	int32_t p,
	int32_t q)
{
	leaf->tail = tail;
	leaf->blk = blk;

	/* calc ridx */
	int64_t mask_idx = p & (BLK - 1);
	int64_t filled_count = mask_idx + 1;
	int32_t bcnt = _dir_bcnt(_dir_load(blk, mask_idx));
	int32_t acnt = filled_count - bcnt;
	v2i32_t ridx = _add_v2i32(
		_load_v2i32(&(blk - 1)->aridx),
		_seta_v2i32((BLK - 1 - q) - bcnt, q - acnt));
	_store_v2i32(&leaf->aridx, ridx);
	debug("idx(%lld), fcnt(%lld), p(%d), q(%d), cnt(%d, %d), ridx(%u, %u)",
		mask_idx, filled_count, p, q, bcnt, acnt, _hi32(ridx), _lo32(ridx));
	
	/* store p and q */
	leaf->p = p;
	leaf->q = q;
	return;
}

/**
 * @fn leaf_save_phantom_coordinates
 */
static _force_inline
void leaf_save_phantom_coordinates(
	struct gaba_dp_context_s *this,
	struct gaba_joint_tail_s const *tail,
	struct gaba_leaf_s *leaf,
	struct gaba_block_s const *blk,
	uint32_t mask_max)
{
	leaf->tail = tail;
	leaf->blk = blk;
	_store_v2i32(&leaf->aridx, _zero_v2i32());
	leaf->p = -1;
	leaf->q = tzcnt(mask_max);
	return;
}

/**
 * @fn leaf_search
 */
static _force_inline
void leaf_search(
	struct gaba_dp_context_s *this,
	struct gaba_joint_tail_s const *tail,
	struct gaba_leaf_s *leaf)
{
	/* load max vector and create mask */
	struct leaf_max_mask_s m = leaf_load_max_mask(this, tail);

	/* search block */
	struct leaf_max_block_s b = leaf_detect_max_block(
		this, tail, m.offset, m.mask_max, m.max);
	debug("check p(%d)", b.p);
	if(b.p == -1) {
		leaf_save_phantom_coordinates(this, tail, leaf, b.blk, b.mask_max);
		return;
	}

	/* refill detected block */
	int64_t len = MIN2(tail->p - b.p, BLK);
	vec_masku_t mask_max_arr[BLK];
	leaf_refill_block(this, mask_max_arr, len, b.blk, b.max);

	/* detect pos */
	struct leaf_max_pos_s l = leaf_detect_max_pos(
		this, mask_max_arr, len, b.mask_max);

	/* calculate ij-coordinates */
	leaf_save_coordinates(this, tail, leaf, b.blk, b.p + l.p, l.q);
	return;
}

/**
 * @fn gaba_dp_search_max
 */
struct gaba_pos_pair_s suffix(gaba_dp_search_max)(
	struct gaba_dp_context_s *this,
	struct gaba_fill_s const *tail)
{
	struct gaba_leaf_s leaf;
	leaf_search(this, _tail(tail), &leaf);

	struct gaba_joint_tail_s const *atail = _tail(tail), *btail = _tail(tail);
	int32_t alen = atail->alen, blen = btail->blen;
	int32_t aidx = alen - leaf.aridx, bidx = blen - leaf.bridx;

	while(aidx <= 0) {
		for(atail = atail->tail; (atail->stat & GABA_STATUS_UPDATE_A) == 0; atail = atail->tail) {}
		aidx += (alen = atail->alen);
	}
	while(bidx <= 0) {
		for(btail = btail->tail; (btail->stat & GABA_STATUS_UPDATE_B) == 0; btail = btail->tail) {}
		bidx += (blen = btail->blen);
	}
	return((struct gaba_pos_pair_s){
		.apos = aidx - 1,
		.bpos = bidx - 1
	});
}

/* path trace functions */
/**
 * @fn trace_load_section_a, trace_load_section_b
 */
static _force_inline
void trace_load_section_a(
	struct gaba_dp_context_s *this)
{
	debug("load section a, idx(%d), len(%d)", this->w.l.aidx, this->w.l.atail->alen);

	/* load tail pointer (must be inited with leaf tail) */
	struct gaba_joint_tail_s const *tail = this->w.l.atail;
	int32_t len = tail->alen;
	int32_t sum = len;
	int32_t idx = this->w.l.aidx + len;

	while(idx <= 0) {
		for(tail = tail->tail; (tail->stat & GABA_STATUS_UPDATE_A) == 0; tail = tail->tail) {}
		len = tail->alen; sum += len; idx += len;
	}

	/* reload finished, store section info */
	// this->w.l.atail = tail->tail;
	this->w.l.atail = tail;		/* fixme: is this correct?? */
	this->w.l.alen = len;
	this->w.l.aid = tail->aid;
	this->w.l.asum = sum;

	this->w.l.aidx = idx;
	this->w.l.asidx = idx;
	return;
}
static _force_inline
void trace_load_section_b(
	struct gaba_dp_context_s *this)
{
	debug("load section b, idx(%d), len(%d)", this->w.l.bidx, this->w.l.btail->blen);

	/* load tail pointer (must be inited with leaf tail) */
	struct gaba_joint_tail_s const *tail = this->w.l.btail;
	int32_t len = tail->blen;
	int32_t sum = len;
	int32_t idx = this->w.l.bidx + len;

	while(idx <= 0) {
		for(tail = tail->tail; (tail->stat & GABA_STATUS_UPDATE_B) == 0; tail = tail->tail) {}
		len = tail->blen; sum += len; idx += len;
	}

	/* reload finished, store section info */
	// this->w.l.btail = tail->tail;
	this->w.l.btail = tail;		/* fixme: is this correct?? */
	this->w.l.blen = len;
	this->w.l.bid = tail->bid;
	this->w.l.bsum = sum;

	this->w.l.bidx = idx;
	this->w.l.bsidx = idx;
	return;
}

/**
 * @macro _trace_*_load_context
 * @brief load context onto registers
 */
#define _trace_load_context(t) \
	register v2i32_t idx = _load_v2i32(&(t)->w.l.aidx), gc = _zero_v2i32(); \
	static uint32_t const hc[2] __attribute__(( aligned(32) )) = { 0, -1 }; \
	static uint32_t const dc[2] __attribute__(( aligned(32) )) = { 0, 0 }; \
	static uint32_t const vc[2] __attribute__(( aligned(32) )) = { -1, 0 }; \
	v2i32_t const hterm = _load_v2i32(hc); \
	v2i32_t const dterm = _load_v2i32(dc); \
	v2i32_t const vterm = _load_v2i32(vc); \
	_print_v2i32(hterm); \
	_print_v2i32(dterm); \
	_print_v2i32(vterm); \
	struct gaba_block_s const *blk = (t)->w.l.blk; \
	int64_t p = (t)->w.l.p, q = (t)->w.l.q; \
	union gaba_dir_u dir = _dir_load(blk, p & (BLK - 1)); \
	union gaba_mask_pair_u const *ptr = &blk->mask[p & (BLK - 1)]; \
	_print_v2i32(idx); \
	debug("p(%lld), q(%lld), mask_h(%x), mask_v(%x), path_array(%llx)", \
		p, q, ptr->pair.h.all, ptr->pair.v.all, path_array);

#define _trace_forward_load_context(t) \
	uint32_t *path = (t)->w.l.path.phead; \
	int64_t ofs = (t)->w.l.path.phofs; \
	uint64_t path_array = _loadu_u64(path - 1)>>(2*BLK - ofs); \
	_trace_load_context(t);

#define _trace_reverse_load_context(t) \
	uint32_t *path = (t)->w.l.path.ptail; \
	int64_t ofs = (t)->w.l.path.ptofs; \
	uint64_t path_array = _loadu_u64(path)<<(2*BLK - ofs); \
	_trace_load_context(t);


/**
 * @macro _trace_reload_ptr
 */
#define _trace_reload_ptr(_idx) { \
	ptr = &(--blk)->mask[(_idx)]; \
	dir = _dir_load(blk, (_idx)); \
}

/**
 * @macro _trace_reload_tail
 */
#define _trace_reload_tail(t) { \
	debug("tail(%p), next tail(%p), p(%d), psum(%lld), ssum(%d)", \
		(t)->w.l.tail, (t)->w.l.tail->tail, (t)->w.l.tail->tail->p, \
		(t)->w.l.tail->tail->psum, (t)->w.l.tail->tail->ssum); \
	/* update psum */ \
	(t)->w.l.psum -= (t)->w.l.p; \
	/* load section lengths */ \
	struct gaba_joint_tail_s const *tail = (t)->w.l.tail; \
	register v2i32_t len = _load_v2i32(&tail->alen); \
	/* reload tail */ \
	tail = (t)->w.l.tail = tail->tail; \
	blk = _last_block(tail) + 1; \
	p = ((t)->w.l.p = tail->p) - 1; \
	debug("updated psum(%lld), w.l.p(%d), p(%lld)", (t)->w.l.psum, (t)->w.l.p, p); \
	/* adjust sum lengths */ \
	static uint32_t const mc[2] __attribute__(( aligned(32) )) = { \
		GABA_STATUS_UPDATE_A, GABA_STATUS_UPDATE_B \
	}; \
	v2i32_t const mask = _load_v2i32(mc); \
	register v2i32_t stat = _set_v2i32(tail->stat); \
	register v2i32_t sum = _load_v2i32(&(t)->w.l.asum); \
	sum = _sub_v2i32(sum, _and_v2i32(_eq_v2i32(_and_v2i32(stat, mask), mask), len)); \
	_store_v2i32(&(t)->w.l.asum, sum); \
	debug("adjusted sum(%u, %u), len(%u, %u), stat(%u)", _hi32(sum), _lo32(sum), _hi32(len), _lo32(len), tail->stat); \
	/* reload dir and mask pointer */ \
	_trace_reload_ptr(p & (BLK - 1)); \
}

/**
 * @macro _trace_*_cap_update_path
 * @brief store path array and update path pointer
 */
#define _trace_forward_cap_update_path() { \
	int64_t _cnt = (p & (BLK - 1)) + 1 - (ptr - blk->mask + 1); \
	debug("path_array(%llx), cnt(%lld), ptr(%p), base(%p)", path_array, _cnt, ptr, blk->mask - 1); \
	_storeu_u64(path - 1, path_array<<(2*BLK - (ofs + _cnt))); \
	path -= ((ofs + (_cnt)) & BLK) != 0; \
	ofs = (ofs + (_cnt)) & (BLK - 1); \
	p -= _cnt; \
}
#define _trace_reverse_cap_update_path() { \
	int64_t _cnt = (p & (BLK - 1)) + 1 - (ptr - blk->mask + 1); \
	debug("path_array(%llx), cnt(%lld), ptr(%p), base(%p)", path_array, _cnt, ptr, blk->mask - 1); \
	_storeu_u64(path, path_array>>(2*BLK - (ofs + _cnt))); \
	path += ((ofs + (_cnt)) & BLK) != 0; \
	ofs = (ofs + (_cnt)) & (BLK - 1); \
	p -= _cnt; \
}

/**
 * @macro _trace_bulk_*_update_path
 * @brief store path array
 */
#define _trace_forward_bulk_update_path() { \
	debug("path_array(%llx)", path_array); \
	_storeu_u64(path - 1, path_array<<(BLK - ofs)); \
	path--; \
	p -= BLK; \
}
#define _trace_reverse_bulk_update_path() { \
	debug("path_array(%llx)", path_array); \
	_storeu_u64(path, path_array>>(BLK - ofs)); \
	path++; \
	p -= BLK; \
}

/**
 * @macro _trace_*_calc_index
 * @brief restore index from ridx and q
 */
#define _trace_forward_calc_index(t) { \
	_print_v2i32(idx); \
	/* calc idx of the head of the block from ridx */ \
	v2i32_t ridx = _load_v2i32(&blk->aridx); \
	v2i32_t sum = _load_v2i32(&(t)->w.l.asum); \
	idx = _sub_v2i32(_sub_v2i32(sum, ridx), _seta_v2i32((BW - 1) - q, q)); \
	debug("calc_index p(%lld), q(%lld)", p, q); \
	_print_v2i32(ridx); \
	_print_v2i32(sum); \
	_print_v2i32(idx); \
}
#define _trace_reverse_calc_index(t) { \
	_print_v2i32(idx); \
	/* calc idx of the head of the block from ridx */ \
	v2i32_t ridx = _load_v2i32(&blk->aridx); \
	v2i32_t sum = _load_v2i32(&(t)->w.l.asum); \
	idx = _sub_v2i32(_sub_v2i32(sum, ridx), _seta_v2i32((BW - 1) - q, q)); \
	debug("calc_index p(%lld), q(%lld)", p, q); \
	_print_v2i32(ridx); \
	_print_v2i32(sum); \
	_print_v2i32(idx); \
}

/**
 * @macro _trace_forward_*_load
 */
#define _trace_forward_head_load(t, _jump_to) { \
	if(_unlikely(ptr == blk->mask - 1)) { \
		_trace_forward_cap_update_path(); \
		_trace_reload_ptr(BLK - 1); \
		debug("jump to %s", #_jump_to); \
		goto _jump_to; \
	} \
}
#define _trace_forward_bulk_load(t, _jump_to) { \
	if(_unlikely(ptr == blk->mask - 1)) { \
		_trace_forward_bulk_update_path(); \
		_trace_reload_ptr(BLK - 1); \
		if(_unlikely(p < BLK)) { \
			_trace_forward_calc_index(t); \
			debug("jump to %s", #_jump_to); \
			goto _jump_to; \
		} \
	} \
}
#define _trace_forward_tail_load(t, _jump_to) { \
	if(_unlikely(ptr == blk->mask - 1)) { \
		debug("load block, blk(%p), next_blk(%p), p(%lld)", \
			blk, blk-1, p); \
		/* store path_array */ \
		_trace_forward_cap_update_path(); \
		if(p < 0) { \
			debug("w.l.psum(%lld), w.l.p(%d), p(%lld)", (t)->w.l.psum, (t)->w.l.p, p); \
			if((t)->w.l.psum < (t)->w.l.p - p) { \
				goto _trace_forward_index_break; \
			} \
			_trace_reload_tail(t); \
			debug("jump to %s", #_jump_to); \
			goto _jump_to; \
		} \
		/* load dir and update mask pointer */ \
		_trace_reload_ptr(BLK - 1); \
	} \
}

/**
 * @macro _trace_reverse_*_load
 */
#define _trace_reverse_head_load(t, _jump_to) { \
	if(_unlikely(ptr == blk->mask - 1)) { \
		_trace_reverse_cap_update_path(); \
		_trace_reload_ptr(BLK - 1); \
		debug("jump to %s", #_jump_to); \
		goto _jump_to; \
	} \
}
#define _trace_reverse_bulk_load(t, _jump_to) { \
	if(_unlikely(ptr == blk->mask - 1)) { \
		_trace_reverse_bulk_update_path(); \
		_trace_reload_ptr(BLK - 1); \
		if(_unlikely(p < BLK)) { \
			_trace_reverse_calc_index(t); \
			debug("jump to %s", #_jump_to); \
			goto _jump_to; \
		} \
	} \
}
#define _trace_reverse_tail_load(t, _jump_to) { \
	if(_unlikely(ptr == blk->mask - 1)) { \
		debug("load block, blk(%p), next_blk(%p), p(%lld)", \
			blk, blk-1, p); \
		/* store path_array */ \
		_trace_reverse_cap_update_path(); \
		if(p < 0) { \
			debug("w.l.psum(%lld), w.l.p(%d), p(%lld)", (t)->w.l.psum, (t)->w.l.p, p); \
			if((t)->w.l.psum < (t)->w.l.p - p) { \
				goto _trace_reverse_index_break; \
			} \
			_trace_reload_tail(t); \
			debug("jump to %s", #_jump_to); \
			goto _jump_to; \
		} \
		/* load dir and update mask pointer */ \
		_trace_reload_ptr(BLK - 1); \
	} \
}

/**
 * @macro _trace_inc_*
 * @brief increment gap counters
 */
#define _trace_inc_gi()	{ \
	gc = _sub_v2i32(gc, vterm); \
}
#define _trace_inc_ge()	{ \
	gc = _sub_v2i32(gc, hterm); \
}

/**
 * @macro _trace_head_*_test_index
 * @brief test and update indices
 */
#define _trace_head_v_test_index()		( 0 )
#define _trace_head_d_test_index()		( 0 )
#define _trace_head_h_test_index()		( 0 )

#define _trace_bulk_v_test_index()		( 0 )
#define _trace_bulk_d_test_index()		( 0 )
#define _trace_bulk_h_test_index()		( 0 )

#define _trace_tail_v_test_index() ( \
	_mask_v2i32(_eq_v2i32(idx, vterm)) \
)
#define _trace_tail_d_test_index() ( \
	_mask_v2i32(_eq_v2i32(idx, dterm)) \
)
#define _trace_tail_h_test_index() ( \
	_mask_v2i32(_eq_v2i32(idx, hterm)) \
)

/**
 * @macro _trace_bulk_*_update_index
 * @brief test and update indices
 */
#define _trace_head_v_update_index()	;
#define _trace_head_h_update_index()	;

#define _trace_bulk_v_update_index()	;
#define _trace_bulk_h_update_index()	;

#define _trace_tail_v_update_index() { \
	idx = _add_v2i32(idx, hterm); \
	_print_v2i32(hterm); \
	_print_v2i32(idx); \
}
#define _trace_tail_h_update_index() { \
	idx = _add_v2i32(idx, vterm); \
	_print_v2i32(vterm); \
	_print_v2i32(idx); \
}

/**
 * @macro _trace_test_*
 * @brief test mask
 */
#if MODEL == LINEAR
#define _trace_test_diag_h()			( (ptr->pair.h.all>>q) & 0x01 )
#define _trace_test_diag_v()			( (ptr->pair.v.all>>q) & 0x01 )
#define _trace_test_gap_h()				( (ptr->pair.h.all>>q) & 0x01 )
#define _trace_test_gap_v()				( (ptr->pair.v.all>>q) & 0x01 )
#else /* MODEL == AFFINE */
#define _trace_test_diag_h()			( (ptr->pair.h.all>>q) & 0x01 )
#define _trace_test_diag_v()			( (ptr->pair.v.all>>q) & 0x01 )
#define _trace_test_gap_h()				( (ptr->pair.e.all>>q) & 0x01 )
#define _trace_test_gap_v()				( (ptr->pair.f.all>>q) & 0x01 )
#endif

/**
 * @macro _trace_*_*_update_path_q
 */
#define _trace_forward_h_update_path_q() { \
	path_array = path_array<<1; \
	ptr--; \
	q += _dir_is_down(dir); \
	_dir_windback(dir); \
}
#define _trace_forward_v_update_path_q() { \
	path_array = (path_array<<1) | 0x01; \
	ptr--; \
	q += _dir_is_down(dir) - 1; \
	_dir_windback(dir); \
}
#define _trace_reverse_h_update_path_q() { \
	path_array = path_array>>1; \
	ptr--; \
	q += _dir_is_down(dir); \
	_dir_windback(dir); \
}
#define _trace_reverse_v_update_path_q() { \
	path_array = (path_array>>1) | 0x8000000000000000; \
	ptr--; \
	q += _dir_is_down(dir) - 1; \
	_dir_windback(dir); \
}

/**
 * @macro _trace_*_save_context
 * @brief save context to working buffer
 */
#define _trace_save_context(t) { \
	(t)->w.l.blk = blk; \
	_store_v2i32(&(t)->w.l.aidx, idx); \
	_print_v2i32(idx); \
	(t)->w.l.psum -= (t)->w.l.p - p; \
	(t)->w.l.p = p; \
	(t)->w.l.q = q; \
	_store_v2i32(&(t)->w.l.path.gic, _add_v2i32(gc, _load_v2i32(&(t)->w.l.path.gic))); \
	debug("p(%lld), psum(%lld), q(%llu)", p, (t)->w.l.psum, q); \
}
#define _trace_forward_save_context(t) { \
	(t)->w.l.path.phead = path; \
	(t)->w.l.path.phofs = ofs; \
	_trace_save_context(t); \
}
#define _trace_reverse_save_context(t) { \
	(t)->w.l.path.ptail = path; \
	(t)->w.l.path.ptofs = ofs; \
	_trace_save_context(t); \
}

/**
 * @fn trace_forward_body
 */
static _force_inline
void trace_forward_body(
	struct gaba_dp_context_s *this)
{
	#define _trace_forward_gap_loop(t, _type, _next, _label) { \
		while(1) { \
		_trace_forward_##_type##_##_label##_head: \
			if(_trace_test_gap_##_label() == 0) { \
				goto _trace_forward_##_type##_d_head; \
			} \
			if(_trace_##_type##_##_label##_test_index()) { \
				goto _trace_forward_index_break; \
			} \
			_trace_inc_ge(); \
			debug("go %s (%s), dir(%llx), mask_h(%x), mask_v(%x), p(%lld), q(%lld), ptr(%p), path_array(%llx)", \
				#_label, #_type, ((uint64_t)dir.dynamic.array), ptr->pair.h.all, ptr->pair.v.all, p, q, ptr, path_array); \
			_trace_##_type##_##_label##_update_index(); \
			_trace_forward_##_label##_update_path_q(); \
			_trace_forward_##_type##_load(t, _trace_forward_##_next##_##_label##_head); \
		} \
	}

	#define _trace_forward_diag_loop(t, _type, _next) { \
		while(1) { \
		_trace_forward_##_type##_d_head: \
			if(_trace_test_diag_h() != 0) { \
				_trace_inc_gi(); goto _trace_forward_##_type##_h_head; \
			} \
			if(_trace_##_type##_d_test_index()) { \
				_trace_forward_cap_update_path(); \
				goto _trace_forward_index_break; \
			} \
			debug("go d (%s), dir(%llx), mask_h(%x), mask_v(%x), p(%lld), q(%lld), ptr(%p), path_array(%llx)", \
				#_type, ((uint64_t)dir.dynamic.array), ptr->pair.h.all, ptr->pair.v.all, p, q, ptr, path_array); \
			_trace_##_type##_h_update_index(); \
			_trace_forward_h_update_path_q(); \
			_trace_forward_##_type##_load(t, _trace_forward_##_next##_d_mid); \
		_trace_forward_##_type##_d_mid: \
			_trace_##_type##_v_update_index(); \
			_trace_forward_v_update_path_q(); \
			_trace_forward_##_type##_load(t, _trace_forward_##_next##_d_tail); \
		_trace_forward_##_type##_d_tail: \
			if(_trace_test_diag_v() != 0) { \
				_trace_inc_gi(); goto _trace_forward_##_type##_v_head; \
			} \
		} \
	}

	_trace_forward_load_context(this);

	/* v loop */
	_trace_forward_loop_v_head: {
		if(p < 2 * BLK) {
			goto _trace_forward_tail_v_head;
		}
		_trace_forward_gap_loop(this, head, bulk, v);
		_trace_forward_gap_loop(this, bulk, tail, v);
		_trace_forward_gap_loop(this, tail, loop, v);
	}

	/* d dispatchers */
	_trace_forward_loop_d_mid: {
		if(p < 2 * BLK) {
			goto _trace_forward_tail_d_mid;
		} else {
			goto _trace_forward_head_d_mid;
		}
	}
	_trace_forward_loop_d_tail: {
		if(p < 2 * BLK) {
			goto _trace_forward_tail_d_tail;
		} else {
			goto _trace_forward_head_d_tail;
		}
	}
	/* d loop */ {
		_trace_forward_diag_loop(this, head, bulk);
		_trace_forward_diag_loop(this, bulk, tail);
		_trace_forward_diag_loop(this, tail, loop);
	}

	/* h loop */
	_trace_forward_loop_h_head: {
		if(p < 2 * BLK) {
			goto _trace_forward_tail_h_head;
		}
		_trace_forward_gap_loop(this, head, bulk, h);
		_trace_forward_gap_loop(this, bulk, tail, h);
		_trace_forward_gap_loop(this, tail, loop, h);
	}

_trace_forward_index_break:;
	_trace_forward_save_context(this);
	return;
}

/**
 * @fn trace_reverse_body
 */
static _force_inline
void trace_reverse_body(
	struct gaba_dp_context_s *this)
{
	#define _trace_reverse_gap_loop(t, _type, _next, _label) { \
		while(1) { \
		_trace_reverse_##_type##_##_label##_head: \
			if(_trace_test_gap_##_label() == 0) { \
				goto _trace_reverse_##_type##_d_head; \
			} \
			if(_trace_##_type##_##_label##_test_index()) { \
				goto _trace_reverse_index_break; \
			} \
			_trace_inc_ge(); \
			debug("go %s (%s), dir(%llx), mask_h(%x), mask_v(%x), p(%lld), q(%lld), ptr(%p), path_array(%llx)", \
				#_label, #_type, ((uint64_t)dir.dynamic.array), ptr->pair.h.all, ptr->pair.v.all, p, q, ptr, path_array); \
			_trace_##_type##_##_label##_update_index(); \
			_trace_reverse_##_label##_update_path_q(); \
			_trace_reverse_##_type##_load(t, _trace_reverse_##_next##_##_label##_head); \
		} \
	}

	#define _trace_reverse_diag_loop(t, _type, _next) { \
		while(1) { \
		_trace_reverse_##_type##_d_head: \
			if(_trace_test_diag_v() != 0) { \
				_trace_inc_gi(); goto _trace_reverse_##_type##_v_head; \
			} \
			if(_trace_##_type##_d_test_index()) { \
				_trace_reverse_cap_update_path(); \
				goto _trace_reverse_index_break; \
			} \
			debug("go d (%s), dir(%llx), mask_h(%x), mask_v(%x), p(%lld), q(%lld), ptr(%p), path_array(%llx)", \
				#_type, ((uint64_t)dir.dynamic.array), ptr->pair.h.all, ptr->pair.v.all, p, q, ptr, path_array); \
			_trace_##_type##_v_update_index(); \
			_trace_reverse_v_update_path_q(); \
			_trace_reverse_##_type##_load(t, _trace_reverse_##_next##_d_mid); \
		_trace_reverse_##_type##_d_mid: \
			_trace_##_type##_h_update_index(); \
			_trace_reverse_h_update_path_q(); \
			_trace_reverse_##_type##_load(t, _trace_reverse_##_next##_d_tail); \
		_trace_reverse_##_type##_d_tail: \
			if(_trace_test_diag_h() != 0) { \
				_trace_inc_gi(); goto _trace_reverse_##_type##_h_head; \
			} \
		} \
	}

	_trace_reverse_load_context(this);

	/* h loop */
	_trace_reverse_loop_h_head: {
		if(p < 2 * BLK) {
			goto _trace_reverse_tail_h_head;
		}
		_trace_reverse_gap_loop(this, head, bulk, h);
		_trace_reverse_gap_loop(this, bulk, tail, h);
		_trace_reverse_gap_loop(this, tail, loop, h);
	}

	/* d dispatchers */
	_trace_reverse_loop_d_mid: {
		if(p < 2 * BLK) {
			goto _trace_reverse_tail_d_mid;
		} else {
			goto _trace_reverse_head_d_mid;
		}
	}
	_trace_reverse_loop_d_tail: {
		if(p < 2 * BLK) {
			goto _trace_reverse_tail_d_tail;
		} else {
			goto _trace_reverse_head_d_tail;
		}
	}
	/* d loop */ {
		_trace_reverse_diag_loop(this, head, bulk);
		_trace_reverse_diag_loop(this, bulk, tail);
		_trace_reverse_diag_loop(this, tail, loop);
	}

	/* v loop */
	_trace_reverse_loop_v_head: {
		if(p < 2 * BLK) {
			goto _trace_reverse_tail_v_head;
		}
		_trace_reverse_gap_loop(this, head, bulk, v);
		_trace_reverse_gap_loop(this, bulk, tail, v);
		_trace_reverse_gap_loop(this, tail, loop, v);
	}

_trace_reverse_index_break:;
	_trace_reverse_save_context(this);
	return;
}

/**
 * @fn trace_forward_push
 */
static _force_inline
void trace_forward_push(
	struct gaba_dp_context_s *this)
{
	/* windback pointer */
	this->w.l.path.shead--;

	/* load section info */
	v2i32_t id = _load_v2i32(&this->w.l.aid);
	v2i32_t idx = _load_v2i32(&this->w.l.aidx);
	v2i32_t sidx = _load_v2i32(&this->w.l.asidx);

	/* adjust breakpoint */
	uint64_t const path_array = *((uint64_t *)this->w.l.path.phead)>>(32 - this->w.l.path.phofs);
	v2i32_t mask = _eq_v2i32(idx, _zero_v2i32());
	v2i32_t adj = _and_v2i32(
		_andn_v2i32(mask, _swap_v2i32(mask)),
		_seta_v2i32(tzcnt(~path_array) - 1, tzcnt(path_array)));
	idx = _min_v2i32(_add_v2i32(idx, adj), sidx);

	debug("path_array(%llx), adj(%u, %u), mask(%u, %u), idx(%u, %u), sidx(%u, %u)",
		path_array,
		_hi32(adj), _lo32(adj),
		_hi32(mask), _lo32(mask),
		_hi32(idx), _lo32(idx),
		_hi32(sidx), _lo32(sidx));

	/* calc path length */
	v2i32_t tlen = _sub_v2i32(sidx, idx);
	// int64_t plen = _lo32(tlen) + _hi32(tlen);

	/* store section info */
	_store_v2i32(&this->w.l.path.shead->aid, id);
	_store_v2i32(&this->w.l.path.shead->apos, idx);
	_store_v2i32(&this->w.l.path.shead->alen, tlen);

	/* store path length */
	// this->w.l.path.shead->plen = plen;
	this->w.l.path.shead->ppos = 0;

	debug("push current section forward a(%u, %u, %u), b(%u, %u, %u), len(%u)",
		this->w.l.path.shead->aid,
		this->w.l.path.shead->apos,
		this->w.l.path.shead->alen,
		this->w.l.path.shead->bid,
		this->w.l.path.shead->bpos,
		this->w.l.path.shead->blen,
		_plen(this->w.l.path.shead));

	/* update rsidx */
	_store_v2i32(&this->w.l.asidx, idx);
	return;
}

/**
 * @fn trace_reverse_push
 */
static _force_inline
void trace_reverse_push(
	struct gaba_dp_context_s *this)
{
	/* push section info to section array */
	v2i32_t const mask = _set_v2i32(0x01);
	v2i32_t len = _load_v2i32(&this->w.l.alen);
	v2i32_t id = _load_v2i32(&this->w.l.aid);
	v2i32_t idx = _load_v2i32(&this->w.l.aidx);
	v2i32_t sidx = _load_v2i32(&this->w.l.asidx);

	/* calc path pos and len */
	v2i32_t tlen = _sub_v2i32(sidx, idx);
	int64_t plen = _lo32(tlen) + _hi32(tlen);
	int64_t ppos = this->w.l.pspos;

	/* store revcomped section */
	_store_v2i32(&this->w.l.path.stail->aid, _xor_v2i32(id, mask));
	_store_v2i32(&this->w.l.path.stail->apos, _sub_v2i32(len, sidx));
	_store_v2i32(&this->w.l.path.stail->alen, tlen);

	/* store path length */
	// this->w.l.path.stail->plen = plen;
	this->w.l.path.stail->ppos = ppos;

	debug("push current section reverse a(%u, %u, %u), b(%u, %u, %u), pos(%lld), len(%u)",
		this->w.l.path.stail->aid,
		this->w.l.path.stail->apos,
		this->w.l.path.stail->alen,
		this->w.l.path.stail->bid,
		this->w.l.path.stail->bpos,
		this->w.l.path.stail->blen,
		this->w.l.path.stail->ppos,
		_plen(this->w.l.path.stail));

	/* update rsidx */
	_store_v2i32(&this->w.l.asidx, idx);

	/* update pspos */
	this->w.l.pspos = ppos + plen;

	/* windback pointer */
	this->w.l.path.stail++;
	return;
}

/**
 * @fn trace_init_work
 */
static _force_inline
void trace_init_work(
	struct gaba_dp_context_s *this,
	struct gaba_leaf_s const *leaf,
	struct gaba_path_intl_s const *path)
{
	/* store tail pointers */
	struct gaba_joint_tail_s const *tail = leaf->tail;
	this->w.l.tail = tail;
	this->w.l.atail = tail;
	this->w.l.btail = tail;

	/* store path object and section array object */
	this->w.l.path = *path;
	struct gaba_block_s const *blk = leaf->blk;
	this->w.l.blk = blk;

	#if 0
	/* copy ids and lengths */
	v2i32_t len = _load_v2i32(&tail->alen);
	v2i32_t id = _load_v2i32(&tail->aid);
	_store_v2i32(&this->w.l.alen, len);
	_store_v2i32(&this->w.l.aid, id);

	v2i32_t ridx = _load_v2i32(&leaf->aridx);
	_store_v2i32(&this->w.l.aidx, _sub_v2i32(len, ridx));
	_store_v2i32(&this->w.l.asidx, _sub_v2i32(len, ridx));
	#else
	/* init ids and lengths with invalid */
	_store_v2i32(&this->w.l.alen, _zero_v2i32());
	_store_v2i32(&this->w.l.aid, _set_v2i32(-1));

	v2i32_t idx = _sub_v2i32(_zero_v2i32(), _load_v2i32(&leaf->aridx));
	_store_v2i32(&this->w.l.aidx, idx);
	_store_v2i32(&this->w.l.asidx, idx);
	#endif

	/* i \in [0 .. BLK) */
	/* adjust global coordinates with local coordinate */
	this->w.l.p = leaf->p;
	this->w.l.q = leaf->q;

	this->w.l.psum = tail->psum - tail->p + leaf->p;
	this->w.l.pspos = 0;

	/* save section info */
	// this->w.l.path = *sec;
	return;
}

/**
 * @fn trace_clean_work
 */
static _force_inline
void trace_clean_work(
	struct gaba_dp_context_s *this,
	struct gaba_leaf_s const *leaf,
	struct gaba_path_intl_s *path)
{
	*path = this->w.l.path;
	return;
}

/**
 * @fn trace_forward_generate_alignment, trace_reverse_generate_alignment
 */
static _force_inline
int64_t trace_forward_generate_alignment(
	struct gaba_dp_context_s *this,
	struct gaba_leaf_s const *leaf,
	struct gaba_path_intl_s *path)
{
	trace_init_work(this, leaf, path);

	while(this->w.l.psum >= 0) {
		/* update section info */
		if(this->w.l.aidx <= 0) {
			trace_load_section_a(this);
		}
		if(this->w.l.bidx <= 0) {
			trace_load_section_b(this);
		}

		/* fragment trace */
		trace_forward_body(this);
		debug("p(%d), psum(%lld), q(%d)", this->w.l.p, this->w.l.psum, this->w.l.q);

		/* check sanity of the q-coordinate */
		if((uint32_t)this->w.l.q >= 32) { return(-1); }		/* abort */

		/* push section info to section array */
		trace_forward_push(this);
	}

	trace_clean_work(this, leaf, path);
	return(0);
}
static _force_inline
int64_t trace_reverse_generate_alignment(
	struct gaba_dp_context_s *this,
	struct gaba_leaf_s const *leaf,
	struct gaba_path_intl_s *path)
{
	trace_init_work(this, leaf, path);

	while(this->w.l.psum >= 0) {
		/* update section info */
		if(this->w.l.aidx <= 0) {
			trace_load_section_a(this);
		}
		if(this->w.l.bidx <= 0) {
			trace_load_section_b(this);
		}

		/* fragment trace */
		trace_reverse_body(this);
		debug("p(%d), psum(%lld), q(%d)", this->w.l.p, this->w.l.psum, this->w.l.q);

		/* check sanity of the q-coordinate */
		if((uint32_t)this->w.l.q >= 32) { return(-1); }

		/* push section info to section array */
		trace_reverse_push(this);
	}

	trace_clean_work(this, leaf, path);
	return(0);
}

/**
 * @fn trace_finalize_path
 */
static _force_inline
struct gaba_path_s *trace_finalize_path(
	struct gaba_path_intl_s p)
{
	uint32_t *head = p.phead, *tail = p.ptail;
	int64_t len = 32 * (p.ptail - p.phead) + p.ptofs;	/* p.phofs must be zero */

	debug("head(%p), tail(%p), hofs(%u), tofs(%u)", p.phead, p.ptail, p.phofs, p.ptofs);


	/* add terminator at the end of path array */
	tail[0] |= 0x55555555<<p.ptofs;
	tail[1] = 0x55555555;

	#if 0
	/* make the pointer 8byte aligned */
	uint32_t ofs = 0;
	if(((uint64_t)head & sizeof(uint32_t)) != 0) {
		debug("fix pointer, %p, %p", head, head - 1);

		head--;
		ofs += 32;
	}
	#endif

	/* create path object */
	struct gaba_path_s *path = (struct gaba_path_s *)(
		(uint8_t *)head - sizeof(struct gaba_path_s));
	path->len = len;
	return(path);
}

/**
 * @fn trace_cat_path
 */
static _force_inline
void trace_cat_path(
	struct gaba_dp_context_s *this,
	struct gaba_path_intl_s *dst,
	struct gaba_path_intl_s const *src)
{
	/* load pointers */
	uint32_t const *sh = src->phead;
	uint32_t const *st = src->ptail;
	uint32_t *dt = dst->ptail;

	/* concatenate heads */
	uint64_t prev_array = ((uint64_t)*dt) & ((0x01ULL<<dst->ptofs) - 1);
	uint64_t curr_array = ((uint64_t)*sh++)>>(32 - src->phofs);

	dt[0] = prev_array | (curr_array<<dst->ptofs);
	dt[1] = curr_array>>(32 - dst->ptofs);

	debug("dt[0](%x), dt[1](%x), dh[0](%x), dh[1](%x)", dt[0], dt[1], dst->phead[0], dst->phead[1]);

	dt += (((dst->ptofs + src->phofs) & 32) != 0) ? 1 : 0;
	uint64_t ofs = (dst->ptofs + src->phofs) & (32 - 1);

	/* cat */
	if(sh <= st) {
		/* load the first element */
		prev_array = *dt;

		while(sh <= st) {
			curr_array = *sh++;
			*dt++ = prev_array | (curr_array<<ofs);
			prev_array = curr_array>>(32 - ofs);

			debug("sh(%p), dt(%p), curr_array(%llx), prev_array(%llx)", sh, dt, curr_array, prev_array);
		}

		/* write back last element */
		*dt = prev_array;

		/* fix tail offset */
		dt -= (((ofs + src->ptofs) & 32) != 0) ? 1 : 0;
		ofs = (ofs + src->ptofs) & (32 - 1);
	}

	/* fix dst path object */
	dt[1] = 0;
	dt[2] = 0;
	dst->ptail = dt;
	dst->ptofs = ofs;

	debug("head(%p), tail(%p), ofs(%u)", dst->phead, dst->ptail, dst->ptofs);
	return;
}

/**
 * @fn trace_cat_section
 */
struct trace_boundary_s {
	struct gaba_path_section_s const *ptr;
	int64_t ppos;
};
static _force_inline
struct trace_boundary_s trace_cat_section(
	struct gaba_dp_context_s *this,
	struct gaba_path_intl_s *dst,
	struct gaba_path_intl_s const *src)
{
	/* load pointers */
	struct gaba_path_section_s *dh = dst->shead;
	struct gaba_path_section_s *dt = dst->stail;
	struct gaba_path_section_s const *sh = src->shead;
	struct gaba_path_section_s const *st = src->stail;

	debug("dh(%p), dt(%p), sh(%p), st(%p)", dh, dt, sh, st);

	/*
	 * load tail ppos, supposing dt[-1] is all zeros when dh == dt.
	 */
	// int64_t ppos = (dh == dt) ? 0 : dt[-1].ppos + dt[-1].plen;
	int64_t ppos = dt[-1].ppos + _plen(&dt[-1]);

	/* check if two sections can be merged */
	// int merged = 0;
	struct trace_boundary_s b = {
		.ptr = dt - (sh == st),
		.ppos = (sh == st) ? _plen(&dt[-1]) : 0
	};
	debug("init, ptr(%p), ppos(%lld)", b.ptr, b.ppos);

	/*
	 * fixme: apos != 0 and bpos != 0 suppose sections are
	 * aligned from their heads.
	 * (sh->aid == dt->aid && sh->apos == dt->apos + dt->alen)
	 * && ...
	 * might work well...
	 */
	if(dh != dt && sh != st && sh->apos != 0 && sh->bpos != 0) {
		// merged = 1;
		b.ptr--;
		b.ppos = _plen(b.ptr);

		debug("dst: id(%u, %u), pos(%u, %u), len(%u, %u)", dt[-1].bid, dt[-1].aid, dt[-1].bpos, dt[-1].apos, dt[-1].blen, dt[-1].alen);
		debug("src: id(%u, %u), pos(%u, %u), len(%u, %u)", sh[0].bid, sh[0].aid, sh[0].bpos, sh[0].apos, sh[0].blen, sh[0].alen);

		dt[-1].alen += sh->alen;
		dt[-1].blen += sh->blen;
		// dt[-1].plen += _plen(sh);
		ppos += _plen(sh); sh++;
	}
	debug("adjust, ptr(%p), ppos(%lld)", b.ptr, b.ppos);

	/* copy sections */
	while(sh < st) {
		debug("dt(%p), sh(%p), dt->plen(%u), sh->plen(%u), a(%u, %u, %u), b(%u, %u, %u), ppos(%llu)",
			dt, sh, _plen(dt), _plen(sh),
			sh->aid, sh->apos, sh->alen, sh->bid, sh->bpos, sh->blen, ppos);
		*dt = *sh;
		dt++->ppos = ppos;
		ppos += _plen(sh); sh++;
	}

	/* write back pointer */
	dst->stail = dt;
	return(b);
}

/**
 * @struct gaba_result_s
 */
struct gaba_result_s {
	struct gaba_alignment_s *aln;
	struct gaba_path_intl_s rv, fw;
};

/**
 * @fn trace_init_alignment
 */
static _force_inline
struct gaba_result_s trace_init_alignment(
	struct gaba_dp_context_s *this,
	struct gaba_joint_tail_s const *fw_tail,
	struct gaba_joint_tail_s const *rv_tail,
	struct gaba_trace_params_s const *params)
{
	/* calculate array lengths */
	uint64_t ssum = fw_tail->ssum + rv_tail->ssum;
	uint64_t psum = _roundup(MAX2(fw_tail->psum, 1), 32)
				  + _roundup(MAX2(rv_tail->psum, 1), 32);

	/* malloc trace working area */
	uint64_t sec_len = 2 * ssum;
	uint64_t path_len = _roundup(psum / 32, sizeof(uint32_t)) + 2;
	debug("psum(%lld), path_len(%llu), sec_len(%llu)", psum, path_len, sec_len);

	/* malloc pointer */
	uint64_t sec_size = sizeof(struct gaba_path_section_s) * (sec_len + 1);
	uint64_t path_size = sizeof(uint32_t) * (path_len + 4);
	uint64_t size = sizeof(struct gaba_alignment_s) + path_size + sec_size
			+ this->head_margin + this->tail_margin;

	lmm_t *lmm = (lmm_t *)params->lmm;
	struct gaba_alignment_s *aln = (struct gaba_alignment_s *)(this->head_margin
		+ (uintptr_t)((lmm == NULL) ? gaba_dp_malloc(this, size) : lmm_malloc(lmm, size)));
	// struct gaba_alignment_s *aln = (struct gaba_alignment_s *)(this->head_margin + lmm_malloc(lmm, size));

	debug("malloc trace mem(%p), lmm(%p), lim(%p)", aln, lmm, (lmm != NULL) ? lmm->lim : NULL);

	aln->lmm = (void *)lmm;
	aln->score = fw_tail->max + rv_tail->max + this->m * params->k;
	// aln->reserved1 = sec_size;
	// aln->reserved2 = path_size;
	aln->reserved3 = this->head_margin;

	/* set pointers */
	struct gaba_path_section_s *msec = (struct gaba_path_section_s *)(aln + 1);
	struct gaba_path_intl_s rv = {
		/* path arrays */
		.phead = (uint32_t *)(msec + sec_len + 1) + 2,
		.ptail = (uint32_t *)(msec + sec_len + 1) + 2,
		.phofs = 0,
		.ptofs = 0,

		/* section arrays */
		.shead = msec + 1,
		.stail = msec + 1,

		/* gap counters */
		.gic = 0,
		.gec = 0
	};
	struct gaba_path_intl_s fw = {
		/* path arrays */
		.phead = rv.phead + path_len,
		.ptail = rv.phead + path_len,
		.phofs = 0,
		.ptofs = 0,

		/* section arrays */
		.shead = rv.stail + sec_len,
		.stail = rv.stail + sec_len,

		/* gap counters */
		.gic = 0,
		.gec = 0
	};

	/* clear array */
	msec[0] = (struct gaba_path_section_s){ 0 };
	rv.phead[-1] = 0;
	rv.phead[0] = 0;
	fw.ptail[0] = 0;
	fw.ptail[1] = 0;

	return((struct gaba_result_s){
		.aln = aln,
		.rv = rv,
		.fw = fw
	});
}

/**
 * @fn trace_refine_alignment
 */
static _force_inline
struct gaba_alignment_s *trace_refine_alignment(
	struct gaba_dp_context_s *this,
	struct gaba_alignment_s *aln,
	struct gaba_path_intl_s rv,
	struct gaba_path_intl_s fw,
	struct gaba_trace_params_s const *params)
{
	struct trace_boundary_s b;
	if(params->sec != NULL) {
		/* set root section info */
		aln->rapos = params->sec[0].apos;
		aln->rbpos = params->sec[0].bpos;

		/* cat reverse and seed sections */
		uint32_t seed[2] = { 0x55555555, 0x55555555 };
		struct gaba_path_intl_s ss = {
			/* path info */
			.phead = seed,
			.ptail = seed,
			.phofs = 2 * params->k,
			.ptofs = 0,

			/* section info */
			.shead = (struct gaba_path_section_s *)params->sec,
			.stail = (struct gaba_path_section_s *)params->sec + params->slen
		};
		b = trace_cat_section(this, &rv, &ss);
		trace_cat_path(this, &rv, &ss);

		/* cat forward section */
		trace_cat_section(this, &rv, &fw);
		trace_cat_path(this, &rv, &fw);

	} else {
		/* set root section info */
		if(fw.shead == fw.stail) {
			aln->rapos = rv.stail[-1].apos + rv.stail[-1].alen;
			aln->rbpos = rv.stail[-1].bpos + rv.stail[-1].blen;
		} else {
			aln->rapos = fw.shead->apos;
			aln->rbpos = fw.shead->bpos;
		}

		/* append forward section */
		b = trace_cat_section(this, &rv, &fw);
		trace_cat_path(this, &rv, &fw);
	}
	debug("ridx(%x), rppos(%x), rapos(%x), rbpos(%x)",
		aln->rsidx, aln->rppos, aln->rapos, aln->rbpos);

	/* store root info */
	aln->rppos = b.ppos;
	aln->rsidx = b.ptr - rv.shead;

	/* store section and path info */
	aln->slen = rv.stail - rv.shead;
	aln->sec = rv.shead;
	aln->path = trace_finalize_path(rv);

	/* calc mismatche count */
	int64_t m = this->m, x = this->x, gi = this->gi, ge = this->ge;
	int64_t gic = fw.gic + rv.gic;
	int64_t gec = fw.gec + rv.gec;
	aln->xcnt =
		(m * ((aln->path->len - gec)>>1) + gi * gic + ge * gec - aln->score)
		/ (m - x);
	aln->gicnt = gic;
	aln->gecnt = gec;
	debug("plen(%lld), gic(%lld), gec(%lld), dcnt(%lld), xcnt(%lld)",
		aln->path->len, gic, gec, (aln->path->len - gec)>>1, aln->xcnt);

	return(aln);
}

/**
 * @fn gaba_dp_trace
 */
struct gaba_alignment_s *suffix(gaba_dp_trace)(
	struct gaba_dp_context_s *this,
	struct gaba_fill_s const *fw_tail,
	struct gaba_fill_s const *rv_tail,
	struct gaba_trace_params_s const *params)
{
	/* substitute tail if NULL */
	fw_tail = (fw_tail == NULL) ? _fill(&this->tail) : fw_tail;
	rv_tail = (rv_tail == NULL) ? _fill(&this->tail) : rv_tail;

	/* restore default params if NULL */
	struct gaba_trace_params_s const default_params = {
		.lmm = NULL,
		.sec = NULL,
		.slen = 0,
		.k = 0
	};
	params = (params == NULL) ? &default_params : params;

	/* search */
	struct gaba_leaf_s fw_leaf, rv_leaf;
	leaf_search(this, _tail(fw_tail), &fw_leaf);
	leaf_search(this, _tail(rv_tail), &rv_leaf);

	/* create alignment object */
	struct gaba_result_s res = trace_init_alignment(this,
		_tail(fw_tail), _tail(rv_tail), params);

	/* generate paths, may fail when path got lost out of the band */
	if(trace_forward_generate_alignment(this, &fw_leaf, &res.fw) < 0
	|| trace_reverse_generate_alignment(this, &rv_leaf, &res.rv) < 0) {
		lmm_t *lmm = (lmm_t *)params->lmm;
		lmm_free(lmm, (void *)((uint8_t *)res.aln - this->head_margin));
		return(NULL);
	}

	/* concatenate paths */
	return(trace_refine_alignment(this, res.aln, res.rv, res.fw, params));
}

/**
 * @fn gaba_dp_recombine
 */
struct gaba_alignment_s *suffix(gaba_dp_recombine)(
	struct gaba_dp_context_s *this,
	struct gaba_alignment_s *x,
	uint32_t xsid,
	struct gaba_alignment_s *y,
	uint32_t ysid)
{
	debug("recombine called, x(%p, %lld), a(%u, %u), b(%u, %u), y(%p, %lld), a(%u, %u), b(%u, %u)",
		x, x->score, x->sec[0].apos, x->sec[0].alen, x->sec[0].bpos, x->sec[0].blen,
		y, y->score, y->sec[0].apos, y->sec[0].alen, y->sec[0].bpos, y->sec[0].blen);

	gaba_dp_res_free(y);
	return(x);
}

/**
 * @fn gaba_dp_res_free
 */
void suffix(gaba_dp_res_free)(
	struct gaba_alignment_s *aln)
{
	if(aln->lmm != NULL) {
		lmm_t *lmm = (lmm_t *)aln->lmm;
		debug("free mem, ptr(%p), lmm(%p)", (void *)aln - aln->reserved3, lmm);
		lmm_free(lmm, (void *)((uint8_t *)aln - aln->reserved3));
	}
	return;
}

/**
 * @fn parse_load_uint64
 */
static inline
uint64_t parse_load_uint64(
	uint64_t const *ptr,
	int64_t pos)
{
	int64_t rem = pos & 63;
	uint64_t a = (ptr[pos>>6]>>rem) | ((ptr[(pos>>6) + 1]<<(63 - rem))<<1);
	debug("load arr(%llx)", a);
	return(a);
}

/**
 * @fn parse_dump_match_string
 */
static _force_inline
int64_t parse_dump_match_string(
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
 * @fn parse_dump_gap_string
 */
static _force_inline
int64_t parse_dump_gap_string(
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
 * @macro _parse_count_match_forward, _parse_count_gap_forward
 */
#define _parse_count_match_forward(_arr) ({ \
	tzcnt((_arr) ^ 0x5555555555555555); \
})
#define _parse_count_gap_forward(_arr) ({ \
	uint64_t _a = (_arr); \
	uint64_t mask = 0ULL - (_a & 0x01); \
	uint64_t gc = tzcnt(_a ^ mask) + (uint64_t)mask; \
	debug("arr(%llx), mask(%llx), gc(%lld)", _a, mask, gc); \
	gc; \
})

/**
 * @fn gaba_dp_print_cigar_forward
 * @brief parse path string and print cigar to file
 */
uint64_t suffix(gaba_dp_print_cigar_forward)(
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

	debug("path(%p), lim(%lld), ridx(%lld), mod(%lld)", p, lim, ridx, ridx % 64);

	while(1) {
		uint64_t rsidx = ridx;
		while(1) {
			uint64_t m = _parse_count_match_forward(parse_load_uint64(p, lim - ridx));
			uint64_t a = MIN2(m, ridx) & ~0x01;
			ridx -= a;
			volatile uint64_t c = a;
			if(c < 64) { break; }

			debug("bulk match, a(%llu)", a);
		}
		uint64_t m = (rsidx - ridx)>>1;
		if(m > 0) {
			clen += printer(fp, m, 'M');
			debug("match m(%lld)", m);
		}
		if(ridx == 0) { break; }

		uint64_t arr;
		uint64_t g = MIN2(
			_parse_count_gap_forward(arr = parse_load_uint64(p, lim - ridx)),
			ridx);
		if(g > 0) {
			clen += printer(fp, g, 'D' + ((char)(0ULL - (arr & 0x01)) & ('I' - 'D')));
			debug("gap g(%lld)", g);
		}
		if((ridx -= g) <= 1) { break; }
	}
	return(clen);
}

/**
 * @fn gaba_dp_dump_cigar_forward
 * @brief parse path string and store cigar to buffer
 */
uint64_t suffix(gaba_dp_dump_cigar_forward)(
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

	debug("path(%p), lim(%lld), ridx(%lld)", p, lim, ridx);

	while(1) {
		uint64_t rsidx = ridx;
		while(1) {
			uint64_t m = _parse_count_match_forward(parse_load_uint64(p, lim - ridx));
			uint64_t a = MIN2(m, ridx) & ~0x01;
			debug("a(%lld), ridx(%lld), ridx&~0x01(%lld)", a, ridx, ridx & ~0x01);
			ridx -= a;
			volatile uint64_t c = a;
			if(c < 64) { break; }

			debug("bulk match");
		}
		uint64_t m = (rsidx - ridx)>>1;
		if(m > 0) {
			b += parse_dump_match_string(b, m);
			debug("match m(%lld)", m);
		}
		if(ridx == 0 || b > blim) { break; }

		uint64_t arr;
		uint64_t g = MIN2(
			_parse_count_gap_forward(arr = parse_load_uint64(p, lim - ridx)),
			ridx);
		if(g > 0) {
			b += parse_dump_gap_string(b, g, 'D' + ((char)(0ULL - (arr & 0x01)) & ('I' - 'D')));
			debug("gap g(%lld)", g);
		}
		if((ridx -= g) <= 1 || b > blim) { break; }
	}
	*b = '\0';
	return(b - buf);
}

/**
 * @macro _parse_count_match_reverse, _parse_count_gap_reverse
 */
#define _parse_count_match_reverse(_arr) ({ \
	lzcnt((_arr) ^ 0x5555555555555555); \
})
#define _parse_count_gap_reverse(_arr) ({ \
	uint64_t _a = (_arr); \
	uint64_t mask = (uint64_t)(((int64_t)_a)>>63); \
	uint64_t gc = lzcnt(_a ^ mask) - ((int64_t)mask + 1); \
	debug("arr(%llx), mask(%llx), gc(%lld)", _a, mask, gc); \
	gc; \
})

/**
 * @fn gaba_dp_print_cigar_reverse
 * @brief parse path string and print cigar to file
 */
uint64_t suffix(gaba_dp_print_cigar_reverse)(
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

	debug("path(%p, %x), p(%p, %llx), idx(%lld), mod(%lld)",
		path, *path, p, *p, idx, idx % 64);

	while(1) {
		uint64_t sidx = idx;
		while(1) {
			uint64_t m = _parse_count_match_reverse(parse_load_uint64(p, idx + ofs));
			uint64_t a = MIN2(m, idx) & ~0x01;
			idx -= a;
			if(a < 64) { break; }

			debug("bulk match");
		}
		uint64_t m = (sidx - idx)>>1;
		if(m > 0) {
			clen += printer(fp, m, 'M');
			debug("match m(%lld)", m);
		}
		if(idx == 0) { break; }

		uint64_t arr;
		uint64_t g = MIN2(
			_parse_count_gap_reverse(arr = parse_load_uint64(p, idx + ofs)),
			idx);
		if(g > 0) {
			clen += printer(fp, g, 'D' + ((char)(((int64_t)arr)>>63) & ('I' - 'D')));
			debug("gap g(%lld)", g);
		}
		if((idx -= g) <= 1) { break; }
	}
	return(clen);
}

/**
 * @fn gaba_dp_dump_cigar_reverse
 * @brief parse path string and store cigar to buffer
 */
uint64_t suffix(gaba_dp_dump_cigar_reverse)(
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
			uint64_t m = _parse_count_match_reverse(parse_load_uint64(p, idx + ofs));
			uint64_t a = MIN2(m, idx) & ~0x01;
			idx -= a;
			if(a < 64) { break; }

			debug("bulk match");
		}
		uint64_t m = (sidx - idx)>>1;
		if(m > 0) {
			b += parse_dump_match_string(b, m);
			debug("match m(%lld)", m);
		}
		if(idx == 0 || b > blim) { break; }

		uint64_t arr;
		uint64_t g = MIN2(
			_parse_count_gap_reverse(arr = parse_load_uint64(p, idx + ofs)),
			idx);
		if(g > 0) {
			b += parse_dump_gap_string(b, g, 'D' + ((char)(((int64_t)arr)>>63) & ('I' - 'D')));
			debug("gap g(%lld)", g);
		}
		if((idx -= g) <= 1 || b > blim) { break; }
	}
	*b = '\0';
	return(b - buf);
}

/**
 * @fn gaba_init_restore_default_params
 */
static _force_inline
void gaba_init_restore_default_params(
	struct gaba_params_s *params)
{
	#define restore(_name, _default) { \
		params->_name = ((uint64_t)(params->_name) == 0) ? (_default) : (params->_name); \
	}
	if(params->m == 0 && params->x == 0 && params->gi == 0 && params->ge == 0) {
		params->m = 1;
		params->x = 1;
		params->gi = 1;
		params->ge = 1;
	}
	restore(xdrop, 				50);
	restore(filter_thresh,		0);
	restore(head_margin, 		0);
	restore(tail_margin, 		0);
	return;
}

/**
 * @fn gaba_init_check_score
 * @brief return non-zero if the score is not applicable
 */
static _force_inline
int gaba_init_check_score(
	struct gaba_params_s const *params)
{
	int8_t m = params->m, ge = -params->ge, gi = -params->gi;

	#if MODEL == LINEAR
		if(m - 2 * (ge + gi) > 255) { return(-1); }
		if((ge + gi) > 0) { return(-1); }
	#else
		if(m - 2 * (ge + gi) > 31) { return(-1); }
		if(ge < gi) { return(-1); }
		if((ge + gi) < -7) { return(-1); }
	#endif
	return(0);
}

/**
 * @fn gaba_init_create_score_vector
 */
static _force_inline
struct gaba_score_vec_s gaba_init_create_score_vector(
	struct gaba_params_s const *params)
{
	int8_t m = params->m, x = -params->x, ge = -params->ge, gi = -params->gi;
	int8_t sb[16] __attribute__(( aligned(16) ));
	struct gaba_score_vec_s sc __attribute__(( aligned(BW) ));

	sb[0] = x - 2 * (ge + gi);
	for(int i = 1; i < 16; i++) {
		sb[i] = m - 2 * (ge + gi);
	}
	_store_sb(sc, _load_v16i8(sb));

	#if MODEL == LINEAR
		_store_adjh(sc, 0, 0, ge + gi, ge + gi);
		_store_adjv(sc, 0, 0, ge + gi, ge + gi);
		_store_ofsh(sc, 0, 0, ge + gi, ge + gi);
		_store_ofsv(sc, 0, 0, ge + gi, ge + gi);
	#else
		_store_adjh(sc, -gi, -gi, -(ge + gi), ge + gi);
		_store_adjv(sc, -gi, -gi, -(ge + gi), ge + gi);
		_store_ofsh(sc, -gi, -gi, -(ge + gi), ge + gi);
		_store_ofsv(sc, -gi, -gi, -(ge + gi), ge + gi);
	#endif
	return(sc);
}

/**
 * @fn gaba_init_create_dir_dynamic
 */
static _force_inline
union gaba_dir_u gaba_init_create_dir_dynamic(
	struct gaba_params_s const *params)
{
	int8_t m = params->m, ge = -params->ge, gi = -params->gi;

	#if MODEL == LINEAR
		int16_t coef = -m + 2 * (ge + gi);
		int16_t ofs = 0;
	#else
		int16_t coef = -m + 2 * ge;
		int16_t ofs = gi;
	#endif

	int64_t acc = (ofs + coef * BW/2) - (ofs + coef * (BW/2 - 1));
	return((union gaba_dir_u) {
		.dynamic = {
			.acc = acc,
			.array = 0x80000000	/* (0, 0) -> (0, 1) (down) */
		}
	});
}

/**
 * @fn gaba_init_create_small_delta
 */
static _force_inline
struct gaba_small_delta_s gaba_init_create_small_delta(
	struct gaba_params_s const *params)
{
	int8_t relax = -128 / BW;

	struct gaba_small_delta_s sd;
	#if MODEL == LINEAR
		for(int i = 0; i < BW/2; i++) {
			sd.delta[i] = relax * (BW/2 - i);
			sd.delta[BW/2 + i] = relax * i;
			sd.max[i] = relax * (BW/2 - i);
			sd.max[BW/2 + i] = relax * i;
		}
	#else
		for(int i = 0; i < BW/2; i++) {
			sd.delta[i] = relax * (BW/2 - i);
			sd.delta[BW/2 + i] = relax * i;
			sd.max[i] = relax * (BW/2 - i);
			sd.max[BW/2 + i] = relax * i;
		}
	#endif

	return(sd);
}

/**
 * @fn gaba_init_create_middle_delta
 */
static _force_inline
struct gaba_middle_delta_s gaba_init_create_middle_delta(
	struct gaba_params_s const *params)
{
	int8_t m = params->m, ge = -params->ge, gi = -params->gi;
	int8_t relax = 128 / BW;

	#if MODEL == LINEAR
		int16_t coef = -m + 2 * (ge + gi) + relax;
		int16_t ofs = 0;
	#else
		int16_t coef = -m + 2*ge + relax;
		int16_t ofs = gi;
	#endif

	struct gaba_middle_delta_s md;
	for(int i = 0; i < BW/2; i++) {
		md.delta[i] = ofs + coef * (BW/2 - i);
		md.delta[BW/2 + i] = ofs + coef * i;
	}
	md.delta[BW/2] = 0;
	return(md);
}

/**
 * @fn gaba_init_create_diff_vectors
 *
 * @detail
 * dH[i, j] = S[i, j] - S[i - 1, j]
 * dV[i, j] = S[i, j] - S[i, j - 1]
 * dE[i, j] = E[i, j] - S[i, j]
 * dF[i, j] = F[i, j] - S[i, j]
 */
#if MODEL == LINEAR
static _force_inline
struct gaba_diff_vec_s gaba_init_create_diff_vectors(
	struct gaba_params_s const *params)
{
	int8_t m = params->m, ge = -params->ge, gi = -params->gi;
	int8_t drop = 0;
	int8_t raise = m - 2 * (ge + gi);

	int8_t dh[BW] __attribute__(( aligned(BW) ));
	int8_t dv[BW] __attribute__(( aligned(BW) ));

	struct gaba_diff_vec_s diff __attribute__(( aligned(BW) ));
	/**
	 * dh: dH[i, j] - gh
	 * dv: dV[i, j] - gv
	 */
	/* calc dh and dv */
	for(int i = 0; i < BW/2; i++) {
		dh[i] = drop;
		dh[BW/2 + i] = raise;
		dv[i] = raise;
		dv[BW/2 + i] = drop;
	}
 	dh[BW/2] = raise;
 	dv[BW/2] = raise;

	_store(&diff.dh, _load(dh));
	_store(&diff.dv, _load(dv));
	return(diff);
}
#else
static _force_inline
struct gaba_diff_vec_s gaba_init_create_diff_vectors(
	struct gaba_params_s const *params)
{
	int8_t m = params->m, ge = -params->ge, gi = -params->gi;

	int8_t ofs_dh = -(ge + gi);
	int8_t ofs_dv = -(ge + gi);
	int8_t ofs_de = -gi;
	int8_t ofs_df = -gi;
	
	int8_t drop_dh = ge + ofs_dh;
	int8_t raise_dh = m - ge + ofs_dh;
	int8_t drop_dv = ge + ofs_dv;
	int8_t raise_dv = m - ge + ofs_dv;
	int8_t drop_de = gi + ofs_de;
	int8_t raise_de = ofs_de;
	int8_t drop_df = gi + ofs_df;
	int8_t raise_df = ofs_df;

	int8_t dh[BW] __attribute__(( aligned(BW) ));
	int8_t dv[BW] __attribute__(( aligned(BW) ));
	int8_t de[BW] __attribute__(( aligned(BW) ));
	int8_t df[BW] __attribute__(( aligned(BW) ));

	struct gaba_diff_vec_s diff __attribute__(( aligned(BW) ));
	/**
	 * dh: dH[i, j] - ge
	 * dv: dV[i, j] - gev
	 * de: dE[i, j] + gi + dV[i, j] - gev
	 * df: dF[i, j] + gi + dH[i, j] - ge
	 */
	/* calc dh and dv */
	for(int i = 0; i < BW/2; i++) {
		dh[i] = drop_dh;
		dh[BW/2 + i] = raise_dh;
		dv[i] = raise_dv;
		dv[BW/2 + i] = drop_dv;
	}
 	dh[BW/2] = raise_dh - gi;
 	dv[BW/2] = raise_dv - gi;

	/* calc de and df */
 	for(int i = 0; i < BW/2; i++) {
 		de[i] = raise_de;
 		de[BW/2 + i] = drop_de;
 		df[i] = drop_df;
 		df[BW/2 + i] = raise_df;
 	}
 	de[BW/2] = drop_de;
 	df[BW/2] = drop_df;

 	vec_t _dh = _load(dh);
 	vec_t _dv = _load(dv);
 	vec_t _de = _load(de);
 	vec_t _df = _load(df);

 	_print(_dh);
 	_print(_dv);
 	_print(_de);
 	_print(_df);

	_dh = _shl(_dh, 3);
	_dv = _shl(_dv, 3);
	_store(&diff.dh, _add(_dh, _de));
	_store(&diff.dv, _add(_dv, _df));
	_print(_add(_dh, _de));
	_print(_add(_dv, _df));

	return(diff);
}
#endif

/**
 * @fn gaba_init_create_char_vector
 */
static _force_inline
struct gaba_char_vec_s gaba_init_create_char_vector(
	void)
{
	struct gaba_char_vec_s ch;

	for(int i = 0; i < BW; i++) {
		ch.w[i] = 0;
	}
	return(ch);
}

/**
 * @fn gaba_init
 */
gaba_t *suffix(gaba_init)(
	struct gaba_params_s const *params)
{
	if(params == NULL) {
		debug("params must not be NULL");
		return(NULL);
	}

	/* copy params to local stack */
	struct gaba_params_s params_intl = *params;

	/* restore defaults */
	gaba_init_restore_default_params(&params_intl);

	/* check the scores are applicable */
	if(gaba_init_check_score(&params_intl) != 0) {
		return(NULL);
	}

	/* malloc gaba_context_s */
	struct gaba_context_s *ctx = (struct gaba_context_s *)gaba_aligned_malloc(
		sizeof(struct gaba_context_s),
		MEM_ALIGN_SIZE);
	if(ctx == NULL) {
		return(NULL);
	}

	/* fill context */
	*ctx = (struct gaba_context_s) {
		/* template */
		.k = (struct gaba_dp_context_s) {
			/* memory management */
			.mem = { 0 },
			.curr_mem = NULL,
			.stack_top = NULL,						/* stored on init */
			.stack_end = NULL,						/* stored on init */

			/* score vectors */
			.scv = gaba_init_create_score_vector(&params_intl),
			.m = params_intl.m,
			.x = -params_intl.x,
			.gi = (MODEL == LINEAR)
				? 0
				: -params_intl.gi,
			.ge = (MODEL == LINEAR)
				? -(params_intl.gi + params_intl.ge)
				: -params_intl.ge,
			.tx = params_intl.xdrop,
			.tf = params_intl.filter_thresh,

			/* input and output options */
			.head_margin = _roundup(params_intl.head_margin, MEM_ALIGN_SIZE),
			.tail_margin = _roundup(params_intl.tail_margin, MEM_ALIGN_SIZE),

			/* phantom block at root */
			.blk = (struct gaba_phantom_block_s) {
				/* direction array */
				.dir = gaba_init_create_dir_dynamic(&params_intl),
				
				/* offset, diffs and deltas */
				.offset = 0,
				.diff = gaba_init_create_diff_vectors(&params_intl),
				.sd = gaba_init_create_small_delta(&params_intl),
				.md = &ctx->md,

				/* char vectors */
				.ch = gaba_init_create_char_vector(),

				/* indices */
				.aridx = 0,
				.bridx = 0
			},

			/* tail of root section */
			.tail = (struct gaba_joint_tail_s) {
				/* coordinates */
				.psum = PSUM_BASE - BW,
				.p = 0,
				.ssum = 0,

				/* max */
				.max = 0,

				/* status */
				.stat = CONT,

				/* internals */
				.rem_len = 0,
				.tail = NULL,
				.apos = 0,
				.bpos = 0,
				.alen = 0,
				.blen = 0,
				.aid = 0xfffc,
				.bid = 0xfffd
			},
		},

		/* default middle delta vector */
		.md = gaba_init_create_middle_delta(&params_intl)
	};

	return((gaba_t *)ctx);
}

/**
 * @fn gaba_clean
 */
void suffix(gaba_clean)(
	struct gaba_context_s *ctx)
{
	gaba_aligned_free(ctx);
	return;
}

/**
 * @fn gaba_dp_init
 */
struct gaba_dp_context_s *suffix(gaba_dp_init)(
	struct gaba_context_s const *ctx,
	uint8_t const *alim,
	uint8_t const *blim)
{
	/* malloc stack memory */
	struct gaba_dp_context_s *this = (struct gaba_dp_context_s *)gaba_aligned_malloc(
		MEM_INIT_SIZE, MEM_ALIGN_SIZE);
	if(this == NULL) {
		debug("failed to malloc memory");
		return(NULL);
	}

	/* init stack pointers */
	this->stack_top = (uint8_t *)(this + 1);
	this->stack_end = (uint8_t *)this + MEM_INIT_SIZE - MEM_MARGIN_SIZE;

	/* init seq lims */
	this->w.r.alim = alim;
	this->w.r.blim = blim;

	/* copy template */
	_memcpy_blk_aa(
		(uint8_t *)this + GABA_DP_CONTEXT_LOAD_OFFSET,
		(uint8_t *)&ctx->k + GABA_DP_CONTEXT_LOAD_OFFSET,
		GABA_DP_CONTEXT_LOAD_SIZE);

	/* init mem object */
	this->curr_mem = &this->mem;
	this->mem = (struct gaba_mem_block_s){
		.next = NULL,
		.prev = NULL,
		.size = MEM_INIT_SIZE
	};
	return(this);
}

/**
 * @fn gaba_dp_add_stack
 */
static _force_inline
int32_t gaba_dp_add_stack(
	struct gaba_dp_context_s *this,
	uint64_t size)
{
	if(this->curr_mem->next == NULL) {
		/* add new block */
		uint64_t next_size = this->curr_mem->size * 2;
		struct gaba_mem_block_s *mem = this->curr_mem->next =
			(struct gaba_mem_block_s *)gaba_aligned_malloc(
				next_size, MEM_ALIGN_SIZE);
		if(mem == NULL) { return(GABA_ERROR_OUT_OF_MEM); }

		mem->next = NULL;
		mem->prev = this->curr_mem;
		mem->size = next_size;
	}

	/* follow the forward link */
	this->curr_mem = this->curr_mem->next;

	/* init stack pointers */
	this->stack_top = (uint8_t *)(this->curr_mem + 1);
	this->stack_end = (uint8_t *)this->curr_mem + this->curr_mem->size - MEM_MARGIN_SIZE;
	return(GABA_SUCCESS);
}

/**
 * @fn gaba_dp_flush
 */
void suffix(gaba_dp_flush)(
	struct gaba_dp_context_s *this,
	uint8_t const *alim,
	uint8_t const *blim)
{
	/* init seq lims */
	this->w.r.alim = alim;
	this->w.r.blim = blim;

	this->curr_mem = &this->mem;
	this->stack_top = (uint8_t *)(this + 1);
	this->stack_end = (uint8_t *)this + this->mem.size - MEM_MARGIN_SIZE;
	return;
}

/**
 * @fn gaba_dp_save_stack
 */
struct gaba_stack_s const *suffix(gaba_dp_save_stack)(
	struct gaba_dp_context_s *this)
{
	struct gaba_mem_block_s *mem = this->curr_mem;
	uint8_t *stack_top = this->stack_top;
	uint8_t *stack_end = this->stack_end;

	/* save */
	struct gaba_stack_s *stack = (struct gaba_stack_s *)gaba_dp_malloc(
		this, sizeof(struct gaba_stack_s));
	stack->mem = mem;
	stack->stack_top = stack_top;
	stack->stack_end = stack_end;
	debug("save stack(%p, %p, %p)", stack->mem, stack->stack_top, stack->stack_end);
	return(stack);
}

/**
 * @fn gaba_dp_flush_stack
 */
void suffix(gaba_dp_flush_stack)(
	struct gaba_dp_context_s *this,
	gaba_stack_t const *stack)
{
	if(stack == NULL) {
		return;
	}

	this->curr_mem = stack->mem;
	this->stack_top = stack->stack_top;
	this->stack_end = stack->stack_end;
	debug("restore stack(%p, %p, %p)", stack->mem, stack->stack_top, stack->stack_end);
	return;
}

/**
 * @fn gaba_dp_malloc
 */
static _force_inline
void *gaba_dp_malloc(
	struct gaba_dp_context_s *this,
	uint64_t size)
{
	/* roundup */
	size = _roundup(size, MEM_ALIGN_SIZE);

	/* malloc */
	debug("this(%p), stack_top(%p), size(%llu)", this, this->stack_top, size);
	if((uint64_t)(this->stack_end - this->stack_top) < size) {
		if(gaba_dp_add_stack(this, size) != GABA_SUCCESS) {
			return(NULL);
		}
		debug("stack_top(%p)", this->stack_top);
	}
	this->stack_top += size;
	return((void *)(this->stack_top - size));
}

/**
 * @fn gaba_dp_clean
 */
void suffix(gaba_dp_clean)(
	struct gaba_dp_context_s *this)
{
	if(this == NULL) {
		return;
	}

	struct gaba_mem_block_s *m = this->mem.next;
	while(m != NULL) {
		struct gaba_mem_block_s *mnext = m->next;
		free(m); m = mnext;
	}
	gaba_aligned_free(this);
	return;
}

/* unittests */
#if UNITTEST == 1

#include <string.h>
#include <sys/types.h>
#include <unistd.h>

/**
 * @fn unittest_build_context
 * @brief build context for unittest
 */
#if MODEL == LINEAR
static struct gaba_params_s const *unittest_default_params = GABA_PARAMS(
	.xdrop = 100,
	GABA_SCORE_SIMPLE(2, 3, 0, 6));
#else
static struct gaba_params_s const *unittest_default_params = GABA_PARAMS(
	.xdrop = 100,
	GABA_SCORE_SIMPLE(2, 3, 5, 1));
#endif
static
void *unittest_build_context(void *params)
{
	/* build context */
	gaba_t *ctx = gaba_init(unittest_default_params);
	return((void *)ctx);
}

/**
 * @fn unittest_clean_context
 */
static
void unittest_clean_context(void *ctx)
{
	gaba_clean((struct gaba_context_s *)ctx);
	return;
}

/**
 * @struct unittest_seqs_s
 * @brief sequence container
 */
struct unittest_seqs_s {
	char const *a;
	char const *b;
};
#define unittest_seq_pair(_a, _b) ( \
	(void *)&((struct unittest_seqs_s const){ \
		.a = _a "GGGGGGGGGGGGGGGGGGGG", \
		.b = _b "CCCCCCCCCCCCCCCCCCCC" \
	}) \
)

/**
 * @struct unittest_sections_s
 * @brief section container
 */
struct unittest_sections_s {
	uint8_t const *a, *b;
	uint8_t const *alim, *blim;
	uint64_t alen, blen;
	struct gaba_section_s afsec, aftail, bfsec, bftail;
	struct gaba_section_s arsec, artail, brsec, brtail;
	struct gaba_section_s assec, bssec;
};

/**
 * @fn unittest_encode_base
 * @brief mapping IUPAC amb. to 2bit / 4bit encoding
 */
static _force_inline
uint8_t unittest_encode_base(
	char c)
{
	/* convert to upper case and subtract offset by 0x40 */
	#define _b(x)	( (x) & 0x1f )

	/* conversion tables */
	enum bases { A = 0x01, C = 0x02, G = 0x04, T = 0x08 };
	static uint8_t const table[] = {
		[_b('A')] = A,
		[_b('C')] = C,
		[_b('G')] = G,
		[_b('T')] = T,
		[_b('U')] = T,
		[_b('R')] = A | G,
		[_b('Y')] = C | T,
		[_b('S')] = G | C,
		[_b('W')] = A | T,
		[_b('K')] = G | T,
		[_b('M')] = A | C,
		[_b('B')] = C | G | T,
		[_b('D')] = A | G | T,
		[_b('H')] = A | C | T,
		[_b('V')] = A | C | G,
		[_b('N')] = 0,		/* treat 'N' as a gap */
		[_b('_')] = 0		/* sentinel */
	};
	return(table[_b((uint8_t)c)]);

	#undef _b
}

/**
 * @fn unittest_build_seqs
 * @brief build seq_pair and sections
 */
static
void *unittest_build_seqs(void *params)
{
	struct unittest_seqs_s const *seqs = (struct unittest_seqs_s const *)params;
	char const *a = seqs->a;
	char const *b = seqs->b;

	uint64_t atot = strlen(a);
	uint64_t btot = strlen(b);
	uint32_t alen = atot - 20;
	uint32_t blen = btot - 20;

	uint64_t margin = 64;
	struct unittest_sections_s *sec = malloc(
		sizeof(struct unittest_sections_s) + (atot + 1) + (btot + 1) + margin);

	/* copy sequences */
	uint8_t *ca = (uint8_t *)(sec + 1);
	uint8_t *cb = ca + atot + 1;
	uint8_t *cm = cb + btot + 1;

	for(uint64_t i = 0; i < atot; i++) {
		ca[i] = unittest_encode_base(a[i]);
	}
	for(uint64_t i = 0; i < btot; i++) {
		cb[i] = unittest_encode_base(b[i]);
	}
	ca[atot] = cb[btot] = '\0';
	memset(cm, 0, margin);

	/* build context */
	uint8_t const *const alim = (void const *)0x800000000000;
	uint8_t const *const blim = (void const *)0x800000000000;
	*sec = (struct unittest_sections_s){
		.a = ca,
		.b = cb,
		.alim = alim,
		.blim = blim,
		.alen = atot,
		.blen = btot,
		
		/* forward */
		.afsec = gaba_build_section(0, ca, alen),
		.aftail = gaba_build_section(2, ca + alen, 20),
		.bfsec = gaba_build_section(4, cb, blen),
		.bftail = gaba_build_section(6, cb + blen, 20),

		/* reverse */
		.arsec = gaba_build_section(1, _rev(ca + alen - 1, alim), alen),
		.artail = gaba_build_section(3, _rev(ca + atot - 1, alim), 20),
		.brsec = gaba_build_section(5, _rev(cb + blen - 1, blim), blen),
		.brtail = gaba_build_section(7, _rev(cb + btot - 1, blim), 20),

		/* seed section */
		.assec = gaba_build_section(100, ca, 14),
		.bssec = gaba_build_section(102, cb, 14)
	};
	return((void *)sec);
}

/**
 * @fn unittest_clean_seqs
 */
static
void unittest_clean_seqs(void *ctx)
{
	free(ctx);
}

/**
 * @macro with_seq_pair
 * @brief constructs .param, .init, and .clean parameters
 */
#define with_seq_pair(_a, _b) \
	.params = unittest_seq_pair(_a, _b), \
	.init = unittest_build_seqs, \
	.clean = unittest_clean_seqs

/**
 * misc macros and functions for assertion
 */
#define check_tail(_t, _max, _p, _psum, _ssum) ( \
	   (_t) != NULL \
	&& (_t)->max == (_max) \
	&& (_t)->p == (_p) \
	&& (_t)->psum == (_psum) \
	&& (_t)->ssum == (_ssum) \
)
#define print_tail(_t) \
	"tail(%p), max(%lld), p(%d), psum(%lld), ssum(%u)", \
	(_t), (_t)->max, (_t)->p, (_t)->psum, (_t)->ssum
#define check_result(_r, _score, _xcnt, _plen, _slen, _rsidx, _rppos, _rapos, _rbpos) ( \
	   (_r) != NULL \
	&& (_r)->sec != NULL \
	&& (_r)->path != NULL \
	&& (_r)->path->len == (_plen) \
	&& (_r)->slen == (_slen) \
	&& (_r)->score == (_score) \
	&& (_r)->xcnt == (_xcnt) \
	&& (_r)->rsidx == (_rsidx) \
	&& (_r)->rppos == (_rppos) \
	&& (_r)->rapos == (_rapos) \
	&& (_r)->rbpos == (_rbpos) \
)
#define print_result(_r) \
	"res(%p), score(%lld), xcnt(%lld), plen(%u), slen(%u), rsid(%u), rppos(%u), rapos(%u), rbpos(%u)", \
	(_r), (_r)->score, (_r)->xcnt, (_r)->path->len, (_r)->slen, (_r)->rsidx, (_r)->rppos, (_r )->rapos, (_r)->rbpos

static
int check_path(
	struct gaba_alignment_s const *aln,
	char const *str)
{
	int64_t plen = aln->path->len, slen = strlen(str);
	uint32_t const *p = aln->path->array;
	debug("%s", str);

	/* first check length */
	if(plen != slen) {
		debug("plen(%lld), slen(%lld)", plen, slen);
		return(0);
	}

	/* next compare encoded string (bit string) */
	while(plen > 0) {
		uint32_t array = 0;
		for(int64_t i = 0; i < 32; i++) {
			if(plen-- == 0) {
				array = (array>>(32 - i)) | ((uint64_t)0x55555555<<i);
				break;
			}
			array = (array>>1) | ((*str++ == 'D') ? 0x80000000 : 0);
			debug("%c, %x", str[-1], array);
		}
		debug("path(%x), array(%x)", *p, array);
		if(*p++ != array) {
			return(0);
		}
	}
	return(1);
}

static
int check_cigar(
	struct gaba_alignment_s const *aln,
	char const *cigar)
{
	char buf[1024];

	debug("path(%x), len(%lld)", aln->path->array[0], aln->path->len);

	uint64_t l = gaba_dp_dump_cigar_forward(
		buf, 1024, aln->path->array, 0, aln->path->len);

	debug("cigar(%s)", buf);

	/* first check length */
	if(strlen(cigar) != l) { return(0); }

	/* next compare cigar string */
	return((strcmp(buf, cigar) == 0) ? 1 : 0);
}

#define decode_path(_r) ({ \
	uint64_t plen = (_r)->path->len, cnt = 0; \
	uint32_t const *path = (_r)->path->array; \
	uint32_t path_array = *path; \
	char *ptr = alloca(plen); \
	char *p = ptr; \
 	while(plen-- > 0) { \
		*p++ = (path_array & 0x01) ? 'D' : 'R'; \
		path_array >>= 1; \
		if(++cnt == 32) { \
			path_array = *++path; \
			cnt = 0; \
		} \
	} \
	*p = '\0'; \
	ptr; \
})
#define print_path(_r)			"%s", decode_path(_r)
#define check_section(_s, _a, _apos, _alen, _b, _bpos, _blen, _ppos, _pl) ( \
	   (_s).aid == (_a).id \
	&& (_s).apos == (_apos) \
	&& (_s).alen == (_alen) \
	&& (_s).bid == (_b).id \
	&& (_s).bpos == (_bpos) \
	&& (_s).blen == (_blen) \
	&& (_s).ppos == (_ppos) \
	&& _plen(&(_s)) == (_pl) \
)
#define print_section(_s) \
	"a(%u), apos(%u), alen(%u), b(%u), bpos(%u), blen(%u), ppos(%u), plen(%u)", \
	(_s).aid, (_s).apos, (_s).alen, \
	(_s).bid, (_s).bpos, (_s).blen, \
	(_s).ppos, _plen(&(_s))

/* global configuration of the tests */
unittest_config(
	.name = "gaba",
	.init = unittest_build_context,
	.clean = unittest_clean_context
);

/**
 * check if gaba_init returns a valid pointer to a context
 */
unittest()
{
	struct gaba_context_s const *c = (struct gaba_context_s const *)gctx;
	assert(c != NULL, "%p", c);
}

/**
 * check if unittest_build_seqs returns a valid seq_pair and sections
 */
unittest(with_seq_pair("A", "A"))
{
	struct unittest_sections_s const *s = (struct unittest_sections_s const *)ctx;

	/* check pointer */
	assert(s != NULL, "%p", s);

	/* check sequences */
	assert(strncmp((char const *)s->a,
		"\x01\x04\x04\x04\x04\x04\x04\x04\x04\x04\x04\x04\x04\x04\x04\x04\x04\x04\x04\x04\x04\0",
		22) == 0, "%s", (char const *)s->a);
	assert(strncmp((char const *)s->b,
		"\x01\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\0",
		22) == 0, "%s", (char const *)s->b);
	assert(s->alen == 21, "%llu", s->alen);
	assert(s->blen == 21, "%llu", s->blen);

	/* check forward sections */
	assert(s->afsec.id == 0, "%d", s->afsec.id);
	assert((uintptr_t)s->afsec.base == (uintptr_t)s->a + 0, "%p", s->afsec.base);
	assert(s->afsec.len == 1, "%u", s->afsec.len);

	assert(s->aftail.id == 2, "%d", s->aftail.id);
	assert((uintptr_t)s->aftail.base == (uintptr_t)s->a + 1, "%p", s->aftail.base);
	assert(s->aftail.len == 20, "%u", s->aftail.len);

	assert(s->bfsec.id == 4, "%d", s->bfsec.id);
	assert((uintptr_t)s->bfsec.base == (uintptr_t)s->b + 0, "%p", s->bfsec.base);
	assert(s->bfsec.len == 1, "%u", s->bfsec.len);

	assert(s->bftail.id == 6, "%d", s->bftail.id);
	assert((uintptr_t)s->bftail.base == (uintptr_t)s->b + 1, "%p", s->bftail.base);
	assert(s->bftail.len == 20, "%u", s->bftail.len);

	/* check reverse sections */
	assert(s->arsec.id == 1, "%d", s->arsec.id);
	assert((uintptr_t)s->arsec.base == (uintptr_t)0x1000000000000 - (uintptr_t)s->a - 1, "%p", s->arsec.base);
	assert(s->arsec.len == 1, "%u", s->arsec.len);

	assert(s->artail.id == 3, "%d", s->artail.id);
	assert((uintptr_t)s->artail.base == (uintptr_t)0x1000000000000 - (uintptr_t)s->a - 21, "%p", s->artail.base);
	assert(s->artail.len == 20, "%u", s->artail.len);

	assert(s->brsec.id == 5, "%d", s->brsec.id);
	assert((uintptr_t)s->brsec.base == (uintptr_t)0x1000000000000 - (uintptr_t)s->b - 1, "%p", s->brsec.base);
	assert(s->brsec.len == 1, "%u", s->brsec.len);

	assert(s->brtail.id == 7, "%d", s->brtail.id);
	assert((uintptr_t)s->brtail.base == (uintptr_t)0x1000000000000 - (uintptr_t)s->b - 21, "%p", s->brtail.base);
	assert(s->brtail.len == 20, "%u", s->brtail.len);
}

/**
 * check if gaba_dp_init returns a vaild pointer to a dp context
 */
#define omajinai() \
	struct gaba_context_s const *c = (struct gaba_context_s const *)gctx; \
	struct unittest_sections_s const *s = (struct unittest_sections_s const *)ctx; \
	struct gaba_dp_context_s *d = gaba_dp_init(c, s->alim, s->blim);

unittest(with_seq_pair("A", "A"))
{
	omajinai();

	assert(d != NULL, "%p", d);
	gaba_dp_clean(d);
}

/**
 * check if gaba_dp_fill_root and gaba_dp_fill returns a correct score
 */
unittest(with_seq_pair("A", "A"))
{
	omajinai();

	/* fill root section */
	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 0, 0, -29, 1), print_tail(f));

	/* fill again */
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 0, 0, -27, 2), print_tail(f));

	/* fill tail section */
	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 4, 13, 13, 3), print_tail(f));

	/* fill tail section again */
	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);
	#if MODEL == LINEAR
		assert(f->status == 0x1ff, "%x", f->status);
		assert(check_tail(f, 4, 40, 53, 4), print_tail(f));
	#else
		assert(f->status == 0x10f, "%x", f->status);
		assert(check_tail(f, 4, 31, 44, 4), print_tail(f));
	#endif

	gaba_dp_clean(d);
}

/* reverse fetch */
unittest(with_seq_pair("A", "A"))
{
	omajinai();

	/* fill root section */
	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->arsec, 0, &s->brsec, 0);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 0, 0, -29, 1), print_tail(f));

	/* fill again */
	f = gaba_dp_fill(d, f, &s->arsec, &s->brsec);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 0, 0, -27, 2), print_tail(f));

	/* fill tail section */
	f = gaba_dp_fill(d, f, &s->artail, &s->brtail);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 4, 13, 13, 3), print_tail(f));

	/* fill tail section again */
	f = gaba_dp_fill(d, f, &s->artail, &s->brtail);
	#if MODEL == LINEAR
		assert(f->status == 0x1ff, "%x", f->status);
		assert(check_tail(f, 4, 40, 53, 4), print_tail(f));
	#else
		assert(f->status == 0x10f, "%x", f->status);
		assert(check_tail(f, 4, 31, 44, 4), print_tail(f));
	#endif

	gaba_dp_clean(d);
}

/* with longer sequences */
unittest(with_seq_pair("ACGTACGTACGT", "ACGTACGTACGT"))
{
	omajinai();

	/* fill root section */
	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 0, 0, -7, 1), print_tail(f));

	/* fill again */
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 16, 17, 17, 2), print_tail(f));

	/* fill tail section */
	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 48, 40, 57, 3), print_tail(f));

	/* fill tail section again */
	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);
	#if MODEL == LINEAR
		assert(f->status == 0x1ff, "%x", f->status);
		assert(check_tail(f, 48, 40, 97, 4), print_tail(f));
	#else
		assert(f->status == 0x10f, "%x", f->status);
		assert(check_tail(f, 48, 31, 88, 4), print_tail(f));
	#endif

	gaba_dp_clean(d);
}

/* reverse fetch */
unittest(with_seq_pair("ACGTACGTACGT", "ACGTACGTACGT"))
{
	omajinai();

	/* fill root section */
	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->arsec, 0, &s->brsec, 0);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 0, 0, -7, 1), print_tail(f));

	/* fill again */
	f = gaba_dp_fill(d, f, &s->arsec, &s->brsec);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 16, 17, 17, 2), print_tail(f));

	/* fill tail section */
	f = gaba_dp_fill(d, f, &s->artail, &s->brtail);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 48, 40, 57, 3), print_tail(f));

	/* fill tail section again */
	f = gaba_dp_fill(d, f, &s->artail, &s->brtail);
	#if MODEL == LINEAR
		assert(f->status == 0x1ff, "%x", f->status);
		assert(check_tail(f, 48, 40, 97, 4), print_tail(f));
	#else
		assert(f->status == 0x10f, "%x", f->status);
		assert(check_tail(f, 48, 31, 88, 4), print_tail(f));
	#endif

	gaba_dp_clean(d);
}

/* sequences with different lengths (consumed as mismatches) */
unittest(with_seq_pair("GAAAAAAAA", "AAAAAAAA"))
{
	omajinai();

	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	assert(f->status == 0x01ff, "%x", f->status);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x01f0, "%x", f->status);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bftail);
	assert(f->status == 0x010f, "%x", f->status);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);
	assert(f->status == 0x01f0, "%x", f->status);
	assert(check_tail(f, 22, 37, 42, 4), print_tail(f));

	gaba_dp_clean(d);
}

/* another pair with different lengths */
unittest(with_seq_pair("TTTTTTTT", "CTTTTTTTT"))
{
	omajinai();

	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	assert(f->status == 0x010f, "%x", f->status);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x01f0, "%x", f->status);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x010f, "%x", f->status);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bfsec);
	assert(f->status == 0x01f0, "%x", f->status);

	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);
	#if MODEL == LINEAR
		assert(f->status == 0x010f, "%x", f->status);
		assert(check_tail(f, 22, 36, 42, 5), print_tail(f));
	#else
		assert(f->status == 0x010f, "%x", f->status);
		assert(check_tail(f, 22, 35, 41, 5), print_tail(f));
	#endif

	gaba_dp_clean(d);
}

/* with deletions */
unittest(with_seq_pair("GACGTACGT", "ACGTACGT"))
{
	omajinai();

	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	assert(f->status == 0x01ff, "%x", f->status);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x01f0, "%x", f->status);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bftail);
	assert(f->status == 0x010f, "%x", f->status);

	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);
	#if MODEL == LINEAR
		assert(f->status == 0x01f0, "%x", f->status);
		assert(check_tail(f, 20, 37, 42, 4), print_tail(f));
	#else
		assert(f->status == 0x01ff, "%x", f->status);
		assert(check_tail(f, 20, 38, 43, 4), print_tail(f));
	#endif

	gaba_dp_clean(d);
}

/* with insertions */
unittest(with_seq_pair("ACGTACGT", "GACGTACGT"))
{
	omajinai();

	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	assert(f->status == 0x010f, "%x", f->status);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x01f0, "%x", f->status);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x010f, "%x", f->status);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bfsec);
	assert(f->status == 0x01f0, "%x", f->status);


	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);
	#if MODEL == LINEAR
		assert(f->status == 0x010f, "%x", f->status);
		assert(check_tail(f, 20, 35, 41, 5), print_tail(f));
	#else
		assert(f->status == 0x010f, "%x", f->status);
		assert(check_tail(f, 20, 36, 42, 5), print_tail(f));
	#endif

	gaba_dp_clean(d);
}

/* print_cigar test */
static
int ut_printer(
	void *pbuf,
	int64_t len,
	char c)
{
	char *buf = *((char **)pbuf);
	len = sprintf(buf, "%" PRId64 "%c", len, c);
	*((char **)pbuf) += len;
	return(len);
}

unittest()
{
	char *buf = (char *)malloc(16384);
	char *p = buf;

	#define _arr(...)		( (uint32_t const []){ 0, 0, __VA_ARGS__, 0, 0 } + 2 )
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55555555), 0, 32);
	assert(strcmp(buf, "16M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55555555, 0x55555555), 0, 64);
	assert(strcmp(buf, "32M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55555555, 0x55555555, 0x55555555, 0x55555555), 0, 128);
	assert(strcmp(buf, "64M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55550000, 0x55555555, 0x55555555, 0x55555555), 16, 112);
	assert(strcmp(buf, "56M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55555000, 0x55555555, 0x55555555, 0x55555555), 12, 116);
	assert(strcmp(buf, "58M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55), 0, 8);
	assert(strcmp(buf, "4M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55555000, 0x55555555, 0x55555555, 0x55), 12, 92);
	assert(strcmp(buf, "46M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0x55550555), 0, 32);
	assert(strcmp(buf, "6M4D8M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0x5555f555), 0, 32);
	assert(strcmp(buf, "6M4I8M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0xaaaa0555), 0, 33);
	assert(strcmp(buf, "6M5D8M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0xaaabf555), 0, 33);
	assert(strcmp(buf, "6M5I8M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0xaaabf555, 0xaaaa0556), 0, 65);
	assert(strcmp(buf, "6M5I8M1I5M5D8M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0xaaabf555, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "6M5I8M1I5M5D8M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_forward(ut_printer, (void *)&p, _arr(0xaaabf554, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "2D5M5I8M1I5M5D8M") == 0, "%s", buf);

	#undef _arr
	free(buf);
}

unittest()
{
	char *buf = (char *)malloc(16384);
	char *p = buf;

	#define _arr(...)		( (uint32_t const []){ 0, 0, __VA_ARGS__, 0, 0 } + 2 )
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55555555), 0, 32);
	assert(strcmp(buf, "16M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55555555, 0x55555555), 0, 64);
	assert(strcmp(buf, "32M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55555555, 0x55555555, 0x55555555, 0x55555555), 0, 128);
	assert(strcmp(buf, "64M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55550000, 0x55555555, 0x55555555, 0x55555555), 16, 112);
	assert(strcmp(buf, "56M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55555000, 0x55555555, 0x55555555, 0x55555555), 12, 116);
	assert(strcmp(buf, "58M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55), 0, 8);
	assert(strcmp(buf, "4M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55555000, 0x55555555, 0x55555555, 0x55), 12, 92);
	assert(strcmp(buf, "46M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x55550555), 0, 32);
	assert(strcmp(buf, "8M4D6M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0x5555f555), 0, 32);
	assert(strcmp(buf, "8M4I6M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0xaaaa0555), 0, 33);
	assert(strcmp(buf, "8M5D6M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0xaaabf555), 0, 33);
	assert(strcmp(buf, "8M5I6M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0xaaabf555, 0xaaaa0556), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I6M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0xaaabf555, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I6M") == 0, "%s", buf);

	p = buf;
	gaba_dp_print_cigar_reverse(ut_printer, (void *)&p, _arr(0xaaabf554, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I5M2D") == 0, "%s", buf);

	#undef _arr
	free(buf);
}

/* dump_cigar test */
unittest()
{
	uint64_t const len = 16384;
	char *buf = (char *)malloc(len);

	#define _arr(...)		( (uint32_t const []){ 0, 0, __VA_ARGS__, 0, 0 } + 2 )
	gaba_dp_dump_cigar_forward(buf, len, _arr(0x55555555), 0, 32);
	assert(strcmp(buf, "16M") == 0, "%s", buf);

	gaba_dp_dump_cigar_forward(buf, len, _arr(0x55555555, 0x55555555), 0, 64);
	assert(strcmp(buf, "32M") == 0, "%s", buf);

	gaba_dp_dump_cigar_forward(buf, len, _arr(0x55555555, 0x55555555, 0x55555555, 0x55555555), 0, 128);
	assert(strcmp(buf, "64M") == 0, "%s", buf);

	gaba_dp_dump_cigar_forward(buf, len, _arr(0x55550000, 0x55555555, 0x55555555, 0x55555555), 16, 112);
	assert(strcmp(buf, "56M") == 0, "%s", buf);

	gaba_dp_dump_cigar_forward(buf, len, _arr(0x55555000, 0x55555555, 0x55555555, 0x55555555), 12, 116);
	assert(strcmp(buf, "58M") == 0, "%s", buf);

	gaba_dp_dump_cigar_forward(buf, len, _arr(0x55), 0, 8);
	assert(strcmp(buf, "4M") == 0, "%s", buf);

	gaba_dp_dump_cigar_forward(buf, len, _arr(0x55555000, 0x55555555, 0x55555555, 0x55), 12, 92);
	assert(strcmp(buf, "46M") == 0, "%s", buf);

	gaba_dp_dump_cigar_forward(buf, len, _arr(0x55550555), 0, 32);
	assert(strcmp(buf, "6M4D8M") == 0, "%s", buf);

	gaba_dp_dump_cigar_forward(buf, len, _arr(0x5555f555), 0, 32);
	assert(strcmp(buf, "6M4I8M") == 0, "%s", buf);

	gaba_dp_dump_cigar_forward(buf, len, _arr(0xaaaa0555), 0, 33);
	assert(strcmp(buf, "6M5D8M") == 0, "%s", buf);

	gaba_dp_dump_cigar_forward(buf, len, _arr(0xaaabf555), 0, 33);
	assert(strcmp(buf, "6M5I8M") == 0, "%s", buf);

	gaba_dp_dump_cigar_forward(buf, len, _arr(0xaaabf555, 0xaaaa0556), 0, 65);
	assert(strcmp(buf, "6M5I8M1I5M5D8M") == 0, "%s", buf);

	gaba_dp_dump_cigar_forward(buf, len, _arr(0xaaabf555, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "6M5I8M1I5M5D8M") == 0, "%s", buf);

	gaba_dp_dump_cigar_forward(buf, len, _arr(0xaaabf554, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "2D5M5I8M1I5M5D8M") == 0, "%s", buf);

	#undef _arr
	free(buf);
}

unittest()
{
	uint64_t const len = 16384;
	char *buf = (char *)malloc(len);

	#define _arr(...)		( (uint32_t const []){ 0, 0, __VA_ARGS__, 0, 0 } + 2 )
	gaba_dp_dump_cigar_reverse(buf, len, _arr(0x55555555), 0, 32);
	assert(strcmp(buf, "16M") == 0, "%s", buf);

	gaba_dp_dump_cigar_reverse(buf, len, _arr(0x55555555, 0x55555555), 0, 64);
	assert(strcmp(buf, "32M") == 0, "%s", buf);

	gaba_dp_dump_cigar_reverse(buf, len, _arr(0x55555555, 0x55555555, 0x55555555, 0x55555555), 0, 128);
	assert(strcmp(buf, "64M") == 0, "%s", buf);

	gaba_dp_dump_cigar_reverse(buf, len, _arr(0x55550000, 0x55555555, 0x55555555, 0x55555555), 16, 112);
	assert(strcmp(buf, "56M") == 0, "%s", buf);

	gaba_dp_dump_cigar_reverse(buf, len, _arr(0x55555000, 0x55555555, 0x55555555, 0x55555555), 12, 116);
	assert(strcmp(buf, "58M") == 0, "%s", buf);

	gaba_dp_dump_cigar_reverse(buf, len, _arr(0x55), 0, 8);
	assert(strcmp(buf, "4M") == 0, "%s", buf);

	gaba_dp_dump_cigar_reverse(buf, len, _arr(0x55555000, 0x55555555, 0x55555555, 0x55), 12, 92);
	assert(strcmp(buf, "46M") == 0, "%s", buf);

	gaba_dp_dump_cigar_reverse(buf, len, _arr(0x55550555), 0, 32);
	assert(strcmp(buf, "8M4D6M") == 0, "%s", buf);

	gaba_dp_dump_cigar_reverse(buf, len, _arr(0x5555f555), 0, 32);
	assert(strcmp(buf, "8M4I6M") == 0, "%s", buf);

	gaba_dp_dump_cigar_reverse(buf, len, _arr(0xaaaa0555), 0, 33);
	assert(strcmp(buf, "8M5D6M") == 0, "%s", buf);

	gaba_dp_dump_cigar_reverse(buf, len, _arr(0xaaabf555), 0, 33);
	assert(strcmp(buf, "8M5I6M") == 0, "%s", buf);

	gaba_dp_dump_cigar_reverse(buf, len, _arr(0xaaabf555, 0xaaaa0556), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I6M") == 0, "%s", buf);

	gaba_dp_dump_cigar_reverse(buf, len, _arr(0xaaabf555, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I6M") == 0, "%s", buf);

	gaba_dp_dump_cigar_reverse(buf, len, _arr(0xaaabf554, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I5M2D") == 0, "%s", buf);

	#undef _arr
	free(buf);
}

/**
 * check if gaba_dp_trace returns a correct path
 */
/* with empty sequences */
unittest(with_seq_pair("A", "A"))
{
	omajinai();

	/* fill sections */
	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);

	/* forward-only traceback */
	struct gaba_alignment_s *r = gaba_dp_trace(d, f, NULL, NULL);
	assert(check_result(r, 0, 0, 0, 0, (uint32_t)-1, 0, 0, 0), print_result(r));

	/* forward-reverse traceback */
	r = gaba_dp_trace(d, f, f, NULL);
	assert(check_result(r, 0, 0, 0, 0, (uint32_t)-1, 0, 0, 0), print_result(r));

	/* section added */
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);

	/* forward-only traceback */
	r = gaba_dp_trace(d, f, NULL, NULL);
	assert(check_result(r, 0, 0, 0, 0, (uint32_t)-1, 0, 0, 0), print_result(r));

	/* forward-reverse traceback */
	r = gaba_dp_trace(d, f, f, NULL);
	assert(check_result(r, 0, 0, 0, 0, (uint32_t)-1, 0, 0, 0), print_result(r));

	gaba_dp_clean(d);
}

/* with short sequences */
unittest(with_seq_pair("A", "A"))
{
	omajinai();

	/* fill sections */
	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);

	/* forward-only traceback */
	struct gaba_alignment_s *r = gaba_dp_trace(d, f, NULL, NULL);
	assert(check_result(r, 4, 0, 4, 2, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "DRDR"), print_path(r));
	assert(check_cigar(r, "2M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 1, s->bfsec, 0, 1, 0, 2), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 0, 1, s->bfsec, 0, 1, 2, 2), print_section(r->sec[1]));

	/* reverse-only traceback */
	r = gaba_dp_trace(d, NULL, f, NULL);
	assert(check_result(r, 4, 0, 4, 2, 1, 2, 1, 1), print_result(r));
	assert(check_path(r, "DRDR"), print_path(r));
	assert(check_cigar(r, "2M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 1, s->brsec, 0, 1, 0, 2), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 1, s->brsec, 0, 1, 2, 2), print_section(r->sec[1]));

	/* forward-reverse traceback */
	r = gaba_dp_trace(d, f, f, NULL);
	assert(check_result(r, 8, 0, 8, 4, 2, 0, 0, 0), print_result(r));
	assert(check_path(r, "DRDRDRDR"), print_path(r));
	assert(check_cigar(r, "4M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 1, s->brsec, 0, 1, 0, 2), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 1, s->brsec, 0, 1, 2, 2), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 0, 1, s->bfsec, 0, 1, 4, 2), print_section(r->sec[2]));
	assert(check_section(r->sec[3], s->afsec, 0, 1, s->bfsec, 0, 1, 6, 2), print_section(r->sec[3]));

	/* forward-reverse traceback with seed */
	r = gaba_dp_trace(d, f, f, GABA_TRACE_PARAMS(
		.sec = &((struct gaba_path_section_s){
			.aid = 100, .bid = 102,
			.apos = 0, .bpos = 0,
			.alen = 14, .blen = 14,
			.ppos = 0, /*.plen = 14*/
		}),
		.slen = 1,
		.k = 14));
	assert(check_result(r, 36, 0, 36, 5, 2, 0, 0, 0), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "18M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 1, s->brsec, 0, 1, 0, 2), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 1, s->brsec, 0, 1, 2, 2), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->assec, 0, 14, s->bssec, 0, 14, 4, 28), print_section(r->sec[2]));
	assert(check_section(r->sec[3], s->afsec, 0, 1, s->bfsec, 0, 1, 32, 2), print_section(r->sec[3]));
	assert(check_section(r->sec[4], s->afsec, 0, 1, s->bfsec, 0, 1, 34, 2), print_section(r->sec[4]));

	gaba_dp_clean(d);
}

/* with longer sequences */
unittest(with_seq_pair("ACGTACGTACGT", "ACGTACGTACGT"))
{
	omajinai();

	/* fill sections */
	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);

	/* fw */
	struct gaba_alignment_s *r = gaba_dp_trace(d, f, NULL, NULL);
	assert(check_result(r, 48, 0, 48, 2, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "24M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 12, s->bfsec, 0, 12, 0, 24), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 0, 12, s->bfsec, 0, 12, 24, 24), print_section(r->sec[1]));

	/* rv */
	r = gaba_dp_trace(d, NULL, f, NULL);
	assert(check_result(r, 48, 0, 48, 2, 1, 24, 12, 12), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "24M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 12, s->brsec, 0, 12, 0, 24), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 12, s->brsec, 0, 12, 24, 24), print_section(r->sec[1]));

	/* fw-rv */
	r = gaba_dp_trace(d, f, f, NULL);
	assert(check_result(r, 96, 0, 96, 4, 2, 0, 0, 0), print_result(r));
	assert(check_path(r,
		"DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"
		"DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"),
		print_path(r));
	assert(check_cigar(r, "48M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 12, s->brsec, 0, 12, 0, 24), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 12, s->brsec, 0, 12, 24, 24), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 0, 12, s->bfsec, 0, 12, 48, 24), print_section(r->sec[2]));
	assert(check_section(r->sec[3], s->afsec, 0, 12, s->bfsec, 0, 12, 72, 24), print_section(r->sec[3]));

	gaba_dp_clean(d);
}

/* concatenate */
unittest(with_seq_pair("ACGTACGTACGT", "ACGTACGTACGT"))
{
	omajinai();

	/* fill forward root section */
	struct gaba_fill_s *f1 = gaba_dp_fill_root(d, &s->afsec, 6, &s->bfsec, 6);
	assert(f1->status == 0x1ff, "%x", f1->status);
	assert(check_tail(f1, 0, 0, -19, 1), print_tail(f1));

	/* fill tail section */
	f1 = gaba_dp_fill(d, f1, &s->aftail, &s->bftail);
	assert(f1->status == 0x1ff, "%x", f1->status);
	assert(check_tail(f1, 12, 21, 21, 2), print_tail(f1));

	/* fill forward root section */
	struct gaba_fill_s *f2 = gaba_dp_fill_root(d, &s->arsec, 6, &s->brsec, 6);
	assert(f2->status == 0x1ff, "%x", f2->status);
	assert(check_tail(f2, 0, 0, -19, 1), print_tail(f2));

	/* fill tail section */
	f2 = gaba_dp_fill(d, f2, &s->artail, &s->brtail);
	assert(f2->status == 0x1ff, "%x", f2->status);
	assert(check_tail(f2, 12, 21, 21, 2), print_tail(f2));

	/* fw-rv */
	struct gaba_alignment_s *r = gaba_dp_trace(d, f1, f2, NULL);
	assert(check_result(r, 24, 0, 24, 1, 0, 12, 6, 6), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "12M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 12, s->bfsec, 0, 12, 0, 24), print_section(r->sec[0]));

	gaba_dp_clean(d);
}

/* sequences with different lengths (consumed as mismatches) */
unittest(with_seq_pair("GAAAAAAAA", "AAAAAAAA"))
{
	omajinai();

	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bftail);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);

	/* fw */	
	struct gaba_alignment_s *r = gaba_dp_trace(d, f, NULL, NULL);
	assert(check_result(r, 22, 2, 32, 3, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "16M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 8, s->bfsec, 0, 8, 0, 16), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 8, 1, s->bfsec, 0, 1, 16, 2), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 0, 7, s->bfsec, 1, 7, 18, 14), print_section(r->sec[2]));

	/* rv */
	r = gaba_dp_trace(d, NULL, f, NULL);
	assert(check_result(r, 22, 2, 32, 3, 2, 16, 9, 8), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "16M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 2, 7, s->brsec, 0, 7, 0, 14), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 1, s->brsec, 7, 1, 14, 2), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->arsec, 1, 8, s->brsec, 0, 8, 16, 16), print_section(r->sec[2]));

	/* fw-rv */
	r = gaba_dp_trace(d, f, f, NULL);
	assert(check_result(r, 44, 4, 64, 6, 3, 0, 0, 0), print_result(r));
	assert(check_path(r,
		"DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"
		"DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"),
		print_path(r));
	assert(check_cigar(r, "32M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 2, 7, s->brsec, 0, 7, 0, 14), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 1, s->brsec, 7, 1, 14, 2), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->arsec, 1, 8, s->brsec, 0, 8, 16, 16), print_section(r->sec[2]));
	assert(check_section(r->sec[3], s->afsec, 0, 8, s->bfsec, 0, 8, 32, 16), print_section(r->sec[3]));
	assert(check_section(r->sec[4], s->afsec, 8, 1, s->bfsec, 0, 1, 48, 2), print_section(r->sec[4]));
	assert(check_section(r->sec[5], s->afsec, 0, 7, s->bfsec, 1, 7, 50, 14), print_section(r->sec[5]));

	gaba_dp_clean(d);
}

/* another pair with different lengths */
unittest(with_seq_pair("TTTTTTTT", "CTTTTTTTT"))
{
	omajinai();

	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bfsec);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);

	/* fw */
	struct gaba_alignment_s *r = gaba_dp_trace(d, f, NULL, NULL);
	assert(check_result(r, 22, 2, 32, 3, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "16M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 8, s->bfsec, 0, 8, 0, 16), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 0, 1, s->bfsec, 8, 1, 16, 2), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 1, 7, s->bfsec, 0, 7, 18, 14), print_section(r->sec[2]));

	/* rv */
	r = gaba_dp_trace(d, NULL, f, NULL);
	assert(check_result(r, 22, 2, 32, 3, 2, 16, 8, 9), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "16M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 7, s->brsec, 2, 7, 0, 14), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 7, 1, s->brsec, 0, 1, 14, 2), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->arsec, 0, 8, s->brsec, 1, 8, 16, 16), print_section(r->sec[2]));

	/* fw-rv */
	r = gaba_dp_trace(d, f, f, NULL);
	assert(check_result(r, 44, 4, 64, 6, 3, 0, 0, 0), print_result(r));
	assert(check_path(r,
		"DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"
		"DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "32M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 7, s->brsec, 2, 7, 0, 14), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 7, 1, s->brsec, 0, 1, 14, 2), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->arsec, 0, 8, s->brsec, 1, 8, 16, 16), print_section(r->sec[2]));
	assert(check_section(r->sec[3], s->afsec, 0, 8, s->bfsec, 0, 8, 32, 16), print_section(r->sec[3]));
	assert(check_section(r->sec[4], s->afsec, 0, 1, s->bfsec, 8, 1, 48, 2), print_section(r->sec[4]));
	assert(check_section(r->sec[5], s->afsec, 1, 7, s->bfsec, 0, 7, 50, 14), print_section(r->sec[5]));

	gaba_dp_clean(d);
}

/* with deletions */
unittest(with_seq_pair("GACGTACGT", "ACGTACGT"))
{
	omajinai();

	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bftail);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);

	/* fw */
	struct gaba_alignment_s *r = gaba_dp_trace(d, f, NULL, NULL);
	assert(check_result(r, 20, 0, 34, 2, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "RDRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "1D8M1D8M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 9, s->bfsec, 0, 8, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 0, 9, s->bfsec, 0, 8, 17, 17), print_section(r->sec[1]));

	/* rv */
	r = gaba_dp_trace(d, NULL, f, NULL);
	assert(check_result(r, 20, 0, 34, 2, 1, 17, 9, 8), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDRR"), print_path(r));
	assert(check_cigar(r, "8M1D8M1D"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 9, s->brsec, 0, 8, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 9, s->brsec, 0, 8, 17, 17), print_section(r->sec[1]));

	/* fw-rv */
	r = gaba_dp_trace(d, f, f, NULL);
	/* fixme!! continuous gaps at the root must be concatenated! */
	assert(check_result(r, 40, 0, 68, 4, 2, 0, 0, 0), print_result(r));
	assert(check_path(r,
		"DRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDRR"
		"RDRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "8M1D8M2D8M1D8M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 9, s->brsec, 0, 8, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 9, s->brsec, 0, 8, 17, 17), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 0, 9, s->bfsec, 0, 8, 34, 17), print_section(r->sec[2]));
	assert(check_section(r->sec[3], s->afsec, 0, 9, s->bfsec, 0, 8, 51, 17), print_section(r->sec[3]));

	gaba_dp_clean(d);
}

/* with insertions */
unittest(with_seq_pair("ACGTACGT", "GACGTACGT"))
{
	omajinai();

	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bfsec);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);

	/* fw */
	struct gaba_alignment_s *r = gaba_dp_trace(d, f, NULL, NULL);
	assert(check_result(r, 20, 0, 34, 2, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "DDRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "1I8M1I8M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 8, s->bfsec, 0, 9, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 0, 8, s->bfsec, 0, 9, 17, 17), print_section(r->sec[1]));

	/* rv */
	r = gaba_dp_trace(d, NULL, f, NULL);
	assert(check_result(r, 20, 0, 34, 2, 1, 17, 8, 9), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDRD"), print_path(r));
	assert(check_cigar(r, "8M1I8M1I"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 8, s->brsec, 0, 9, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 8, s->brsec, 0, 9, 17, 17), print_section(r->sec[1]));

	/* fw-rv */
	r = gaba_dp_trace(d, f, f, NULL);
	assert(check_result(r, 40, 0, 68, 4, 2, 0, 0, 0), print_result(r));
	assert(check_path(r,
		"DRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDRD"
		"DDRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "8M1I8M2I8M1I8M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 8, s->brsec, 0, 9, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 8, s->brsec, 0, 9, 17, 17), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 0, 8, s->bfsec, 0, 9, 34, 17), print_section(r->sec[2]));
	assert(check_section(r->sec[3], s->afsec, 0, 8, s->bfsec, 0, 9, 51, 17), print_section(r->sec[3]));

	gaba_dp_clean(d);
}

/* breakpoint adjustment */
unittest(with_seq_pair("GACGTACGTGACGTACGT", "ACGTACGT"))
{
	omajinai();

	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bftail);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);

	/* fw */
	struct gaba_alignment_s *r = gaba_dp_trace(d, f, NULL, NULL);
	assert(check_result(r, 20, 0, 34, 2, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "RDRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "1D8M1D8M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 10, s->bfsec, 0, 8, 0, 18), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 10, 8, s->bfsec, 0, 8, 18, 16), print_section(r->sec[1]));

	/* rv */
	r = gaba_dp_trace(d, NULL, f, NULL);
	assert(check_result(r, 20, 0, 34, 2, 1, 17, 18, 8), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDRR"), print_path(r));
	assert(check_cigar(r, "8M1D8M1D"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 9, s->brsec, 0, 8, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 9, 9, s->brsec, 0, 8, 17, 17), print_section(r->sec[1]));

	/* fw-rv */
	r = gaba_dp_trace(d, f, f, NULL);
	assert(check_result(r, 40, 0, 68, 4, 2, 0, 0, 0), print_result(r));
	assert(check_path(r,
		"DRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDRR"
		"RDRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "8M1D8M2D8M1D8M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 9, s->brsec, 0, 8, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 9, 9, s->brsec, 0, 8, 17, 17), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 0, 10, s->bfsec, 0, 8, 34, 18), print_section(r->sec[2]));
	assert(check_section(r->sec[3], s->afsec, 10, 8, s->bfsec, 0, 8, 52, 16), print_section(r->sec[3]));

	gaba_dp_clean(d);
}

unittest(with_seq_pair("ACGTACGT", "GACGTACGTGACGTACGT"))
{
	omajinai();

	struct gaba_fill_s *f = gaba_dp_fill_root(d, &s->afsec, 0, &s->bfsec, 0);
	f = gaba_dp_fill(d, f, &s->afsec, &s->bfsec);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bfsec);
	f = gaba_dp_fill(d, f, &s->aftail, &s->bftail);

	/* fw */
	struct gaba_alignment_s *r = gaba_dp_trace(d, f, NULL, NULL);
	assert(check_result(r, 20, 0, 34, 2, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "DDRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "1I8M1I8M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 8, s->bfsec, 0, 10, 0, 18), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 0, 8, s->bfsec, 10, 8, 18, 16), print_section(r->sec[1]));

	/* rv */
	r = gaba_dp_trace(d, NULL, f, NULL);
	assert(check_result(r, 20, 0, 34, 2, 1, 17, 8, 18), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDRD"), print_path(r));
	assert(check_cigar(r, "8M1I8M1I"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 8, s->brsec, 0, 9, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 8, s->brsec, 9, 9, 17, 17), print_section(r->sec[1]));

	/* fw-rv */
	r = gaba_dp_trace(d, f, f, NULL);
	assert(check_result(r, 40, 0, 68, 4, 2, 0, 0, 0), print_result(r));
	assert(check_path(r,
		"DRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDRD"
		"DDRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "8M1I8M2I8M1I8M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 8, s->brsec, 0, 9, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 8, s->brsec, 9, 9, 17, 17), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 0, 8, s->bfsec, 0, 10, 34, 18), print_section(r->sec[2]));
	assert(check_section(r->sec[3], s->afsec, 0, 8, s->bfsec, 10, 8, 52, 16), print_section(r->sec[3]));

	gaba_dp_clean(d);
}


/* cross tests */

/**
 * @struct unittest_naive_result_s
 * @brief result container
 */
struct unittest_naive_result_s {
	int32_t score;
	uint32_t path_length;
	int64_t apos, bpos;
	int64_t alen, blen;
	char *path;
};

/**
 * @fn unittest_naive_encode_a
 */
static inline
int8_t unittest_naive_encode(char a)
{
	return(0x03 & ((a>>1) ^ (a>>2)));
}

/**
 * @fn unittest_naive
 *
 * @brief naive implementation of the forward semi-global alignment algorithm
 * left-aligned gap and left-aligned deletion
 */
#define UNITTEST_SEQ_MARGIN			( 8 )			/* add margin to avoid warnings in the glibc strlen */
#define UNITTEST_NAIVE_FORWARD 		( 0 )
#define UNITTEST_NAIVE_REVERSE 		( 1 )
#if MODEL == LINEAR
static
struct unittest_naive_result_s unittest_naive(
	struct gaba_params_s const *sc,
	char const *a,
	char const *b,
	int dir)
{
	/* utils */
	#define _a(p, q, plen)	( (q) * ((plen) + 1) + (p) )
	#define s(p, q)			_a(p, (q), alen)
	#define m(p, q)			( a[(p) - 1] == b[(q) - 1] ? m : x )

	/* load gap penalties */
	int8_t m = sc->m;
	int8_t x = -sc->x;
	int8_t g = -(sc->gi + sc->ge);

	/* calc lengths */
	int64_t alen = strlen(a);
	int64_t blen = strlen(b);

	/* calc min */
	int64_t min = INT16_MIN + sc->x - 2 * g;

	/* malloc matrix */
	int16_t *mat = (int16_t *)malloc(
		(alen + 1) * (blen + 1) * sizeof(int16_t));

	/* init */
	struct unittest_naive_maxpos_s {
		int16_t score;
		int64_t apos;
		int64_t bpos;
	};

	struct unittest_naive_maxpos_s max = { 0, 0, 0 };

	mat[s(0, 0)] = 0;
	for(int64_t i = 1; i < alen+1; i++) {
		mat[s(i, 0)] = MAX2(min, i * g);
	}
	for(int64_t j = 1; j < blen+1; j++) {
		mat[s(0, j)] = MAX2(min, j * g);
	}

	for(int64_t j = 1; j < blen+1; j++) {
		for(int64_t i = 1; i < alen+1; i++) {
			int16_t score = mat[s(i, j)] = MAX4(min,
				mat[s(i - 1, j - 1)] + m(i, j),
				mat[s(i - 1, j)] + g,
				mat[s(i, j - 1)] + g);
			if(score > max.score
			|| (score == max.score && (i + j) < (max.apos + max.bpos))) {
				max = (struct unittest_naive_maxpos_s){
					score, i, j
				};
			}
		}
	}
	if(max.score == 0) {
		max = (struct unittest_naive_maxpos_s){ 0, 0, 0 };
	}

	debug("max(%d), apos(%lld), bpos(%lld)", max.score, max.apos, max.bpos);

	struct unittest_naive_result_s result = {
		.score = max.score,
		.apos = max.apos,
		.bpos = max.bpos,
		.path_length = max.apos + max.bpos + 1,
		.path = (char *)malloc(max.apos + max.bpos + UNITTEST_SEQ_MARGIN)
	};
	if(dir == UNITTEST_NAIVE_FORWARD) {
		/* forward trace */

		int64_t path_index = max.apos + max.bpos + 1;
		while(max.apos > 0 || max.bpos > 0) {
			debug("path_index(%llu), apos(%lld), bpos(%lld)", path_index, max.apos, max.bpos);

			/* M > I > D > X */
			if(max.bpos > 0
			&& mat[s(max.apos, max.bpos)] == mat[s(max.apos, max.bpos - 1)] + g) {
				max.bpos--;
				result.path[--path_index] = 'D';
			} else if(max.apos > 0
			&& mat[s(max.apos, max.bpos)] == mat[s(max.apos - 1, max.bpos)] + g) {
				max.apos--;
				result.path[--path_index] = 'R';
			} else {
				result.path[--path_index] = 'R';
				result.path[--path_index] = 'D';
				max.apos--;
				max.bpos--;
			}
		}
		result.alen = result.apos - max.apos;
		result.blen = result.bpos - max.bpos;
		result.apos = max.apos;
		result.bpos = max.bpos;

		result.path_length -= path_index;
		for(uint64_t i = 0; i < result.path_length; i++) {
			result.path[i] = result.path[path_index++];
		}
		result.path[result.path_length] = '\0';

	} else {
		int64_t path_index = 0;
		while(max.apos > 0 || max.bpos > 0) {
			/* M > I > D > X */
			if(max.apos > 0
			&& mat[s(max.apos, max.bpos)] == mat[s(max.apos - 1, max.bpos)] + g) {
				max.apos--;
				result.path[path_index++] = 'R';
			} else if(max.bpos > 0
			&& mat[s(max.apos, max.bpos)] == mat[s(max.apos, max.bpos - 1)] + g) {
				max.bpos--;
				result.path[path_index++] = 'D';
			} else {
				result.path[path_index++] = 'D';
				result.path[path_index++] = 'R';
				max.apos--;
				max.bpos--;
			}
		}
		result.alen = result.apos - max.apos;
		result.blen = result.bpos - max.bpos;
		result.apos = alen - result.apos;
		result.bpos = blen - result.bpos;

		result.path_length = path_index;
		result.path[result.path_length] = '\0';
	}
	free(mat);

	#undef _a
	#undef s
	#undef m
	return(result);
}
#else /* MODEL == AFFINE */
static
struct unittest_naive_result_s unittest_naive(
	struct gaba_params_s const *sc,
	char const *a,
	char const *b,
	int dir)
{
	/* utils */
	#define _a(p, q, plen)	( (q) * ((plen) + 1) + (p) )
	#define s(p, q)			_a(p, 3*(q), alen)
	#define e(p, q)			_a(p, 3*(q)+1, alen)
	#define f(p, q)			_a(p, 3*(q)+2, alen)
	#define m(p, q)			( a[(p) - 1] == b[(q) - 1] ? m : x )

	/* load gap penalties */
	int8_t m = sc->m;
	int8_t x = -sc->x;
	int8_t gi = -sc->gi;
	int8_t ge = -sc->ge;

	/* calc lengths */
	int64_t alen = strlen(a);
	int64_t blen = strlen(b);

	/* calc min */
	int64_t min = INT16_MIN + sc->x - 2*gi;

	/* malloc matrix */
	int16_t *mat = (int16_t *)malloc(
		3 * (alen + 1) * (blen + 1) * sizeof(int16_t));

	/* init */
	struct unittest_naive_maxpos_s {
		int16_t score;
		int64_t apos;
		int64_t bpos;
	};

	struct unittest_naive_maxpos_s max = { 0, 0, 0 };

	mat[s(0, 0)] = mat[e(0, 0)] = mat[f(0, 0)] = 0;
	for(int64_t i = 1; i < alen+1; i++) {
		mat[s(i, 0)] = mat[e(i, 0)] = MAX2(min, gi + i * ge);
		mat[f(i, 0)] = MAX2(min, gi + i * ge + gi - 1);
	}
	for(int64_t j = 1; j < blen+1; j++) {
		mat[s(0, j)] = mat[f(0, j)] = MAX2(min, gi + j * ge);
		mat[e(0, j)] = MAX2(min, gi + j * ge + gi - 1);
	}

	for(int64_t j = 1; j < blen+1; j++) {
		for(int64_t i = 1; i < alen+1; i++) {
			int16_t score_e = mat[e(i, j)] = MAX2(
				mat[s(i - 1, j)] + gi + ge,
				mat[e(i - 1, j)] + ge);
			int16_t score_f = mat[f(i, j)] = MAX2(
				mat[s(i, j - 1)] + gi + ge,
				mat[f(i, j - 1)] + ge);
			int16_t score = mat[s(i, j)] = MAX4(min,
				mat[s(i - 1, j - 1)] + m(i, j),
				score_e, score_f);
			if(score > max.score
			|| (score == max.score && (i + j) < (max.apos + max.bpos))) {
				max = (struct unittest_naive_maxpos_s){
					score, i, j
				};
			}
		}
	}
	if(max.score == 0) {
		max = (struct unittest_naive_maxpos_s){ 0, 0, 0 };
	}

	struct unittest_naive_result_s result = {
		.score = max.score,
		.apos = max.apos,
		.bpos = max.bpos,
		.path_length = max.apos + max.bpos + 1,
		.path = (char *)malloc(max.apos + max.bpos + UNITTEST_SEQ_MARGIN)
	};
	if(dir == UNITTEST_NAIVE_FORWARD) {
		int64_t path_index = max.apos + max.bpos + 1;
		while(max.apos > 0 || max.bpos > 0) {
			/* M > I > D > X */
			if(mat[s(max.apos, max.bpos)] == mat[f(max.apos, max.bpos)]) {
				while(mat[f(max.apos, max.bpos)] == mat[f(max.apos, max.bpos - 1)] + ge) {
					max.bpos--;
					result.path[--path_index] = 'D';
				}
				max.bpos--;
				result.path[--path_index] = 'D';
			} else if(mat[s(max.apos, max.bpos)] == mat[e(max.apos, max.bpos)]) {
				while(mat[e(max.apos, max.bpos)] == mat[e(max.apos - 1, max.bpos)] + ge) {
					max.apos--;
					result.path[--path_index] = 'R';
				}
				max.apos--;
				result.path[--path_index] = 'R';
			} else {
				result.path[--path_index] = 'R';
				result.path[--path_index] = 'D';
				max.apos--;
				max.bpos--;
			}
		}

		result.alen = result.apos - max.apos;
		result.blen = result.bpos - max.bpos;
		result.apos = max.apos;
		result.bpos = max.bpos;

		result.path_length -= path_index;
		for(uint64_t i = 0; i < result.path_length; i++) {
			result.path[i] = result.path[path_index++];
		}
		result.path[result.path_length] = '\0';

	} else {
		int64_t path_index = 0;
		while(max.apos > 0 || max.bpos > 0) {
			/* M > I > D > X */
			if(mat[s(max.apos, max.bpos)] == mat[e(max.apos, max.bpos)]) {
				while(mat[e(max.apos, max.bpos)] == mat[e(max.apos - 1, max.bpos)] + ge) {
					max.apos--;
					result.path[path_index++] = 'R';
				}
				max.apos--;
				result.path[path_index++] = 'R';
			} else if(mat[s(max.apos, max.bpos)] == mat[f(max.apos, max.bpos)]) {
				while(mat[f(max.apos, max.bpos)] == mat[f(max.apos, max.bpos - 1)] + ge) {
					max.bpos--;
					result.path[path_index++] = 'D';
				}
				max.bpos--;
				result.path[path_index++] = 'D';
			} else {
				result.path[path_index++] = 'D';
				result.path[path_index++] = 'R';
				max.apos--;
				max.bpos--;
			}
		}

		result.alen = result.apos - max.apos;
		result.blen = result.bpos - max.bpos;
		result.apos = alen - result.apos;
		result.bpos = blen - result.bpos;

		result.path_length = path_index;
		result.path[result.path_length] = '\0';
	}
	free(mat);

	#undef _a
	#undef s
	#undef e
	#undef f
	#undef m
	return(result);
}
#endif /* MODEL */

/**
 * @fn unittest_random_base
 */
static _force_inline
char unittest_random_base(void)
{
	char const table[4] = {'A', 'C', 'G', 'T'};
	return(table[rand() % 4]);
}

/**
 * @fn unittest_generate_random_sequence
 */
static _force_inline
char *unittest_generate_random_sequence(
	int64_t len)
{
	char *seq;		/** a pointer to sequence */
	seq = (char *)malloc(sizeof(char) * (len + UNITTEST_SEQ_MARGIN));

	if(seq == NULL) { return NULL; }
	for(int64_t i = 0; i < len; i++) {
		seq[i] = unittest_random_base();
	}
	seq[len] = '\0';
	return seq;
}

/**
 * @fn unittest_generate_mutated_sequence
 */
static _force_inline
char *unittest_generate_mutated_sequence(
	char const *seq,
	double x,
	double d,
	int bw)
{
	if(seq == NULL) { return NULL; }

	int64_t wave = 0;			/** wave is q-coordinate of the alignment path */
	int64_t len = strlen(seq);
	char *mutated_seq = (char *)malloc(sizeof(char) * (len + UNITTEST_SEQ_MARGIN));
	if(mutated_seq == NULL) { return NULL; }
	for(int64_t i = 0, j = 0; i < len; i++) {
		if(((double)rand() / (double)RAND_MAX) < x) {
			mutated_seq[i] = unittest_random_base();	j++;	/** mismatch */
		} else if(((double)rand() / (double)RAND_MAX) < d) {
			if(rand() & 0x01 && wave > -bw+1) {
				mutated_seq[i] = (j < len) ? seq[j++] : unittest_random_base();
				j++; wave--;						/** deletion */
			} else if(wave < bw-2) {
				mutated_seq[i] = unittest_random_base();
				wave++;								/** insertion */
			} else {
				mutated_seq[i] = (j < len) ? seq[j++] : unittest_random_base();
			}
		} else {
			mutated_seq[i] = (j < len) ? seq[j++] : unittest_random_base();
		}
	}
	mutated_seq[len] = '\0';
	return mutated_seq;
}

/**
 * @fn unittest_add_tail
 */
static _force_inline
char *unittest_add_tail(
	char *seq,
	char c,
	int64_t tail_len)
{
	int64_t len = strlen(seq);
	seq = realloc(seq, len + tail_len + UNITTEST_SEQ_MARGIN);

	for(int64_t i = 0; i < tail_len; i++) {
		seq[len + i] = (c == 0) ? unittest_random_base() : c;
	}
	seq[len + tail_len] = '\0';
	return(seq);
}

/* test if the naive implementation is sane */
#define check_naive_result(_r, _score, _path) ( \
	   (_r).score == (_score) \
	&& strcmp((_r).path, (_path)) == 0 \
	&& (_r).path_length == strlen(_path) \
)
#define print_naive_result(_r) \
	"score(%d), path(%s), len(%d)", \
	(_r).score, (_r).path, (_r).path_length

static
char *string_pair_diff(
	char const *a,
	char const *b)
{
	uint64_t len = 2 * (strlen(a) + strlen(b));
	char *base = malloc(len);
	char *ptr = base, *tail = base + len - 1;
	uint64_t state = 0;

	#define push(ch) { \
		*ptr++ = (ch); \
		if(ptr == tail) { \
			base = realloc(base, 2 * len); \
			ptr = base + len; \
			tail = base + 2 * len; \
			len *= 2; \
		} \
	}
	#define push_str(str) { \
		for(uint64_t i = 0; i < strlen(str); i++) { \
			push(str[i]); \
		} \
	}

	uint64_t i;
	for(i = 0; i < MIN2(strlen(a), strlen(b)); i++) {
		if(state == 0 && a[i] != b[i]) {
			push_str("\x1b[31m"); state = 1;
		} else if(state == 1 && a[i] == b[i]) {
			push_str("\x1b[39m"); state = 0;
		}
		push(a[i]);
	}
	if(state == 1) { push_str("\x1b[39m"); state = 0; }
	for(; i < strlen(a); i++) { push(a[i]); }

	push('\n');
	for(uint64_t i = 0; i < strlen(b); i++) {
		push(b[i]);
	}

	push('\0');
	return(base);
}
#define format_string_pair_diff(_a, _b) ({ \
	char *str = string_pair_diff(_a, _b); \
	char *copy = alloca(strlen(str) + 1); \
	strcpy(copy, str); \
	free(str); \
	copy; \
})
#define print_string_pair_diff(_a, _b)		"\n%s", format_string_pair_diff(_a, _b)

#if MODEL == LINEAR
unittest()
{
	struct gaba_params_s const *p = unittest_default_params;
	struct unittest_naive_result_s n;

	/* all matches */
	n = unittest_naive(p, "AAAA", "AAAA", UNITTEST_NAIVE_FORWARD);
	assert(check_naive_result(n, 8, "DRDRDRDR"), print_naive_result(n));
	free(n.path);

	/* with deletions */
	n = unittest_naive(p, "TTTTACGTACGT", "TTACGTACGT", UNITTEST_NAIVE_FORWARD);
	assert(check_naive_result(n, 8, "DRDRRRDRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path);

	/* with insertions */
	n = unittest_naive(p, "TTACGTACGT", "TTTTACGTACGT", UNITTEST_NAIVE_FORWARD);
	assert(check_naive_result(n, 8, "DRDRDDDRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path);
}
#else /* MODEL == AFFINE */
unittest()
{
	struct gaba_params_s const *p = unittest_default_params;
	struct unittest_naive_result_s n;

	/* all matches */
	n = unittest_naive(p, "AAAA", "AAAA", UNITTEST_NAIVE_FORWARD);
	assert(check_naive_result(n, 8, "DRDRDRDR"), print_naive_result(n));
	free(n.path);

	/* with deletions */
	n = unittest_naive(p, "TTTTACGTACGT", "TTACGTACGT", UNITTEST_NAIVE_FORWARD);
	assert(check_naive_result(n, 13, "DRDRRRDRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path);

	/* with insertions */
	n = unittest_naive(p, "TTACGTACGT", "TTTTACGTACGT", UNITTEST_NAIVE_FORWARD);
	assert(check_naive_result(n, 13, "DRDRDDDRDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path);

	/* ins-match-del */
	n = unittest_naive(p, "ATGAAGCTGCGAGGC", "TGATGGCTTGCGAGGC", UNITTEST_NAIVE_FORWARD);
	assert(check_naive_result(n, 6, "DDDRDRDRRRDRDRDRDDRDRDRDRDRDRDR"), print_naive_result(n));
	free(n.path);
}
#endif /* MODEL */


#if 1
/* cross test */
unittest()
{
	struct gaba_context_s const *c = (struct gaba_context_s const *)gctx;
	struct gaba_params_s const *p = unittest_default_params;

	/* seed rand */
	#ifndef SEED
	int32_t seed = getpid();
	#else
	int32_t seed = SEED;
	#endif
	srand(seed);

	// int64_t cross_test_count = 10000000;
	int64_t cross_test_count = 1000;
	for(int64_t i = 0; i < cross_test_count; i++) {
		/* generate sequences */
		char *a = unittest_generate_random_sequence(1000);
		char *b = unittest_generate_mutated_sequence(a, 0.1, 0.1, 500);

		/* add random sequences at the tail */
		a = unittest_add_tail(a, 0, 64);
		b = unittest_add_tail(b, 0, 64);

		/* add tail margin */
		int64_t const mlen = 20;
		a = unittest_add_tail(a, 'C', mlen);
		b = unittest_add_tail(b, 'G', mlen);


		/* naive */
		struct unittest_naive_result_s nf = unittest_naive(p, a, b, UNITTEST_NAIVE_FORWARD);
		struct unittest_naive_result_s nr = unittest_naive(p, a, b, UNITTEST_NAIVE_REVERSE);
		assert(nf.score == nr.score, "(%d, %d)", nf.score, nr.score);
		assert(nf.alen == nr.alen, "(%lld, %lld)", nf.alen, nr.alen);
		assert(nf.blen == nr.blen, "(%lld, %lld)", nf.blen, nr.blen);
		assert(nf.path_length == nr.path_length, "(%u, %u)", nf.path_length, nr.path_length);

		/* build section */
		struct unittest_sections_s *sec = unittest_build_seqs(
			&((struct unittest_seqs_s){ .a = a, .b = b }));

		debug("seed(%d)\n%s", seed, format_string_pair_diff(a, b));

		/* generate dp context */
		struct gaba_dp_context_s *d = gaba_dp_init(c, sec->alim, sec->blim);

		/* fill section */
		struct gaba_section_s const *as = &sec->afsec;
		struct gaba_section_s const *bs = &sec->bfsec;
		struct gaba_fill_s *f = gaba_dp_fill_root(d, as, 0, bs, 0);
		struct gaba_fill_s *m = f;
		
		/* fill tail (1) */
		as = (f->status & GABA_STATUS_UPDATE_A) ? &sec->aftail : as;
		bs = (f->status & GABA_STATUS_UPDATE_B) ? &sec->bftail : bs;
		struct gaba_fill_s *t1 = gaba_dp_fill(d, f, as, bs);
		m = (t1->max > m->max) ? t1 : m;

		/* fill tail (2) */
		as = (t1->status & GABA_STATUS_UPDATE_A) ? &sec->aftail : as;
		bs = (t1->status & GABA_STATUS_UPDATE_B) ? &sec->bftail : bs;
		struct gaba_fill_s *t2 = gaba_dp_fill(d, t1, as, bs);
		m = (t2->max > m->max) ? t2 : m;

		/* check max */
		assert(m->max == nf.score, "m->max(%lld), f(%lld, %u), t1->max(%lld, %u), t2->max(%lld, %u), n.score(%d)",
			m->max, f->max, f->status, t1->max, t1->status, t2->max, t2->status, nf.score);
		if(m->max != nf.score) {
			struct gaba_fill_s *f2 = gaba_dp_fill_root(d, &sec->afsec, 0, &sec->bfsec, 0);
			(void)f2;
			debug("refill f2(%lld, %u)", f2->max, f2->status);
		}

		/* forward trace */
		struct gaba_alignment_s *rf = gaba_dp_trace(d, m, NULL, NULL);

		/* check results */
		assert(rf->score == nf.score, "m->max(%lld), rf->score(%lld), nf.score(%d)",
			m->max, rf->score, nf.score);
		assert(rf->sec[0].apos == nf.apos, "apos(%u, %lld)", rf->sec[0].apos, nf.apos);
		assert(rf->sec[0].bpos == nf.bpos, "bpos(%u, %lld)", rf->sec[0].bpos, nf.bpos);
		assert(rf->sec[0].alen == nf.alen, "alen(%u, %lld)", rf->sec[0].alen, nf.alen);
		assert(rf->sec[0].blen == nf.blen, "blen(%u, %lld)", rf->sec[0].blen, nf.blen);
		assert(check_path(rf, nf.path), "\n%s\n%s\n%s",
			a, b, format_string_pair_diff(decode_path(rf), nf.path));

		int64_t acnt = 0, bcnt = 0;
		for(int64_t i = 0; i < rf->path->len; i++) {
			if(((rf->path->array[i / 32]>>(i & 31)) & 0x01) == 0) {
				acnt++;
			} else {
				bcnt++;
			}
		}
		assert(acnt == rf->sec[0].alen, "acnt(%lld), alen(%u)", acnt, rf->sec[0].alen);
		assert(bcnt == rf->sec[0].blen, "bcnt(%lld), blen(%u)", bcnt, rf->sec[0].blen);

		debug("score(%lld, %d), alen(%lld), blen(%lld)\n%s",
			rf->score, nf.score, nf.alen, nf.blen,
			format_string_pair_diff(decode_path(rf), nf.path));

		/* reverse trace */
		struct gaba_alignment_s *rr = gaba_dp_trace(d, NULL, m, NULL);

		/* check results */
		assert(rr->score == nr.score, "m->max(%lld), rr->score(%lld), nr.score(%d)",
			m->max, rr->score, nr.score);
		assert(rr->sec[0].apos == nr.apos - mlen, "apos(%u, %lld)", rr->sec[0].apos, nr.apos);
		assert(rr->sec[0].bpos == nr.bpos - mlen, "bpos(%u, %lld)", rr->sec[0].bpos, nr.bpos);
		assert(rr->sec[0].alen == nr.alen, "alen(%u, %lld)", rr->sec[0].alen, nr.alen);
		assert(rr->sec[0].blen == nr.blen, "blen(%u, %lld)", rr->sec[0].blen, nr.blen);
		assert(check_path(rr, nr.path), "\n%s\n%s\n%s",
			a, b, format_string_pair_diff(decode_path(rr), nr.path));

		acnt = 0, bcnt = 0;
		for(int64_t i = 0; i < rr->path->len; i++) {
			if(((rr->path->array[i / 32]>>(i & 31)) & 0x01) == 0) {
				acnt++;
			} else {
				bcnt++;
			}
		}
		assert(acnt == rr->sec[0].alen, "acnt(%lld), alen(%u)", acnt, rr->sec[0].alen);
		assert(bcnt == rr->sec[0].blen, "bcnt(%lld), blen(%u)", bcnt, rr->sec[0].blen);

		debug("score(%lld, %d), alen(%lld), blen(%lld)\n%s",
			rr->score, nr.score, nr.alen, nr.blen,
			format_string_pair_diff(decode_path(rr), nr.path));


		/* cleanup */
		gaba_dp_clean(d);
		free(sec);
		free(nf.path);
		free(nr.path);
		free(a);
		free(b);
	}
}
#else
/* for debugging */
unittest(with_seq_pair(
"CTGCGCGAGTCTGCCATGAAATCGAGCTTACAATCCCGATCTTCTCAGCCCTATTGCGGATAGTAGTATATTCA",
"ACGTGCGCGGTGGTTGCTCTTCTGGACGCGTTCGACACGTATTACGAAGTCCTTACCGCTATAAATCACAACGC"))
{
	omajinai();
	struct gaba_params_s const *p = unittest_default_params;

	/* fill section */
	struct gaba_section_s const *as = &s->afsec;
	struct gaba_section_s const *bs = &s->bfsec;
	struct gaba_fill_s *f = gaba_dp_fill_root(d, as, 0, bs, 0);
	
	/* fill tail (1) */
	as = (f->status & GABA_STATUS_UPDATE_A) ? &s->aftail : as;
	bs = (f->status & GABA_STATUS_UPDATE_B) ? &s->bftail : bs;
	struct gaba_fill_s *t1 = gaba_dp_fill(d, f, as, bs);
	f = (t1->max > f->max) ? t1 : f;

	/* fill tail (2) */
	as = (f->status & GABA_STATUS_UPDATE_A) ? &s->aftail : as;
	bs = (f->status & GABA_STATUS_UPDATE_B) ? &s->bftail : bs;
	struct gaba_fill_s *t2 = gaba_dp_fill(d, t1, as, bs);
	f = (t2->max > f->max) ? t2 : f;
	struct gaba_result_s *r = gaba_dp_trace(d, f, NULL, NULL);

	/* naive */
	char const *a = (char const *)s->a;
	char const *b = (char const *)s->b;
	struct unittest_naive_result_s n = unittest_naive(p, a, b);

	/* check scores */
	assert(r->score == n.score, "f->max(%lld), r->score(%lld), n.score(%d)",
		f->max, r->score, n.score);
	assert(check_path(r, n.path), "\n%s\n%s\n%s",
		a, b, format_string_pair_diff(decode_path(r), n.path));

	debug("score(%lld, %d), alen(%lld), blen(%lld)\n%s",
		r->score, n.score, n.alen, n.blen,
		format_string_pair_diff(decode_path(r), n.path));

	/* cleanup */
	gaba_dp_clean(d);
}
#endif

#endif /* UNITTEST */

/**
 * end of gaba.c
 */
