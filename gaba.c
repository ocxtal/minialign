
/**
 * @file gaba.c
 *
 * @brief libgaba (libsea3) DP routine implementation
 *
 * @author Hajime Suzuki
 * @date 2016/1/11
 * @license Apache v2
 */
// #define DEBUG
/* make sure POSIX APIs are properly activated */
#if defined(__linux__) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE		200112L
#endif

#if defined(__darwin__) && !defined(_BSD_SOURCE)
#  define _BSD_SOURCE
#endif

/* import general headers */
#include <stdio.h>				/* sprintf in dump_path */
#include <stdint.h>				/* uint32_t, uint64_t, ... */
#include <stddef.h>				/* offsetof */
#include <string.h>				/* memset, memcpy */
#include "gaba.h"
#include "log.h"
#include "sassert.h"
#include "arch/arch.h"			/* architecture dependents */


/*
 * gap penalty model configuration: choose one of the following three by a define macro
 *   linear-gap penalty function: g(k) = ge * k
 *   affine-gap penalty function: g(k) = gi + ge * k
 *   combined function: g(k) = min(gf * k, gi + ge * k) where ge < gf
 */
#define LINEAR 						1
#define AFFINE						2
#define COMBINED					3

#ifdef MODEL
#  if !(MODEL == LINEAR || MODEL == AFFINE || MODEL == COMBINED)
#    error "MODEL must be LINEAR (1), AFFINE (2), or COMBINED (3)."
#  endif
#else
#  define MODEL 					AFFINE
#endif

#if MODEL == LINEAR
#  define MODEL_LABEL				linear
#elif MODEL == AFFINE
#  define MODEL_LABEL				affine
#else
#  define MODEL_LABEL				combined
#endif


/* define ENABLE_FILTER to enable gapless alignment filter */
// #define ENABLE_FILTER


/* bandwidth-specific configurations aliasing vector macros */
#define BW_MAX						64
#ifndef BW
#  define BW						64
#endif

#if BW == 16
#  define _NVEC_ALIAS_PREFIX		v16i8
#  define _WVEC_ALIAS_PREFIX		v16i16
#  define DP_CTX_INDEX				2
#elif BW == 32
#  define _NVEC_ALIAS_PREFIX		v32i8
#  define _WVEC_ALIAS_PREFIX		v32i16
#  define DP_CTX_INDEX				1
#elif BW == 64
#  define _NVEC_ALIAS_PREFIX		v64i8
#  define _WVEC_ALIAS_PREFIX		v64i16
#  define DP_CTX_INDEX				0
#else
#  error "BW must be either 16, 32, or 64."
#endif
#include "arch/vector_alias.h"

#define DP_CTX_MAX					( 3 )
// #define _dp_ctx_index(_bw)			( DP_CTX_MAX - ((_bw)>>4) )
#define _dp_ctx_index(_bw)			( ((_bw) == 64) ? 0 : (((_bw) == 32) ? 1 : 2) )
_static_assert(_dp_ctx_index(BW) == DP_CTX_INDEX);


/* add suffix for gap-model- and bandwidth-wrapper (see gaba_wrap.h) */
#ifdef SUFFIX
#  define _suffix_cat3_2(a, b, c)	a##_##b##_##c
#  define _suffix_cat3(a, b, c)		_suffix_cat3_2(a, b, c)
#  define _suffix(_base)			_suffix_cat3(_base, MODEL_LABEL, BW)
#else
#  define _suffix(_base)			_base
#endif


/* add namespace for arch wrapper (see main.c) */
#ifdef NAMESPACE
#  define _export_cat(x, y)			x##_##y
#  define _export_cat2(x, y)		_export_cat(x, y)
#  define _export(_base)			_export_cat2(NAMESPACE, _suffix(_base))
#else
#  define _export(_base)			_suffix(_base)
#endif


/* import unittest */
#ifndef UNITTEST_UNIQUE_ID
#  if MODEL == LINEAR
#    if BW == 16
#      define UNITTEST_UNIQUE_ID	31
#    elif BW == 32
#      define UNITTEST_UNIQUE_ID	32
#    else
#      define UNITTEST_UNIQUE_ID	33
#    endif
#  elif MODEL == AFFINE
#    if BW == 16
#      define UNITTEST_UNIQUE_ID	34
#    elif BW == 32
#      define UNITTEST_UNIQUE_ID	35
#    else
#      define UNITTEST_UNIQUE_ID	36
#    endif
#  else
#    if BW == 16
#      define UNITTEST_UNIQUE_ID	37
#    elif BW == 32
#      define UNITTEST_UNIQUE_ID	38
#    else
#      define UNITTEST_UNIQUE_ID	39
#    endif
#  endif
#endif
#include  "unittest.h"


/* internal constants */
#define BLK_BASE					( 5 )
#define BLK 						( 0x01<<BLK_BASE )

#define MIN_BULK_BLOCKS				( 32 )
#define MEM_ALIGN_SIZE				( 32 )		/* 32byte aligned for AVX2 environments */
#define MEM_INIT_SIZE				( (uint64_t)256 * 1024 * 1024 )
#define MEM_MARGIN_SIZE				( 2048 )	/* tail margin of internal memory blocks */
#define GP_INIT						( 1 )		/* global p coordinate where total fetched sequence length becomes BW */
#define GP_ROOT						( -1 )		/* global p coordinate for the first vector of the first block */

/* test consistency of exported macros */
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
_static_assert(sizeof(struct gaba_params_s) == 24);
_static_assert(sizeof(struct gaba_section_s) == 16);
_static_assert(sizeof(struct gaba_fill_s) == 24);
_static_assert(sizeof(struct gaba_path_section_s) == 32);
_static_assert(sizeof(struct gaba_alignment_s) == 64);
_static_assert(sizeof(nvec_masku_t) == BW / 8);

/**
 * @macro _plen
 * @brief extract plen from path_section_s
 */
#define _plen(sec)		( (sec)->alen + (sec)->blen )

/* forward declarations */
static int32_t gaba_dp_add_stack(struct gaba_dp_context_s *self, uint64_t size);
static void *gaba_dp_malloc(struct gaba_dp_context_s *self, uint64_t size);
static void gaba_dp_free(struct gaba_dp_context_s *self, void *ptr);				/* do nothing */
struct gaba_dp_context_s;


/**
 * @struct gaba_small_delta_s
 */
struct gaba_small_delta_s {
	int8_t delta[BW];			/** (32) small delta */
};
_static_assert(sizeof(struct gaba_small_delta_s) == BW);

/**
 * @struct gaba_drop_s
 */
struct gaba_drop_s {
	int8_t drop[BW];			/** (32) max */
};
_static_assert(sizeof(struct gaba_small_delta_s) == BW);

/**
 * @struct gaba_middle_delta_s
 */
struct gaba_middle_delta_s {
	int16_t delta[BW];			/** (64) middle delta */
};
_static_assert(sizeof(struct gaba_middle_delta_s) == sizeof(int16_t) * BW);

/**
 * @struct gaba_mask_pair_u
 */
#if MODEL == LINEAR
struct gaba_mask_pair_s {
	nvec_masku_t h;				/** (4) horizontal mask vector */
	nvec_masku_t v;				/** (4) vertical mask vector */
};
_static_assert(sizeof(struct gaba_mask_pair_s) == BW / 4);
#else	/* affine and combined */
struct gaba_mask_pair_s {
	nvec_masku_t h;				/** (4) horizontal mask vector */
	nvec_masku_t v;				/** (4) vertical mask vector */
	nvec_masku_t e;				/** (4) e mask vector */
	nvec_masku_t f;				/** (4) f mask vector */
};
_static_assert(sizeof(struct gaba_mask_pair_s) == BW / 2);
#endif

/**
 * @struct gaba_diff_vec_s
 */
#if MODEL == LINEAR
struct gaba_diff_vec_s {
	uint8_t dh[BW];				/** (32) dh */
	uint8_t dv[BW];				/** (32) dv */
};
_static_assert(sizeof(struct gaba_diff_vec_s) == 2 * BW);
#else	/* affine and combined gap penalty */
struct gaba_diff_vec_s {
	uint8_t dh[BW];				/** (32) dh */
	uint8_t dv[BW];				/** (32) dv */
	uint8_t de[BW];				/** (32) de */
	uint8_t df[BW];				/** (32) df */
};
_static_assert(sizeof(struct gaba_diff_vec_s) == 4 * BW);
#endif

/**
 * @struct gaba_char_vec_s
 */
struct gaba_char_vec_s {
	uint8_t w[BW];				/** (32) a in the lower 4bit, b in the higher 4bit */
};
_static_assert(sizeof(struct gaba_char_vec_s) == BW);

/**
 * @struct gaba_block_s
 * @brief a unit of banded matrix, 32 vector updates will be recorded in a single block object.
 * phantom is an alias of the block struct as a head cap of contiguous blocks.
 */
struct gaba_block_s {
	struct gaba_mask_pair_s mask[BLK];	/** (256 / 512) traceback capability flag vectors (set if transition to the ajdacent cell is possible) */
	struct gaba_small_delta_s sd;		/** (16, 32, 64) difference from the previous block */
	uint32_t dir_mask;					/** (4) extension direction bit array */
	int16_t acnt, bcnt;					/** (4) forwarded length */
	uint64_t max_update_mask;			/** (8) lanewise update mask (set if the lane contains the current max) */
	struct gaba_diff_vec_s diff; 		/** (64, 128, 256) diff variables of the last vector */
};
struct gaba_phantom_s {
	uint32_t root;						/** (4) root: -1, head: 0 (FIXME: overlaps with dir_mask) */
	int16_t acnt, bcnt;					/** (4) prefetched sequence lengths (only effective at the root, otherwise zero) */
	struct gaba_joint_tail_s const *tail;/** (8) link to the previous tail */
	struct gaba_diff_vec_s diff; 		/** (64, 128, 256) diff variables of the last (just before the head) vector */
};
_static_assert(sizeof(struct gaba_block_s) % 16 == 0);
_static_assert(sizeof(struct gaba_phantom_s) % 16 == 0);
#define _last_block(x)				( (struct gaba_block_s *)(x) - 1 )
#define _last_phantom(x)			( (struct gaba_phantom_s *)(x) - 1 )

/**
 * @struct gaba_section_pair_s
 */
struct gaba_section_pair_s {
	uint8_t const *atail, *btail;		/** (16) tail of the current section */
	uint32_t alen, blen;				/** (8) lengths of the current section */
	uint32_t aid, bid;					/** (8) ids */
};

/**
 * @struct gaba_joint_tail_s
 * @brief (internal) tail cap of a contiguous matrix blocks, contains a context of the blocks
 * (band) and can be connected to the next blocks.
 */
struct gaba_joint_tail_s {
	/* char vector and delta vectors */
	struct gaba_char_vec_s ch;		/** (16, 32, 64) char vector */
	struct gaba_drop_s xd;			/** (16, 32, 64) */
	struct gaba_middle_delta_s md;	/** (32, 64, 128) */

	uint32_t aridx, bridx;			/** (8) reverse indices for the tails */
	uint32_t asridx, bsridx;		/** (8) start reverse indices (for internal use) */
	uint32_t pridx;					/** (4) p reverse index */
	int32_t acc;					/** (4) accumulator */
	int64_t offset;					/** (8) large offset */
	struct gaba_fill_s f;			/** (24) */

	/* tail pointer */
	struct gaba_joint_tail_s const *tail;/** (8) the previous tail */

	/* section info */
	struct gaba_section_pair_s s;	/** (32) */
};
_static_assert((sizeof(struct gaba_joint_tail_s) % 32) == 0);
#define TAIL_BASE				( offsetof(struct gaba_joint_tail_s, f) )
#define _tail(x)				( (struct gaba_joint_tail_s *)((uint8_t *)(x) - TAIL_BASE) )
#define _fill(x)				( (struct gaba_fill_s *)((uint8_t *)(x) + TAIL_BASE) )

/**
 * @struct gaba_joint_block_s
 * @brief (internal) result container for block fill functions
 */
struct gaba_joint_block_s {
	struct gaba_block_s *blk;
	uint32_t stat;
};

/**
 * @struct gaba_root_block_s
 */
struct gaba_root_block_s {
	struct gaba_phantom_s blk;
	struct gaba_joint_tail_s tail;
	uint8_t _pad[640 - sizeof(struct gaba_phantom_s) - sizeof(struct gaba_joint_tail_s)];
};
_static_assert(sizeof(struct gaba_root_block_s) == 640);

/**
 * @struct gaba_reader_work_s
 * @brief (internal) working buffer for fill-in functions, contains sequence prefetch buffer
 * and middle/max-small deltas.
 */
struct gaba_reader_work_s {
	/** 64byte aligned */
	uint8_t bufa[BW_MAX + BLK];			/** (128) */
	uint8_t bufb[BW_MAX + BLK];			/** (128) */
	/** 256 */

	/** 64byte alidned */
	struct gaba_section_pair_s s;		/** (32) section pair */
	int32_t acc;						/** (4) direction deteminer accumulator */
	uint32_t ofsd;						/** (4) delta of large offset */
	uint32_t aridx, bridx;				/** (8) current ridx */
	uint32_t asridx, bsridx;			/** (8) start ridx */
	uint32_t pridx, psridx;				/** (8) p-coordinate (currently unused) */
	/** 64 */

	/** 64byte aligned */
	struct gaba_drop_s xd;				/** (16, 32, 64) current drop from max */
#if BW != 64
	uint8_t _pad[BW == 16 ? 16 : 32];	/** padding to align to 64-byte boundary */
#endif
	struct gaba_middle_delta_s md;		/** (32, 64, 128) */
};
_static_assert((sizeof(struct gaba_reader_work_s) % 64) == 0);

/**
 * @struct gaba_aln_intl_s
 * @brief internal alias of gaba_alignment_t, allocated in the local context and used as a working buffer.
 */
struct gaba_aln_intl_s {
	void *opaque;
	gaba_free_t free;

	uint32_t head_margin;
	uint32_t sec_len;
	struct gaba_path_section_s const *sec;

	uint64_t path_len;

	int64_t score;				/** (8) score */
	uint32_t mcnt, xcnt;		/** (8) #mismatchs */
	uint32_t gicnt, gecnt;		/** (8) #gap opens, #gap bases */
};
_static_assert(sizeof(struct gaba_alignment_s) == sizeof(struct gaba_aln_intl_s));
_static_assert(offsetof(struct gaba_alignment_s, sec_len) == offsetof(struct gaba_aln_intl_s, sec_len));
_static_assert(offsetof(struct gaba_alignment_s, sec) == offsetof(struct gaba_aln_intl_s, sec));
_static_assert(offsetof(struct gaba_alignment_s, path_len) == offsetof(struct gaba_aln_intl_s, path_len));
_static_assert(offsetof(struct gaba_alignment_s, mcnt) == offsetof(struct gaba_aln_intl_s, mcnt));
_static_assert(offsetof(struct gaba_alignment_s, xcnt) == offsetof(struct gaba_aln_intl_s, xcnt));
_static_assert(offsetof(struct gaba_alignment_s, gicnt) == offsetof(struct gaba_aln_intl_s, gicnt));
_static_assert(offsetof(struct gaba_alignment_s, gecnt) == offsetof(struct gaba_aln_intl_s, gecnt));

/**
 * @struct gaba_leaf_s
 * @brief working buffer for max score search (currently allocated on stack)
 */
struct gaba_leaf_s {
	struct gaba_joint_tail_s const *tail;
	struct gaba_block_s const *blk;
	uint32_t p, q;						/** (8) local p (to restore mask pointer), local q */
	uint64_t ppos;						/** (8) global p (== resulting path length) */
	uint32_t aridx, bridx;
};

/**
 * @struct gaba_writer_work_s
 * @brief working buffer for traceback (allocated in the thread-local context)
 */
struct gaba_writer_work_s {
	/** local context */
	struct gaba_aln_intl_s a;			/** (64) working buffer, copied to the result object */

	/** work */
	uint32_t *path;						/** (8) path array pointer */
	uint32_t ofs, _pad;					/** (8) path array offset */
	struct gaba_block_s const *blk;		/** (8) current block */
	uint32_t p, q;						/** (8) local p, q-coordinate, [0, BW) */

	/** save */
	uint32_t agidx, bgidx;				/** (8) grid indices of the current trace */
	uint32_t asgidx, bsgidx;			/** (8) base indices of the current trace */
	uint32_t aid, bid;					/** (8) section ids */

	/** section info */
	struct gaba_joint_tail_s const *atail;/** (8) */
	struct gaba_joint_tail_s const *btail;/** (8) */
	struct gaba_alignment_s *aln;		/** (8) */

	uint64_t _pad1[6];
	/** 64, 192 */
};
_static_assert((sizeof(struct gaba_writer_work_s) % 64) == 0);

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
};
_static_assert(sizeof(struct gaba_mem_block_s) == 24);

/**
 * @struct gaba_stack_s
 * @brief save stack pointer
 */
struct gaba_stack_s {
	struct gaba_mem_block_s *mem;
	uint8_t *top;
	uint8_t *end;
};
_static_assert(sizeof(struct gaba_stack_s) == 24);

/**
 * @struct gaba_dp_context_s
 *
 * @brief (internal) container for dp implementations
 */
struct gaba_dp_context_s {
	/* working buffers */
	union gaba_work_s {
		struct gaba_reader_work_s r;	/** (192) */
		struct gaba_writer_work_s l;	/** (192) */
	} w;
	/** 64byte aligned */

	/* tail-of-userland pointers */
	uint8_t const *alim, *blim;			/** (16) max index of seq array */

	/* memory management */
	struct gaba_mem_block_s mem;		/** (24) root memory block */
	struct gaba_stack_s stack;			/** (24) current stack */
	/** 64byte aligned */

	/** loaded on init */
	struct gaba_score_vec_s scv;		/** (80) substitution matrix and gaps */

	int8_t m;							/** (1) match award */
	int8_t x;							/** (1) mismatch penalty (neg.int) */
	int8_t gi;							/** (1) gap open penalty */
	int8_t ge;							/** (1) gap extension penalty */
	int8_t gf;							/** (1) linear-gap extension penalty for short indels (gf > ge) */
	int8_t tx;							/** (1) xdrop threshold */
	uint8_t tf;							/** (1) ungapped alignment filter threshold */
	uint8_t _pad1;

	/** output options */
	uint32_t head_margin;				/** (1) margin at the head of gaba_res_t */
	uint32_t tail_margin;				/** (1) margin at the tail of gaba_res_t */
	void *_pad2[4];
	/** 128; 64byte aligned */

	/** phantom vectors */
	struct gaba_root_block_s ph[3];		/** (768) [0] for 16-cell, [1] for 32-cell, [2] for 64-cell */
};
_static_assert((sizeof(struct gaba_dp_context_s) % 64) == 0);
#define GABA_DP_CONTEXT_LOAD_OFFSET	( offsetof(struct gaba_dp_context_s, scv) )
#define GABA_DP_CONTEXT_LOAD_SIZE	( sizeof(struct gaba_dp_context_s) - GABA_DP_CONTEXT_LOAD_OFFSET )
_static_assert((GABA_DP_CONTEXT_LOAD_OFFSET % 64) == 0);
_static_assert((GABA_DP_CONTEXT_LOAD_SIZE % 64) == 0);
#define _proot(_t)					( &(_t)->ph[_dp_ctx_index(BW)] )
#define _ptail(_t)					( &_proot(_t)->tail )

/**
 * @struct gaba_opaque_s
 */
struct gaba_opaque_s {
	void *api[4];
};
#define _export_dp_context(_t) ( \
	(struct gaba_dp_context_s *)(((struct gaba_opaque_s *)(_t)) - DP_CTX_MAX + _dp_ctx_index(BW)) \
)
#define _restore_dp_context(_t) ( \
	(struct gaba_dp_context_s *)(((struct gaba_opaque_s *)(_t)) - _dp_ctx_index(BW) + DP_CTX_MAX) \
)
#define _export_dp_context_global(_t) ( \
	(struct gaba_dp_context_s *)(((struct gaba_opaque_s *)(_t)) - DP_CTX_MAX + _dp_ctx_index(BW)) \
)
#define _restore_dp_context_global(_t) ( \
	(struct gaba_dp_context_s *)(((struct gaba_opaque_s *)(_t)) - _dp_ctx_index(BW) + DP_CTX_MAX) \
)

/**
 * @struct gaba_context_s
 *
 * @brief (API) an algorithmic context.
 *
 * @sa gaba_init, gaba_close
 */
struct gaba_context_s {
	/** opaque pointers for function dispatch */
	struct gaba_opaque_s api[4];		/** (128) */
	/** 64byte aligned */

	/** templates */
	struct gaba_dp_context_s k;			/** (896) */
	/** 64byte aligned */
};
_static_assert((sizeof(struct gaba_dp_context_s) % 32) == 0);
#define _pmd(_c)					( &(_c)->md[_dp_ctx_index(BW)] )

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
void *gaba_malloc(
	size_t size)
{
	void *ptr = NULL;
	if(posix_memalign(&ptr, MEM_ALIGN_SIZE, size) != 0) {
		debug("posix_memalign failed");
		return(NULL);
	}
	debug("posix_memalign(%p)", ptr);
	return(ptr);
}
static inline
void gaba_free(
	void *ptr)
{
	free(ptr);
	return;
}


/* matrix fill functions */

/* direction macros */
/**
 * @struct gaba_dir_s
 */
struct gaba_dir_s {
	uint32_t mask;
	int32_t acc;
};

/**
 * @macro _dir_init
 */
#define _dir_init(_self) ((struct gaba_dir_s){ .mask = 0, .acc = (_self)->w.r.acc })
/**
 * @macro _dir_fetch
 */
#define _dir_fetch(_d) { \
	(_d).mask <<= 1; (_d).mask |= (uint32_t)((_d).acc < 0);  \
	debug("fetched dir(%x), %s", _d, _dir_is_down(_d) ? "go down" : "go right"); \
}
/**
 * @macro _dir_update
 * @brief update direction determiner for the next band
 */
#define _dir_update(_d, _vector, _sign) { \
	(_d).acc += (_sign) * (_ext_n(_vector, 0) - _ext_n(_vector, BW-1)); \
	/*debug("acc(%d), (%d, %d)", _dir_acc, _ext_n(_vector, 0), _ext_n(_vector, BW-1));*/ \
}
/**
 * @macro _dir_adjust_remainder
 * @brief adjust direction array when termination is detected in the middle of the block
 */
#define _dir_adjust_remainder(_d, _filled_count) { \
	debug("adjust remainder, array(%x), shifted array(%x)", (_d).mask, (_d).mask<<(BLK - (_filled_count))); \
	(_d).mask <<= (BLK - (_filled_count)); \
}
/**
 * @macro _dir_is_down, _dir_is_right
 * @brief direction indicator (_dir_is_down returns true if dir == down)
 */
#define _dir_is_down(_d)			( ((uint64_t)(_d).mask) & 0x01 )
#define _dir_is_right(_d)			( ~((uint64_t)(_d).mask) & 0x01 )
/**
 * @macro _dir_save
 */
#define _dir_save(_self, _blk, _d) { \
	(_blk)->dir_mask = (_d).mask;	/* store mask */ \
	(_self)->w.r.acc = (_d).acc; 	/* store accumulator */ \
}
#if 1
/**
 * @macro _dir_load
 */
#define _dir_load(_blk, _filled_count) ({ \
	struct gaba_dir_s _d = (struct gaba_dir_s){ \
		.mask = (_blk)->dir_mask, \
		.acc = 0 \
	}; \
	debug("load dir idx(%d), mask(%x), shifted mask(%x)", (int32_t)_filled_count, _d.mask, _d.mask>>(BLK - (_filled_count))); \
	_d.mask >>= (BLK - (_filled_count)); \
	_d; \
})
#else
/**
 * @macro _dir_load
 */
#define _dir_load(_blk, _local_idx) ({ \
	struct gaba_dir_s _d = (struct gaba_dir_s){ \
		.mask = (_blk)->dir_mask, \
		.acc = 0 \
	}; \
	debug("load dir idx(%d), mask(%x), shifted mask(%x)", (int32_t)_local_idx, _d.mask, _d.mask>>(BLK - (_local_idx) - 1)); \
	_d.mask >>= (BLK - (_local_idx) - 1); \
	_d; \
})
#endif
/**
 * @macro _dir_bcnt
 */
#define _dir_bcnt(_d) ( \
	popcnt((_d).mask) \
)
/**
 * @macro _dir_windback
 */
#define _dir_windback(_d) { \
	(_d).mask >>= 1; \
}

/**
 * @macro _match_n
 * @brief alias to sequence matcher macro
 */
#ifndef _match_n
#  define _match_n				_and_n		/* 4bit encoded */
#endif /* _match */

/**
 * seqreader macros
 */
#define _rd_bufa_base(k)		( (k)->w.r.bufa + BLK + BW )
#define _rd_bufb_base(k)		( (k)->w.r.bufb )
#define _rd_bufa(k, pos, len)	( _rd_bufa_base(k) - (pos) - (len) )
#define _rd_bufb(k, pos, len)	( _rd_bufb_base(k) + (pos) )
#define _lo64(v)				_ext_v2i64(v, 0)
#define _hi64(v)				_ext_v2i64(v, 1)
#define _lo32(v)				_ext_v2i32(v, 0)
#define _hi32(v)				_ext_v2i32(v, 1)


/**
 * @fn fill_load_seq_a
 */
static _force_inline
void fill_load_seq_a(
	struct gaba_dp_context_s *self,
	uint8_t const *pos,
	uint64_t len)
{
	if(pos < self->alim) {
		debug("reverse fetch a: pos(%p), len(%llu)", pos, len);
		/* reverse fetch: 2 * alen - (2 * alen - pos) + (len - 32) */
		v32i8_t a = _loadu_v32i8(pos + (len - BLK));
		_storeu_v32i8(_rd_bufa(self, BW, len), _swap_v32i8(a));
	} else {
		debug("forward fetch a: pos(%p), len(%llu)", pos, len);
		/* take complement */
		static uint8_t const comp[16] __attribute__(( aligned(16) )) = {
			0x00, 0x08, 0x04, 0x0c, 0x02, 0x0a, 0x06, 0x0e,
			0x01, 0x09, 0x05, 0x0d, 0x03, 0x0b, 0x07, 0x0f
		};
		v32i8_t const cv = _from_v16i8_v32i8(_load_v16i8(comp));

		/* forward fetch: 2 * alen - pos */
		v32i8_t a = _loadu_v32i8(_rev(pos, self->alim) - (len - 1));
		_storeu_v32i8(_rd_bufa(self, BW, len), _shuf_v32i8(cv, a));
	}
	return;
}

/**
 * @fn fill_load_seq_b
 */
static _force_inline
void fill_load_seq_b(
	struct gaba_dp_context_s *self,
	uint8_t const *pos,
	uint64_t len)
{
	if(pos < self->blim) {
		debug("forward fetch b: pos(%p), len(%llu)", pos, len);
		/* forward fetch: pos */
		v32i8_t b = _loadu_v32i8(pos);
		_storeu_v32i8(_rd_bufb(self, BW, len), b);
	} else {
		debug("reverse fetch b: pos(%p), len(%llu)", pos, len);
		/* take complement */
		static uint8_t const comp[16] __attribute__(( aligned(16) )) = {
			0x00, 0x08, 0x04, 0x0c, 0x02, 0x0a, 0x06, 0x0e,
			0x01, 0x09, 0x05, 0x0d, 0x03, 0x0b, 0x07, 0x0f
		};
		v32i8_t const cv = _from_v16i8_v32i8(_load_v16i8(comp));

		/* reverse fetch: 2 * blen - pos + (len - 32) */
		v32i8_t b = _loadu_v32i8(_rev(pos, self->blim) - (BLK - 1));
		_storeu_v32i8(_rd_bufb(self, BW, len), _shuf_v32i8(cv, _swap_v32i8(b)));
	}
	return;
}

/**
 * @fn fill_fetch_intl
 * @brief fetch sequences from current sections, lengths must be shorter than 32 (= BLK)
 */
static _force_inline
void fill_fetch_intl(
	struct gaba_dp_context_s *self,
	struct gaba_block_s const *blk,
	v2i32_t const len)
{
	/* fetch seq a */
	nvec_t a = _load_n(_rd_bufa(self, blk->acnt, BW));
	fill_load_seq_a(self, self->w.r.s.atail - self->w.r.aridx, _lo32(len));
	_store_n(_rd_bufa(self, 0, BW), a);

	/* fetch seq b */
	nvec_t b = _load_n(_rd_bufb(self, blk->bcnt, BW));
	_store_n(_rd_bufb(self, 0, BW), b);
	fill_load_seq_b(self, self->w.r.s.btail - self->w.r.bridx, _hi32(len));
	return;
}

/**
 * @fn fill_bulk_fetch
 * @brief fetch 32bases from current section
 */
static _force_inline
void fill_bulk_fetch(
	struct gaba_dp_context_s *self,
	struct gaba_block_s const *blk)
{
	/* bulk fetch */
	v2i32_t const len = _set_v2i32(BLK);
	fill_fetch_intl(self, blk, len);
	return;
}

/**
 * @fn fill_cap_fetch
 */
static _force_inline
void fill_cap_fetch(
	struct gaba_dp_context_s *self,
	struct gaba_block_s const *blk)
{
	/* fetch: len might be clipped by ridx */
	v2i32_t const ridx = _load_v2i32(&self->w.r.aridx);
	v2i32_t const lim = _set_v2i32(BLK);
	fill_fetch_intl(self, blk, _min_v2i32(ridx, lim));
	return;
}

/**
 * @fn fill_init_fetch
 * @brief similar to cap fetch, updating ridx and rem
 */
static _force_inline
void fill_init_fetch(
	struct gaba_dp_context_s *self,
	struct gaba_phantom_s *blk,
	struct gaba_joint_tail_s const *prev_tail)
{
	/* restore remaining head margin lengths */
	#if 0
	v2i32_t const adj = _seta_v2i32(0, 1);
	int64_t prem = -prev_tail->f.ppos;
	v2i32_t rem = _sar_v2i32(_add_v2i32(_set_v2i32(prem), adj), 1);
	#else
	v2i32_t const adj = _seta_v2i32(1, 0);
	int64_t prem = -prev_tail->f.ppos;
	v2i32_t rem = _sar_v2i32(_sub_v2i32(_set_v2i32(prem), adj), 1);
	#endif

	/* load ridx: no need to update ridx since advance counters are always zero */
	v2i32_t ridx = _load_v2i32(&self->w.r.aridx);

	/* fetch, bounded by remaining sequence lengths and remaining head margins */
	v2i32_t len = _min_v2i32(
		rem,									/* if remaining head margin is the shortest */
		_min_v2i32(
			ridx,								/* if remaining sequence length is the shortest */
			_add_v2i32(_swap_v2i32(ridx), adj)	/* if the opposite sequence length is the shortest */
		)
	);
	/* fetch sequence and store at (0, 0), then save newly loaded sequence lengths */
	fill_fetch_intl(self, (struct gaba_block_s *)(blk + 1) - 1, len);
	_store_v2i32(&blk->acnt, _cvt_v2i32_v2i8(len));
	return;
}

/**
 * @fn fill_restore_fetch
 * @brief fetch sequence from an existing block
 */
static _force_inline
void fill_restore_fetch(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *tail,
	struct gaba_block_s const *blk)
{

	#if 0
	nvec_t const mask = _set_n(0x0f);

	/* calc cnt */
	v2i32_t curr_len = _load_v2i32(&blk->aridx);
	v2i32_t prev_len = _load_v2i32(&(blk - 1)->aridx);
	v2i32_t cnt = _sub_v2i32(prev_len, curr_len);

	/* from the current block */
	nvec_t cw = _load_n(&blk->ch.w);
	nvec_t ca = _and_n(mask, cw);
	nvec_t cb = _and_n(mask, _shr_n(cw, 4));
	_storeu_n(_rd_bufa(self, _lo32(cnt), BW), ca);
	_storeu_n(_rd_bufb(self, _hi32(cnt), BW), cb);

	/* from the previous block */
	nvec_t pw = _load_n(&(blk - 1)->ch.w);
	nvec_t pa = _and_n(mask, pw);
	nvec_t pb = _and_n(mask, _shr_n(pw, 4));
	_store_n(_rd_bufa(self, 0, BW), pa);
	_store_n(_rd_bufb(self, 0, BW), pb);

	_print_n(pa);
	_print_n(pb);
	#endif
	return;
}

/*
 * @fn fill_load_section
 */
static _force_inline
void fill_load_section(
	struct gaba_dp_context_s *self,
	struct gaba_section_s const *a,
	struct gaba_section_s const *b,
	v2i32_t ridx,
	uint32_t pridx)
{
	/* load current section lengths */
	v2i64_t asec = _loadu_v2i64(a);					/* tuple of (64bit ptr, 32-bit id, 32-bit len) */
	v2i64_t bsec = _loadu_v2i64(b);

	/* transpose sections */
	v2i32_t aid_alen = _cast_v2i64_v2i32(asec);		/* extract lower 64bit */
	v2i32_t bid_blen = _cast_v2i64_v2i32(bsec);

	v2i32_t id = _lo_v2i32(aid_alen, bid_blen);
	v2i32_t len = _hi_v2i32(aid_alen, bid_blen);
	v2i64_t base = _hi_v2i64(asec, bsec);

	/* calc tail pointer */
	v2i64_t tail = _add_v2i64(base, _cvt_v2i32_v2i64(len));

	_print_v2i32(id);
	_print_v2i32(len);
	_print_v2i64(tail);

	/* store sections */
	_store_v2i64(&self->w.r.s.atail, tail);
	_store_v2i32(&self->w.r.s.alen, len);
	_store_v2i32(&self->w.r.s.aid, id);

	/* calc ridx */
	v2i32_t const z = _zero_v2i32();
	v2i32_t len = _load_v2i32(&self->w.r.s.alen);
	ridx = _sel_v2i32(_eq_v2i32(ridx, z),		/* if ridx is zero (occurs when section is updated) */
		len,									/* load the next section */
		ridx									/* otherwise use the same section as the previous */
	);
	_store_v2i32(&self->w.r.aridx, ridx);
	_store_v2i32(&self->w.r.asridx, ridx);
	self->w.r.pridx = self->w.r.psridx = pridx;
	_print_v2i32(ridx);
	return;
}

/**
 * @fn fill_load_work
 */
static _force_inline
void fill_load_work(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *tail)
{
	/* copy max and middle delta vectors */
	nvec_t xd = _load_n(&prev_tail->xd);
	wvec_t md = _load_w(&prev_tail->md);
	_store_n(&self->w.r.xd, xd);
	_store_w(&self->w.r.md, md);

	/* clear offset */
	self->w.r.acc = prev_tail->acc;
	self->w.r.ofsd = 0;
	return;
}

/**
 * @fn fill_save_work
 */
static _force_inline
void fill_save_work(
	struct gaba_dp_context_s *self)
{
	return;
}

/**
 * @fn fill_load_tail
 * @brief load sequences and indices from the previous tail
 */
static _force_inline
struct gaba_block_s *fill_load_tail(
	struct gaba_dp_context_s *self,
	struct gaba_fill_s const *fill,
	struct gaba_section_s const *a,
	uint32_t aridx,
	struct gaba_section_s const *b,
	uint32_t bridx,
	uint32_t pridx)
{
	debug("create head prev_tail(%p), p(%d), ppos(%lld), scnt(%d)",
		prev_tail, prev_tail->p, prev_tail->f.ppos, prev_tail->scnt);

	/* load sequences */
	nvec_t const mask = _set_n(0x0f);
	nvec_t ch = _load_n(&prev_tail->ch.w);
	nvec_t a = _and_n(mask, ch);
	nvec_t b = _and_n(mask, _shr_n(ch, 4));		/* bit 7 must be cleared not to affect shuffle mask */
	_store_n(_rd_bufa(self, 0, BW), a);
	_store_n(_rd_bufb(self, 0, BW), b);

	/* check if init fetch is needed */
	uint32_t stat = CONT;
	if(prev_tail->f.ppos < GP_ROOT) {
		fill_init_fetch(self, blk, prev_tail);
		if(prev_tail->f.ppos + blk->acnt + blk->bcnt < GP_ROOT) {
			return(NULL);
		}
		/* init fetch done, issue ungapped here filter if needed */
	}

	/* extract block pointer, pass to fill-in loop */
	return((struct gaba_block_s *)prev_tail);
}

/**
 * @fn fill_create_tail
 * @brief create joint_tail at the end of the blocks
 */
static _force_inline
struct gaba_joint_tail_s *fill_create_tail(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *prev_tail,
	struct gaba_block_s *blk,
	uint32_t stat)
{
	debug("create tail, p(%lld)", p);

	/* create joint_tail */
	struct gaba_joint_tail_s *tail = (struct gaba_joint_tail_s *)blk;
	self->stack.top = (void *)(tail + 1);	/* write back stack_top */
	debug("end stack_top(%p), stack_end(%p), blocks(%lld)", self->stack.top, self->stack.end, (p + BLK - 1) / BLK);

	/* store char vector */
	nvec_t a = _loadu_n(_rd_bufa(self, blk->acnt, BW));
	nvec_t b = _loadu_n(_rd_bufb(self, blk->bcnt, BW));
	_store_n(&tail->ch, _or_n(a, _shl_n(b, 4)));

	/* copy delta vectors */
	nvec_t xd = _load_n(&self->w.r.xd);
	wvec_t md = _load_w(&self->w.r.md);
	_store_n(&tail->xd, xd);
	_store_w(&tail->md, md);

	/* search max section */
	md = _sub_w(md, _cvt_n_w(xd));					/* xd holds drop from max */
	uint64_t offset = tail->offset + self->w.r.ofsd;
	uint64_t max = ((int16_t)_hmax_w(md)) + offset;
	_print_w(md);
	debug("offset(%lld), max(%d)", tail->offset, _hmax_w(md));

	/* reverse indices */
	v2i32_t ridx = _load_v2i32(&self->w.r.aridx);
	v2i32_t sridx = _load_v2i32(&self->w.r.asridx);
	_store_v2i32(&tail->aridx, ridx);
	_store_v2i32(&tail->aridx, sridx);
	tail->pridx = self->w.r.pridx;					/* is this really needed? */

	/* store scores */
	tail->acc = self->w.r.acc;
	tail->offset = offset;
	tail->f.max = max;

	/* status flags, coordinates, link pointer, and sections */
	v2i32_t update = _eq_v2i32(ridx, _zero_v2i32());
	v2i32_t adv = _sub_v2i32(sridx, ridx);
	tail->f.stat = stat | _mask_v2i32(update);
	tail->f.scnt = prev_tail->f.scnt - _hi32(update) - _lo32(update);
	tail->f.ppos = prev_tail->f.ppos + _hi32(adv) + _lo32(adv);
	tail->tail = prev_tail;
	_memcpy_blk_aa(&tail->s.atail, &self->w.r.s.atail, sizeof(struct gaba_section_pair_s));
	return(tail);
}

/**
 * @fn fill_create_phantom
 * @brief initialize phantom block
 */
static _force_inline
struct gaba_block_s *fill_create_phantom(
	struct gaba_dp_context_s *self,
	struct gaba_block_s const *prev_blk)
{
	struct gaba_phantom_s *ph = (struct gaba_phantom_s *)self->stack.top;
	debug("start stack_top(%p), stack_end(%p)", self->stack.top, self->stack.end);

	/* head sequence buffers are already filled, continue to body fill-in (first copy phantom block) */
	ph->reserved = 0;							/* overlaps with mask */
	ph->acnt = ph->bcnt = 0;					/* clear counter (a and b sequences are aligned at the head of the buffer) */
	ph->xstat = HEAD;							/* FIXME: mark head */
	ph->blk = prev_blk;
	_memcpy_blk_aa(&ph->diff, &prev_blk->diff, sizeof(struct gaba_diff_vec_s));
	return((struct gaba_block_s *)(ph + 1));
}

/**
 * @macro _fill_load_context
 * @brief load vectors onto registers
 */
#if MODEL == LINEAR
#define _fill_load_context(_blk) \
	debug("blk(%p)", (_blk)); \
	/* load sequence buffer offset */ \
	uint8_t const *aptr = _rd_bufa(self, 0, BW); \
	uint8_t const *bptr = _rd_bufb(self, 0, BW); \
	/* load mask pointer */ \
	struct gaba_mask_pair_s *ptr = (_blk)->mask; \
	/* load vector registers */ \
	register nvec_t dh = _load_n(((_blk) - 1)->diff.dh); \
	register nvec_t dv = _load_n(((_blk) - 1)->diff.dv); \
	_print_n(_add_n(dh, _load_ofsh(self->scv))); \
	_print_n(_add_n(dv, _load_ofsv(self->scv))); \
	/* load delta vectors */ \
	register nvec_t delta = _zero_n(); \
	register nvec_t drop = _load_n(self->w.r.xd.drop); \
	_print_n(drop); \
	_print_w(_add_w(_cvt_n_w(delta), _load_w(_last_block(_ptail(self))->md))); \
	/* load direction determiner */ \
	struct gaba_dir_s dir = _dir_init(self);
#else
#define _fill_load_context(_blk) \
	debug("blk(%p)", (_blk)); \
	/* load sequence buffer offset */ \
	uint8_t const *aptr = _rd_bufa(self, 0, BW); \
	uint8_t const *bptr = _rd_bufb(self, 0, BW); \
	/* load mask pointer */ \
	struct gaba_mask_pair_s *ptr = (_blk)->mask; \
	/* load vector registers */ \
	register nvec_t dh = _load_n(((_blk) - 1)->diff.dh); \
	register nvec_t dv = _load_n(((_blk) - 1)->diff.dv); \
	register nvec_t de = _load_n(((_blk) - 1)->diff.de); \
	register nvec_t df = _load_n(((_blk) - 1)->diff.df); \
	_print_n(_add_n(dh, _load_ofsh(self->scv))); \
	_print_n(_add_n(dv, _load_ofsv(self->scv))); \
	_print_n(_sub_n(_sub_n(de, dv), _load_adjh(self->scv))); \
	_print_n(_sub_n(_add_n(df, dh), _load_adjv(self->scv))); \
	/* load delta vectors */ \
	register nvec_t delta = _zero_n(); \
	register nvec_t drop = _load_n(self->w.r.xd.drop); \
	_print_n(drop); \
	_print_w(_add_w(_cvt_n_w(delta), _load_w(_last_block(_ptail(self))->md))); \
	/* load direction determiner */ \
	struct gaba_dir_s dir = _dir_init(self);
#endif

/**
 * @macro _fill_body
 * @brief update vectors
 */
#if MODEL == LINEAR
#define _fill_body() { \
	register nvec_t t = _match_n(_loadu_n(aptr), _loadu_n(bptr)); \
	t = _shuf_n(_load_sb(self->scv), t); \
	_print_n(t); \
	t = _max_n(dh, t); \
	t = _max_n(dv, t); \
	ptr->h.mask = _mask_n(_eq_n(t, dv)); \
	ptr->v.mask = _mask_n(_eq_n(t, dh)); \
	debug("mask(%x, %x)", ptr->h.all, ptr->v.all); \
	ptr++; \
	nvec_t _dv = _sub_n(t, dh); \
	dh = _sub_n(t, dv); \
	dv = _dv; \
	_print_n(_add_n(dh, _load_ofsh(self->scv))); \
	_print_n(_add_n(dv, _load_ofsv(self->scv))); \
}
#else /* MODEL == AFFINE */
#define _fill_body() { \
	register nvec_t t = _match_n(_loadu_n(aptr), _loadu_n(bptr)); \
	_print_n(_loadu_n(aptr)); \
	_print_n(_loadu_n(bptr)); \
	t = _shuf_n(_load_sb(self->scv), t); \
	_print_n(t); \
	t = _max_n(de, t); \
	t = _max_n(df, t); \
	ptr->h.mask = _mask_n(_eq_n(t, de)); \
	ptr->v.mask = _mask_n(_eq_n(t, df)); \
	debug("mask(%x, %x)", ptr->h.all, ptr->v.all); \
	/* update de and dh */ \
	de = _add_n(de, _load_adjh(self->scv)); \
	nvec_t te = _max_n(de, t); \
	ptr->e.mask = _mask_n(_eq_n(te, de)); \
	de = _add_n(te, dh); \
	dh = _add_n(dh, t); \
	/* update df and dv */ \
	df = _add_n(df, _load_adjv(self->scv)); \
	nvec_t tf = _max_n(df, t); \
	ptr->f.mask = _mask_n(_eq_n(tf, df)); \
	df = _sub_n(tf, dv); \
	t = _sub_n(dv, t); \
	ptr++; \
	dv = dh; dh = t; \
	_print_n(_add_n(dh, _load_ofsh(self->scv))); \
	_print_n(_add_n(dv, _load_ofsv(self->scv))); \
	_print_n(_sub_n(_sub_n(de, dv), _load_adjh(self->scv))); \
	_print_n(_sub_n(_add_n(df, dh), _load_adjv(self->scv))); \
}
#endif /* MODEL */

/**
 * @macro _fill_update_delta
 * @brief update small delta vector and max vector
 */
#define _fill_update_delta(_op_add, _op_subs, _vector, _offset, _sign) ({ \
	nvec_t _t = _add_n(_vector, _offset); \
	nvec_t _m = _gt_n(drop, _t); \
	delta = _op_add(delta, _t); \
	drop = _op_subs(drop, _t); \
	_dir_update(dir, _vector, _sign); \
	_print_w(_add_w(_set_w(offset), _add_w(_cvt_n_w(delta), _load_w(&Wself->w.r.md)))); \
	_print_w(_add_w(_set_w(offset), _add_w(_cvt_n_w(drop), _load_w(&self->w.r.md)))); \
	_m;		/* update-flag vector: each element is set only when the new cell value *exceeded* the previous max */ \
})

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
#define _fill_right() ({ \
	dh = _bsl_n(dh, 1);	/* shift left dh */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(_add_n, _subs_n, dh, _load_ofsh(self->scv), 1); \
})
#else
#define _fill_right() ({ \
	dh = _bsl_n(dh, 1);	/* shift left dh */ \
	df = _bsl_n(df, 1);	/* shift left df */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(_sub_n, _adds_n, dh, _load_ofsh(self->scv), -1); \
})
#endif /* MODEL */
#define _fill_down_update_ptr() { \
	bptr++;				/* increment sequence buffer pointer */ \
}
#define _fill_down_windback_ptr() { \
	bptr--; \
}
#if MODEL == LINEAR
#define _fill_down() ({ \
	dv = _bsr_n(dv, 1);	/* shift right dv */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(_add_n, _subs_n, dv, _load_ofsv(self->scv), 1); \
})
#else
#define _fill_down() ({ \
	dv = _bsr_n(dv, 1);	/* shift right dv */ \
	de = _bsr_n(de, 1);	/* shift right de */ \
	_fill_body();		/* update vectors */ \
	_fill_update_delta(_add_n, _subs_n, dv, _load_ofsv(self->scv), 1); \
})
#endif /* MODEL */

/**
 * @macro _fill_store_context
 * @brief store vectors at the end of the block
 */
#define _fill_store_context_head(_blk) { \
	/* store direction array */ \
	_dir_save(self, _blk, dir); \
	/* update offset */ \
	self->w.r.ofsd += _ext_n(delta, BW/2); \
	/* calc cnt and update ridxs */ \
	int32_t acnt = _rd_bufa(self, 0, BW) - aptr; \
	int32_t bcnt = bptr - _rd_bufb(self, 0, BW); \
	(_blk)->acnt = acnt; (_blk)->bcnt = bcnt; \
	self->w.r.aridx -= acnt; self->w.r.bridx -= bcnt; \
	/* update max and middle vectors in the working buffer */ \
	nvec_t prev_drop = _load_n(&self->w.r.xd); \
	prev_drop = _add_n(prev_drop, delta); \
	(_blk)->max_update_mask = ((nvec_masku_t){ \
		.mask = _mask_n(_gt_n(drop, prev_drop)) \
	}).all; \
}
#define _fill_store_context_tail(_blk) { \
	/* update max delta vector */ \
	_store_n(&self->w.r.xd, drop); \
	/* update middle delta vector */ \
	wvec_t md = _load_w(&self->w.r.md); \
	md = _add_w(md, _cvt_n_w(delta)); \
	_store_w(&self->w.r.md, md); \
}
#if MODEL == LINEAR
#define _fill_store_context(_blk) ({ \
	_fill_store_context_head(_blk); \
	_store_n((_blk)->diff.dh, dh); _print_n(dh); \
	_store_n((_blk)->diff.dv, dv); _print_n(dv); \
	_fill_store_context_tail(_blk); \
})
#else
#define _fill_store_context(_blk) ({ \
	_fill_store_context_head(_blk); \
	_store_n((_blk)->diff.dh, dh); _print_n(dh); \
	_store_n((_blk)->diff.dv, dv); _print_n(dv); \
	_store_n((_blk)->diff.dh, de); _print_n(de); \
	_store_n((_blk)->diff.dv, df); _print_n(df); \
	_fill_store_context_tail(_blk); \
})
#endif

/**
 * @fn fill_test_xdrop
 * @brief returns negative if terminate-condition detected
 */
static _force_inline
int64_t fill_test_xdrop(
	struct gaba_dp_context_s const *self)
{
	return(self->tx - self->w.r.xd.drop[BW/2]);
}

/**
 * @fn fill_bulk_test_idx
 * @brief returns negative if ij-bound (for the bulk fill) is invaded
 */
static _force_inline
int64_t fill_bulk_test_idx(
	struct gaba_dp_context_s const *self)
{
	debug("test(%lld, %lld), len(%d, %d)",
		(int64_t)self->w.r.aridx - BW,
		(int64_t)self->w.r.bridx - BW,
		self->w.r.aridx, self->w.r.bridx);
	return(((int64_t)self->w.r.aridx - BW) | ((int64_t)self->w.r.bridx - BW));
}

/**
 * @fn fill_cap_test_idx
 * @brief returns negative if ij-bound (for the cap fill) is invaded
 */
#define _fill_cap_test_idx_init() \
	uint8_t *alim = _rd_bufa(self, self->w.r.aridx, BW); \
	uint8_t *blim = _rd_bufb(self, self->w.r.bridx, BW);
#define _fill_cap_test_idx() ( \
	((int64_t)aptr - (int64_t)alim) | ((int64_t)blim - (int64_t)bptr) \
)

/**
 * @fn fill_bulk_block
 * @brief fill a block
 */
static _force_inline
void fill_bulk_block(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk)
{
	/* fetch sequence */
	fill_bulk_fetch(self, blk);

	/* load vectors onto registers */
	debug("blk(%p)", blk);
	_fill_load_context(blk);
	/**
	 * @macro _fill_block
	 * @brief unit unrolled fill-in loop
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

	/* store vectors */
	_fill_store_context(blk);
	return;
}

/**
 * @fn fill_bulk_predetd_blocks
 * @brief fill <cnt> contiguous blocks without ij-bound test
 */
static _force_inline
struct gaba_block_s *fill_bulk_predetd_blocks(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk,
	uint64_t cnt)
{
	/* decrement blk to expose phantom block */
	blk--;

	/* bulk fill loop, first check termination */
	struct gaba_block_s *tblk = blk + cnt - 1;
	while((blk->xstat | (ptrdiff_t)(tblk - blk)) >= 0) {
		/* bulk fill */
		debug("blk(%p)", blk);
		fill_bulk_block(self, ++blk);
	}

	/* fix status flag */
	blk->xstat = blk->xstat < 0 ? TERM : CONT;
	return(blk);
}

/**
 * @fn fill_bulk_seq_bounded
 * @brief fill blocks with ij-bound test
 */
static _force_inline
struct gaba_block_s *fill_bulk_seq_bounded(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk)
{
	/* decrement blk to expose phantom block */
	blk--;

	/* bulk fill loop, first check termination */
	while((blk->xstat | fill_bulk_test_idx(self)) >= 0) {
		debug("blk(%p)", blk + 1);
		fill_bulk_block(self, ++blk);
	}

	/* fix status flag */
	blk->xstat = blk->xstat < 0 ? TERM : CONT;
	return(blk);
}

/**
 * @fn fill_cap_seq_bounded
 * @brief fill blocks with cap test
 */
static _force_inline
struct gaba_block_s *fill_cap_seq_bounded(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk)
{
	#define _fill_cap_seq_bounded_core(_dir) { \
		/* update sequence coordinate and then check term */ \
		_fill_##_dir##_update_ptr(); \
		if(_fill_cap_test_idx() < 0) { \
			_fill_##_dir##_windback_ptr(); \
			_dir_windback(dir); \
			break; \
		} \
		_fill_##_dir();		/* update band */ \
	}

	blk--;									/* decrement blk to expose phantom block */
	while(blk->xstat >= 0) {
		/* fetch sequence */
		fill_cap_fetch(self, ++blk);
		_fill_cap_test_idx_init();
		_fill_load_context(blk);			/* contains ptr as struct gaba_mask_pair_s *ptr = blk->mask; */

		/* update diff vectors */
		struct gaba_mask_pair_s *tptr = &blk->mask[BLK];
		while(ptr < tptr) {					/* ptr is automatically incremented in _fill_right() or _fill_down() */
			_dir_fetch(dir);				/* determine direction */
			if(_dir_is_right(dir)) {
				_fill_cap_seq_bounded_core(right);
			} else {
				_fill_cap_seq_bounded_core(down);
			}
		}

		uint64_t i = ptr - blk->mask;		/* calc filled count */
		_dir_adjust_remainder(dir, i);		/* adjust dir remainder */
		_fill_store_context(blk);			/* store mask and vectors */
		if(i != BLK) { blk -= i == 0; break; }/* break if not filled full length */
	}

	/* fix status flag */
	blk->xstat = blk->xstat < 0 ? TERM : UPDATE;
	return(blk);
}

/**
 * @fn max_blocks_mem
 * @brief calculate maximum number of blocks (limited by stack size)
 */
static _force_inline
uint64_t max_blocks_mem(
	struct gaba_dp_context_s const *self)
{
	uint64_t mem_size = self->stack.end - self->stack.top;
	uint64_t blk_cnt = mem_size / sizeof(struct gaba_block_s);
	debug("calc_max_block_mem, stack_top(%p), stack_end(%p), mem_size(%llu), cnt(%llu)",
		self->stack.top, self->stack.end, mem_size, (blk_cnt > 3) ? (blk_cnt - 3) : 0);
	return(((blk_cnt > 3) ? blk_cnt : 3) - 3);
}

/**
 * @fn max_blocks_idx
 * @brief calc max #expected blocks from block,
 * used to determine #blocks which can be filled without mem boundary check.
 */
static _force_inline
uint64_t max_blocks_idx(
	struct gaba_dp_context_s const *self)
{
	uint64_t p = MIN2(self->w.r.aridx, self->w.r.bridx);
	return((2*p + p/2) / BLK);
}

/**
 * @fn min_blocks_idx
 * @brief calc min #expected blocks from block,
 * #blocks filled without seq boundary check.
 */
static _force_inline
uint64_t min_blocks_idx(
	struct gaba_dp_context_s const *self)
{
	uint64_t p = MIN2(self->w.r.aridx, self->w.r.bridx);
	return((p + p/2) / BLK);
}

/**
 * @fn fill_seq_bounded
 * @brief fill blocks with seq bound tests (without mem test), adding head and tail
 */
static _force_inline
struct gaba_block_s *fill_seq_bounded(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk)
{
	uint64_t cnt;						/* #blocks filled in bulk */
	while((cnt = min_blocks_idx(self)) > MIN_BULK_BLOCKS) {
		/* bulk fill without ij-bound test */
		if((blk = fill_bulk_predetd_blocks(self, blk, cnt))->xstat != CONT) {
			return(blk);				/* xdrop termination detected, skip cap */
		}
	}

	/* bulk fill with ij-bound test */
	if((blk = fill_bulk_seq_bounded(self, blk))->xstat != CONT) {
		return(blk);					/* xdrop termination detected, skip cap */
	}

	/* cap fill (without p-bound test) */
	return(fill_cap_seq_bounded(self, blk));
}

/**
 * @fn fill_section_seq_bounded
 * @brief fill dp matrix inside section pairs
 */
static _force_inline
struct gaba_block_s *fill_section_seq_bounded(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *prev_blk)
{
	/* create phantom head */
	struct gaba_block_s *blk = fill_create_phantom(self, prev_blk);

	/* extra large bulk fill (with stack allocation) */
	uint64_t mem_cnt, seq_cnt;			/* #blocks filled in bulk */
	while((mem_cnt = max_blocks_mem(self)) < (seq_cnt = max_blocks_idx(self))) {
		debug("mem_cnt(%llu), seq_cnt(%llu)", mem_cnt, seq_cnt);

		/* here mem limit reaches first (seq sections too long), so fill without seq test */
		if((mem_cnt = MIN2(mem_cnt, min_blocks_idx(self))) > MIN_BULK_BLOCKS) {
			debug("mem bounded fill");
			if((blk = fill_bulk_predetd_blocks(self, blk, mem_cnt))->xstat != CONT) {
				return(blk);
			}
		}

		/* memory ran out: malloc a next stack and create a new phantom head */
		debug("add stack");
		if(gaba_dp_add_stack(self, 0) != 0) { return(NULL); }
		blk = fill_create_phantom(self, blk);
	}
	debug("v(%p), ppos(%lld), p(%d)", _last_block(tail)->md, tail->f.ppos, tail->p);

	/* bulk fill with seq bound check */
	return(fill_seq_bounded(self, blk));
}

/**
 * @fn gaba_dp_fill_root
 *
 * @brief build_root API
 */
struct gaba_fill_s *_export(gaba_dp_fill_root)(
	struct gaba_dp_context_s *self,
	struct gaba_section_s const *a,
	uint32_t apos,
	struct gaba_section_s const *b,
	uint32_t bpos)
{
	/* restore dp context pointer by adding offset */
	self = _restore_dp_context(self);

	/* load sections and extract the last block pointer */
	struct gaba_block_s *blk = fill_load_tail(
		self, _pfill(self),
		a, a->len - apos,
		b, b->len - bpos,
		UINT32_MAX
	);

	/* fill blocks then create a tail cap */
	return(fill_create_tail(self, fill_section_seq_bounded(self, blk)));
}

/**
 * @fn gaba_dp_fill
 *
 * @brief fill API
 */
struct gaba_fill_s *_export(gaba_dp_fill)(
	struct gaba_dp_context_s *self,
	struct gaba_fill_s const *fill,
	struct gaba_section_s const *a,
	struct gaba_section_s const *b)
{
	self = _restore_dp_context(self);

	/* load sections and extract the last block pointer */
	struct gaba_block_s *blk = fill_load_tail(
		self, _pfill(self),
		a, _tail(fill)->aridx,
		b, _tail(fill)->bridx,
		_tail(fill)->pridx
	);

	/* fill blocks then create a tail cap */
	return(fill_create_tail(self, fill_section_seq_bounded(self, blk)));
}


/* trace leaf search functions */
/**
 * @fn leaf_load_max_mask
 */
static _force_inline
uint64_t leaf_load_max_mask(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *tail)
{
	debug("p(%d), ppos(%lld), scnt(%d), offset(%lld)",
		tail->p, tail->f.ppos, tail->scnt, tail->offset);

	/* load max vector, create mask */
	nvec_t drop = _load_n(tail->xd.drop);
	wvec_t md = _load_w(tail->md.delta);
	uint64_t max_mask = ((nvec_masku_t){
		.mask = _mask_w(_eq_w(
			_set_w(tail->max - tail->offset),
			_add_w(md, _cvt_n_w(drop))
		))
	}).all;
	debug("max_mask(%x)", max_mask);
	_print_w(_set_w(tail->max - tail->offset));
	_print_w(_add_w(md, _cvt_n_w(drop)));
	return(max_mask);
}

/**
 * @fn leaf_refill_block
 */
static _force_inline
void leaf_refill_block(
	struct gaba_dp_context_s *self,
	struct gaba_block_s *blk,
	nvec_masku_t *max_mask_arr)
{
	/* vectors on registers inside self block */
	#define _fill_block_leaf(_m) { \
		_dir_fetch(dir); \
		nvec_t update; \
		if(_dir_is_right(dir)) { \
			_fill_right_update_ptr(); \
			update = _fill_right(); \
		} else { \
			_fill_down_update_ptr(); \
			update = _fill_down(); \
		} \
		_m++->mask = _mask_n(update); \
		debug("mask(%x)", ((nvec_masku_t){ .mask = (_m - 1)->mask }).all); \
	}

	/* load contexts and overwrite max vector */
	_fill_load_context(blk);
	for(int64_t i = 0; i < blk->acnt + blk->bcnt; i++) {
		_fill_block_leaf(max_mask_arr);
	}
}

/**
 * @fn leaf_search
 */
static _force_inline
void leaf_search(
	struct gaba_dp_context_s *self,
	struct gaba_joint_tail_s const *tail,
	struct gaba_leaf_s *leaf)
{
	/* load mask and block pointers */
	uint64_t max_mask = leaf_load_max_mask(self, tail);
	struct gaba_block_s *b = _last_block(tail) + 1;

	/*
	 * iteratively clear lanes with longer paths;
	 * max_mask will be zero if the block contains the maximum scoring cell with the shortest path
	 */
	v2i32_t ridx = _load_v2i32(&tail->aridx);
	while(_last_phantom(b)->root != -1 && (max_mask &= ~(--b->max_update_mask)) != 0) {
		v2i8_t cnt = _load_v2i8(&b->acnt);
		ridx = _add_v2i32(ridx, _cvt_v2i8_v2i32(cnt));
	}

	/* fetch from existing blocks */
	fill_restore_fetch(self, tail, b, ridx);

	/* refill block to obtain cell-wise update-masks */
	nvec_masku_t max_mask_arr[BLK];		/* cell-wise update-mask array */
	leaf_refill_block(self, b, max_mask_arr);

	/* search max cell, probe from the tail to the head */
	nvec_masku_t *m = &max_mask_arr[b->acnt + b->bcnt];
	while(m > max_mask_arr && (--m->all & b->max_update_mask) != 0) {}

	/* calc (p, q) coordinates */
	leaf->blk = b;
	leaf->p = m - max_mask_arr;
	leaf->q = tzcnt(m->all & b->max_update_mask);

	/* restore reverse indices */
	int64_t filled_count = (m - max_mask_arr) + 1;		/* p + 1 */
	int32_t bcnt = _dir_bcnt(_dir_load(b, filled_count));
	int32_t acnt = filled_count - bcnt;
	ridx = _add_v2i32(ridx,
		_seta_v2i32(
			(b->bcnt - bcnt) + (BW - leaf->q),
			(b->acnt - acnt) + (1 + leaf->q)
		)
	);

	v2i32_t eridx = _load_v2i32(&tail->aridx);
	v2i32_t rem = _sub_v2i32(ridx, eridx);
	leaf->ppos = tail->f.ppos - _hi32(rem) - _lo32(rem);
	_store_v2i32(&leaf->aridx, ridx);
	return;
}


/**
 * @fn gaba_dp_search_max
 */
struct gaba_pos_pair_s _export(gaba_dp_search_max)(
	struct gaba_dp_context_s *self,
	struct gaba_fill_s const *fill)
{
	self = _restore_dp_context(self);

	struct gaba_leaf_s leaf;
	leaf_search(self, _tail(fill), &leaf);

	struct gaba_joint_tail_s const *atail = _tail(fill), *btail = _tail(fill);
	int32_t aidx = atail->alen - leaf.aridx, bidx = btail->blen - leaf.bridx;

	while(aidx < 0) {
		for(atail = atail->tail; (atail->stat & GABA_STATUS_UPDATE_A) == 0; atail = atail->tail) {}
		aidx += atail->alen;
	}
	while(bidx < 0) {
		for(btail = btail->tail; (btail->stat & GABA_STATUS_UPDATE_B) == 0; btail = btail->tail) {}
		bidx += btail->blen;
	}
	return((struct gaba_pos_pair_s){
		.apos = aidx,
		.bpos = bidx
	});
}


/* path trace functions */
/**
 * @fn trace_reload_section
 * @brief reload section and adjust gidx, pass 0 for section a or 1 for section b
 */
static _force_inline
void trace_reload_section(
	struct gaba_dp_context_s *self,
	uint64_t i)
{
	#define _r(_x, _idx)		( (&(_x))[(_idx)] )
	static uint32_t const mask[2] = {
		GABA_STATUS_UPDATE_A,
		GABA_STATUS_UPDATE_B
	};

	debug("load section a, idx(%d), len(%d)", self->w.l.aidx, self->w.l.atail->alen);

	/* load tail pointer (must be inited with leaf tail) */
	struct gaba_joint_tail_s const *tail = _r(self->w.l.atail, i);
	int32_t gidx = _r(self->w.l.agidx, i) + _r(tail->alen, i);

	while(gidx <= 0) {
		for(tail = tail->tail; (tail->stat & mask[i]) == 0; tail = tail->tail) {}
		gidx += _r(tail->alen, i);
	}

	/* reload finished, store section info */
	_r(self->w.l.atail, i) = tail;		/* fixme: is self correct?? */
	_r(self->w.l.aid, i) = _r(tail->aid, i);

	_r(self->w.l.agidx, i) = gidx;
	_r(self->w.l.asgidx, i) = gidx;
	return;

	#undef _r
}

/**
 * @macro _trace_inc_*
 * @brief increment gap counters
 */
#define _trace_inc_gi()					{ gc = _sub_v2i32(gc, m01); }
#define _trace_dec_gi()					{ gc = _add_v2i32(gc, m01); }
#define _trace_inc_ge()					{ gc = _sub_v2i32(gc, m10); }
#define _trace_dec_ge()					{ gc = _add_v2i32(gc, m10); }

/**
 * @macro _trace_*_*_test_index
 * @brief test if path reached a section boundary
 */
#define _trace_bulk_v_test_index()		( 0 )
#define _trace_bulk_d_test_index()		( 0 )
#define _trace_bulk_h_test_index()		( 0 )
#define _trace_tail_v_test_index()		( _test_v2i32(gidx, m10) /* bgidx == 0 */ )
#define _trace_tail_d_test_index()		( _test_v2i32(_eq_v2i32(gidx, m00), m11) /* agidx == 0 || bgidx == 0 */ )
#define _trace_tail_h_test_index()		( _test_v2i32(gidx, m01) /* agidx == 0 */ )

/**
 * @macro _trace_*_*_update_index
 * @brief update indices by one
 */
#define _trace_bulk_v_update_index()	;
#define _trace_bulk_h_update_index()	;
#define _trace_tail_v_update_index()	{ gidx = _add_v2i32(gidx, m10); _print_v2i32(gidx); }
#define _trace_tail_h_update_index()	{ gidx = _add_v2i32(gidx, m01); _print_v2i32(gidx); }

/**
 * @macro _trace_test_*
 * @brief test mask
 */
#if MODEL == LINEAR
#define _trace_test_diag_h()			( (mask->h.all>>q) & 0x01 )
#define _trace_test_diag_v()			( (mask->v.all>>q) & 0x01 )
#define _trace_test_gap_h()				( (mask->h.all>>q) & 0x01 )
#define _trace_test_gap_v()				( (mask->v.all>>q) & 0x01 )
#else /* MODEL == AFFINE */
#define _trace_test_diag_h()			( (mask->h.all>>q) & 0x01 )
#define _trace_test_diag_v()			( (mask->v.all>>q) & 0x01 )
#define _trace_test_gap_h()				( (mask->e.all>>q) & 0x01 )
#define _trace_test_gap_v()				( (mask->f.all>>q) & 0x01 )
#endif

/**
 * @macro _trace_*_update_path_q
 */
#define _trace_h_update_path_q() { \
	path_array = path_array<<1; \
	mask--; \
	q += _dir_is_down(dir); \
	_dir_windback(dir); \
}
#define _trace_v_update_path_q() { \
	path_array = (path_array<<1) | 0x01; \
	mask--; \
	q += _dir_is_down(dir) - 1; \
	_dir_windback(dir); \
}

/**
 * @macro _trace_load_block_rem
 * @brief reload mask pointer and direction mask, and adjust path offset
 */
#define _trace_load_block_rem(_idx) ({ \
	path -= ofs < (_idx) + 1; \
	ofs = (ofs - ((_idx) + 1)) & (BLK - 1); \
	_dir_load(blk, (_idx) + 1).mask; \
})

/**
 * @macro _trace_reload_tail
 * @brief reload tail, issued at each band-segment boundaries
 */
#define _trace_reload_tail(t) { \
	/* reload block pointer */ \
	blk = _last_phantom(blk)->blk; \
	/* reload dir and mask pointer, adjust path offset */ \
	uint64_t _idx = blk->acnt + blk->bcnt; \
	mask = &blk->mask[_idx]; \
	dir_mask = _trace_load_block_rem(_idx); \
}

/**
 * @macro _trace_reload_block
 * @brief reload block for bulk trace
 */
#define _trace_reload_block() { \
	mask = &(--blk)->mask[BLK - 1]; \
	dir_mask = _dir_load(blk, BLK).mask; \
}

/**
 * @macro _trace_store_path
 * @brief store path array
 */
#define _trace_store_path() { \
	debug("path_array(%llx)", path_array); \
	_storeu_u64(path, path_array<<ofs); \
	path--; \
}

/**
 * @macro _trace_test_bulk
 * @brief update gidx and test if section boundary can be found in the next block.
 * adjust gidx to the next base and return 0 if bulk loop can be continued,
 * otherwise gidx will not updated. must be called just after reload_ptr or reload_tail.
 */
#define _trace_test_bulk() { \
	v2i8_t cnt = _load_v2i8(&blk->acnt); \
	v2i32_t _gidx = _sub_v2i32(gidx, _cvt_v2i8_v2i32(cnt)); \
	_test_v2i32(_gt_v2i32(_gidx, m00), m11) ? (gidx = _gidx, 1) : 0; \
}

/**
 * @macro _trace_*_load
 */
#define _trace_bulk_load_n(t, _jump_to) { \
	if(_unlikely(mask < blk->mask)) { \
		_trace_store_path();					/* store path_array */ \
		_trace_reload_block(); \
		if(_unlikely(_trace_test_bulk())) {		/* adjust gidx */ \
			debug("jump to %s", #_jump_to); \
			goto _jump_to;						/* jump to tail loop */ \
		} \
	} \
}
#define _trace_tail_load_n(t, _jump_to) { \
	if(_unlikely(mask < blk->mask)) { \
		_trace_store_path();					/* store path_array */ \
		if(_last_phantom(blk)->root + 1 < 2) {	/* head (phantom) block is marked mask == -1 or 0 */ \
			_trace_reload_tail(t);				/* fetch the previous tail */ \
			if(!_trace_test_bulk()) {			/* adjust gidx */ \
				goto _jump_to;					/* dispatch again (bulk loop) */ \
			} \
		} else { \
			/* load dir and update mask pointer */ \
			_trace_reload_block();				/* not reached the head yet */ \
		} \
	} \
}

/**
 * @fn trace_body
 */
static _force_inline
void trace_body(
	struct gaba_dp_context_s *self)
{
	#define _trace_gap_loop(t, _type, _next, _label) { \
		_trace_forward_##_type##_##_label##_head: \
		_trace_inc_gi();		/* increment #gap regions at the head */ \
		while(1) { \
			if(_trace_test_gap_##_label() == 0) { \
				goto _trace_##_type##_d_head; \
			} \
			if(_trace_##_type##_##_label##_test_index()) { \
				goto _trace_index_break; \
			} \
			_trace_inc_ge();	/* increment #gap bases on every iter */ \
			debug("go %s (%s), dir(%llx), mask_h(%x), mask_v(%x), p(%lld), q(%lld), ptr(%p), path_array(%llx)", \
				#_label, #_type, ((uint64_t)dir.dynamic.array), ptr->pair.h.all, ptr->pair.v.all, p, q, ptr, path_array); \
			_trace_##_type##_##_label##_update_index(); \
			_trace_##_label##_update_path_q(); \
			_trace_##_type##_load_n(t, _trace_##_next##_##_label##_head); \
		} \
	}

	#define _trace_diag_loop(t, _type, _next) { \
		while(1) { \
		_trace_##_type##_d_head: \
			if(_trace_test_diag_h() != 0) { \
				goto _trace_##_type##_h_head; \
			} \
			if(_trace_##_type##_d_test_index()) { \
				goto _trace_index_break; \
			} \
			debug("go d (%s), dir(%llx), mask_h(%x), mask_v(%x), p(%lld), q(%lld), ptr(%p), path_array(%llx)", \
				#_type, ((uint64_t)dir.dynamic.array), ptr->pair.h.all, ptr->pair.v.all, p, q, ptr, path_array); \
			_trace_##_type##_h_update_index(); \
			_trace_h_update_path_q(); \
			_trace_##_type##_load_n(t, _trace_##_next##_d_mid); \
		_trace_##_type##_d_mid: \
			_trace_##_type##_v_update_index(); \
			_trace_v_update_path_q(); \
			_trace_##_type##_load_n(t, _trace_##_next##_d_tail); \
		_trace_##_type##_d_tail: \
			if(_trace_test_diag_v() != 0) { \
				goto _trace_##_type##_v_head; \
			} \
		} \
	}

	/* constants for agidx and bgidx end detection */
	static uint32_t const c00[2] __attribute__(( aligned(16) )) = { 0, 0 };
	static uint32_t const c01[2] __attribute__(( aligned(16) )) = { -1, 0 };
	static uint32_t const c10[2] __attribute__(( aligned(16) )) = { 0, -1 };
	static uint32_t const c11[2] __attribute__(( aligned(16) )) = { -1, -1 };
	v2i32_t const m00 = _load_v2i32(c00);
	v2i32_t const m01 = _load_v2i32(c01);
	v2i32_t const m10 = _load_v2i32(c10);
	v2i32_t const m11 = _load_v2i32(c11);

	/* gap counts; #gap regions in lower and #gap bases in higher */
	register v2i32_t gc = _load_v2i32(&self->w.l.a.gicnt);

	/* load path array, adjust path offset to align the head of the current block */
	uint32_t *path = self->w.l.path;
	uint64_t ofs = self->w.l.ofs;	/* global p-coordinate */
	uint64_t path_array = _loadu_u64(path)>>ofs;

	/* load pointers and coordinates */
	struct gaba_block_s const *blk = self->w.l.blk;
	struct gaba_mask_pair_s const *mask = &blk->mask[self->w.l.p];
	uint64_t q = self->w.l.q;
	uint32_t dir_mask = _trace_load_block_rem(self->w.l.p);

	/* load grid indices; assigned for each of N+1 boundaries of length-N sequence */
	register v2i32_t gidx = _load_v2i32(&self->w.l.agidx);

	_trace_dec_gi();				/* compensate gap count for the first dispatch */
	if(_trace_test_bulk()) {		/* test if the first loop contains section boundary */
		/* bulk loops; v is the beginning state */
		_trace_gap_loop(self, bulk, tail, v);
		_trace_diag_loop(self, bulk, tail);
		_trace_gap_loop(self, bulk, tail, h);
	} else {
		/* tail loops */
		_trace_gap_loop(self, tail, bulk, v);
		_trace_diag_loop(self, tail, bulk);
		_trace_gap_loop(self, tail, bulk, h);
	}

_trace_index_break:;
	/* reached a boundary of sections, compensate ofs */
	uint64_t rem = mask - blk->mask + 1;
	path += ofs + rem >= BLK;
	ofs = (ofs + rem) & (BLK - 1);
	_storeu_u64(path, path_array<<ofs);

	/* save gap counts */
	_store_v2i32(&self->w.l.path.gic, gc);

	/* store pointers and coordinates */
	self->w.l.path = path;
	self->w.l.ofs = ofs;			/* global p-coordinate */
	self->w.l.blk = blk;
	self->w.l.p = mask - blk->mask;	/* local p-coordinate */
	self->w.l.q = q;				/* q coordinate */
	_store_v2i32(&self->w.l.agidx, gidx);
	debug("p(%lld), ppos(%lld), q(%llu)", p, self->w.l.ppos, q);
	return;
}

/**
 * @fn trace_push_section
 */
static _force_inline
void trace_push_section(
	struct gaba_dp_context_s *self)
{
	/* windback pointer */
	self->w.l.sec_len++;
	self->w.l.sec--;

	/* calc ppos */
	uint64_t ppos = (self->w.l.path - self->w.l.aln->path) * 32 + self->w.l.ofs;

	/* load section info */
	v2i32_t gidx = _load_v2i32(&self->w.l.agidx);
	v2i32_t sgidx = _load_v2i32(&self->w.l.asgidx);
	v2i32_t id = _load_v2i32(&self->w.l.aid);

	/* store section info */
	_store_v2i32(&self->w.l.sec->apos, gidx);
	_store_v2i32(&self->w.l.sec->alen, _sub_v2i32(sgidx, gidx));
	_store_v2i32(&self->w.l.sec->aid, id);
	self->w.l.sec->ppos = ppos;

	/* update rsgidx */
	_store_v2i32(&self->w.l.asgidx, gidx);
	return;
}

/**
 * @fn trace_init_alignment
 */
static _force_inline
void trace_init_alignment(
	struct gaba_dp_context_s *self,
	struct gaba_leaf_s const *leaf,
	struct gaba_allocator_s const *alloc)
{
	struct gaba_joint_tail_s const *tail = leaf->tail;

	/* malloc container */
	uint64_t slen = 2 * tail->scnt, plen = _roundup((tail->f.ppos + 31) / 32 + 1, 8);
	uint64_t size = (
		  sizeof(struct gaba_alignment_s)				/* base */
		+ sizeof(uint32_t) * plen						/* path array and its margin */
		+ sizeof(struct gaba_path_section_s) * slen		/* section array */
		+ self->head_margin + self->tail_margin			/* head and tail margins */
	);

	/* save aln pointer to working buffer */
	self->w.l.aln = alloc->malloc(alloc->opaque, size) + self->head_margin;
	debug("malloc trace mem(%p), opaque(%p)", aln, opaque);

	/* use gaba_alignment_s buffer instead in the traceback loop */
	self->w.l.a = (struct gaba_alignment_s){
		/* internals */
		.opaque = alloc->opaque,						/* save opaque pointer */
		.free = alloc->free,
		.head_margin = self->head_margin,				/* required when free the object */

		/* score and path element counts */
		.score = tail->max,								/* just copy */
		.mcnt = 0,
		.xcnt = 0,
		.gicnt = 0,
		.gecnt = 0,

		/* section */
		.sec_len = 0,
		.sec = (struct gaba_path_section_s *)(self->w.l.aln->path + plen, 8) + slen - 1,

		/* path */
		.path_len = leaf->ppos
	};

	/* clear array */
	self->w.l.a.path[0] = 0;
	self->w.l.a.path[1] = 0;

	/* store block and coordinates */
	self->w.l.path = aln->path + plen - 1;
	self->w.l.ofs = self->w.l.a.path_len & (BLK - 1);
	self->w.l.blk = leaf->blk;
	self->w.l.p = leaf->p;
	self->w.l.q = leaf->q;

	/* increment by one to convert char index to grid index */
	v2i32_t gidx = _sub_v2i32(
		_set_v2i32(1),
		_load_v2i32(&leaf->aridx)
	);
	_store_v2i32(&self->w.l.agidx, gidx);				/* current gidx */
	_store_v2i32(&self->w.l.asgidx, gidx);				/* tail gidx (pos) */
	_store_v2i32(&self->w.l.aid, _zero_v2i32());		/* invalid ids, must be loaded afterward */

	/* store tail pointers for sequence segment tracking */
	self->w.l.atail = leaf->tail;
	self->w.l.btail = leaf->tail;
	return;
}

/**
 * @fn trace_generate_alignment
 */
static _force_inline
struct gaba_alignment_s *trace_generate_alignment(
	struct gaba_dp_context_s *self,
	struct gaba_leaf_s const *leaf)
{
	/* blockwise traceback loop, until ppos reaches the root */
	while( _last_phantom(self->w.l.blk)->root != -1
		|| self->w.l.mask != &self->w.l.blk->mask[1]
	) {
		/* update section info */
		if(self->w.l.agidx <= 0) { trace_reload_section(self, 0); }
		if(self->w.l.bgidx <= 0) { trace_reload_section(self, 1); }

		/* fragment trace: q must be inside [0, BW) */
		trace_body(self);
		debug("p(%d), ppos(%lld), q(%d)", self->w.l.p, self->w.l.ppos, self->w.l.q);
		if(_unlikely(self->w.l.q >= BW)) {
			/* out of band: abort */
			self->w.l.a.free(
				self->w.l.a.opaque,
				(void *)((uint8_t *)self->w.l.aln - self->w.l.a.head_margin)
			);
			return(NULL);
		}

		/* push section info to section array */
		trace_push_section(self);
	}

	/* calc mismatch and match counts */
	self->w.l.a.xcnt = (
		(
			  self->m * ((self->w.l.a.path_len - self->w.l.a.gecnt)>>1)
			- (self->w.l.a.score - self->gi * self->w.l.a.gicnt - self->ge * gecnt)
		) / (m - x)
	);
	self->w.l.a.mcnt = ((self->w.l.a.path_len - self->w.l.a.gecnt)>>1) - self->w.l.a.xcnt;
	debug("plen(%lld), gic(%lld), gec(%lld), dcnt(%lld), xcnt(%lld)",
		self->w.l.a.path->len, gic, gec, (self->w.l.a.path->len - gec)>>1, self->w.l.a.xcnt);

	/* copy */
	_memcpy_blk_aa(self->w.l.aln, &self->w.l.a, sizeof(struct gaba_alignment_s));
	return(self->w.l.aln);
}

/**
 * @fn gaba_dp_trace
 */
struct gaba_alignment_s *_export(gaba_dp_trace)(
	struct gaba_dp_context_s *self,
	struct gaba_fill_s const *fill,
	struct gaba_allocator_s const *alloc)
{
	/* restore dp context pointer by adding offset */
	self = _restore_dp_context(self);

	/* restore default alloc if NULL */
	struct gaba_allocator_s const default_alloc = {
		.opaque = (void *)self,
		.malloc = (gaba_malloc_t)gaba_dp_malloc,
		.free = (gaba_free_t)gaba_dp_free
	};
	alloc = (alloc == NULL) ? &default_alloc : alloc;

	/* search */
	struct gaba_leaf_s leaf;
	leaf_search(self, _tail(fill), &leaf);

	/* create alignment object */
	trace_init_alignment(self, &leaf, alloc);

	/* generate paths, may fail when path got lost out of the band */
	return(trace_generate_alignment(self));
}


/**
 * @fn gaba_dp_res_free
 */
void _export(gaba_dp_res_free)(
	struct gaba_alignment_s *aln)
{
	debug("free mem, ptr(%p), lmm(%p)", (void *)aln - aln->head_margin, lmm);
	aln->free(aln->opaque, (void *)((uint8_t *)aln - aln->head_margin));
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
uint64_t _export(gaba_dp_print_cigar_forward)(
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
			ZCNT_RESULT uint64_t c = a;
			if(c < 64) { break; }

			debug("bulk match, a(%llu), pos(%llu), ridx(%llu)", a, lim - ridx, ridx);
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
			debug("gap g(%lld), pos(%llu), ridx(%llu)", g, lim - ridx, ridx);
		}
		if((ridx -= g) <= 1) { break; }
	}
	return(clen);
}

/**
 * @fn gaba_dp_dump_cigar_forward
 * @brief parse path string and store cigar to buffer
 */
uint64_t _export(gaba_dp_dump_cigar_forward)(
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
			ZCNT_RESULT uint64_t c = a;
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
uint64_t _export(gaba_dp_print_cigar_reverse)(
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
uint64_t _export(gaba_dp_dump_cigar_reverse)(
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
	struct gaba_score_vec_s sc __attribute__(( aligned(32) ));

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
void gaba_init_fill_middle_delta(
	struct gaba_middle_delta_s *md,
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

	for(int i = 0; i < BW/2; i++) {
		md->delta[i] = ofs + coef * (BW/2 - i);
		md->delta[BW/2 + i] = ofs + coef * i;
	}
	md->delta[BW/2] = 0;
	return;
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

	_store_n(&diff.dh, _load_n(dh));
	_store_n(&diff.dv, _load_n(dv));
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

 	nvec_t _dh = _load_n(dh);
 	nvec_t _dv = _load_n(dv);
 	nvec_t _de = _load_n(de);
 	nvec_t _df = _load_n(df);

 	_print_n(_dh);
 	_print_n(_dv);
 	_print_n(_de);
 	_print_n(_df);

	_dh = _shl_n(_dh, 3);
	_dv = _shl_n(_dv, 3);
	_store_n(&diff.dh, _add_n(_dh, _de));
	_store_n(&diff.dv, _add_n(_dv, _df));
	_print_n(_add_n(_dh, _de));
	_print_n(_add_n(_dv, _df));

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
 * @fn gaba_init_fill_phantom
 * @brief phantom block at root
 */
static _force_inline
void gaba_init_fill_phantom(
	struct gaba_root_block_s *ph,
	struct gaba_params_s *params_intl,
	struct gaba_middle_delta_box_s *md)
{
	/* 192 bytes are reserved for phantom block */
	*_last_phantom(&ph->tail) = (struct gaba_phantom_s){
		.root = -1,
		.acnt = 0,
		.bcnt = 0,
		.tail = &ph->tail
	};

	/* fill root tail object */
	ph->tail = (struct gaba_joint_tail_s){
		/* coordinates */
		.ppos = GP_INIT - BW,
		.p = 0,
		.scnt = 0,

		/* max */
		.max = 0,

		/* status */
		.stat = CONT,

		/* internals */
		.tail = NULL,
		.apos = 0,
		.bpos = 0,
		.alen = 0,
		.blen = 0,
		.aid = 0xfffc,
		.bid = 0xfffd
	};
	return;
}

/**
 * @fn gaba_init_fill_dp_context
 */
static _force_inline
void gaba_init_fill_dp_context(
	struct gaba_dp_context_s *dp,
	struct gaba_params_s *params_intl)
{
	*dp = (struct gaba_dp_context_s){
		/* memory management */
		.mem = { 0 },
		.stack = { 0 },

		/* score vectors */
		.scv = gaba_init_create_score_vector(params_intl),
		.m = params_intl->m,
		.x = -params_intl->x,
		.gi = (MODEL == LINEAR)
			? 0
			: -params_intl->gi,
		.ge = (MODEL == LINEAR)
			? -(params_intl->gi + params_intl->ge)
			: -params_intl->ge,
		.gf = -params_intl->gf,
		.tx = params_intl->xdrop,
		.tf = params_intl->filter_thresh,

		/* input and output options */
		.head_margin = _roundup(params_intl->head_margin, MEM_ALIGN_SIZE),
		.tail_margin = _roundup(params_intl->tail_margin, MEM_ALIGN_SIZE)
	};
	return;
}

/**
 * @fn gaba_init
 */
gaba_t *_export(gaba_init)(
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
	struct gaba_context_s *ctx = NULL;
	if(params_intl.reserved == NULL) {
		if((ctx = gaba_malloc(sizeof(struct gaba_context_s))) == NULL) {
			return(NULL);
		}
		gaba_init_fill_dp_context(&ctx->k, &params_intl);
	} else {
		/* fill phantom objects of existing dp context */
		ctx = (struct gaba_context_s *)params_intl.reserved;
	}

	/* load default phantom objects */
	gaba_init_fill_phantom(_proot(&ctx->k), &params_intl, _pmd(ctx));
	gaba_init_fill_middle_delta((struct gaba_middle_delta_s *)_pmd(ctx), &params_intl);
	return((gaba_t *)ctx);
}

/**
 * @fn gaba_clean
 */
void _export(gaba_clean)(
	struct gaba_context_s *ctx)
{
	gaba_free(ctx);
	return;
}

/**
 * @fn gaba_dp_init
 */
struct gaba_dp_context_s *_export(gaba_dp_init)(
	struct gaba_context_s const *ctx,
	uint8_t const *alim,
	uint8_t const *blim)
{
	/* malloc stack memory */
	struct gaba_dp_context_s *self = gaba_malloc(sizeof(gaba_dp_context_s) + MEM_INIT_SIZE);
	if(self == NULL) {
		debug("failed to malloc memory");
		return(NULL);
	}

	/* add offset */
	self = _restore_dp_context_global(self);

	/* init stack pointers */
	self->stack.mem = &self->mem;
	self->stack.top = (uint8_t *)(self + 1);
	self->stack.end = (uint8_t *)self + MEM_INIT_SIZE - MEM_MARGIN_SIZE;

	/* init mem object */
	self->mem = (struct gaba_mem_block_s){
		.next = NULL,
		.prev = NULL,
		.size = MEM_INIT_SIZE
	};

	/* init seq lims */
	self->w.alim = alim;
	self->w.blim = blim;

	/* copy template */
	_memcpy_blk_aa(
		(uint8_t *)self + GABA_DP_CONTEXT_LOAD_OFFSET,
		(uint8_t *)&ctx->k + GABA_DP_CONTEXT_LOAD_OFFSET,
		GABA_DP_CONTEXT_LOAD_SIZE
	);

	/* return offsetted pointer */
	return(_export_dp_context(self));
}

/**
 * @fn gaba_dp_add_stack
 * @brief returns zero when succeeded
 */
static _force_inline
int32_t gaba_dp_add_stack(
	struct gaba_dp_context_s *self,
	uint64_t size)
{
	if(self->curr_mem->next == NULL) {
		/* add new block */
		uint64_t next_size = self->curr_mem->size * 2;
		struct gaba_mem_block_s *mem = gaba_malloc(next_size);
		if(mem == NULL) { return(GABA_ERROR_OUT_OF_MEM); }
		self->curr_mem->next = mem;

		mem->next = NULL;
		mem->prev = self->curr_mem;
		mem->size = next_size;
	}

	/* follow the forward link */
	self->curr_mem = self->curr_mem->next;

	/* init stack pointers */
	self->stack_top = (uint8_t *)_roundup((uintptr_t)(self->curr_mem + 1), MEM_ALIGN_SIZE);
	self->stack_end = (uint8_t *)self->curr_mem + self->curr_mem->size - MEM_MARGIN_SIZE;
	return(0);
}

/**
 * @fn gaba_dp_flush
 */
void _export(gaba_dp_flush)(
	struct gaba_dp_context_s *self,
	uint8_t const *alim,
	uint8_t const *blim)
{
	/* restore dp context pointer by adding offset */
	self = _restore_dp_context(self);

	/* init seq lims */
	self->w.alim = alim;
	self->w.blim = blim;

	self->curr_mem = &self->mem;
	self->stack_top = (uint8_t *)_roundup((uintptr_t)(self->curr_mem + 1), MEM_ALIGN_SIZE);
	self->stack_end = (uint8_t *)self + self->mem.size - MEM_MARGIN_SIZE;
	return;
}

/**
 * @fn gaba_dp_save_stack
 */
struct gaba_stack_s const *_export(gaba_dp_save_stack)(
	struct gaba_dp_context_s *self)
{
	self = _restore_dp_context(self);

	struct gaba_mem_block_s *mem = self->curr_mem;
	uint8_t *stack_top = self->stack_top;
	uint8_t *stack_end = self->stack_end;

	debug("save stack, self(%p), mem(%p, %p, %p)", self, self->curr_mem, self->stack_top, self->stack_end);

	/* save */
	struct gaba_stack_s *stack = gaba_dp_malloc(self, sizeof(struct gaba_stack_s));
	stack->mem = mem;
	stack->stack_top = stack_top;
	stack->stack_end = stack_end;
	debug("save stack(%p), (%p, %p, %p)", stack, stack->mem, stack->stack_top, stack->stack_end);
	return(stack);
}

/**
 * @fn gaba_dp_flush_stack
 */
void _export(gaba_dp_flush_stack)(
	struct gaba_dp_context_s *self,
	gaba_stack_t const *stack)
{
	if(stack == NULL) {
		return;
	}

	self = _restore_dp_context(self);
	self->curr_mem = stack->mem;
	self->stack_top = stack->stack_top;
	self->stack_end = stack->stack_end;
	debug("restore stack, self(%p), mem(%p, %p, %p)", self, stack->mem, stack->stack_top, stack->stack_end);
	return;
}

/**
 * @fn gaba_dp_malloc
 */
static _force_inline
void *gaba_dp_malloc(
	struct gaba_dp_context_s *self,
	uint64_t size)
{
	/* roundup */
	size = _roundup(size, MEM_ALIGN_SIZE);

	/* malloc */
	debug("self(%p), stack_top(%p), size(%llu)", self, self->stack_top, size);
	if((uint64_t)(self->stack_end - self->stack_top) < size) {
		if(gaba_dp_add_stack(self, size) != 0) {
			return(NULL);
		}
		debug("stack_top(%p)", self->stack_top);
	}
	self->stack_top += size;
	return((void *)(self->stack_top - size));
}

/**
 * @fn gaba_dp_free
 * @brief do nothing
 */
static _force_inline
void gaba_dp_free(
	struct gaba_dp_context_s *self,
	void *ptr)
{
	return;
}

/**
 * @fn gaba_dp_clean
 */
void _export(gaba_dp_clean)(
	struct gaba_dp_context_s *self)
{
	if(self == NULL) {
		return;
	}

	/* restore dp context pointer by adding offset */
	self = _restore_dp_context(self);

	struct gaba_mem_block_s *m = self->mem.next;
	while(m != NULL) {
		struct gaba_mem_block_s *mnext = m->next;
		free(m); m = mnext;
	}
	gaba_free(_export_dp_context_global(self));
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
	gaba_t *ctx = _export(gaba_init)(unittest_default_params);
	return((void *)ctx);
}

/**
 * @fn unittest_clean_context
 */
static
void unittest_clean_context(void *ctx)
{
	_export(gaba_clean)((struct gaba_context_s *)ctx);
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
#define check_tail(_t, _max, _p, _psum, _scnt) ( \
	   (_t) != NULL \
	&& (_t)->max == (_max) \
	&& (_t)->p == (_p) \
	&& (_t)->ppos == (_psum) \
	&& (_t)->scnt == (_scnt) \
)
#define print_tail(_t) \
	"tail(%p), max(%lld), p(%d), ppos(%lld), scnt(%u)", \
	(_t), (_t)->max, (_t)->p, (_t)->ppos, (_t)->scnt
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

	uint64_t l = _export(gaba_dp_dump_cigar_forward)(
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
	struct gaba_dp_context_s *d = _export(gaba_dp_init)(c, s->alim, s->blim);

unittest(with_seq_pair("A", "A"))
{
	omajinai();

	assert(d != NULL, "%p", d);
	_export(gaba_dp_clean)(d);
}

/**
 * check if gaba_dp_fill_root and gaba_dp_fill returns a correct score
 */
unittest(with_seq_pair("A", "A"))
{
	omajinai();

	/* fill root section */
	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 0, 0, -29, 1), print_tail(f));

	/* fill again */
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 0, 0, -27, 2), print_tail(f));

	/* fill tail section */
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 4, 13, 13, 3), print_tail(f));

	/* fill tail section again */
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);
	#if MODEL == LINEAR
		assert(f->status == 0x1ff, "%x", f->status);
		assert(check_tail(f, 4, 40, 53, 4), print_tail(f));
	#else
		assert(f->status == 0x10f, "%x", f->status);
		assert(check_tail(f, 4, 31, 44, 4), print_tail(f));
	#endif

	_export(gaba_dp_clean)(d);
}

/* reverse fetch */
unittest(with_seq_pair("A", "A"))
{
	omajinai();

	/* fill root section */
	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->arsec, 0, &s->brsec, 0);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 0, 0, -29, 1), print_tail(f));

	/* fill again */
	f = _export(gaba_dp_fill)(d, f, &s->arsec, &s->brsec);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 0, 0, -27, 2), print_tail(f));

	/* fill tail section */
	f = _export(gaba_dp_fill)(d, f, &s->artail, &s->brtail);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 4, 13, 13, 3), print_tail(f));

	/* fill tail section again */
	f = _export(gaba_dp_fill)(d, f, &s->artail, &s->brtail);
	#if MODEL == LINEAR
		assert(f->status == 0x1ff, "%x", f->status);
		assert(check_tail(f, 4, 40, 53, 4), print_tail(f));
	#else
		assert(f->status == 0x10f, "%x", f->status);
		assert(check_tail(f, 4, 31, 44, 4), print_tail(f));
	#endif

	_export(gaba_dp_clean)(d);
}

/* with longer sequences */
unittest(with_seq_pair("ACGTACGTACGT", "ACGTACGTACGT"))
{
	omajinai();

	/* fill root section */
	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 0, 0, -7, 1), print_tail(f));

	/* fill again */
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 16, 17, 17, 2), print_tail(f));

	/* fill tail section */
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 48, 40, 57, 3), print_tail(f));

	/* fill tail section again */
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);
	#if MODEL == LINEAR
		assert(f->status == 0x1ff, "%x", f->status);
		assert(check_tail(f, 48, 40, 97, 4), print_tail(f));
	#else
		assert(f->status == 0x10f, "%x", f->status);
		assert(check_tail(f, 48, 31, 88, 4), print_tail(f));
	#endif

	_export(gaba_dp_clean)(d);
}

/* reverse fetch */
unittest(with_seq_pair("ACGTACGTACGT", "ACGTACGTACGT"))
{
	omajinai();

	/* fill root section */
	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->arsec, 0, &s->brsec, 0);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 0, 0, -7, 1), print_tail(f));

	/* fill again */
	f = _export(gaba_dp_fill)(d, f, &s->arsec, &s->brsec);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 16, 17, 17, 2), print_tail(f));

	/* fill tail section */
	f = _export(gaba_dp_fill)(d, f, &s->artail, &s->brtail);
	assert(f->status == 0x1ff, "%x", f->status);
	assert(check_tail(f, 48, 40, 57, 3), print_tail(f));

	/* fill tail section again */
	f = _export(gaba_dp_fill)(d, f, &s->artail, &s->brtail);
	#if MODEL == LINEAR
		assert(f->status == 0x1ff, "%x", f->status);
		assert(check_tail(f, 48, 40, 97, 4), print_tail(f));
	#else
		assert(f->status == 0x10f, "%x", f->status);
		assert(check_tail(f, 48, 31, 88, 4), print_tail(f));
	#endif

	_export(gaba_dp_clean)(d);
}

/* sequences with different lengths (consumed as mismatches) */
unittest(with_seq_pair("GAAAAAAAA", "AAAAAAAA"))
{
	omajinai();

	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	assert(f->status == 0x01ff, "%x", f->status);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x01f0, "%x", f->status);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bftail);
	assert(f->status == 0x010f, "%x", f->status);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);
	assert(f->status == 0x01f0, "%x", f->status);
	assert(check_tail(f, 22, 37, 42, 4), print_tail(f));

	_export(gaba_dp_clean)(d);
}

/* another pair with different lengths */
unittest(with_seq_pair("TTTTTTTT", "CTTTTTTTT"))
{
	omajinai();

	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	assert(f->status == 0x010f, "%x", f->status);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x01f0, "%x", f->status);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x010f, "%x", f->status);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bfsec);
	assert(f->status == 0x01f0, "%x", f->status);

	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);
	#if MODEL == LINEAR
		assert(f->status == 0x010f, "%x", f->status);
		assert(check_tail(f, 22, 36, 42, 5), print_tail(f));
	#else
		assert(f->status == 0x010f, "%x", f->status);
		assert(check_tail(f, 22, 35, 41, 5), print_tail(f));
	#endif

	_export(gaba_dp_clean)(d);
}

/* with deletions */
unittest(with_seq_pair("GACGTACGT", "ACGTACGT"))
{
	omajinai();

	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	assert(f->status == 0x01ff, "%x", f->status);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x01f0, "%x", f->status);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bftail);
	assert(f->status == 0x010f, "%x", f->status);

	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);
	#if MODEL == LINEAR
		assert(f->status == 0x01f0, "%x", f->status);
		assert(check_tail(f, 20, 37, 42, 4), print_tail(f));
	#else
		assert(f->status == 0x01ff, "%x", f->status);
		assert(check_tail(f, 20, 38, 43, 4), print_tail(f));
	#endif

	_export(gaba_dp_clean)(d);
}

/* with insertions */
unittest(with_seq_pair("ACGTACGT", "GACGTACGT"))
{
	omajinai();

	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	assert(f->status == 0x010f, "%x", f->status);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x01f0, "%x", f->status);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	assert(f->status == 0x010f, "%x", f->status);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bfsec);
	assert(f->status == 0x01f0, "%x", f->status);


	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);
	#if MODEL == LINEAR
		assert(f->status == 0x010f, "%x", f->status);
		assert(check_tail(f, 20, 35, 41, 5), print_tail(f));
	#else
		assert(f->status == 0x010f, "%x", f->status);
		assert(check_tail(f, 20, 36, 42, 5), print_tail(f));
	#endif

	_export(gaba_dp_clean)(d);
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
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0x55555555), 0, 32);
	assert(strcmp(buf, "16M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0x55555555, 0x55555555), 0, 64);
	assert(strcmp(buf, "32M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0x55555555, 0x55555555, 0x55555555, 0x55555555), 0, 128);
	assert(strcmp(buf, "64M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0x55550000, 0x55555555, 0x55555555, 0x55555555), 16, 112);
	assert(strcmp(buf, "56M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0x55555000, 0x55555555, 0x55555555, 0x55555555), 12, 116);
	assert(strcmp(buf, "58M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0x55), 0, 8);
	assert(strcmp(buf, "4M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0x55555000, 0x55555555, 0x55555555, 0x55), 12, 92);
	assert(strcmp(buf, "46M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0x55550555), 0, 32);
	assert(strcmp(buf, "6M4D8M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0x5555f555), 0, 32);
	assert(strcmp(buf, "6M4I8M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0xaaaa0555), 0, 33);
	assert(strcmp(buf, "6M5D8M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0xaaabf555), 0, 33);
	assert(strcmp(buf, "6M5I8M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0xaaabf555, 0xaaaa0556), 0, 65);
	assert(strcmp(buf, "6M5I8M1I5M5D8M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0xaaabf555, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "6M5I8M1I5M5D8M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_forward)(ut_printer, (void *)&p, _arr(0xaaabf554, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "2D5M5I8M1I5M5D8M") == 0, "%s", buf);

	#undef _arr
	free(buf);
}

unittest()
{
	char *buf = (char *)malloc(16384);
	char *p = buf;

	#define _arr(...)		( (uint32_t const []){ 0, 0, __VA_ARGS__, 0, 0 } + 2 )
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0x55555555), 0, 32);
	assert(strcmp(buf, "16M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0x55555555, 0x55555555), 0, 64);
	assert(strcmp(buf, "32M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0x55555555, 0x55555555, 0x55555555, 0x55555555), 0, 128);
	assert(strcmp(buf, "64M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0x55550000, 0x55555555, 0x55555555, 0x55555555), 16, 112);
	assert(strcmp(buf, "56M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0x55555000, 0x55555555, 0x55555555, 0x55555555), 12, 116);
	assert(strcmp(buf, "58M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0x55), 0, 8);
	assert(strcmp(buf, "4M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0x55555000, 0x55555555, 0x55555555, 0x55), 12, 92);
	assert(strcmp(buf, "46M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0x55550555), 0, 32);
	assert(strcmp(buf, "8M4D6M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0x5555f555), 0, 32);
	assert(strcmp(buf, "8M4I6M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0xaaaa0555), 0, 33);
	assert(strcmp(buf, "8M5D6M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0xaaabf555), 0, 33);
	assert(strcmp(buf, "8M5I6M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0xaaabf555, 0xaaaa0556), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I6M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0xaaabf555, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I6M") == 0, "%s", buf);

	p = buf;
	_export(gaba_dp_print_cigar_reverse)(ut_printer, (void *)&p, _arr(0xaaabf554, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
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
	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0x55555555), 0, 32);
	assert(strcmp(buf, "16M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0x55555555, 0x55555555), 0, 64);
	assert(strcmp(buf, "32M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0x55555555, 0x55555555, 0x55555555, 0x55555555), 0, 128);
	assert(strcmp(buf, "64M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0x55550000, 0x55555555, 0x55555555, 0x55555555), 16, 112);
	assert(strcmp(buf, "56M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0x55555000, 0x55555555, 0x55555555, 0x55555555), 12, 116);
	assert(strcmp(buf, "58M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0x55), 0, 8);
	assert(strcmp(buf, "4M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0x55555000, 0x55555555, 0x55555555, 0x55), 12, 92);
	assert(strcmp(buf, "46M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0x55550555), 0, 32);
	assert(strcmp(buf, "6M4D8M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0x5555f555), 0, 32);
	assert(strcmp(buf, "6M4I8M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0xaaaa0555), 0, 33);
	assert(strcmp(buf, "6M5D8M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0xaaabf555), 0, 33);
	assert(strcmp(buf, "6M5I8M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0xaaabf555, 0xaaaa0556), 0, 65);
	assert(strcmp(buf, "6M5I8M1I5M5D8M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0xaaabf555, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "6M5I8M1I5M5D8M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_forward)(buf, len, _arr(0xaaabf554, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "2D5M5I8M1I5M5D8M") == 0, "%s", buf);

	#undef _arr
	free(buf);
}

unittest()
{
	uint64_t const len = 16384;
	char *buf = (char *)malloc(len);

	#define _arr(...)		( (uint32_t const []){ 0, 0, __VA_ARGS__, 0, 0 } + 2 )
	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0x55555555), 0, 32);
	assert(strcmp(buf, "16M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0x55555555, 0x55555555), 0, 64);
	assert(strcmp(buf, "32M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0x55555555, 0x55555555, 0x55555555, 0x55555555), 0, 128);
	assert(strcmp(buf, "64M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0x55550000, 0x55555555, 0x55555555, 0x55555555), 16, 112);
	assert(strcmp(buf, "56M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0x55555000, 0x55555555, 0x55555555, 0x55555555), 12, 116);
	assert(strcmp(buf, "58M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0x55), 0, 8);
	assert(strcmp(buf, "4M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0x55555000, 0x55555555, 0x55555555, 0x55), 12, 92);
	assert(strcmp(buf, "46M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0x55550555), 0, 32);
	assert(strcmp(buf, "8M4D6M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0x5555f555), 0, 32);
	assert(strcmp(buf, "8M4I6M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0xaaaa0555), 0, 33);
	assert(strcmp(buf, "8M5D6M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0xaaabf555), 0, 33);
	assert(strcmp(buf, "8M5I6M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0xaaabf555, 0xaaaa0556), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I6M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0xaaabf555, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
	assert(strcmp(buf, "8M5D5M1I8M5I6M") == 0, "%s", buf);

	_export(gaba_dp_dump_cigar_reverse)(buf, len, _arr(0xaaabf554, 0xaaaa0556, 0xaaaaaaaa), 0, 65);
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
	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);

	/* forward-only traceback */
	struct gaba_alignment_s *r = _export(gaba_dp_trace)(d, f, NULL, NULL);
	assert(check_result(r, 0, 0, 0, 0, (uint32_t)-1, 0, 0, 0), print_result(r));

	/* forward-reverse traceback */
	r = _export(gaba_dp_trace)(d, f, f, NULL);
	assert(check_result(r, 0, 0, 0, 0, (uint32_t)-1, 0, 0, 0), print_result(r));

	/* section added */
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);

	/* forward-only traceback */
	r = _export(gaba_dp_trace)(d, f, NULL, NULL);
	assert(check_result(r, 0, 0, 0, 0, (uint32_t)-1, 0, 0, 0), print_result(r));

	/* forward-reverse traceback */
	r = _export(gaba_dp_trace)(d, f, f, NULL);
	assert(check_result(r, 0, 0, 0, 0, (uint32_t)-1, 0, 0, 0), print_result(r));

	_export(gaba_dp_clean)(d);
}

/* with short sequences */
unittest(with_seq_pair("A", "A"))
{
	omajinai();

	/* fill sections */
	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);

	/* forward-only traceback */
	struct gaba_alignment_s *r = _export(gaba_dp_trace)(d, f, NULL, NULL);
	assert(check_result(r, 4, 0, 4, 2, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "DRDR"), print_path(r));
	assert(check_cigar(r, "2M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 1, s->bfsec, 0, 1, 0, 2), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 0, 1, s->bfsec, 0, 1, 2, 2), print_section(r->sec[1]));

	/* reverse-only traceback */
	r = _export(gaba_dp_trace)(d, NULL, f, NULL);
	assert(check_result(r, 4, 0, 4, 2, 1, 2, 1, 1), print_result(r));
	assert(check_path(r, "DRDR"), print_path(r));
	assert(check_cigar(r, "2M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 1, s->brsec, 0, 1, 0, 2), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 1, s->brsec, 0, 1, 2, 2), print_section(r->sec[1]));

	/* forward-reverse traceback */
	r = _export(gaba_dp_trace)(d, f, f, NULL);
	assert(check_result(r, 8, 0, 8, 4, 2, 0, 0, 0), print_result(r));
	assert(check_path(r, "DRDRDRDR"), print_path(r));
	assert(check_cigar(r, "4M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 1, s->brsec, 0, 1, 0, 2), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 1, s->brsec, 0, 1, 2, 2), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 0, 1, s->bfsec, 0, 1, 4, 2), print_section(r->sec[2]));
	assert(check_section(r->sec[3], s->afsec, 0, 1, s->bfsec, 0, 1, 6, 2), print_section(r->sec[3]));

	/* forward-reverse traceback with seed */
	r = _export(gaba_dp_trace)(d, f, f, GABA_TRACE_PARAMS(
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

	_export(gaba_dp_clean)(d);
}

/* with longer sequences */
unittest(with_seq_pair("ACGTACGTACGT", "ACGTACGTACGT"))
{
	omajinai();

	/* fill sections */
	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);

	/* fw */
	struct gaba_alignment_s *r = _export(gaba_dp_trace)(d, f, NULL, NULL);
	assert(check_result(r, 48, 0, 48, 2, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "24M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 12, s->bfsec, 0, 12, 0, 24), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 0, 12, s->bfsec, 0, 12, 24, 24), print_section(r->sec[1]));

	/* rv */
	r = _export(gaba_dp_trace)(d, NULL, f, NULL);
	assert(check_result(r, 48, 0, 48, 2, 1, 24, 12, 12), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "24M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 12, s->brsec, 0, 12, 0, 24), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 12, s->brsec, 0, 12, 24, 24), print_section(r->sec[1]));

	/* fw-rv */
	r = _export(gaba_dp_trace)(d, f, f, NULL);
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

	_export(gaba_dp_clean)(d);
}

/* concatenate */
unittest(with_seq_pair("ACGTACGTACGT", "ACGTACGTACGT"))
{
	omajinai();

	/* fill forward root section */
	struct gaba_fill_s *f1 = _export(gaba_dp_fill_root)(d, &s->afsec, 6, &s->bfsec, 6);
	assert(f1->status == 0x1ff, "%x", f1->status);
	assert(check_tail(f1, 0, 0, -19, 1), print_tail(f1));

	/* fill tail section */
	f1 = _export(gaba_dp_fill)(d, f1, &s->aftail, &s->bftail);
	assert(f1->status == 0x1ff, "%x", f1->status);
	assert(check_tail(f1, 12, 21, 21, 2), print_tail(f1));

	/* fill forward root section */
	struct gaba_fill_s *f2 = _export(gaba_dp_fill_root)(d, &s->arsec, 6, &s->brsec, 6);
	assert(f2->status == 0x1ff, "%x", f2->status);
	assert(check_tail(f2, 0, 0, -19, 1), print_tail(f2));

	/* fill tail section */
	f2 = _export(gaba_dp_fill)(d, f2, &s->artail, &s->brtail);
	assert(f2->status == 0x1ff, "%x", f2->status);
	assert(check_tail(f2, 12, 21, 21, 2), print_tail(f2));

	/* fw-rv */
	struct gaba_alignment_s *r = _export(gaba_dp_trace)(d, f1, f2, NULL);
	assert(check_result(r, 24, 0, 24, 1, 0, 12, 6, 6), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "12M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 12, s->bfsec, 0, 12, 0, 24), print_section(r->sec[0]));

	_export(gaba_dp_clean)(d);
}

/* sequences with different lengths (consumed as mismatches) */
unittest(with_seq_pair("GAAAAAAAA", "AAAAAAAA"))
{
	omajinai();

	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bftail);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);

	/* fw */
	struct gaba_alignment_s *r = _export(gaba_dp_trace)(d, f, NULL, NULL);
	assert(check_result(r, 22, 2, 32, 3, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "16M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 8, s->bfsec, 0, 8, 0, 16), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 8, 1, s->bfsec, 0, 1, 16, 2), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 0, 7, s->bfsec, 1, 7, 18, 14), print_section(r->sec[2]));

	/* rv */
	r = _export(gaba_dp_trace)(d, NULL, f, NULL);
	assert(check_result(r, 22, 2, 32, 3, 2, 16, 9, 8), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "16M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 2, 7, s->brsec, 0, 7, 0, 14), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 1, s->brsec, 7, 1, 14, 2), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->arsec, 1, 8, s->brsec, 0, 8, 16, 16), print_section(r->sec[2]));

	/* fw-rv */
	r = _export(gaba_dp_trace)(d, f, f, NULL);
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

	_export(gaba_dp_clean)(d);
}

/* another pair with different lengths */
unittest(with_seq_pair("TTTTTTTT", "CTTTTTTTT"))
{
	omajinai();

	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bfsec);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);

	/* fw */
	struct gaba_alignment_s *r = _export(gaba_dp_trace)(d, f, NULL, NULL);
	assert(check_result(r, 22, 2, 32, 3, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "16M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 8, s->bfsec, 0, 8, 0, 16), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 0, 1, s->bfsec, 8, 1, 16, 2), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 1, 7, s->bfsec, 0, 7, 18, 14), print_section(r->sec[2]));

	/* rv */
	r = _export(gaba_dp_trace)(d, NULL, f, NULL);
	assert(check_result(r, 22, 2, 32, 3, 2, 16, 8, 9), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "16M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 7, s->brsec, 2, 7, 0, 14), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 7, 1, s->brsec, 0, 1, 14, 2), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->arsec, 0, 8, s->brsec, 1, 8, 16, 16), print_section(r->sec[2]));

	/* fw-rv */
	r = _export(gaba_dp_trace)(d, f, f, NULL);
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

	_export(gaba_dp_clean)(d);
}

/* with deletions */
unittest(with_seq_pair("GACGTACGT", "ACGTACGT"))
{
	omajinai();

	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bftail);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);

	/* fw */
	struct gaba_alignment_s *r = _export(gaba_dp_trace)(d, f, NULL, NULL);
	assert(check_result(r, 20, 0, 34, 2, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "RDRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "1D8M1D8M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 9, s->bfsec, 0, 8, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 0, 9, s->bfsec, 0, 8, 17, 17), print_section(r->sec[1]));

	/* rv */
	r = _export(gaba_dp_trace)(d, NULL, f, NULL);
	assert(check_result(r, 20, 0, 34, 2, 1, 17, 9, 8), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDRR"), print_path(r));
	assert(check_cigar(r, "8M1D8M1D"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 9, s->brsec, 0, 8, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 9, s->brsec, 0, 8, 17, 17), print_section(r->sec[1]));

	/* fw-rv */
	r = _export(gaba_dp_trace)(d, f, f, NULL);
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

	_export(gaba_dp_clean)(d);
}

/* with insertions */
unittest(with_seq_pair("ACGTACGT", "GACGTACGT"))
{
	omajinai();

	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bfsec);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);

	/* fw */
	struct gaba_alignment_s *r = _export(gaba_dp_trace)(d, f, NULL, NULL);
	assert(check_result(r, 20, 0, 34, 2, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "DDRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "1I8M1I8M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 8, s->bfsec, 0, 9, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 0, 8, s->bfsec, 0, 9, 17, 17), print_section(r->sec[1]));

	/* rv */
	r = _export(gaba_dp_trace)(d, NULL, f, NULL);
	assert(check_result(r, 20, 0, 34, 2, 1, 17, 8, 9), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDRD"), print_path(r));
	assert(check_cigar(r, "8M1I8M1I"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 8, s->brsec, 0, 9, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 8, s->brsec, 0, 9, 17, 17), print_section(r->sec[1]));

	/* fw-rv */
	r = _export(gaba_dp_trace)(d, f, f, NULL);
	assert(check_result(r, 40, 0, 68, 4, 2, 0, 0, 0), print_result(r));
	assert(check_path(r,
		"DRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDRD"
		"DDRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "8M1I8M2I8M1I8M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 8, s->brsec, 0, 9, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 8, s->brsec, 0, 9, 17, 17), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 0, 8, s->bfsec, 0, 9, 34, 17), print_section(r->sec[2]));
	assert(check_section(r->sec[3], s->afsec, 0, 8, s->bfsec, 0, 9, 51, 17), print_section(r->sec[3]));

	_export(gaba_dp_clean)(d);
}

/* breakpoint adjustment */
unittest(with_seq_pair("GACGTACGTGACGTACGT", "ACGTACGT"))
{
	omajinai();

	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bftail);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);

	/* fw */
	struct gaba_alignment_s *r = _export(gaba_dp_trace)(d, f, NULL, NULL);
	assert(check_result(r, 20, 0, 34, 2, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "RDRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "1D8M1D8M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 10, s->bfsec, 0, 8, 0, 18), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 10, 8, s->bfsec, 0, 8, 18, 16), print_section(r->sec[1]));

	/* rv */
	r = _export(gaba_dp_trace)(d, NULL, f, NULL);
	assert(check_result(r, 20, 0, 34, 2, 1, 17, 18, 8), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDRR"), print_path(r));
	assert(check_cigar(r, "8M1D8M1D"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 9, s->brsec, 0, 8, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 9, 9, s->brsec, 0, 8, 17, 17), print_section(r->sec[1]));

	/* fw-rv */
	r = _export(gaba_dp_trace)(d, f, f, NULL);
	assert(check_result(r, 40, 0, 68, 4, 2, 0, 0, 0), print_result(r));
	assert(check_path(r,
		"DRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDRR"
		"RDRDRDRDRDRDRDRDRRDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "8M1D8M2D8M1D8M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 9, s->brsec, 0, 8, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 9, 9, s->brsec, 0, 8, 17, 17), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 0, 10, s->bfsec, 0, 8, 34, 18), print_section(r->sec[2]));
	assert(check_section(r->sec[3], s->afsec, 10, 8, s->bfsec, 0, 8, 52, 16), print_section(r->sec[3]));

	_export(gaba_dp_clean)(d);
}

unittest(with_seq_pair("ACGTACGT", "GACGTACGTGACGTACGT"))
{
	omajinai();

	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, &s->afsec, 0, &s->bfsec, 0);
	f = _export(gaba_dp_fill)(d, f, &s->afsec, &s->bfsec);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bfsec);
	f = _export(gaba_dp_fill)(d, f, &s->aftail, &s->bftail);

	/* fw */
	struct gaba_alignment_s *r = _export(gaba_dp_trace)(d, f, NULL, NULL);
	assert(check_result(r, 20, 0, 34, 2, 0, 0, 0, 0), print_result(r));
	assert(check_path(r, "DDRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "1I8M1I8M"), print_path(r));
	assert(check_section(r->sec[0], s->afsec, 0, 8, s->bfsec, 0, 10, 0, 18), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->afsec, 0, 8, s->bfsec, 10, 8, 18, 16), print_section(r->sec[1]));

	/* rv */
	r = _export(gaba_dp_trace)(d, NULL, f, NULL);
	assert(check_result(r, 20, 0, 34, 2, 1, 17, 8, 18), print_result(r));
	assert(check_path(r, "DRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDRD"), print_path(r));
	assert(check_cigar(r, "8M1I8M1I"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 8, s->brsec, 0, 9, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 8, s->brsec, 9, 9, 17, 17), print_section(r->sec[1]));

	/* fw-rv */
	r = _export(gaba_dp_trace)(d, f, f, NULL);
	assert(check_result(r, 40, 0, 68, 4, 2, 0, 0, 0), print_result(r));
	assert(check_path(r,
		"DRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDRD"
		"DDRDRDRDRDRDRDRDRDDRDRDRDRDRDRDRDR"), print_path(r));
	assert(check_cigar(r, "8M1I8M2I8M1I8M"), print_path(r));
	assert(check_section(r->sec[0], s->arsec, 0, 8, s->brsec, 0, 9, 0, 17), print_section(r->sec[0]));
	assert(check_section(r->sec[1], s->arsec, 0, 8, s->brsec, 9, 9, 17, 17), print_section(r->sec[1]));
	assert(check_section(r->sec[2], s->afsec, 0, 8, s->bfsec, 0, 10, 34, 18), print_section(r->sec[2]));
	assert(check_section(r->sec[3], s->afsec, 0, 8, s->bfsec, 10, 8, 52, 16), print_section(r->sec[3]));

	_export(gaba_dp_clean)(d);
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
		struct gaba_dp_context_s *d = _export(gaba_dp_init)(c, sec->alim, sec->blim);

		/* fill section */
		struct gaba_section_s const *as = &sec->afsec;
		struct gaba_section_s const *bs = &sec->bfsec;
		struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, as, 0, bs, 0);
		struct gaba_fill_s *m = f;

		/* fill tail (1) */
		as = (f->status & GABA_STATUS_UPDATE_A) ? &sec->aftail : as;
		bs = (f->status & GABA_STATUS_UPDATE_B) ? &sec->bftail : bs;
		struct gaba_fill_s *t1 = _export(gaba_dp_fill)(d, f, as, bs);
		m = (t1->max > m->max) ? t1 : m;

		/* fill tail (2) */
		as = (t1->status & GABA_STATUS_UPDATE_A) ? &sec->aftail : as;
		bs = (t1->status & GABA_STATUS_UPDATE_B) ? &sec->bftail : bs;
		struct gaba_fill_s *t2 = _export(gaba_dp_fill)(d, t1, as, bs);
		m = (t2->max > m->max) ? t2 : m;

		/* check max */
		assert(m->max == nf.score, "m->max(%lld), f(%lld, %u), t1->max(%lld, %u), t2->max(%lld, %u), n.score(%d)",
			m->max, f->max, f->status, t1->max, t1->status, t2->max, t2->status, nf.score);
		if(m->max != nf.score) {
			struct gaba_fill_s *f2 = _export(gaba_dp_fill_root)(d, &sec->afsec, 0, &sec->bfsec, 0);
			(void)f2;
			debug("refill f2(%lld, %u)", f2->max, f2->status);
		}

		/* forward trace */
		struct gaba_alignment_s *rf = _export(gaba_dp_trace)(d, m, NULL, NULL);

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
		struct gaba_alignment_s *rr = _export(gaba_dp_trace)(d, NULL, m, NULL);

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
		_export(gaba_dp_clean)(d);
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
	struct gaba_fill_s *f = _export(gaba_dp_fill_root)(d, as, 0, bs, 0);

	/* fill tail (1) */
	as = (f->status & GABA_STATUS_UPDATE_A) ? &s->aftail : as;
	bs = (f->status & GABA_STATUS_UPDATE_B) ? &s->bftail : bs;
	struct gaba_fill_s *t1 = _export(gaba_dp_fill)(d, f, as, bs);
	f = (t1->max > f->max) ? t1 : f;

	/* fill tail (2) */
	as = (f->status & GABA_STATUS_UPDATE_A) ? &s->aftail : as;
	bs = (f->status & GABA_STATUS_UPDATE_B) ? &s->bftail : bs;
	struct gaba_fill_s *t2 = _export(gaba_dp_fill)(d, t1, as, bs);
	f = (t2->max > f->max) ? t2 : f;
	struct gaba_result_s *r = _export(gaba_dp_trace)(d, f, NULL, NULL);

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
	_export(gaba_dp_clean)(d);
}
#endif

#endif /* UNITTEST */

/**
 * end of gaba.c
 */
