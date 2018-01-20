
/**
 * @file gaba.h
 *
 * @brief C header of the libgaba (libsea3) API
 *
 * @author Hajime Suzuki
 * @date 2014/12/29
 * @license Apache v2
 *
 * @detail
 * a header for libgaba (libsea3): a fast banded seed-and-extend alignment library.
 */

#ifndef _GABA_H_INCLUDED
#define _GABA_H_INCLUDED

#include <stdlib.h>		/** NULL and size_t */
#include <stdint.h>		/** uint8_t, int32_t, int64_t */

/**
 * @macro GABA_EXPORT_LEVEL
 */
#ifdef _GABA_WRAP_H_INCLUDED
/* included from gaba_wrap.h */
#  define GABA_EXPORT_LEVEL		static inline
#else
/* single, linked to an object compiled without -DSUFFIX */
#  define GABA_EXPORT_LEVEL
#endif

/**
 * @enum gaba_status
 */
enum gaba_status {
	GABA_CONT 		= 0,		/* continue, call again the function with the same args (but rarely occurrs) */
	GABA_UPDATE_A 	= 0x000f,	/* update required on section a (always combined with GABA_UPDATE) */
	GABA_UPDATE_B 	= 0x00f0,	/* update required on section b (always combined with GABA_UPDATE) */
	GABA_TERM		= 0x8000,	/* extension terminated by X-drop */
	GABA_OOM		= 0x0400	/* out of memory (indicates malloc returned NULL) */
};

/**
 * @type gaba_lmalloc_t, gaba_free_t
 * @brief external malloc can be passed, otherwise system malloc will be used
 */
typedef void *(*gaba_lmalloc_t)(void *opaque, size_t size);
typedef void (*gaba_lfree_t)(void *opaque, void *ptr);

/**
 * @struct gaba_alloc_s
 * @brief optional memory allocator, malloc and free pair must not be NULL.
 * memory block must be freed when size == 0:
 *
 *	void *alloc(void *opaque, void *ptr, size_t size) {
 *		if(size == 0) { free(ptr); }
 *		return(realloc(ptr, size));
 *	}
 */
struct gaba_alloc_s {
	void *opaque;				/** local memory arena */
	gaba_lmalloc_t lmalloc;		/** local malloc; dedicated for alignment path generation */
	gaba_lfree_t lfree;			/** local free */
};
typedef struct gaba_alloc_s gaba_alloc_t;

/**
 * @struct gaba_params_s
 * @brief input parameters of gaba_init
 */
struct gaba_params_s {
	/** scoring parameters */
	int8_t score_matrix[16];	/** score matrix (substitution matrix) max must not exceed 7 */
	int8_t gi;					/** gap open penalty (0 for the linear-gap penalty; positive integer) */
	int8_t ge;					/** gap extension penalty (positive integer) */
	int8_t gfa, gfb;			/** linear-gap extension penalty for short indels (combined-gap penalty; gf > ge) */

	/** score parameters */
	int8_t xdrop;				/** X-drop threshold, positive, less than 128 */

	/** filtering parameters */
	uint8_t filter_thresh;		/** popcnt filter threshold, set zero if you want to disable it */

	/** output options */
	uint32_t head_margin;		/** margin at the head of gaba_res_t */
	uint32_t tail_margin;		/** margin at the tail of gaba_res_t */

	/* internal */
	void *reserved;
	uint64_t _pad;
};
typedef struct gaba_params_s gaba_params_t;

/**
 * @macro GABA_PARAMS
 * @brief utility macro for gaba_init, see example on header.
 */
#define GABA_PARAMS(...)			( &((struct gaba_params_s const) { __VA_ARGS__ }) )

/**
 * @macro GABA_SCORE_SIMPLE
 * @brief utility macro for constructing score parameters.
 */
#define GABA_SCORE_SIMPLE(_m, _x, _gi, _ge)	.score_matrix = { (_m),-(_x),-(_x),-(_x),-(_x),(_m),-(_x),-(_x),-(_x),-(_x),(_m),-(_x),-(_x),-(_x),-(_x),(_m) }, .gi = (_gi), .ge = (_ge)

/**
 * @type gaba_t
 *
 * @brief (API) an alias to `struct gaba_context_s'.
 */
typedef struct gaba_context_s gaba_t;

/**
 * @type gaba_stack_t
 *
 * @brief stack context container
 */
typedef struct gaba_stack_s gaba_stack_t;

/**
 * @struct gaba_section_s
 *
 * @brief section container, a tuple of (id, length, head position).
 */
struct gaba_section_s {
	uint32_t id;				/** (4) section id */
	uint32_t len;				/** (4) length of the seq */
	uint8_t const *base;		/** (8) pointer to the head of the sequence */
};
typedef struct gaba_section_s gaba_section_t;
#define gaba_build_section(_id, _base, _len) ( \
	(struct gaba_section_s){ \
		.id = (_id), \
		.base = (uint8_t const *)(_base), \
		.len = (_len) \
	} \
)
#define gaba_rev(pos, len)		( (len) + (uint64_t)(len) - (uint64_t)(pos) - 1 )

/**
 * @type gaba_dp_t
 *
 * @brief an alias to `struct gaba_dp_context_s`.
 */
#ifndef _GABA_WRAP_H_INCLUDED
typedef struct gaba_dp_context_s gaba_dp_t;
#endif

/**
 * @struct gaba_fill_s
 */
struct gaba_fill_s {
	uint32_t aid, bid;			/** (8) the last-filled section ids */
	uint32_t ascnt, bscnt;		/** (8) aligned section counts */
	uint64_t apos, bpos;		/** (16) #fetched bases from the head (ppos = apos + bpos) */
	int64_t max;				/** (8) max score in the entire band */
	uint32_t stat;				/** (4) status (section update flags) */
	// int32_t ppos;				/** (8) #vectors from the head (FIXME: should be 64bit int) */
	uint32_t reserved[5];
};
typedef struct gaba_fill_s gaba_fill_t;

/**
 * @struct gaba_pos_pair_s
 */
struct gaba_pos_pair_s {
	uint32_t aid, bid;
	uint32_t apos, bpos;
};
typedef struct gaba_pos_pair_s gaba_pos_pair_t;

/**
 * @struct gaba_segment_s
 */
struct gaba_segment_s {
	uint32_t aid, bid;			/** (8) id of the sections */
	uint32_t apos, bpos;		/** (8) pos in the sections */
	uint32_t alen, blen;		/** (8) lengths of the segments */
	uint64_t ppos;				/** (8) path string position (offset) */
};
typedef struct gaba_segment_s gaba_path_section_t;
#define gaba_plen(seg)		( (seg)->alen + (seg)->blen )

/**
 * @struct gaba_alignment_s
 */
struct gaba_alignment_s {
	/* reserved for internal use */
	void *reserved1[2];
	uint32_t reserved2;

	uint32_t slen;				/* segment length */
	struct gaba_segment_s const *seg;

	uint64_t plen;				/* path length */
	int64_t score;				/** score */
	// uint32_t mcnt, xcnt;		/** #matches, #mismatches */
	// uint32_t gicnt, gecnt;		/** #gap opens, #gap bases */
	uint32_t reserved3[4];
	uint32_t path[];
};
typedef struct gaba_alignment_s gaba_alignment_t;

/**
 * @fn gaba_init
 * @brief (API) gaba_init new API
 */
GABA_EXPORT_LEVEL
gaba_t *gaba_init(gaba_params_t const *params);

/**
 * @fn gaba_clean
 * @brief (API) clean up the alignment context structure.
 */
GABA_EXPORT_LEVEL
void gaba_clean(gaba_t *ctx);

/**
 * @fn gaba_dp_init
 * @brief create thread-local context deriving the global context (ctx)
 * with local memory arena and working buffers. alim and blim are respectively
 * the tails of sequence arrays.
 */
GABA_EXPORT_LEVEL
gaba_dp_t *gaba_dp_init(
	gaba_t const *ctx,
	uint8_t const *alim,
	uint8_t const *blim);

/**
 * @fn gaba_dp_flush
 * @brief flush stack (flush all if NULL) 
 */
GABA_EXPORT_LEVEL
void gaba_dp_flush(
	gaba_dp_t *dp);

/**
 * @fn gaba_dp_save_stack
 */
GABA_EXPORT_LEVEL
gaba_stack_t const *gaba_dp_save_stack(
	gaba_dp_t *dp);

/**
 * @fn gaba_dp_flush_stack
 */
GABA_EXPORT_LEVEL
void gaba_dp_flush_stack(
	gaba_dp_t *dp,
	gaba_stack_t const *stack);

/**
 * @fn gaba_dp_clean
 */
GABA_EXPORT_LEVEL
void gaba_dp_clean(
	gaba_dp_t *dp);

/**
 * @fn gaba_dp_fill_root
 */
GABA_EXPORT_LEVEL
gaba_fill_t *gaba_dp_fill_root(
	gaba_dp_t *dp,
	gaba_section_t const *a,
	uint32_t apos,
	gaba_section_t const *b,
	uint32_t bpos,
	uint32_t pridx);

/**
 * @fn gaba_dp_fill
 * @brief fill dp matrix inside section pairs
 */
GABA_EXPORT_LEVEL
gaba_fill_t *gaba_dp_fill(
	gaba_dp_t *dp,
	gaba_fill_t const *prev_sec,
	gaba_section_t const *a,
	gaba_section_t const *b,
	uint32_t pridx);

/**
 * @fn gaba_dp_merge
 * @brief merge multiple sections. all the vectors (tail objects) must be aligned on the same ppos,
 * and qofs are the q-distance of the two fill objects.
 */
#define MAX_MERGE_COUNT				( 14 )
GABA_EXPORT_LEVEL
gaba_fill_t *gaba_dp_merge(
	gaba_dp_t *dp,
	gaba_fill_t const *const *sec,
	uint8_t const *qofs,
	uint32_t cnt);

/**
 * @fn gaba_dp_search_max
 */
GABA_EXPORT_LEVEL
gaba_pos_pair_t *gaba_dp_search_max(
	gaba_dp_t *dp,
	gaba_fill_t const *sec);

/**
 * @fn gaba_dp_trace
 * @brief generate alignment result string, alloc->malloc and alloc->free must not be NULL if alloc is not NULL.
 */
GABA_EXPORT_LEVEL
gaba_alignment_t *gaba_dp_trace(
	gaba_dp_t *dp,
	gaba_fill_t const *tail,
	gaba_alloc_t const *alloc);

/**
 * @fn gaba_dp_res_free
 */
GABA_EXPORT_LEVEL
void gaba_dp_res_free(
	gaba_alignment_t *aln);

#endif  /* #ifndef _GABA_H_INCLUDED */

/*
 * end of gaba.h
 */
