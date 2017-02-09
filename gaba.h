

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
 * 
 * from C:
 * Include this header file as #include <gaba.h>. This will enable you to
 * use all the APIs in the gaba_init, and gaba_align form.
 *
 * from C++:
 * Include this header as #include <gaba.h>. The C APIs are wrapped with
 * namespace sea and the C++ class AlignmentContext and AlignmentResult
 * are added. See example.cpp for the detail of the usage in C++.
 */

#ifndef _GABA_H_INCLUDED
#define _GABA_H_INCLUDED

#include <stdlib.h>		/** NULL and size_t */
#include <stdint.h>		/** uint8_t, int32_t, int64_t */

/**
 * @enum gaba_error
 *
 * @brief (API) error flags. see gaba_init function and status member in the gaba_alignment structure for more details.
 */
enum gaba_error {
	GABA_SUCCESS 				=  0,	/*!< success!! */
	GABA_TERMINATED				=  1,	/*!< (internal code) success */
	GABA_ERROR 					= -1,	/*!< unknown error */
	/** errors which occur in an alignment function */
	GABA_ERROR_INVALID_MEM 		= -2,	/*!< invalid pointer to memory */
	GABA_ERROR_INVALID_CONTEXT 	= -3,	/*!< invalid pointer to the alignment context */
	GABA_ERROR_OUT_OF_BAND 		= -4,	/*!< traceback failure. using wider band may resolve this type of error. */
	GABA_ERROR_OUT_OF_MEM 		= -5,	/*!< out of memory error. mostly caused by exessively long queries. */
	GABA_ERROR_OVERFLOW 			= -6, 	/*!< cell overflow error */
	GABA_ERROR_INVALID_ARGS 		= -7,	/*!< inproper input arguments. */
	/** errors which occur in an initialization function */
	GABA_ERROR_UNSUPPORTED_ALG 	= -8,	/*!< unsupported combination of algorithm and processor options. use naive implementations instead. */
	GABA_ERROR_INVALID_COST 		= -9	/*!< invalid alignment cost */
};

/**
 * @enum gaba_clip_type
 */
enum gaba_clip_type {
	GABA_CLIP_SOFT = 'S',
	GABA_CLIP_HARD = 'H'
};

/**
 * @struct gaba_params_s
 * @brief input parameters of gaba_init
 */
struct gaba_params_s {
	/** scoring parameters */
	int8_t m, x, gi, ge;		/** match, mismatch, gap open, and gap extend, all in positive integer */

	/** score parameters */
	int8_t xdrop;

	/** filtering parameters */
	uint8_t filter_thresh;		/** popcnt filter threshold, set zero if you want to disable it */

	/** output options */
	uint8_t head_margin;		/** margin at the head of gaba_res_t */
	uint8_t tail_margin;		/** margin at the tail of gaba_res_t */
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
#define GABA_SCORE_SIMPLE(_m, _x, _gi, _ge)		.m = (_m), .x = (_x), .gi = (_gi), .ge = (_ge)

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
		.base = (_base), \
		.len = (_len) \
	} \
)
#define gaba_rev(pos, len)		( (len) + (uint64_t)(len) - (uint64_t)(pos) - 1 )

/**
 * @type gaba_dp_t
 *
 * @brief an alias to `struct gaba_dp_context_s`.
 */
typedef struct gaba_dp_context_s gaba_dp_t;

/**
 * @struct gaba_fill_s
 */
struct gaba_fill_s {
	/* coordinates */
	int64_t psum;				/** (8) global p-coordinate of the tail of the section */
	int32_t p;					/** (4) local p-coordinate of the tail of the section */
	uint32_t ssum;				/** (4) */

	/* status and max scores */
	int64_t max;				/** (8) max */
	uint32_t status;			/** (4) */

	uint8_t _pad[36];
};
typedef struct gaba_fill_s gaba_fill_t;

/**
 * @enum gaba_status
 */
enum gaba_status {
	GABA_STATUS_CONT 		= 0,
	GABA_STATUS_UPDATE		= 0x100,
	GABA_STATUS_UPDATE_A 	= 0x0f,
	GABA_STATUS_UPDATE_B 	= 0xf0,
	GABA_STATUS_TERM		= 0x200
};

/**
 * @struct gaba_pos_pair_s
 */
struct gaba_pos_pair_s {
	uint32_t apos, bpos;
};
typedef struct gaba_pos_pair_s gaba_pos_pair_t;

/**
 * @struct gaba_path_section_s
 */
struct gaba_path_section_s {
	uint32_t aid, bid;			/** (8) id of the sections */
	uint32_t apos, bpos;		/** (8) pos in the sections */
	uint32_t alen, blen;		/** (8) lengths of the segments */
	int64_t ppos;				/** (8) path string position (offset) */
};
typedef struct gaba_path_section_s gaba_path_section_t;
#define gaba_plen(sec)		( (sec)->alen + (sec)->blen )

/**
 * @struct gaba_path_s
 */
struct gaba_path_s {
	int64_t len;				/** (8) path length (= array bit length) */
	uint32_t array[];			/** () path array */
};
typedef struct gaba_path_s gaba_path_t;

/**
 * @struct gaba_alignment_s
 */
struct gaba_alignment_s {
	void *lmm;
	int64_t score, xcnt;		/** (16) score, #mismatchs */
	int64_t gicnt, gecnt;		/** (16) #gap opens, #gap bases */
	uint32_t rapos, rbpos;
	uint32_t rppos;				/** (4) local path index in the root section */
	uint32_t rsidx;				/** (4) index of the root section */
	uint32_t reserved3;
	uint32_t slen;
	struct gaba_path_section_s const *sec;
	struct gaba_path_s const *path;
};
typedef struct gaba_alignment_s gaba_alignment_t;

/**
 * @fn gaba_init
 * @brief (API) gaba_init new API
 */
gaba_t *gaba_init(gaba_params_t const *params);

/**
 * @fn gaba_clean
 *
 * @brief (API) clean up the alignment context structure.
 *
 * @param[in] ctx : a pointer to the alignment structure.
 *
 * @return none.
 *
 * @sa gaba_init
 */
void gaba_clean(
	gaba_t *ctx);

/**
 * @fn gaba_dp_init
 */
gaba_dp_t *gaba_dp_init(
	gaba_t const *ctx,
	uint8_t const *alim,
	uint8_t const *blim);

/**
 * @fn gaba_dp_flush
 * @brief flush stack (flush all if NULL) 
 */
void gaba_dp_flush(
	gaba_dp_t *dp,
	uint8_t const *alim,
	uint8_t const *blim);

/**
 * @fn gaba_dp_save_stack
 */
gaba_stack_t const *gaba_dp_save_stack(
	gaba_dp_t *dp);

/**
 * @fn gaba_dp_flush_stack
 */
void gaba_dp_flush_stack(
	gaba_dp_t *dp,
	gaba_stack_t const *stack);

/**
 * @fn gaba_dp_clean
 */
void gaba_dp_clean(
	gaba_dp_t *dp);

/**
 * @fn gaba_dp_fill_root
 */
gaba_fill_t *gaba_dp_fill_root(
	gaba_dp_t *dp,
	gaba_section_t const *a,
	uint32_t apos,
	gaba_section_t const *b,
	uint32_t bpos);

/**
 * @fn gaba_dp_fill
 * @brief fill dp matrix inside section pairs
 */
gaba_fill_t *gaba_dp_fill(
	gaba_dp_t *dp,
	gaba_fill_t const *prev_sec,
	gaba_section_t const *a,
	gaba_section_t const *b);

/**
 * @fn gaba_dp_merge
 */
gaba_fill_t *gaba_dp_merge(
	gaba_dp_t *dp,
	gaba_fill_t const *sec_list,
	uint64_t sec_list_len);

/**
 * @fn gaba_dp_search_max
 */
gaba_pos_pair_t gaba_dp_search_max(
	gaba_dp_t *dp,
	gaba_fill_t const *sec);

/**
 * @struct gaba_trace_params_s
 */
struct gaba_trace_params_s {
	void *lmm;
	struct gaba_path_section_s const *sec;
	uint16_t slen;				/* section length */
	uint16_t k;					/* path length (k-mer length) */
	uint16_t xcnt;				/* #mismatches */
	uint16_t _pad;
};
typedef struct gaba_trace_params_s gaba_trace_params_t;

/**
 * @macro GABA_TRACE_PARAMS
 */
#define GABA_TRACE_PARAMS(...)		( &((struct gaba_trace_params_s const) { __VA_ARGS__ }) )
#define GABA_TRACE_NONE				( NULL )

/**
 * @type gaba_alignment_writer
 * @brief pointer to putchar-compatible writer
 */
typedef int (*gaba_alignment_writer)(int c);

/**
 * @fn gaba_dp_trace
 *
 * @brief generate alignment result string
 */
gaba_alignment_t *gaba_dp_trace(
	gaba_dp_t *dp,
	gaba_fill_t const *fw_tail,
	gaba_fill_t const *rv_tail,
	gaba_trace_params_t const *params);

/**
 * @fn gaba_dp_recombine
 *
 * @brief recombine two alignments x and y at xsid and ysid.
 */
gaba_alignment_t *gaba_dp_recombine(
	gaba_dp_t *dp,
	gaba_alignment_t *x,
	uint32_t xsid,
	gaba_alignment_t *y,
	uint32_t ysid);

/**
 * @fn gaba_dp_res_free
 */
void gaba_dp_res_free(
	gaba_alignment_t *aln);

/**
 * @fn gaba_dp_print_cigar_forward
 *
 * @brief convert path string to cigar.
 * @detail
 * printer(void *fp, int64_t len, char c); will be called with a pair;
 * cigar operation (c) and its length (len).
 * void *fp is an opaque pointer to the context of the printer.
 */
typedef int (*gaba_dp_printer_t)(void *, int64_t, char);
uint64_t gaba_dp_print_cigar_forward(
	gaba_dp_printer_t printer,
	void *fp,
	uint32_t const *path,
	uint32_t offset,
	uint32_t len);

/**
 * @fn gaba_dp_print_cigar_reverse
 *
 * @brief convert path string to cigar in reverse direction
 */
uint64_t gaba_dp_print_cigar_reverse(
	gaba_dp_printer_t printer,
	void *fp,
	uint32_t const *path,
	uint32_t offset,
	uint32_t len);

/**
 * @fn gaba_dp_dump_cigar_forward
 */
uint64_t gaba_dp_dump_cigar_forward(
	char *buf,
	uint64_t buf_size,
	uint32_t const *path,
	uint32_t offset,
	uint32_t len);

/**
 * @fn gaba_dp_dump_cigar_reverse
 */
uint64_t gaba_dp_dump_cigar_reverse(
	char *buf,
	uint64_t buf_size,
	uint32_t const *path,
	uint32_t offset,
	uint32_t len);

#endif  /* #ifndef _GABA_H_INCLUDED */

/*
 * end of gaba.h
 */
