
/**
 * @file arch.h
 */
#ifndef _ARCH_H_INCLUDED
#define _ARCH_H_INCLUDED


#ifdef __x86_64__
#  if defined(__AVX2__)
#    include "x86_64_avx2/arch_util.h"
#    include "x86_64_avx2/vector.h"
#  elif defined(__SSE4_1__)
#    include "x86_64_sse41/arch_util.h"
#    include "x86_64_sse41/vector.h"
#  else
#    error "No SIMD instruction set enabled. Check if SSE4.1 or AVX2 instructions are available and add `-msse4.1' or `-mavx2' to CFLAGS."
#  endif

/* map reverse-complement sequence out of the canonical-formed address */
#define GREF_SEQ_LIM			( (uint8_t const *)0x800000000000 )

#endif

#ifdef AARCH64

/* use x86_64 default */
#define GREF_SEQ_LIM			( (uint8_t const *)0x800000000000 )

#endif

#ifdef PPC64

/* use x86_64 default */
#define GREF_SEQ_LIM			( (uint8_t const *)0x800000000000 )

#endif

#if !defined(_ARCH_UTIL_H_INCLUDED) || !defined(_VECTOR_H_INCLUDED)
#  error "No SIMD environment detected. Check CFLAGS."
#endif

#ifndef GREF_SEQ_LIM
#  error "No architecuture detected. Check CFLAGS."
#endif


/* elem_t and move definitions */
#define _rd(p)				( *((elem_t *)p) )
#define _wr(p, k)			{ *((elem_t *)p) = (k); }
#define _ex(k, p)			( ((k)>>((p)*8)) & (WCR_OCC_SIZE-1) )
#define _p(v)				( (elem_t)(v) )
#define _e(v)				( (uint64_t)(v) )


#endif /* #ifndef _ARCH_H_INCLUDED */
/**
 * end of arch.h
 */
