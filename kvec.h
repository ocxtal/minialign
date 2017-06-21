/* The MIT License

   Copyright (c) 2008, by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
  An example:

#include "kvec.h"
int main() {
	kvec_t(int) array;
	kv_init(array);
	kv_push(int, array, 10); // append
	kv_a(int, array, 20) = 5; // dynamic
	kv_A(array, 20) = 4; // static
	kv_destroy(array);
	return 0;
}
*/

/*
  2008-09-22 (0.1.0):

	* The initial version.

*/

#ifndef AC_KVEC_H
#define AC_KVEC_H

#include <stdlib.h>

#define kv_roundup(x, base)			( (((x) + (base) - 1) / (base)) * (base) )
#define kv_roundup32(x)				kv_roundup(x, 32)
#define kv_max2(a, b)				( ((a) < (b)) ? (b) : (a) )
#define kv_min2(a, b)				( ((a) < (b)) ? (a) : (b) )

#define kvec_t(type) struct { size_t n, m; type *a; }
#define kv_init(v) ((v).n = (v).m = 0, (v).a = 0)
#define kv_inits(type) ((kvec_t(type)){ .n = 0, .m = 0, .a = 0 })
#define kv_destroy(v) free((v).a)
#define kv_A(v, i) ((v).a[(i)])
#define kv_pop(v) ((v).a[--(v).n])
#define kv_size(v) ((v).n)
#define kv_max(v) ((v).m)

#define kv_resize(type, v, s) do { \
		if ((v).m < (s)) { \
			(v).m = (s); \
			kv_roundup32((v).m); \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
	} while (0)

#define kv_reserve(type, v, s) ( \
	(v).m > (s) ? 0 : ((v).m = (s), (v).a = realloc((v).a, sizeof(type) * (v).m), 0) )

#define kv_copy(type, v1, v0) do {							\
		if ((v1).m < (v0).n) kv_resize(type, v1, (v0).n);	\
		(v1).n = (v0).n;									\
		memcpy((v1).a, (v0).a, sizeof(type) * (v0).n);		\
	} while (0)												\

#define kv_push(type, v, x) do {									\
		if ((v).n == (v).m) {										\
			(v).m = (v).m? (v).m<<1 : 2;							\
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m);	\
		}															\
		(v).a[(v).n++] = (x);										\
	} while (0)

#define kv_pushp(type, v, p) do { \
		if ((v).n == (v).m) { \
			(v).m = (v).m? (v).m<<1 : 2; \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
		*(p) = &(v).a[(v).n++]; \
	} while (0)

#define kv_pushm(type, v, arr, size) do { \
		if(((v).m - (v).n) < (uint64_t)(size)) { \
			(v).m = kv_max2((v).m * 2, (v).n + (size));		\
			(v).a = (type*)realloc((v).a, sizeof(*(v).a) * (v).m);	\
		} \
		for(uint64_t _i = 0; _i < (uint64_t)(size); _i++) { \
			(v).a[(v).n + _i] = (arr)[_i]; \
		} \
		(v).n += (uint64_t)(size); \
	} while(0)

#define kv_a(type, v, i) ((v).m <= (size_t)(i)?						\
						  ((v).m = (v).n = (i) + 1, kv_roundup32((v).m), \
						   (v).a = (type*)realloc((v).a, sizeof(type) * (v).m), 0) \
						  : (v).n <= (size_t)(i)? (v).n = (i)			\
						  : 0), (v).a[(i)]

#define kv_reverse(type, v, start) do { \
		if ((v).m > 0 && (v).n > (start)) { \
			size_t __i, __end = (v).n - (start); \
			type *__a = (v).a + (start); \
			for (__i = 0; __i < __end>>1; ++__i) { \
				type __t = __a[__end - 1 - __i]; \
				__a[__end - 1 - __i] = __a[__i]; __a[__i] = __t; \
			} \
		} \
	} while (0)

/** heap queue : elements in v must be orderd in heap */
#define kv_hq_init(v)		{ (v).n = (v).m = 1; (v).a = NULL; }
#define kv_hq_inits(type)	((kvec_t(type)){ .n = 0, .m = 0, .a = NULL })
#define kv_hq_destroy(v)	kv_destroy(v)
#define kv_hq_clear(v)		( (v).n = 1 )

#define kv_hq_n(v, i) ( *((int64_t *)&v.a[i]) )
#define kv_hq_push(type, __comp, v, x) { \
	kv_push(type, v, x); \
	uint64_t i = (v).n - 1; \
	while(i > 1 && __comp((v).a[i>>1], (v).a[i]) > 0) { \
		(v).a[0] = (v).a[i>>1]; \
		(v).a[i>>1] = (v).a[i]; \
		(v).a[i] = (v).a[0]; \
		i >>= 1; \
	} \
}
#define kv_hq_pop(type, __comp, v) ({ \
	uint64_t i = 1, j = 2; \
	(v).a[0] = (v).a[i]; \
	(v).a[i] = (v).a[--(v).n]; \
	(v).a[(v).n] = (v).a[0]; \
	while(j < (v).n) { \
		uint64_t k; \
		/*k = (j + 1 < (v).n && kv_hq_n(v, j + 1) < kv_hq_n(v, j)) ? (j + 1) : j; */ \
		k = (j + 1 < (v).n && __comp((v).a[j + 1], (v).a[j]) < 0) ? j + 1 : j; \
		/*k = (kv_hq_n(v, k) < kv_hq_n(v, i)) ? k : 0; */ \
		k = __comp((v).a[k], (v).a[i]) < 0 ? k : 0; \
		if(k == 0) { break; } \
		(v).a[0] = (v).a[k]; \
		(v).a[k] = (v).a[i]; \
		(v).a[i] = (v).a[0]; \
		i = k; j = k<<1; \
	} \
	v.a[v.n]; \
})
#endif
