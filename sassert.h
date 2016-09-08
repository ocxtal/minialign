
/**
 * @file sassert.h
 *
 * @brief static assertion
 */
#ifndef _SASSERT_H_INCLUDED
#define _SASSERT_H_INCLUDED

#include <stddef.h>

/**
 * assert macros
 */
#define _sa_cat_intl(x, y) 		x##y
#define _sa_cat(x, y)			_sa_cat_intl(x, y)
/* static assert */
#define _static_assert(expr)	typedef char _sa_cat(_st, __LINE__)[(expr) ? 1 : -1]
/* check offset equality of elements in two structs */
#define _static_assert_offset(st1, mb1, st2, mb2, ofs) \
	_static_assert(offsetof(st1, mb1) == offsetof(st2, mb2) + ofs)


#endif /* #ifndef _SASSERT_H_INCLUDED */
/**
 * end of sassert.h
 */
