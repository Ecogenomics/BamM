/* Copyright (c) 2009, Thomas Hurst <tom@hur.st>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/* pqsort - Partial quicksort algorithm.
 *
 * This is the sort routine largely as described by Bentley & McIlroy's
 * "Engineering a Sort Function", along with some modifications as found
 * in FreeBSD's libc.
 *
 * "Partial" means pqsort can be used to place only specific items in
 * their final sorted positions, making it useful for pagination, ranking
 * and partitioning.
 *
 * Usage:
 * #include <pqsort.h>
 * define_pqsort(str, char *, strcmp);
 *
 *     str_pqsort(array, nitems, offset, limit);
 *
 * offset=100, limit=100 will place array[100] - array[199] in the correct
 * position, ensure array[0] - array[99] contains only values smaller than
 * that in array[100], and array[200] - array[nitems - 1] only values larger
 * than array[199].
 *
 * Since it's a preprocessor macro, comparison functions can be inlined.
 * This can be quite a win if comparisons are cheap compared with function
 * calls, e.g. when sorting integers.
 *
 * WWW: http://www.aagh.net/projects/pqsort
 *
 * pqsort development was funded by Newzbin Ltd (http://www.newzbin.com/).
 */

#ifndef __MY_PQSORT_H
#define __MY_PQSORT_H

#include <sys/cdefs.h>
#include <sys/types.h>
#include <stdlib.h>

#define PQSORT_MIN(a, b)	(a) < (b) ? a : b
#define define_pqsort(a_name, a_type, a_cmp)				\
static inline size_t							\
a_name##_med3(a_type x[], const size_t a, const size_t b,		\
              const size_t c)						\
{									\
	return a_cmp(x[a], x[b]) < 0 ?					\
	           (a_cmp(x[b], x[c]) < 0 ? b				\
	               : (a_cmp(x[a], x[c]) < 0 ? c : a))		\
	       : (a_cmp(x[b], x[c]) > 0 ? b				\
	           : (a_cmp(x[a], x[c]) < 0 ? a : c));			\
}									\
									\
static inline void							\
a_name##_swap(a_type x[], const size_t a, const size_t b)		\
{									\
	a_type temp = x[a];x[a] = x[b];x[b] = temp;			\
}									\
									\
static void								\
a_name##_insert_sort(a_type x[], const size_t first, const size_t last)	\
{									\
	size_t left,right;						\
	for (left = first + 1; left <= last; left++)			\
	{								\
		for (right = left;					\
		     right > first && a_cmp(x[right - 1], x[right]) > 0;\
		     right--)						\
		{							\
			a_name##_swap(x, right, right - 1);		\
		}							\
	}								\
}									\
									\
static inline size_t							\
a_name##_select_pivot(a_type x[], const size_t first, const size_t last)\
{									\
	size_t length = last -  first;					\
	size_t eigth = length / 8;					\
	size_t middle = first + length / 2;				\
	size_t f = first, m = middle, l = last;				\
	if (length > 40)						\
	{								\
		f = a_name##_med3(x,					\
		    first, first + eigth, first + 2 * eigth);		\
		m = a_name##_med3(x,					\
		    middle - eigth, middle, middle + eigth);		\
		l = a_name##_med3(x,					\
		    last - 2 * eigth, last - eigth, last);		\
	}								\
	return a_name##_med3(x, f, m, l);				\
}									\
									\
static void								\
a_name##_pqsorter(a_type x[], size_t first, const size_t last,		\
                  const size_t offset, const size_t olimit)		\
{									\
	size_t length = (last - first) + 1;				\
	size_t pindex, a, b, c, d;					\
	ssize_t cmp_result;						\
	a_type piv;							\
	int swapped;							\
									\
top:									\
	swapped = 0;							\
									\
	if (length < 7) {						\
		if (length > 1) a_name##_insert_sort(x, first, last);	\
		return;							\
	}								\
									\
	pindex = a_name##_select_pivot(x, first, last);			\
	a_name##_swap(x, first, pindex);				\
	piv = x[first];							\
									\
	a = b = first + 1;						\
	c = d = last;							\
	for (;;)							\
	{								\
		while (b <= c && (cmp_result = a_cmp(x[b], piv)) <= 0)	\
		{							\
			if (cmp_result == 0) {				\
				a_name##_swap(x, a, b);			\
				swapped = 1;				\
				a++;					\
			}						\
			b++;						\
		}							\
		while (b <= c && (cmp_result = a_cmp(x[c], piv)) >= 0)	\
		{							\
			if (cmp_result == 0) {				\
				a_name##_swap(x, c, d);			\
				swapped = 1;				\
				d--;					\
			}						\
			c--;						\
		}							\
		if (b > c) break;					\
		a_name##_swap(x, b, c);					\
		swapped = 1;						\
		b++;							\
		c--;							\
	}								\
	if (swapped == 0) {						\
		a_name##_insert_sort(x, first, last);			\
		return;							\
	}								\
									\
	size_t s, l, h;							\
	s = PQSORT_MIN(a - first, b - a);				\
	for (l = first, h = b - s; s; s--) {				\
		a_name##_swap(x, l, h);					\
		l++, h++;						\
	}								\
	s = PQSORT_MIN(d - c, last - d - 1);				\
	for (l = b, h = last - s; s; s--) {				\
		a_name##_swap(x, l, h);					\
		l++, h++;						\
	}								\
									\
	size_t nl, nr;							\
	nl = b - a;							\
	nr = d - c;							\
									\
	if (nl > 0 && (nl = first + nl) >= offset)			\
	{								\
		a_name##_pqsorter(x, first, nl, offset, olimit);	\
	}								\
	if (nr > 0 && (nr = last - nr) <= olimit)			\
	{								\
		first = nr;						\
		goto top;						\
	}								\
}									\
									\
void a_name##_pqsort(a_type x[], const size_t len, const size_t offset,	\
                     const size_t limit);				\
void a_name##_pqsort(a_type x[], const size_t len, const size_t offset,	\
                     const size_t limit)				\
{									\
	a_name##_pqsorter(x, 0, len-1, offset, offset+limit-1);		\
}
#endif
