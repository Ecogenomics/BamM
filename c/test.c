// compile with: gcc -std=gnu99 test.c

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "pqsort.h"

//===================== Method 1: =============================================
//Algorithm from N. Wirth’s book Algorithms + data structures = programs of 1976    

#ifndef ELEM_SWAP(a,b)
#define ELEM_SWAP(a,b) { register int t=(a);(a)=(b);(b)=t; }

//===================== Method 2: =============================================
//This is the faster median determination method.
//Algorithm from Numerical recipes in C of 1992

int quick_select_median(int arr[], int n)
{
    int low, high ;
    int median;
    int middle, ll, hh;
    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;
        if (high == low + 1) { /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }
        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if (arr[middle] > arr[high])
            ELEM_SWAP(arr[middle], arr[high]) ;
        if (arr[low] > arr[high])
            ELEM_SWAP(arr[low], arr[high]) ;
        if (arr[middle] > arr[low])
            ELEM_SWAP(arr[middle], arr[low]) ;
        /* Swap low item (now in position middle) into position (low+1) */
        ELEM_SWAP(arr[middle], arr[low+1]) ;
        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        for (;;) {
            do ll++; while (arr[low] > arr[ll]) ;
            do hh--; while (arr[hh] > arr[low]) ;
            if (hh < ll)
                break;
            ELEM_SWAP(arr[ll], arr[hh]) ;
        }
        /* Swap middle item (in position low) back into correct position */
        ELEM_SWAP(arr[low], arr[hh]) ;
        /* Re-set active partition */
        if (hh <= median)
            low = ll;
        if (hh >= median)
            high = hh - 1;
    }
    return arr[median] ;
}
#endif

static inline int IntCmp(const int a, const int b) {
        return (a - b);
}

define_pqsort(my_int, int, IntCmp);

#define DUMP()  printf("{ "); \
        for (int i = 0; i < nvalues; i++)  \
                printf("%d, ", values[i]); \
        printf("}\n");

int main(void) {
        int values[] = { 42, 98, 56, 23, 45, 63, 56, 80, 102, 2 };
        int nvalues = sizeof(values) / sizeof(values[0]);

        // partition the list in two
        my_int_pqsort(values, nvalues, 5, 0);
        DUMP(); // => { 23, 2, 42, 45, 56, 80, 102, 98, 63, 56, }

        // sort the first half
        my_int_pqsort(values, nvalues, 0, 5);
        DUMP(); // => { 2, 23, 42, 45, 56, 102, 98, 63, 80, 56, }

        // sort the rest
        my_int_pqsort(values + 5, nvalues - 5, 0, nvalues);
        DUMP(); // => { 2, 23, 42, 45, 56, 56, 63, 80, 98, 102, }

        fprintf(stdout, "-->%d\n", quick_select_median(values, 10));

        return 0;
}
