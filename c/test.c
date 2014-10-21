// compile with: gcc -std=gnu99 test.c

#include <stdio.h>
#include "pqsort.h"

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

        return 0;
}
