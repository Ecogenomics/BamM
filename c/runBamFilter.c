// runBamFilter.c
//
// Tim Lamberton

#include <stdio.h>
#include "bamFilter.h"


int main( int argc, char *argv[] ) {
    if ( argc != 2 ) {
        printf( "USAGE: %s BAM_FILE\n", argv[0]);
        return 1;
    }
    else {
        filterReads(argv[1], "out.bam", 0, 0, 1000, 0.9, 0.9, 0, 0);
        return remove("out.bam");
    }
}
