// runBamProfiler.c
//
// Tim Lamberton

#include <stdio.h>
#include "bamProfiler.h"


int main( int argc, char *argv[] ) {
    if ( argc != 2 ) {
        printf( "USAGE: %s BAM_FILE\n", argv[0]);
        return 1;
    }
    else {
        profileReads(argv[1], 0, 0);
        profileReads(argv[1], 1, 1);
        return 0;
    }
}
