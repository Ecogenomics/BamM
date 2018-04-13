// runBamParser.c
//
// Tim Lamberton

#include <stdio.h>
#include "bamParser.h"
#include "coverageEstimators.h"


int main( int argc, char *argv[] ) {
    int result;
    int doLinks = 1;
    int types_c_array[1];
    BM_coverageType BCT;
    BM_coverageType *pBCT;
    char* bamfiles_c_array[1];
    BM_fileInfo BFI;
    BM_fileInfo *pBFI;
    
    if ( argc != 2 ) {
        printf( "USAGE: %s BAM_FILE\n", argv[0]);
        return 1;
    }
    else {
        types_c_array[0] = 1;
        BCT.type = CT_P_MEAN_TRIMMED;
        BCT.upperCut = 0.0;
        BCT.lowerCut = 0.0;
        pBCT = &BCT;
        bamfiles_c_array[0] = argv[1];
        pBFI = &BFI;
        if (doLinks) {
            result = parseCoverageAndLinks(1, // doLinks
                                           1, // doCov
                                           1, // numBams
                                           0, // baseQuality
                                           0, // mappingQuality
                                           0, // minLength
                                           1000, // maxMismatches
                                           types_c_array,
                                           1, // ignoreSuppAlignments,
                                           1, // ignoreSecondaryAlignments,
                                           pBCT,
                                           bamfiles_c_array,
                                           pBFI
                                          );
        }
        else {
            // types only
            BCT.type = CT_NONE;
            result = parseCoverageAndLinks(0,
                                           0,
                                           1,
                                           0,
                                           0,
                                           0,
                                           0,
                                           types_c_array,
                                           1,
                                           1,
                                           pBCT,
                                           bamfiles_c_array,
                                           pBFI
                                          );
        }
        return result;
    }
}
