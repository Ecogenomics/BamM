//#############################################################################
//
//   bamFilter.h
//
//   Filter reads from BAM files
//
//   Copyright (C) Tim Lamberton
//
//   This library is free software; you can redistribute it and/or
//   modify it under the terms of the GNU Lesser General Public
//   License as published by the Free Software Foundation; either
//   version 3.0 of the License, or (at your option) any later version.
//
//   This library is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//   Lesser General Public License for more details.
//
//   You should have received a copy of the GNU Lesser General Public
//   License along with this library.
//
//#############################################################################

#ifndef BM_BAM_FILTER_H
   #define BM_BAM_FILTER_H


// system includes
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
//#include <libgen.h>

// htslib
#include "sam.h"

typedef BGZF bamFile;

#ifdef __cplusplus
extern "C" {
#endif


/*!
 * @abstract Filter reads in a BAM file
 *
 * @param  outFile              filename of filtered BAM file
 * @param  bamFile              filename of BAM file to filter
 * @param  minMapQual           mapping quality threshold
 * @param  minLen               min query length
 * @param  maxMisMatches        maximum number of mismatches to accept (NM flag)
 * @param  minPcId              minimum percentage identity to accept (int between 0 and 100)
 * @param  minPcAln             minimum percentage alignment to accept (int between 0 and 100)
 * @param  invertMatch == 1 -> select unmapped reads
 * @param  ignoreSuppAlignments == 1 -> ignore supplmentary alignments
 * @param  ignoreSecondaryAlignments  == 1 -> ignore secondary alignments
 * @return 0 for success
 *
 */
void filterReads(char * bamFile,
                 char * outFile,
                 int mapQ,
                 int minLen,
                 int maxMisMatches,
                 float minPcId,
                 float minPcAln,
                 int invertMatch,
                 int ignoreSuppAlignments,
                 int ignoreSecondaryAlignments);

int bam_cigar2matches(int n_cigar, const uint32_t *cigar);

#ifdef __cplusplus
}
#endif

#endif // BM_BAM_FILTER_H
