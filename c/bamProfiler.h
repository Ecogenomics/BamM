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

#ifndef BM_BAM_PROFILER_H
   #define BM_BAM_PROFILER_H


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
 * @abstract Profile reads in a BAM file
 *
 * @param  bamFile              filename of profiled BAM file
 * @param  ignoreSuppAlignments == 1 -> ignore supplmentary alignments
 * @param  ignoreSecondaryAlignments  == 1 -> ignore secondary alignments
 * @return 0 for success
 *
 */
void profileReads(char* bamFile,
				  int ignoreSuppAlignments,
				  int ignoreSecondaryAlignments);

#ifdef __cplusplus
}
#endif

#endif // BM_BAM_PROFILER_H
