//#############################################################################
//
//   bamExtractor.h
//
//   Extract reads from BAM files given a set of targets
//
//   Copyright (C) Michael Imelfort
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

#ifndef BM_BAM_EXTRACTOR_H
#define BM_BAM_EXTRACTOR_H

// system includes
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>

// htslib
#include "htslib/bgzf.h"
#include "htslib/sam.h"

// local includes
#include "bamParser.h"
#include "bamRead.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * @abstract Extract reads from a BAM file based on which contigs they map to
 *
 * @param  bamFile          BAM filename to extract from
 * @param  contigs          array of contigs (targets) to extract reads for
 * @param  numContigs       number of contigs (targets) to extract reads for
 * @param  groups           array of IDs indicating the target group each contig belongs to
 * @param  prettyName       string used for displaying BAM file name in read header
 * @param  headersOnly      == 0 -> extract only headers otherwise do seq + qual
 * @return BM_mappedRead *  start of list of mapped reads
 *
 * @discussion This function expects BFI to be initialised.
 */
BM_mappedRead * extractReads(char * bamFile,
                             char ** contigs,
                             int numContigs,
                             uint16_t * groups,
                             char * prettyName,
                             int headersOnly);

/*!
 * @abstract Strip read pair information from a read Id
 *
 * @param  longId   the Id to try strip
 * @param  idLen    length of the Id
 * @return char *   the stripped Id or 0 if not stripped
 */
char * pairStripper(char * longId, int idLen);

#ifdef __cplusplus
}
#endif

#endif // BM_BAM_EXTRACTOR_H