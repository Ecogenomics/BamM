//#############################################################################
//
//   bamFilter.c
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

// system includes
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
//#include <libgen.h>

// htslib
#include "sam.h"

// local includes
#include "bamFilter.h"

void filterReads(char * inBamFile,
                 char * outBamFile,
                 int minMapQual,
                 int minLen,
                 int maxMisMatches,
                 float minPcId,
                 float minPcAln,
                 int ignoreSuppAlignments,
                 int ignoreSecondaryAlignments) {
    //
    int result = -1;
    int outResult = -1;

    int supp_check = 0x0;
    if (ignoreSuppAlignments) {
        supp_check |= BAM_FSUPPLEMENTARY;
    }
    if (ignoreSecondaryAlignments) {
        supp_check |= BAM_FSECONDARY;
    }

    // helper variables
    bamFile *in = 0;
    bamFile *out = 0;
    bam1_t *b = bam_init1();
    bam_hdr_t *h;

    // open bam
    if ((in = bgzf_open(inBamFile, "r")) == 0) {
        fprintf(stderr,
               "ERROR: Failed to open \"%s\" for reading.\n",
               inBamFile);
    }
    else if ((h = bam_hdr_read(in)) == 0) {
        fprintf(stderr,
                "ERROR: Failed to read BAM header of file \"%s\".\n",
                inBamFile);
    }
    else if ((out = bgzf_open(outBamFile, "w")) == 0) {
        fprintf(stderr,
               "ERROR: Failed to open \"%s\" for writing.\n",
               outBamFile);
    }
    else {
        bam_hdr_write(out, h);
        bam_hdr_destroy(h);
        // fetch alignments
        int line = 0;
        int matches, mismatches, qLen;
        float pcAln, pcId;
        while ((result = bam_read1(in, b)) >= 0) {
            line += 1;

            // only primary mappings
            if ((b->core.flag & supp_check) != 0)
                continue;

            // only high quality
            if (b->core.qual < minMapQual)
                continue;

            // not too many absolute mismatches
            mismatches = bam_aux2i(bam_aux_get(b, "NM"));
            if (mismatches > maxMisMatches)
                continue;

            // not too short
            qLen = bam_cigar2qlen((&b->core)->n_cigar, bam_get_cigar(b));
            if (qLen < minLen)
                continue;

            // only high percent identity
            matches = bam_cigar2matches((&b->core)->n_cigar, bam_get_cigar(b));
            pcId = (matches - mismatches) / (float)matches; // percentage as float between 0 to 1
            if (pcId < minPcId)
                continue;

            // only high percent alignment
            pcAln = matches / (float)qLen; // percentage as float between 0 to 1
            if (pcAln < minPcAln)
                continue;

            if ((outResult = bam_write1(out, b)) < -1) {
                fprintf(stderr,
                        "ERROR: Attempt to write read no. %d to file \"%s\" failed with code %d.\n",
                        line, outBamFile, outResult);
            }

            bam_destroy1(b);
            b = bam_init1();
        }
        if (result < -1) {
            fprintf(stderr,
                    "ERROR: retrieval of read no. %d from file \"%s\" failed with code %d.\n",
                    line, inBamFile, result);
        }
    }
    if (in) bgzf_close(in);
    if (out) bgzf_close(out);
}


int bam_cigar2matches(int n_cigar, const uint32_t *cigar)
{
    int k, l;
    for (k = l = 0; k < n_cigar; ++k)
        if (bam_cigar_type(bam_cigar_op(cigar[k]))==3)
            l += bam_cigar_oplen(cigar[k]);
    return l;
}
