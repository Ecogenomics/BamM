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
#include "bgzf.h"
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
                 int invertMatch,
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
    BGZF* in = 0;
    BGZF* out = 0;
    bam1_t *b = bam_init1();
    bam_hdr_t *h;

    // open bam
    if ((in = bgzf_open(inBamFile, "r")) == 0) {
        fprintf(stderr,
               "ERROR: Failed to open \"%s\" for reading.\n",
               inBamFile);
    }
    else if ((h = bam_hdr_read(in)) == 0) { // read header
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
        // write and destroy header
        bam_hdr_write(out, h);
        bam_hdr_destroy(h);

        int line = 0;
        int matches, mismatches, qLen, unmapped;
        float pcAln, pcId;
        int showStats = 0;
        uint8_t *aux_mismatches;

        // fetch alignments
        while ((result = bam_read1(in, b)) >= 0) {
            line += 1;
            unmapped = 0;
            
            // only primary mappings even if inverting matches
            if ((b->core.flag & supp_check) != 0) { 
                if (showStats)
                    fprintf(stdout, "Rejected %d, non-primary\n", line);
                continue;
            }
            
            // only high quality
            if ((unmapped = b->core.qual < minMapQual)) {
                if (showStats)
                    fprintf(stdout, "Rejected %d, quality: %d\n", line, b->core.qual);
            }
            
            if (unmapped == 0) {
                // bam_aux_get returns 0 if optional NM tag is missing
                if ((aux_mismatches = bam_aux_get(b, "NM")))
                   mismatches = bam_aux2i(aux_mismatches);
                else
                    mismatches = 0;
                // not too many absolute mismatches
                if ((unmapped = mismatches > maxMisMatches)) {
                    if (showStats)
                        fprintf(stdout, "Rejected %d, mismatches: %d\n", line, mismatches);
                }
            }
            if (unmapped == 0) {
                // not too short
                qLen = bam_cigar2qlen((&b->core)->n_cigar, bam_get_cigar(b));
                if ((unmapped = qLen < minLen)) {
                    if (showStats)
                        fprintf(stdout, "Rejected %d, length: %d\n", line, qLen);
                }
            }
            if (unmapped == 0) {
                // only high percent identity
                matches = bam_cigar2matches((&b->core)->n_cigar, bam_get_cigar(b));
                pcId = (matches - mismatches) / (float)matches; // percentage as float between 0 to 1
                if ((unmapped = pcId < minPcId)) {
                    if (showStats)
                        fprintf(stdout, "Rejected %d, identity pc: %.4f\n", line, pcId);
                }
            }
            if (unmapped == 0) {
                // only high percent alignment
                pcAln = matches / (float)qLen; // percentage as float between 0 to 1
                if ((unmapped = pcAln < minPcAln)) {
                    if (showStats)
                        fprintf(stdout, "Rejected %d, alignment pc: %.4f\n", line, pcAln);
                }
            }

            // write mapped reads (unmapped==0) when invertMatch==0
            // write unmapped reads (unmapped==1) when invertMatch==1
            if (unmapped == invertMatch) {
                if ((outResult = bam_write1(out, b)) < -1) {
                    fprintf(stderr,
                            "ERROR: Attempt to write read no. %d to file \"%s\" failed with code %d.\n",
                            line, outBamFile, outResult);
                }
            }
        }
        if (result < -1) {
            fprintf(stderr,
                    "ERROR: retrieval of read no. %d from file \"%s\" failed with code %d.\n",
                    line, inBamFile, result);
        }
    }
    if (in) bgzf_close(in);
    if (out) bgzf_close(out);
    bam_destroy1(b);
}


int bam_cigar2matches(int n_cigar, const uint32_t *cigar)
{
    int k, l;
    for (k = l = 0; k < n_cigar; ++k)
        // bam_cigar_type bit flag == 3 implies a match to both
        // the reference and query
        if (bam_cigar_type(bam_cigar_op(cigar[k]))==3)
            l += bam_cigar_oplen(cigar[k]);
    return l;
}
