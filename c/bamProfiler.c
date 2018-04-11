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
#include "bamProfiler.h"


void profileReads(char* bamFile,
                  int ignoreSuppAlignments,
                  int ignoreSecondaryAlignments) {
    //
    int result = -1;

    int supp_check = 0x0;
    if (ignoreSuppAlignments) {
        supp_check |= BAM_FSUPPLEMENTARY;
    }
    if (ignoreSecondaryAlignments) {
        supp_check |= BAM_FSECONDARY;
    }

    // helper variables
    BGZF* in = 0 ;
    bam1_t *b = bam_init1();
    bam_hdr_t *h;

    // open bam
    if ((in = bgzf_open(bamFile, "r")) == 0) {
        fprintf(stderr,
               "ERROR: Failed to open \"%s\" for reading.\n",
               bamFile);
    }
    else if ((h = bam_hdr_read(in)) == 0) { // read header
        fprintf(stderr,
                "ERROR: Failed to read BAM header of file \"%s\".\n",
                bamFile);
    }
    else {
        // destroy header
        bam_hdr_destroy(h);

        int line = 0;
        
        int supplementary, secondary;
        int mapQual;
        int matches, mismatches, qLen;
        float pcAln, pcId;
        int showStats = 0;
        uint8_t *aux_mismatches;
        
        // print header
        printf("line\tsupp\tsecondary\tmapQ\tmismatches\tmatches\tqLen\tpcId\tpcAln\n");

        // fetch alignments
        while ((result = bam_read1(in, b)) >= 0) {
            line += 1;
            
            
            // only primary mappings
            if ((b->core.flag & supp_check) != 0) { 
                if (showStats)
                    fprintf(stdout, "Rejected %d, non-primary\n", line);
                continue;
            }
            supplementary = (b->core.flag & (1 | BAM_FSUPPLEMENTARY)) != 0;
            secondary = (b->core.flag & (1 | BAM_FSECONDARY)) != 0;
            // quality
            mapQual = b->core.qual;
            // bam_aux_get returns 0 if optional NM tag is missing
            if ((aux_mismatches = bam_aux_get(b, "NM")))
               mismatches = bam_aux2i(aux_mismatches);
            else
                mismatches = 0;
            // length
            qLen = bam_cigar2qlen((&b->core)->n_cigar, bam_get_cigar(b));
            // percent identity
            matches = bam_cigar2matches((&b->core)->n_cigar, bam_get_cigar(b));
            pcId = (matches - mismatches) / (float)matches; // percentage as float between 0 to 1
            // percent alignment
            pcAln = matches / (float)qLen; // percentage as float between 0 to 1
            
            // print read values
            printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\n",
                   line, supplementary, secondary, mapQual, mismatches, matches,
                   qLen, pcId, pcAln);
        }
        if (result < -1) {
            fprintf(stderr,
                    "ERROR: retrieval of read no. %d from file \"%s\" failed with code %d.\n",
                    line, bamFile, result);
        }
    }
    if (in) bgzf_close(in);
    bam_destroy1(b);
}
