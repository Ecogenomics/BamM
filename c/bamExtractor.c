//#############################################################################
//
//   bamExtractor.c
//
//   Engine for extracting reads from BAM files - non-parallel
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

// system includes
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <libgen.h>

// htslib
//#include "htslib/bgzf.h"
#include "htslib/sam.h"

// local includes
#include "bamParser.h"
#include "bamExtractor.h"

    /***********************
    ***      READS       ***
    ***********************/

BM_mappedRead * extractReads(char * bamFile,
                             char ** contigs,
                             int numContigs,
                             uint16_t * bins,
                             char * prettyName,
                             int headersOnly) {
    //-----
    // code uses the pattern outlined in samtools view (sam_view.c)
    // thanks lh3!
    //
    int i = 0;
    int result = -1;
    int hh = 0;

    // helper variables
    samFile *in = 0;
    bam_hdr_t *header = NULL;
    bam1_t *b = bam_init1();

    BM_mappedRead * root = 0;
    BM_mappedRead * prev = 0;

    // open file handlers
    if ((in = sam_open(bamFile, "r")) == 0) {
        fprintf(stderr, "ERROR: Failed to open \"%s\" for reading.\n", bamFile);
    }
    else {
        // retrieve the header
        if ((header = sam_hdr_read(in)) == 0) {
            fprintf(stderr, "ERROR: Failed to read the header from \"%s\".\n", bamFile);
        }
        else {
            // check the index is intact
            hts_idx_t *idx = sam_index_load(in, bamFile); // load index
            if (idx == 0) { // index is unavailable
                fprintf(stderr, "ERROR: Random alignment retrieval only works for indexed BAM or CRAM files.\n");
            }
            else {
                for (hh = 0; hh < numContigs; ++hh) {
                    hts_itr_t *iter = sam_itr_querys(idx, header, contigs[hh]); // parse a region in the format like `chr2:100-200'
                    if (iter == NULL) { // reference name is not found
                        fprintf(stderr, "WARNING: Could not find contig: [%s] in BAM: [%s].\n", contigs[hh], bamFile);
                    }

                    // fetch alignments
                    while ((result = sam_itr_next(in, iter, b)) >= 0) {
                        bam1_core_t core = b->core;
                        char * seq = 0;
                        char * qual = 0;
                        int qual_len = 0;
                        int seq_len = 0;
                        if(0 == headersOnly) {  // no point allocating unused space
                            seq = calloc(core.l_qseq+1, sizeof(char));
                            qual = calloc(core.l_qseq+1, sizeof(char));
                            uint8_t *s = bam_get_seq(b);
                            for (i = 0; i < core.l_qseq; ++i) { seq[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(s, i)]; }
                            seq_len = core.l_qseq;

                            s = bam_get_qual(b);
                            if (s[0] != 0xff) {
                                qual_len = core.l_qseq;
                                for (i = 0; i < core.l_qseq; ++i) { qual[i] = (char)(s[i] + 33); }
                            }
                            else {
                                free(qual);
                                qual = 0;
                            }
                        }

                        #define MAX_SEQ_ID_LEN 80
                        char * seq_id = calloc(MAX_SEQ_ID_LEN, sizeof(char));
                        // allocate the string to the buffer but check to make sure we're not cutting anything off
                        int id_len = snprintf(seq_id, MAX_SEQ_ID_LEN, "f_%s;c_%s;r_%s", prettyName, contigs[hh], bam_get_qname(b));
                        if(id_len >= MAX_SEQ_ID_LEN) {
                            seq_id = calloc(id_len+1, sizeof(char));
                            snprintf(seq_id, id_len+1, "f_%s;c_%s;r_%s", prettyName, contigs[hh], bam_get_qname(b));
                        }

                        uint8_t rpi = RPI_ERROR;
                        if (core.flag&BAM_FPAIRED) {
                            if(core.flag&BAM_FMUNMAP) {
                                rpi = RPI_SNGL;
                            }
                            else {
                                if (core.flag&BAM_FREAD1) {
                                    rpi = RPI_FIR;
                                }
                                else if (core.flag&BAM_FREAD2) {
                                    rpi = RPI_SEC;
                                }
                            }
                        }
                        else {
                            rpi = RPI_SNGL;
                        }
                        //fprintf(stdout, "%d -- %s -- %d \n", rpi, bam_get_qname(b), bins[hh]); fflush(stdout);
                        prev = makeMappedRead(seq_id,
                                              seq,
                                              qual,
                                              id_len,
                                              seq_len,
                                              qual_len,
                                              rpi,
                                              bins[hh],
                                              prev);

                        if (0 == root) { root = prev; }
                    }
                    hts_itr_destroy(iter);
                    if (result < -1) {
                        fprintf(stderr, "ERROR: retrieval of reads from contig: \"%s\" failed due to truncated file or corrupt BAM index file\n", header->target_name[hh]);
                        break;
                    }
                }
            }
            hts_idx_destroy(idx); // destroy the BAM index
        }
    }
    // always do this
    if (in) sam_close(in);
    bam_destroy1(b);
    if ( header ) bam_hdr_destroy(header);

    return root;
}

BM_mappedRead * makeMappedRead(char * seqId,
                                char * seq,
                                char * qual,
                                uint16_t idLen,
                                uint16_t seqLen,
                                uint16_t qualLen,
                                uint8_t rpi,
                                uint16_t bin,
                                BM_mappedRead * prev_MR) {
    BM_mappedRead * MR = calloc(1, sizeof(BM_mappedRead));
    MR->seqId = strdup(seqId);
    if( 0 != seq ) {
    	MR->seq = strdup(seq);
     }
     else {
     	MR->seq = 0;
     }
    MR->rpi = rpi;
    MR->idLen = idLen;
    MR->seqLen = seqLen;
    MR->qualLen = qualLen;
    MR->bin = bin;
    if( 0 != qual )
        MR->qual = strdup(qual);
    else
        MR->qual = 0;

    if(prev_MR)
        prev_MR->nextRead = MR;

    return MR;
}

BM_mappedRead * nextMappedRead(BM_mappedRead * MR) {
    return MR->nextRead;
}


void destroyMappedReads(BM_mappedRead * root_MR) {
    while(root_MR) {
        if(root_MR->seqId) free(root_MR->seqId);
        if(root_MR->seq) free(root_MR->seq);
        if(root_MR->qual) free(root_MR->qual);
        BM_mappedRead * tmp_MR = root_MR->nextRead;
        free(root_MR);
        root_MR = tmp_MR;
    }
}

void printMappedReads(BM_mappedRead * root_MR, FILE * f) {

    if(0 == f) { f = stdout; }

    BM_mappedRead * MR = root_MR;
    if(MR->qual) {
        while(MR) {
            fprintf(f, "@b_%d;%s\n%s\n+\n%s\n", MR->bin, MR->seqId, MR->seq, MR->qual);
            MR = MR->nextRead;
        }
    }
    else {
        while(MR) {
            fprintf(f, ">b_%d;%s\n%s\n", MR->bin, MR->seqId, MR->seq);
            MR = MR->nextRead;
        }
    }
}

