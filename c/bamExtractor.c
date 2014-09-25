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

// cfuhash
#include "cfuhash.h"

// local includes
#include "bamParser.h"
#include "bamExtractor.h"
#include "bamRead.h"

BM_mappedRead * extractReads(char * bamFile,
                             char ** contigs,
                             int numContigs,
                             uint16_t * groups,
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
                cfuhash_table_t *pair_buffer = cfuhash_new_with_initial_size(1000000);
                cfuhash_set_flag(pair_buffer, CFUHASH_FROZEN_UNTIL_GROWS);

                for (hh = 0; hh < numContigs; ++hh) {
                    hts_itr_t *iter = sam_itr_querys(idx, header, contigs[hh]); // parse a region in the format like `chr2:100-200'
                    if (iter == NULL) { // reference name is not found
                        fprintf(stderr, "WARNING: Could not find contig: [%s] in BAM: [%s].\n", contigs[hh], bamFile);
                    }

                    // fetch alignments
                    while ((result = sam_itr_next(in, iter, b)) >= 0) {
                        bam1_core_t core = b->core;
                        char * seqId = bam_get_qname(b);
                        char * seq = 0;
                        char * qual = 0;
                        int qual_len = 0;
                        int seq_len = 0;

                        // get sequence and quality
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

                        // work out pairing information
                        uint8_t rpi = RPI_ERROR;
                        if (core.flag&BAM_FPAIRED) {
                            if(core.flag&BAM_FMUNMAP) {
                                if (core.flag&BAM_FREAD1) {
                                    rpi = RPI_SNGL_FIR;
                                }
                                else if (core.flag&BAM_FREAD2) {
                                    rpi = RPI_SNGL_SEC;
                                }
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

                        // make the funky Id
                        #define MAX_SEQ_ID_LEN 80
                        char * seq_id = calloc(MAX_SEQ_ID_LEN, sizeof(char));
                        // allocate the string to the buffer but check to make sure we're not cutting anything off
                        int id_len = snprintf(seq_id, MAX_SEQ_ID_LEN, "b_%s;c_%s;r_%s", prettyName, contigs[hh], seqId);
                        if(id_len >= MAX_SEQ_ID_LEN) {
                            seq_id = calloc(id_len+1, sizeof(char));
                            snprintf(seq_id, id_len+1, "b_%s;c_%s;r_%s", prettyName, contigs[hh], seqId);
                        }

                        // make the mapped read struct
                        prev = makeMappedRead(seq_id,
                                              seq,
                                              qual,
                                              id_len,
                                              seq_len,
                                              qual_len,
                                              rpi,
                                              groups[hh],
                                              prev);

                        if (0 == root) { root = prev; }

                        if(rpi == RPI_SNGL || rpi == RPI_SNGL_FIR || rpi == RPI_SNGL_SEC) {
                            // we can just add away -> indicate singleton reads by pointing
                            // the parter pointer to itself
                            prev->partnerRead = prev;
                        }
                        else {
                            // RPI_FIR or RPI_SEC
                            // now work out pairing information using the hash table
                            char * stripped_result = pairStripper(seqId, core.l_qname-1);
                            char * stripped = seqId;
                            size_t stripped_len = core.l_qname-1;
                            if(stripped_result)
                                stripped = stripped_result;
                                stripped_len = core.l_qname-3;

                            // now stripped always holds a stripped value
                            // see if it is in the hash already
                            BM_mappedRead * stored_MR = cfuhash_get(pair_buffer, stripped);
                            if (stored_MR != NULL) {
                                // exists in the hash -> We can add the pair info
                                if(rpi == RPI_FIR)
                                    prev->partnerRead = stored_MR;
                                else
                                    stored_MR->partnerRead = prev;

                                // delete the entry from the hash
                                cfuhash_delete_data(pair_buffer, stripped, stripped_len);
                            }
                            else {
                                // we should put it in the hash
                                cfuhash_put(pair_buffer, stripped, prev);
                            }

                            if(stripped_result) // free this!
                                free(stripped_result);
                        }
                    }
                    hts_itr_destroy(iter);
                    if (result < -1) {
                        fprintf(stderr, "ERROR: retrieval of reads from contig: \"%s\" failed due to truncated file or corrupt BAM index file\n", header->target_name[hh]);
                        break;
                    }
                }
                cfuhash_clear(pair_buffer);
                cfuhash_destroy(pair_buffer);
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

char * pairStripper(char * longId, int idLen) {
    // we're only looking at the second last character
    // we only care about '.' or '_'
    --idLen;
    if(longId[idLen] == '1' || longId[idLen] == '2') {
        --idLen;
        if(longId[idLen] == '_' || longId[idLen] == '.') {
            // if we're here then we need to strip
            char * ret_str = strdup(longId); // malloc two extra chars
            ret_str[idLen] = 0;              // chuck the null in early
            return ret_str;
        }
    }
    return 0;
}
