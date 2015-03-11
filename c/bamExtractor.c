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
#include "sam.h"

// cfuhash
#include "cfuhash.h"

// local includes
#include "bamParser.h"
#include "bamExtractor.h"
#include "bamRead.h"

/*

 JUST SOME DEBRI THAT I'M SURE I'LL NEED LATER

if (core.n_cigar) { // cigar
    fprintf(stdout, "-->[[");
    uint32_t *cigar = bam_get_cigar(b);
    for (i = 0; i < core.n_cigar; ++i) {
        fprintf(stdout, "%d", bam_cigar_oplen(cigar[i]));
        fprintf(stdout, "%c", bam_cigar_opchr(cigar[i]));
    }
    fprintf(stdout, "]]<--\n");
}
*/

const char *samFlagToBinary(int x)
{
    static char b[13];
    b[0] = '\0';

    int z;
    for (z = 2048; z > 0; z >>= 1)
    {
        strcat(b, ((x & z) == z) ? "1" : "0");
    }

    return b;
}

void printPairCorruptionWarning(int done) {
    if(!done) {
        fprintf(stderr, "WARNING: read pairings may be corrupted. Most likely "\
                        "due to use of secondary / supplemental mappings\n");
    }
}

BM_mappedRead * extractReads(char * bamFile,
                             char ** contigs,
                             int numContigs,
                             uint16_t * groups,
                             char * prettyName,
                             int headersOnly,
                             int minMapQual,
                             int maxMisMatches,
                             int ignoreSuppAlignments,
                             int ignoreSecondaryAlignments) {
    //-----
    // code uses the pattern outlined in samtools view (sam_view.c)
    // thanks lh3!
    //
    int i = 0;
    int result = -1;
    int hh = 0;

    int supp_check = 0x0; // include supp mappings
    if (ignoreSuppAlignments) {
        supp_check |= BAM_FSUPPLEMENTARY;
    }
    if (ignoreSecondaryAlignments) {
        supp_check |= BAM_FSECONDARY;
    }

    // we need to let the users know if their pairings
    // will be corrupted
    int p_corrupt = 0;

    // helper variables
    samFile *in = 0;
    bam_hdr_t *header = NULL;
    bam1_t *b = bam_init1();

    BM_mappedRead * root = 0;
    BM_mappedRead * prev = 0;

    // open file handlers
    if ((in = sam_open(bamFile, "r")) == 0) {
        fprintf(stderr,
                "ERROR: Failed to open \"%s\" for reading.\n",
                bamFile);
    }
    else {
        // retrieve the header
        if ((header = sam_hdr_read(in)) == 0) {
            fprintf(stderr,
                    "ERROR: Failed to read the header from \"%s\".\n",
                    bamFile);
        }
        else {
            // check the index is intact
            hts_idx_t *idx = sam_index_load(in, bamFile); // load index
            if (idx == 0) { // index is unavailable
                fprintf(stderr,
                        "ERROR: Random retrieval only works "\
                        "for indexed files.\n");
            }
            else {
                cfuhash_table_t *pair_buffer = \
                    cfuhash_new_with_initial_size(1000000);
                cfuhash_set_flag(pair_buffer, CFUHASH_FROZEN_UNTIL_GROWS);

                for (hh = 0; hh < numContigs; ++hh) {
                    // parse a region in the format like `chr2:100-200'
                    hts_itr_t *iter = sam_itr_querys(idx, header, contigs[hh]);
                    if (iter == NULL) { // reference name is not found
                        fprintf(stderr,
                                "WARNING: Could not find contig: "\
                                "[%s] in BAM: [%s].\n",
                                contigs[hh],
                                bamFile);
                    }

                    // fetch alignments
                    int line = 0;
                    while ((result = sam_itr_next(in, iter, b)) >= 0) {
                        bam1_core_t core = b->core;
                        line += 1;
                        // only high quality?, primary? mappings
                        if ( core.qual < minMapQual)
                            continue;
                        if ((core.flag & supp_check) != 0)
                            continue;
                        if(bam_aux2i(bam_aux_get(b, "NM")) > maxMisMatches) {
                            continue;
                        }

                        char * seqId = bam_get_qname(b);
                        char * seq = 0;
                        char * qual = 0;
                        int qual_len = 0;
                        int seq_len = 0;

                        // get sequence and quality
                        if(0 == headersOnly) {
                            // no point allocating unused space
                            seq = calloc(core.l_qseq+1, sizeof(char));
                            qual = calloc(core.l_qseq+1, sizeof(char));
                            uint8_t *s = bam_get_seq(b);
                            if (core.flag&BAM_FREVERSE) {
                                // reverse the read
                                int r = 0;
                                for (i = core.l_qseq-1; i >=0 ; --i) {
                                    seq[r]="=TGKCYSBAWRDMHVN"[bam_seqi(s,
                                                                       i)];
                                    ++r;
                                }
                            }
                            else {
                                for (i = 0; i < core.l_qseq; ++i) {
                                    seq[i]="=ACMGRSVTWYHKDBN"[bam_seqi(s,
                                                                       i)];
                                }
                            }
                            seq_len = core.l_qseq;

                            s = bam_get_qual(b);
                            if (s[0] != 0xff) {
                                qual_len = core.l_qseq;
                                for (i = 0; i < core.l_qseq; ++i) {
                                    qual[i] = (char)(s[i] + 33);
                                }
                            }
                            else if (qual != 0) {
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
                        char * seq_id = calloc(MAX_SEQ_ID_LEN,
                                               sizeof(char));
                        // allocate the string to the buffer but check to
                        // ensure we're not cutting anything off
                        int id_len = snprintf(seq_id,
                                              MAX_SEQ_ID_LEN,
                                              "b_%s;c_%s;r_%s",
                                              prettyName,
                                              contigs[hh],
                                              seqId);
                        if(id_len >= MAX_SEQ_ID_LEN) {
                            seq_id = calloc(id_len+1, sizeof(char));
                            snprintf(seq_id,
                                     id_len+1, // don't forget the NULL!
                                     "b_%s;c_%s;r_%s",
                                     prettyName,
                                     contigs[hh],
                                     seqId);
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

                        if(rpi == RPI_SNGL || \
                           rpi == RPI_SNGL_FIR || \
                           rpi == RPI_SNGL_SEC) {
                            // we can just add away
                            // indicate singleton reads by pointing the
                            // partner pointer to itself
                            prev->partnerRead = prev;
                        }
                        else {
                            // RPI_FIR or RPI_SEC
                            // work out pairing information using the hash
                            // we append a 1 or 2 to the end so that
                            // we don't accidentally pair 1's with 1's etc.
                            char * stripped_result;
                            if(rpi == RPI_FIR) {
                                stripped_result = \
                                    pairStripper(seqId,
                                                 core.l_qname-1,
                                                 '2');
                            }
                            else {
                                stripped_result = \
                                    pairStripper(seqId,
                                                 core.l_qname-1,
                                                 '1');
                            }

                            char * stripped = seqId;
                            if(stripped_result)
                                stripped = stripped_result;

                            //fprintf(stdout, "SEARCH %s\n", stripped);
                            // now stripped always holds a stripped value
                            // see if it is in the hash already
                            BM_mappedRead * stored_MR = \
                                cfuhash_get(pair_buffer,
                                            stripped);

                            if (0 != stored_MR) {
                                // exists in the hash -> Add the pair info
                                if(rpi == RPI_FIR) {
                                    prev->partnerRead = stored_MR;
                                }
                                else {
                                    stored_MR->partnerRead = prev;
                                }

                                // delete the entry from the hash
                                cfuhash_delete(pair_buffer,
                                               stripped);
                            }
                            else {
                                // we should put it in the hash
                                // make sure to change it into something
                                // we will find next time
                                if(rpi == RPI_FIR)
                                    stripped[strlen(stripped)-1] = '1';
                                else
                                    stripped[strlen(stripped)-1] = '2';

                                // check to make sure we're not overwriting
                                // anything important. cfuhash overwrites
                                // duplicate entries, so we need to grab
                                // it and put it to "SNGL_XXX" before we
                                // lose the pointer
                                BM_mappedRead * OWMMR = \
                                    cfuhash_put(pair_buffer,
                                                stripped, prev);
                                if(OWMMR) {
                                    if(OWMMR->rpi == RPI_FIR)
                                        OWMMR->rpi = RPI_SNGL_FIR;
                                    else
                                        OWMMR->rpi = RPI_SNGL_SEC;
                                    OWMMR->partnerRead = OWMMR;
                                    printPairCorruptionWarning(p_corrupt);
                                    p_corrupt = 1;
                                }


                            }

                            if(stripped_result != 0) { // free this!
                                free(stripped_result);
                                stripped_result = 0;
                            }
                        }
                    }
                    hts_itr_destroy(iter);
                    if (result < -1) {
                        fprintf(stderr, "ERROR: retrieval of reads from "\
                                        "contig:  \"%s\" failed due to "\
                                        "truncated file or corrupt BAM index "\
                                        "file\n", header->target_name[hh]);
                        break;
                    }
                }

                // any entries left in the hash are pairs whose mates did
                // not meet quality standards
                size_t key_size = 0;
                char * key;
                BM_mappedRead * LOMMR;
                size_t pr_size = 1;
                if(cfuhash_each_data(pair_buffer,
                                     (void**)&key,
                                     &key_size,
                                     (void**)&LOMMR,
                                     &pr_size)) {
                    do {
                        // get the mapped read
                        // update it's pairing so we know it's really single
                        if (LOMMR->rpi == RPI_FIR)
                            LOMMR->rpi = RPI_SNGL_FIR;
                        else if (LOMMR->rpi == RPI_SEC)
                            LOMMR->rpi = RPI_SNGL_SEC;

                        // indicate singleton reads by pointing the
                        // partner pointer to itself
                        LOMMR->partnerRead = LOMMR;

                    } while(cfuhash_next_data(pair_buffer,
                                              (void**)&key,
                                              &key_size,
                                              (void**)&LOMMR,
                                              &pr_size));
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

char * pairStripper(char * longId, int idLen, char target) {
    // we're only looking at the second last character
    // we only care about '.' or '_'
    // idLen initially points to the null terminator
    --idLen;
    if(longId[idLen] == '1' || longId[idLen] == '2') {
        // looks like the read id may end in a paired read identifier
        --idLen;
        if(longId[idLen] == '_' ||
           longId[idLen] == '.' ||
           longId[idLen] == '/') {
            // it does!
            // we need to strip the old identifier and add the target instead
            char * ret_str = strdup(longId); // malloc two extra chars
            ret_str[idLen] = '_';
            ret_str[idLen+1] = target;
            return ret_str;
        }
        else {
            // the read id just happened to end in a 0 or a 1
            // undo the last decrement
            ++idLen;
        }
    } // else there is no identifier

    // if we get here then there is no indicator on the end of a read to denote
    // it's pairing. use strcpy and then add pairing information
    char *ret_str = calloc(idLen+4, sizeof(char));
    strcpy(ret_str, longId);
    ++idLen;
    ret_str[idLen] = '_';
    ++idLen;
    ret_str[idLen] = target;
    ++idLen;
    ret_str[idLen] = '\0';
    return ret_str;
}
