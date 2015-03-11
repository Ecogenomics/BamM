//#############################################################################
//
//   bamParser.c
//
//   Engine for parsing BAM files - non-parallel
//   Functions include:
//   BAM parsing to produce pileup information
//   Determine "type" of bam file (insert size, orientation of reads etc...)
//   Find linking read pairs
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
#include "sam.h"

// cfuhash
#include "cfuhash.h"

// local includes
#include "bamParser.h"
#include "pairedLink.h"
#include "coverageEstimators.h"
#include "stats.h"

BM_fileInfo * createBFI(void) {
    return (BM_fileInfo*) calloc(1, sizeof(BM_fileInfo));
}

void initBFI(BM_fileInfo * BFI,
             bam_hdr_t * BAM_header,
             int numBams,
             char * bamFiles[],
             int * types,
             int isLinks,
             BM_coverageType * coverageType,
             int ignoreSuppAlignments
            ) {
    int i = 0; int j = 0;
    BFI->numContigs = BAM_header->n_targets;
    BFI->numBams = numBams;
    BFI->isLinks = isLinks;
    BFI->coverageType = coverageType;
    BFI->isIgnoreSupps = ignoreSuppAlignments;
    BFI->links = 0;

    if(BFI->numContigs != 0 && BFI->numBams != 0) {
        if((BFI->coverageType)->type != CT_NONE) {
            // unless this is an extraction run we will
            // need to make room to store read counts
            BFI->coverages = calloc(BFI->numContigs, sizeof(float*));
            for(i = 0; i < BFI->numContigs; ++i) {
                BFI->coverages[i] = calloc(BFI->numBams, sizeof(float));
            }
        }

        // make BAM file structs
        BFI->bamFiles = (BM_bamFile**) calloc(BFI->numBams,
                                              sizeof(BM_bamFile*));
        for (i = 0; i < numBams; ++i) {
            BM_bamFile * BF = (BM_bamFile*) calloc(1, sizeof(BM_bamFile));
            BF->fileName = strdup(bamFiles[i]);
            BF->fileNameLength = strlen(bamFiles[i]);

            // support for mixed insert sizes in one bam (i.e. shadow library)
            BF->numTypes = types[i];
            BF->types = (BM_bamType**) calloc(BF->numTypes,
                                              sizeof(BM_bamType*));
            for(j = 0; j < BF->numTypes; ++j) {
                BM_bamType * BT = \
                    (BM_bamType*) calloc(1, sizeof(BM_bamType));
                BT->orientationType = OT_NONE;
                BT->insertSize = 0;
                BT->insertStdev = 0;
                BT->supporting = 0;
                BF->types[j] = BT;
            }
            BFI->bamFiles[i] = BF;
        }

        // place for contig storage
        BFI->contigNames = calloc(BFI->numContigs, sizeof(char*));
        BFI->contigNameLengths = calloc(BFI->numContigs, sizeof(uint16_t*));
        BFI->contigLengths = calloc(BFI->numContigs, sizeof(uint32_t));
        for(i =0; i < BFI->numContigs; ++i) {
            BFI->contigNames[i] = strdup(BAM_header->target_name[i]);
            BFI->contigNameLengths[i] = strlen(BAM_header->target_name[i]);
            BFI->contigLengths[i] = (uint32_t)BAM_header->target_len[i];
        }

        if(BFI->isLinks) {
            cfuhash_table_t *links = cfuhash_new_with_initial_size(30);
            cfuhash_set_flag(links, CFUHASH_FROZEN_UNTIL_GROWS);
            BFI->links = links;
        }
    }
}

void mergeBFIs(BM_fileInfo * BFI_A, BM_fileInfo * BFI_B) {
    // Check to see that the number of contigs are the same
    int i = 0, j = 0, k = 0;
    if(BFI_A->numContigs != BFI_B->numContigs) {
        char str[80];
        sprintf(str,
                "Unequal number of contigs in BFI structs to be merged " \
                "(Num contigs: %d, %d)",
                BFI_A->numContigs, BFI_B->numContigs);
        printError(str, __LINE__);
        return;
    }

    // now check that the lengths of all the contigs are the same
    // just check every 100th contig
    for(i=0; i < BFI_A->numContigs; ++i) {
        if(i < BFI_A->numContigs) {
            if(BFI_A->contigLengths[i] != BFI_B->contigLengths[i]) {
                char str[80];
                sprintf(str,
                        "Unequal contig lengths in BFI structs to be merged " \
                        "(Contig#: %d, Lengths: %d, %d)",
                        i,
                        BFI_A->numContigs, BFI_B->numContigs);
                printError(str, __LINE__);
                return;
            }
        }
        else {
            break;
        }
    }

    // we can assume that the headers are the same. Merge the data
    // keep a backup of these guys
    uint32_t old_numBams = BFI_A->numBams;
    BM_bamFile ** old_bamFiles = BFI_A->bamFiles;
    float ** old_coverages = BFI_A->coverages;

    //-----
    // BAM files
    BFI_A->numBams += BFI_B->numBams;
    // realloc the memory for BFI_A
    BFI_A->bamFiles = (BM_bamFile**) calloc(BFI_A->numBams,
                                            sizeof(BM_bamFile*));
    for (i = 0; i < old_numBams; ++i) {
        BFI_A->bamFiles[i] = old_bamFiles[i];
    }
    for (j = 0; j < BFI_B->numBams; ++j) {
        BFI_A->bamFiles[i] = BFI_B->bamFiles[j];
        BFI_B->bamFiles[j] = 0; // avoid deleting twice
        ++i;
    }
    // free these original data
    if(old_bamFiles != 0) {
        free(old_bamFiles);
        old_bamFiles = 0;
    }

    //-----
    // Pileups
    BFI_A->coverages = calloc(BFI_A->numContigs, sizeof(float*));
    for(i = 0; i < BFI_A->numContigs; ++i) {
        BFI_A->coverages[i] = calloc(BFI_A->numBams, sizeof(float));
    }
    for(k = 0; k < BFI_A->numContigs; ++k)
    {
        for (i = 0; i < old_numBams; ++i) {
            //printf("A: %d, %d, %d\n", k, i, old_coverages[k][i]);
            BFI_A->coverages[k][i] = old_coverages[k][i];
        }
        for (j = 0; j < BFI_B->numBams; ++j) {
            //printf("B: %d, %d, %d\n", k, i, BFI_B->coverages[k][j]);
            BFI_A->coverages[k][i] = BFI_B->coverages[k][j];
            ++i;
        }
    }
    if(old_coverages != 0) {
        for(i = 0; i < BFI_A->numContigs; ++i) {
            if(old_coverages[i] != 0) {
                free(old_coverages[i]);
                old_coverages[i] = 0;
            }
        }
        free(old_coverages);
        old_coverages = 0;
    }

    //-----
    // Links
    if(BFI_A->isLinks) {
        // break the pattern of overwriting. We can just add the links here
        char **keys = NULL;
        size_t *key_sizes = NULL;
        size_t key_count = 0;
        keys = (char **)cfuhash_keys_data(BFI_B->links,
                                          &key_count,
                                          &key_sizes,
                                          0);
        for (i = 0; i < (int)key_count; i++) {
            BM_linkPair * LP = cfuhash_get(BFI_B->links, keys[i]);
            if (keys[i] != 0) {
                free(keys[i]);
                keys[i] = 0;
            }
            BM_linkInfo* LI = LP->LI;
            do {
                BM_linkInfo * new_LI = cloneLinkInfo(LI);
                // adjust the bid
                new_LI->bid += old_numBams;
                addLink(BFI_A->links,
                        new_LI,
                        LP->cid1,
                        LP->cid2);
            } while(getNextLinkInfo(&LI));
        }
        if (keys != 0) {
            free(keys);
            keys = 0;
        }
        if (key_sizes != 0) {
            free(key_sizes);
            key_sizes = 0;
        }
    }
    destroyBFI(BFI_B);
}

int read_bam(void *data,
             bam1_t *b) // read level filters better go here to avoid pileup
{
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an aux structure
    int ret = \
        aux->iter? \
            hts_itr_next(aux->fp, aux->iter, b, 0) : \
            bam_read1(aux->fp, b);
    if (!(b->core.flag&BAM_FUNMAP)) {
        if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
        else if (aux->min_len &&
                 ( bam_cigar2qlen((&b->core)->n_cigar, bam_get_cigar(b)) <
                 aux->min_len )) {
            b->core.flag |= BAM_FUNMAP;
        }
    }
    return ret;
}

int parseCoverageAndLinks(int doLinks,
                          int doCovs,
                          int numBams,
                          int baseQ,
                          int mapQ,
                          int minLen,
                          int maxMisMatches,
                          int * types,
                          int ignoreSuppAlignments,
                          int ignoreSecondaryAlignments,
                          BM_coverageType * coverageType,
                          char* bamFiles[],
                          BM_fileInfo * BFI
) {

    int supp_check = 0x0; // include supp mappings
    if (ignoreSuppAlignments) {
        supp_check |= BAM_FSUPPLEMENTARY;
    }
    if (ignoreSecondaryAlignments) {
        supp_check |= BAM_FSECONDARY;
    }

    // initialize the auxiliary data structures
    const bam_pileup1_t **plp;
    bam_hdr_t *h = 0; // BAM header of the 1st input
    aux_t **data;
    bam_mplp_t mplp;
    int tid = 0, *n_plp, b = 0;
    // load contig names and BAM index.
    data = (aux_t**) calloc(numBams, sizeof(aux_t*)); // data[b], the i-th input
    int beg = 0, end = 1<<30;  // set the default region

    for (b = 0; b < numBams; ++b) {
        bam_hdr_t *htmp;
        data[b] = (aux_t*) calloc(1, sizeof(aux_t));
        data[b]->fp = bgzf_open(bamFiles[b], "r"); // open BAM
        data[b]->min_mapQ = mapQ;                  // set the mapQ filter
        data[b]->min_len  = minLen;                // set the qlen filter
        // ( I think this must be done for each file for legitness!)
        htmp = bam_hdr_read(data[b]->fp);          // read the BAM header
        if (b == 0) {
            h = htmp; // keep the header of the 1st BAM
        } else { bam_hdr_destroy(htmp); } // if not the 1st BAM, kill the header
    }

    // initialise the BFI struct
    initBFI(BFI,
            h,
            numBams,
            bamFiles,
            types,
            doLinks,
            coverageType,
            ignoreSuppAlignments
           );

    // always type the bam files
    typeBamFiles(BFI);

    // if we're only here to type the file then we're done
    if(0 == (doCovs || doLinks)) {
        bam_hdr_destroy(h);

        for (b = 0; b < numBams; ++b) {
            bgzf_close(data[b]->fp);
            if (data[b]->iter != 0) {
                bam_itr_destroy(data[b]->iter);
            }
            if (data[b] != 0) {
                free(data[b]);
                data[b] = 0;
            }
        }
        if (data != 0) {
            free(data);
            data = 0;
        }
        return 0;
    }

    // the core multi-pileup loop
    mplp = bam_mplp_init(numBams, read_bam, (void**)data); // initialization
    // n_plp[b] is the number of covering reads from the i-th BAM
    n_plp = calloc(numBams, sizeof(int));
    // plp[b] points to the array of covering reads (internal in mplp)
    plp = calloc(numBams, sizeof(void*));

    // initialise
    int prev_tid = -1;  // the id of the previous positions tid
    int r, rejects = 0;
    int pos = 0; // current position in the contig ( 1 indexed )

    // hold the pileup count at each position in the contig
    uint32_t ** pileup_values = 0;
    // container to hold return values of coverage estimators
    float * coverage_values = 0;
    if(doCovs) {
        pileup_values = calloc(numBams, sizeof(uint32_t*));
        coverage_values = calloc(numBams, sizeof(float));
    }

    // go through each of the contigs in the file, from tid == 0 --> end
    while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) {
        // come to the next covered position
        if (pos < beg || pos >= end) continue; // out of range; skip
        if(tid != prev_tid) {  // we've arrived at a new contig
            if(doCovs) {       // only do coverage if forced
                if(prev_tid != -1) {
                    // at the end of a contig
                    estimateCoverages(coverage_values,
                                      pileup_values,
                                      BFI->coverageType,
                                      BFI->contigLengths[prev_tid],
                                      BFI->numBams);

                    for (b = 0; b < numBams; ++b) {
                        // load the coverages into the BFI
                        BFI->coverages[prev_tid][b] = coverage_values[b];
                        if(pileup_values[b] != 0) {
                            free(pileup_values[b]); // free this up
                            pileup_values[b] = 0;
                        }
                    }
                }
                for (b = 0; b < numBams; ++b) {
                    // reset for next contig
                    // this needs to be reset because each contig is a
                    // different length
                    pileup_values[b] = calloc(BFI->contigLengths[tid],
                                                sizeof(uint32_t));
                }
            }
            prev_tid = tid;
        }
        for (b = 0; b < numBams; ++b) {
            rejects = 0;

            // for each read in the pileup
            for (r = 0; r < n_plp[b]; ++r) {
                // DON'T modfity plp[][] unless you really know
                const bam_pileup1_t *p = plp[b] + r;

                // always skip on too many mismatches
                if (bam_aux2i(bam_aux_get(&(p->b[0]), "NM")) > \
                    maxMisMatches) {
                         ++rejects;
                         continue;
                }

                // see if this is the start of a read
                bam1_core_t core = p->b[0].core;
                if (core.pos == pos) {
                    if(doLinks) {
                        // read is first in pair (avoid dupe links)
                        if (0 == (core.flag & BAM_FREAD1)) {}
                        // read is a paired read
                        else if (0 == (core.flag & BAM_FPAIRED)) {}
                        // both ends are mapped
                        else if (0 != (core.flag & BM_BAM_FMAPPED)) {}
                        // is primary mapping (optional)
                        else if (0 != (core.flag & supp_check)) {}
                        // hits different contigs
                        else if(core.tid == core.mtid) {}
                        else {
                            // looks legit
                            // get a link info
                            BM_linkInfo * LI = \
                                makeLinkInfo(core.tid,
                                             core.mtid,
                                             core.pos,
                                             core.mpos,
                                             ((core.flag&BAM_FREVERSE) != 0),
                                             ((core.flag&BAM_FMREVERSE) != 0),
                                             b
                                            );

                            // add the link
                            addLink(BFI->links,
                                    LI,
                                    core.tid,
                                    core.mtid);
                        }
                    }
                    // for read count coverages we just need to add the starts
                    // usurp pileup counts to mean read starts
                    if(coverageType->type == CT_COUNT || \
                       coverageType->type == CT_C_MEAN) {
                        ++pileup_values[b][pos];
                    }
                }

                // for pilup coverages we need to filter based
                // on base-specific metrics
                if ( (p->is_del || p->is_refskip) || \
                     // having dels or refskips at tid:pos
                     (bam_get_qual(p->b)[p->qpos] < baseQ) )
                     // low base quality
                    {++rejects;}
            }

            // add this position's depth
            if(doCovs) {
                if(coverageType->type != CT_COUNT && \
                   coverageType->type != CT_C_MEAN) {
                    pileup_values[b][pos] = n_plp[b] - rejects;
                }
            }
        }
    }

    if(doCovs && prev_tid != -1) {
        // at the end of a contig
        estimateCoverages(coverage_values,
                          pileup_values,
                          BFI->coverageType,
                          BFI->contigLengths[prev_tid],
                          BFI->numBams);
        for (b = 0; b < numBams; ++b) {
            // load the coverages into the BFI
            BFI->coverages[prev_tid][b] = coverage_values[b];
            if(pileup_values[b] != 0) {
                free(pileup_values[b]);
                pileup_values[b] = 0;
            }
        }
    }

    // all done, time to clean up and leave
    if (pileup_values != 0) {
        free(pileup_values);
        pileup_values = 0;
    }
    if (coverage_values != 0) {
        free(coverage_values);
        coverage_values = 0;
    }
    if (n_plp != 0) {
        free(n_plp);
        n_plp = 0;
    }
    if (plp != 0) {
        free(plp);
        plp = 0;
    }
    bam_mplp_destroy(mplp);
    bam_hdr_destroy(h);

    for (b = 0; b < numBams; ++b) {
        bgzf_close(data[b]->fp);
        if (data[b]->iter != 0) bam_itr_destroy(data[b]->iter);
        if (data[b] != 0) {
            free(data[b]);
            data[b] = 0;
        }
    }
    if (data != 0) {
        free(data);
        data = 0;
    }

    /*
    int c = 0;
    for (b = 0; b < numBams; ++b) {
        for (; c < BFI->numContigs; ++c) {
            fprintf(stdout, "%0.2f, ", BFI->coverages[c][b]);
        }
        fprintf(stdout, "\n--------------------------\n");
    }
    */
    return 0;
}

void typeBamFiles(BM_fileInfo * BFI) {
    //-----
    // code uses the pattern outlined in samtools view (sam_view.c)
    // thanks lh3!
    //
    int i = 0, j = 0;

    // always ignoreSuppAlignments
    int supp_check = BAM_FSUPPLEMENTARY | BAM_FSECONDARY;

    // the number of good links found for each bam
    uint32_t * num_found = calloc(BFI->numBams, sizeof(uint32_t));
    // observed insert sizes per bam per OT
    uint32_t *** inserts = calloc(BFI->numBams, sizeof(uint32_t**));
    // hold counts for each of the orientation types we'll see
    int ** orient_counts = calloc(BFI->numBams, sizeof(int*));

    for (i = 0; i < BFI->numBams; ++i) {
        num_found[i] = 0;
        orient_counts[i] = calloc(3, sizeof(uint32_t*));
        inserts[i] = calloc(3, sizeof(int*));
        for (j = 0; j < 3; ++j) {
            inserts[i][j] = calloc(BM_PAIRS_FOR_TYPE, sizeof(uint32_t));
            orient_counts[i][j] = 0;
        }
    }

    samFile *in = 0; j = 0;
    bam_hdr_t *header = NULL;

    for (i = 0; i < BFI->numBams; ++i) {
        // open file handlers
        if ((in = sam_open((BFI->bamFiles[i])->fileName, "r")) == 0) {
            if (in) sam_close(in);
            if ( header ) bam_hdr_destroy(header);
            continue;
        }
        if ((header = sam_hdr_read(in)) == 0) {
            if (in) sam_close(in);
            if ( header ) bam_hdr_destroy(header);
            continue;
        }

        // retrieve alignments in specified regions
        bam1_t *b;
        // load index
        hts_idx_t *idx = sam_index_load(in, (BFI->bamFiles[i])->fileName);
        if (idx == 0) { // index is unavailable
            fprintf(stderr, "ERROR: Random alignment retrieval only works for "\
                            "indexed BAM or CRAM files.\n");
            if (in) sam_close(in);
            if ( header ) bam_hdr_destroy(header);
            continue;
        }
        b = bam_init1();
        int result = -1;

        int hh;
        for (hh = 0; hh < header->n_targets; ++hh) {
            // parse a region in the format like `chr2:100-200'
            hts_itr_t *iter = sam_itr_querys(idx,
                                             header,
                                             header->target_name[hh]);
            if (iter == NULL) { // reference name is not found
                hts_itr_destroy(iter);
                continue;
            }
            int contig_length = header->target_len[hh];
            if(contig_length < (3 * BM_IGNORE_FROM_END))
            {
                hts_itr_destroy(iter);
                continue;
            }
            int upper_bound = contig_length - BM_IGNORE_FROM_END;
            // fetch alignments
            while ((num_found[i] < BM_PAIRS_FOR_TYPE) &&
                   ((result = sam_itr_next(in, iter, b)) >= 0)) {
                bam1_core_t core = b->core;
                // check to see if this is a proper linking paired read
                if ((core.flag & BAM_FPAIRED) &&        // read is a paired read
                    (core.flag & BAM_FREAD1) &&         // read is first in pair
                    ((core.flag & BM_BAM_FMAPPED) == 0) && // both ends mapped
                    ((core.flag & supp_check) == 0) && // is primary mapping
                    core.tid == core.mtid) {           // hits same contig

                    // make sure we're not too close to the end of the contig
                    int isize = abs(core.mpos - core.pos) + core.l_qseq;
                    int cpos = 0;
                    if(core.mpos < core.pos)
                        cpos = (int)(((float)isize)/2 + core.mpos);
                    else
                        cpos = (int)(((float)isize)/2 + core.pos);

                    if((cpos < BM_IGNORE_FROM_END) || (cpos > upper_bound))
                        continue;

                    // get orientation type
                    int ot = OT_NONE;
                    if((core.flag&BAM_FREVERSE) != 0) {      // <-1--
                        if((core.flag&BAM_FMREVERSE) != 0)
                            ot = OT_SAME;                       // <-1-- <-2--
                        else
                            ot = OT_OUT;                        // <-1-- --2->
                    }
                    else {                                   // --1->
                        if((core.flag&BAM_FMREVERSE) != 0)
                            ot = OT_IN;                         // --1-> <-2--
                        else
                            ot = OT_SAME;                       // --1-> --2->
                    }

                    // check to see if it's about face
                    if(core.isize < 0)
                    {
                        if(ot == OT_IN)
                            ot = OT_OUT;
                        else if(ot == OT_OUT)
                            ot = OT_IN;
                    }

                    // increment thse guys
                    inserts[i][ot][orient_counts[i][ot]] += isize;
                    ++orient_counts[i][ot];
                    ++num_found[i];
                }
            }

            hts_itr_destroy(iter);
            if (result < -1) {
                fprintf(stderr, "ERROR: retrieval of reads from: \"%s\" "\
                                "failed due to truncated file or corrupt BAM "\
                                "index file\n", header->target_name[hh]);
                break;
            }

            // have we found enough?
            if(num_found[i] > BM_PAIRS_FOR_TYPE)
                break;
            // else next contig!
        }
        if (in) sam_close(in);
        if ( header ) bam_hdr_destroy(header);
        bam_destroy1(b);
        hts_idx_destroy(idx); // destroy the BAM index
    }

    for (i = 0; i < BFI->numBams; ++i) {
        int num_to_find = BFI->bamFiles[i]->numTypes;
        while(num_to_find > 0) {
            int type = 0;
            int max = orient_counts[i][type];
            for(j = 1; j < 3; ++j) {
                if(orient_counts[i][j] > max) {
                    type = j;
                    max = orient_counts[i][j];
                }
            }
            int type_index = BFI->bamFiles[i]->numTypes - num_to_find;
            ((BFI->bamFiles[i])->types[type_index])->orientationType = type;
            ((BFI->bamFiles[i])->types[type_index])->insertSize =
                BM_mean(inserts[i][type], max);
            ((BFI->bamFiles[i])->types[type_index])->insertStdev =
                BM_stdDev(inserts[i][type],
                          max,
                          ((BFI->bamFiles[i])->types[type_index])->insertSize);
            ((BFI->bamFiles[i])->types[type_index])->supporting = max;
            // reset this so that we find the next best one
            orient_counts[i][type] = 0;
            --num_to_find;
        }
    }

    // clean up
    if (num_found != 0) {
        free(num_found);
        num_found = 0;
    }
    for (i = 0; i < BFI->numBams; ++i) {
        if (orient_counts[i] != 0) {
            free(orient_counts[i]);
            orient_counts[i] = 0;
        }
        for (j = 0; j < 3; ++j) {
            if (inserts[i][j] != 0) {
                free(inserts[i][j]);
                inserts[i][j] = 0;
            }
        }
        if (inserts[i] != 0) {
            free(inserts[i]);
            inserts[i] = 0;
        }
    }
    if (inserts != 0) {
        free(inserts);
        inserts = 0;
    }
    if (orient_counts != 0) {
        free(orient_counts);
        orient_counts = 0;
    }
}

int initLW(BM_LinkWalker * walker, BM_fileInfo * BFI) {
    return initLinkWalker(walker, BFI->links);
}

int stepLW(BM_LinkWalker * walker) {
    return stepLinkWalker(walker);
}

void destroyBFI(BM_fileInfo * BFI) {
    int i = 0;
    if(BFI != 0)
    {
        if(BFI->numContigs != 0 && BFI->numBams != 0) {
            if(BFI->coverages != 0) {
                for(i = 0; i < BFI->numContigs; ++i) {
                    if(BFI->coverages[i] != 0) {
                        free(BFI->coverages[i]);
                        BFI->coverages[i] = 0;
                    }
                }
                free(BFI->coverages);
                BFI->coverages = 0;
            }

            if(BFI->bamFiles != 0) {
                destroyBamFiles(BFI->bamFiles, BFI->numBams);
            }

            if(BFI->contigNames != 0) {
                for(i = 0; i < BFI->numContigs; ++i) {
                    if(BFI->contigNames[i] != 0) {
                        free(BFI->contigNames[i]);
                        BFI->contigNames[i] = 0;
                    }
                }
                free(BFI->contigNames);
                BFI->contigNames = 0;
                if(BFI->contigNameLengths != 0) {
                    free(BFI->contigNameLengths);
                    BFI->contigNameLengths = 0;
                }
            }

            if(BFI->contigLengths != 0)
                free(BFI->contigLengths);
                BFI->contigLengths = 0;
        }

        // destroy paired links
        if(BFI->isLinks) {
            destroyLinks(BFI->links);
            cfuhash_clear(BFI->links);
            cfuhash_destroy(BFI->links);
        }
        // these guys are alloc'd by ctypes, so we shouldn't free them
        // free(BFI);
        // BFI = 0;
    }
}

void destroyLW(BM_LinkWalker * walker) {
    destroyLinkWalker(walker);
}

void destroyBamFiles(BM_bamFile ** BFs, int numBams) {
    int i = 0;
    int j = 0;
    if(BFs != 0) {
        for(i = 0; i < numBams; ++i) {
            if(BFs[i] != 0) {
                if (BFs[i]->fileName != 0) {
                    free(BFs[i]->fileName);
                    BFs[i]->fileName = 0;
                }
                for(j = 0; j < BFs[i]->numTypes; ++j) {
                    if (BFs[i]->types[j] != 0) {
                        free(BFs[i]->types[j]);
                        BFs[i]->types[j] = 0;
                    }
                }
                if (BFs[i]->types != 0) {
                    free(BFs[i]->types);
                    BFs[i]->types = 0;
                }
                free(BFs[i]);
                BFs[i] = 0;
            }
        }
        free(BFs);
        BFs = 0;
    }
}

char * OT2Str(OT type) {
    switch(type) {
        case OT_OUT:
            return "OUT";
            break;
        case OT_SAME:
            return "SAME";
            break;
        case OT_IN:
            return "IN";
            break;
        case OT_NONE:
            return "NONE";
            break;
        case OT_ERROR:
            return "ERROR";
            break;
    }
    return "UNKNOWN";
}

void printError(char* errorMessage, int line) {
    //-----
    // Consistent way to print error messages
    //
    printf("ERROR: At line: %d\n\t%s\n\n", line, errorMessage);
}

void printBFI(BM_fileInfo * BFI) {
    int i = 0, j = 0;
    if(BFI->numContigs != 0 && BFI->numBams != 0) {
        if(BFI->coverages != 0) {
            // get the coverages we want
            // print away!
            printf("#contig\tlength");
            for(j = 0; j < BFI->numBams; ++j) {
                printf("\t%s",(BFI->bamFiles[j])->fileName);
            }
            printf("\n");
            for(i = 0; i < BFI->numContigs; ++i) {
                printf("%s\t%d",
                       BFI->contigNames[i],
                       BFI->contigLengths[i]);
                for(j = 0; j < BFI->numBams; ++j) {
                    printf("\t%0.4f", BFI->coverages[i][j]);
                }
                printf("\n");
            }
            printf("---\n");
            // we're responsible for cleaning up the covs structure
            if(BFI->isLinks) {
                char** bfns = calloc(BFI->numBams, sizeof(char*));
                printf("#file\tinsert\tstdev\torientation\tsupporting\n");
                for(j = 0; j < BFI->numBams; ++j) {
                    bfns[j] = (BFI->bamFiles[j])->fileName;
                    printBamFileType(BFI->bamFiles[j]);
                }
                printLinks(BFI->links, bfns, BFI->contigNames);
                free(bfns);
            }
        }
    }
}

void printBamFileType(BM_bamFile * BF)
{
    int i = 0;
    for(i = 0; i < BF->numTypes; ++i)
    {
        BM_bamType * BT = BF->types[i];
        printf("%s\t%0.4f\t%0.4f\t%s\t%d\n", BF->fileName,
                                             BT->insertSize,
                                             BT->insertStdev,
                                             OT2Str(BT->orientationType),
                                             BT->supporting);
    }
    printf("---\n");
}
