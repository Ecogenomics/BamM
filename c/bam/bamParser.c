//#############################################################################
//
//   bamParser.c
//
//   Determine average coverage values and linking read pairs
//
//   Copyright (C) Michael Imelfort
//
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//#############################################################################

// system includes
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

// htslib
//#include "htslib/bgzf.h"
//#include "htslib/sam.h"

// cfuhash
#include "cfuhash.h"

// local includes
#include "bamParser.h"
#include "pairedLink.h"
#include "stats.h"

// proper linking read is a properly paired, (primary alignment) of the first read in thr pair
#define BM_BAM_FSUPP (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)
#define BM_BAM_FMAPPED (BAM_FMUNMAP | BAM_FUNMAP)

BM_mappingResults * create_MR(void)
{
    //-----
    // allocate space for a new MR struct
    //
    return calloc(1, sizeof(BM_mappingResults));
}

void init_MR(BM_mappingResults * MR,
             bam_hdr_t * BAM_header,
             int numBams,
             char * bamFiles[],
             int * links,
             char * coverageMode,
             int ignoreSuppAlignments
            ) {
    //----
    // initialise the MR object
    //
    // check for a valid coverage mode
    int valid_coverage_mode = strcmp(coverageMode, "vanilla");
    if(valid_coverage_mode != 0)
    { valid_coverage_mode = strcmp(coverageMode, "outlier"); }
    if(valid_coverage_mode != 0)
    {
        char str[80];
        sprintf(str, "Invalid coverage mode '%s'", coverageMode);
        printError(str, __LINE__);
        return;
    }

    int i = 0;
    int j = 0;
    MR->numContigs = BAM_header->n_targets;
    MR->numBams = numBams;
    MR->isLinks = (links != 0);
    MR->coverage_mode = strdup(coverageMode);
    MR->isIgnoreSupps = ignoreSuppAlignments;

    if(MR->numContigs != 0 && MR->numBams != 0) {
        // make room to store read counts
        MR->plpBp = calloc(MR->numContigs, sizeof(uint32_t*));
        for(i = 0; i < MR->numContigs; ++i) {
            MR->plpBp[i] = calloc(MR->numBams, sizeof(uint32_t));
        }

        // make BAM file structs
        MR->bamFiles = calloc(MR->numBams, sizeof(BM_bamFile*));
        for (i = 0; i < numBams; ++i) {
            BM_bamFile * BF = calloc(1, sizeof(BM_bamFile));
            BF->fileName = strdup(bamFiles[i]);
            BF->fileNameLength = strlen(bamFiles[i]);

            // support for mixed insert sizes in one bam (i.e. shadow library)
            BF->numTypes = links[i];
            BF->types = calloc(BF->numTypes, sizeof(BM_bamType*));
            for(j = 0; j < BF->numTypes; ++j) {
                BM_bamType * BT = calloc(1, sizeof(BM_bamType));
                BT->orientationType = OT_NONE;
                BT->insertSize = 0;
                BT->insertStdev = 0;
                BT->supporting = 0;
                BF->types[j] = BT;
            }
            MR->bamFiles[i] = BF;
        }

        // place for contig storage
        MR->contigNames = calloc(MR->numContigs, sizeof(char*));
        MR->contig_name_lengths = calloc(MR->numContigs, sizeof(uint16_t*));
        MR->contigLengths = calloc(MR->numContigs, sizeof(uint32_t));
        for(i =0; i < MR->numContigs; ++i) {
            MR->contigNames[i] = strdup(BAM_header->target_name[i]);
            MR->contig_name_lengths[i] = strlen(BAM_header->target_name[i]);
            MR->contigLengths[i] = (uint32_t)BAM_header->target_len[i];
        }

        //----------------------------
        // only allocate if we NEED to
        //----------------------------
        if (strcmp(MR->coverage_mode, "outlier") == 0) {
            MR->contigLengthCorrectors = calloc(MR->numContigs, sizeof(uint32_t*));
            for(i = 0; i < MR->numContigs; ++i) {
                MR->contigLengthCorrectors[i] = calloc(MR->numBams, sizeof(uint32_t));
            }
        } else {
            MR->contigLengthCorrectors = NULL;
        }

        if(MR->isLinks)
        {
            cfuhash_table_t *links = cfuhash_new_with_initial_size(30);
            cfuhash_set_flag(links, CFUHASH_FROZEN_UNTIL_GROWS);
            MR->links = links;
        }
        else
        {
            MR->links = 0;
        }
    }
}

void merge_MRs(BM_mappingResults * MR_A, BM_mappingResults * MR_B)
{
    //----
    // Merge the contents of MR_B into MR_A
    //
    // Check to see that the number of contigs are the same
    int i = 0, j = 0, k = 0;
    if(MR_A->numContigs != MR_B->numContigs)
    {
        char str[80];
        sprintf(str, "Unequal number of contigs in MR structs to be merged (Num contigs: %d, %d)", MR_A->numContigs, MR_B->numContigs);
        printError(str, __LINE__);
        return;
    }

    // now check that the lengths of all the contigs are the same
    // just check every 100th contig
    for(i=0; i < MR_A->numContigs; ++i)
    {
        if(i < MR_A->numContigs)
        {
            if(MR_A->contigLengths[i] != MR_B->contigLengths[i])
            {
                char str[80];
                sprintf(str, "Unequal contig lengths in MR structs to be merged (Contig#: %d, Lengths: %d, %d)", i, MR_A->numContigs, MR_B->numContigs);
                printError(str, __LINE__);
                return;
            }
        }
        else
        {
            break;
        }
    }

    // we can assume that the headers are the same. So now time to merge the data

    // keep a backup of these guys
    uint32_t old_numBams = MR_A->numBams;
    BM_bamFile ** old_bamFiles = MR_A->bamFiles;
    uint32_t ** old_contigLengthCorrectors = MR_A->contigLengthCorrectors;
    uint32_t ** old_plpBp = MR_A->plpBp;

    //-----
    // Fix the num bams and bam file names
    MR_A->numBams += MR_B->numBams;
    // realloc the memory for MR_A
    MR_A->bamFiles = calloc(MR_A->numBams, sizeof(BM_bamFile*));
    for (i = 0; i < old_numBams; ++i) {
        MR_A->bamFiles[i] = old_bamFiles[i];
    }
    for (j = 0; j < MR_B->numBams; ++j) {
        MR_A->bamFiles[i] = MR_B->bamFiles[j];
        ++i;
    }
    // free these original data
    if(old_bamFiles != 0) {
        free(old_bamFiles);
    }

    //-----
    // Pileups
    MR_A->plpBp = calloc(MR_A->numContigs, sizeof(uint32_t*));
    for(i = 0; i < MR_A->numContigs; ++i) {
        MR_A->plpBp[i] = calloc(MR_A->numBams, sizeof(uint32_t));
    }
    for(k = 0; k < MR_A->numContigs; ++k)
    {
        for (i = 0; i < old_numBams; ++i) {
            //printf("A: %d, %d, %d\n", k, i, old_plpBp[k][i]);
            MR_A->plpBp[k][i] = old_plpBp[k][i];
        }
        for (j = 0; j < MR_B->numBams; ++j) {
            //printf("B: %d, %d, %d\n", k, i, MR_B->plpBp[k][j]);
            MR_A->plpBp[k][i] = MR_B->plpBp[k][j];
            ++i;
        }
    }
    if(old_plpBp != 0) {
        for(i = 0; i < MR_A->numContigs; ++i) {
            if(old_plpBp[i] != 0)
                free(old_plpBp[i]);
        }
        free(old_plpBp);
    }

    //-----
    // Contig length correctors
    //
    if (strcmp(MR_A->coverage_mode, "outlier") == 0) {
        MR_A->contigLengthCorrectors = calloc(MR_A->numContigs, sizeof(uint32_t*));
        for(i = 0; i < MR_A->numContigs; ++i) {
            MR_A->contigLengthCorrectors[i] = calloc(MR_A->numBams, sizeof(uint32_t));
        }
        for(k = 0; k < MR_A->numContigs; ++k)
        {
            for (i = 0; i < old_numBams; ++i) {
                MR_A->contigLengthCorrectors[k][i] = old_contigLengthCorrectors[k][i];
            }
            for (j = 0; j < MR_B->numBams; ++j) {
                MR_A->contigLengthCorrectors[k][i] = MR_B->contigLengthCorrectors[k][j];
                ++i;
            }
        }
        if(old_contigLengthCorrectors != 0) {
            for(i = 0; i < MR_A->numContigs; ++i) {
                if(old_contigLengthCorrectors[i] != 0)
                    free(old_contigLengthCorrectors[i]);
            }
            free(old_contigLengthCorrectors);
        }
    }


    //-----
    // Links
    if(MR_A->isLinks)
    {
        // break the pattern of overwriting. We can just add the links here
        char **keys = NULL;
        size_t *key_sizes = NULL;
        size_t key_count = 0;
        keys = (char **)cfuhash_keys_data(MR_B->links, &key_count, &key_sizes, 0);

        for (i = 0; i < (int)key_count; i++) {
            BM_linkPair * LP = cfuhash_get(MR_B->links, keys[i]);
            free(keys[i]);
            BM_linkInfo* LI = LP->LI;
            do {
                addLink(MR_A->links,
                        LP->cid1,
                        LP->cid2,
                        LI->pos1,
                        LI->pos2,
                        LI->reversed1,
                        LI->reversed2,
                        LI->bid+old_numBams,
                        LT_NONE);
            } while(getNextLinkInfo(&LI));

        }
        free(keys);
        free(key_sizes);
    }
}

void destroy_MR(BM_mappingResults * MR)
{
    //----
    // free the MR object
    //
    int i = 0;
    if(MR->numContigs != 0 && MR->numBams != 0) {
        if(MR->plpBp != 0) {
            for(i = 0; i < MR->numContigs; ++i) {
                if(MR->plpBp[i] != 0)
                    free(MR->plpBp[i]);
            }
            free(MR->plpBp);
        }

        if(MR->bamFiles != 0) {
            destroyBamFiles(MR->bamFiles, MR->numBams);
        }

        if(MR->contigNames != 0) {
            for(i = 0; i < MR->numContigs; ++i) {
                if(MR->contigNames[i] != 0)
                    free(MR->contigNames[i]);
            }
            free(MR->contigNames);
            free(MR->contig_name_lengths);
        }

        if(MR->contigLengths != 0)
            free(MR->contigLengths);

        if(MR->contigLengthCorrectors != 0) {
            for(i = 0; i < MR->numContigs; ++i) {
                if(MR->contigLengthCorrectors[i] != 0)
                    free(MR->contigLengthCorrectors[i]);
            }
            free(MR->contigLengthCorrectors);
        }
    }

    free(MR->coverage_mode);

    // destroy paired links
    if(MR->isLinks)
    {
        destroyLinks(MR->links);
        cfuhash_clear(MR->links);
        cfuhash_destroy(MR->links);
    }
}

// This function reads a BAM alignment from one BAM file.
int read_bam(void *data,
             bam1_t *b) // read level filters better go here to avoid pileup
{
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret = aux->iter? hts_itr_next(aux->fp, aux->iter, b, 0) : bam_read1(aux->fp, b);
    if (!(b->core.flag&BAM_FUNMAP)) {
        if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
        else if (aux->min_len && bam_cigar2qlen((&b->core)->n_cigar, bam_get_cigar(b)) < aux->min_len) b->core.flag |= BAM_FUNMAP;
    }
    return ret;
}

int parseCoverageAndLinks(int numBams,
                          int baseQ,
                          int mapQ,
                          int minLen,
                          int * links,
                          int ignoreSuppAlignments,
                          char* coverageMode,
                          char* bamFiles[],
                          BM_mappingResults * MR
) {
    //-----
    // work out coverage depths and also pairwise linkages if asked to do so
    //

    int supp_check = 0x0; // include supp mappings
    if (ignoreSuppAlignments) {
        supp_check = BM_BAM_FSUPP;
    }

    // initialize the auxiliary data structures
    const bam_pileup1_t **plp;
    bam_hdr_t *h = 0; // BAM header of the 1st input
    aux_t **data;
    bam_mplp_t mplp;
    int tid = 0, *n_plp, i = 0;
    // load contig names and BAM index.
    data = calloc(numBams, sizeof(aux_t*)); // data[i] for the i-th input
    int beg = 0, end = 1<<30;  // set the default region

    for (i = 0; i < numBams; ++i) {
        bam_hdr_t *htmp;
        data[i] = calloc(1, sizeof(aux_t));
        data[i]->fp = bgzf_open(bamFiles[i], "r"); // open BAM
        data[i]->min_mapQ = mapQ;                  // set the mapQ filter
        data[i]->min_len  = minLen;                // set the qlen filter
        htmp = bam_hdr_read(data[i]->fp);          // read the BAM header ( I think this must be done for each file for legitness!)
        if (i == 0) {
            h = htmp; // keep the header of the 1st BAM
        } else { bam_hdr_destroy(htmp); } // if not the 1st BAM, trash the header
    }

    // initialise the mapping results struct
    init_MR(MR,
            h,
            numBams,
            bamFiles,
            links,
            coverageMode,
            ignoreSuppAlignments
           );

    // type the bam files, but only if we're doing links
    if(links != 0) {
        typeBamFiles(MR);
    }

    // the core multi-pileup loop
    mplp = bam_mplp_init(numBams, read_bam, (void**)data); // initialization
    n_plp = calloc(numBams, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp = calloc(numBams, sizeof(void*)); // plp[i] points to the array of covering reads (internal in mplp)

    // initialise
    int prev_tid = -1;  // the id of the previous positions tid
    int j, rejects = 0;
    int pos = 0; // current position in the contig ( 1 indexed )
    uint32_t ** position_holder; // hold the pileup count at each position in the contig
    position_holder = calloc(numBams, sizeof(uint32_t*));
    // go through each of the contigs in the file, from tid == 0 --> end

    while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
        if (pos < beg || pos >= end) continue; // out of range; skip
        if(tid != prev_tid) {  // we've arrived at a new contig
            if(prev_tid != -1) {
                // at the end of a contig
                adjustPlpBp(MR, position_holder, prev_tid);
                for (i = 0; i < numBams; ++i) {
                    free(position_holder[i]); // free this up
                }
            }
            for (i = 0; i < numBams; ++i) {
                // reset for next contig
                position_holder[i] = calloc(MR->contigLengths[tid], sizeof(uint32_t));
            }
            prev_tid = tid;
        }
        for (i = 0; i < numBams; ++i) {
            rejects = 0;
            // for each read in the pileup
            for (j = 0; j < n_plp[i]; ++j) {
                const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
                if (p->is_del || p->is_refskip) {++rejects;} // having dels or refskips at tid:pos
                else if (bam_get_qual(p->b)[p->qpos] < baseQ) {++rejects;} // low base quality
                else if(MR->isLinks) {
                    // now we do links if we've been asked to
                    bam1_core_t core = p->b[0].core;
                    //printf("-->%d %d %d %d %d %d\n", (core.flag & BAM_FPAIRED), (core.flag & BAM_FREAD1), ((core.flag & BM_BAM_FMAPPED) == 0), ((core.flag & supp_check) == 0), core.tid, core.mtid);
                    // check to see if this is a proper linking paired read
                    if ((core.pos == pos) &&                        // make sure it's the first time we've seen this link
                        (core.flag & BAM_FPAIRED) &&                // read is a paired read
                        (core.flag & BAM_FREAD1) &&                 // read is first in pair (avoid dupe links)
                        ((core.flag & BM_BAM_FMAPPED) == 0) &&     // both ends are mapped
                        ((core.flag & supp_check) == 0) &&          // is primary mapping (optional)
                        core.tid != core.mtid) {                    // hits different contigs

                        // looks legit
                        addLink(MR->links,
                                core.tid,                           // contig 1
                                core.mtid,                          // contig 2
                                core.pos,                           // pos 1
                                core.mpos,                          // pos 2
                                ((core.flag&BAM_FREVERSE) != 0),    // 1 == reversed
                                ((core.flag&BAM_FMREVERSE) != 0),   // 0 = agrees
                                i,                                  // bam file ID
                                LT_NONE                             // the type of the link
                                );
                    }
                }
            }
            position_holder[i][pos] = n_plp[i] - rejects; // add this position's depth
        }
    }

    if(prev_tid != -1) {
        // at the end of a contig
        adjustPlpBp(MR, position_holder, prev_tid);
        for (i = 0; i < numBams; ++i) {
            free(position_holder[i]);
        }
    }

    free(position_holder);

    free(n_plp); free(plp);
    bam_mplp_destroy(mplp);
    bam_hdr_destroy(h);

    for (i = 0; i < numBams; ++i) {
        bgzf_close(data[i]->fp);
        if (data[i]->iter) bam_itr_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data);

    return 0;
}

void typeBamFiles(BM_mappingResults * MR) {

    // get a list of all the bam file names
    char** bam_files = calloc(MR->numBams, sizeof(char*));
    int i = 0, tid = 0, *n_plp, j = 0;
    for(j = 0; j < MR->numBams; ++j) {
        bam_files[j] = (MR->bamFiles[j])->fileName;
    }

    // always ignoreSuppAlignments
    int supp_check = BM_BAM_FSUPP;

    // initialize the auxiliary data structures
    const bam_pileup1_t **plp;
    bam_hdr_t *h = 0; // BAM header of the 1st input
    aux_t **data;
    bam_mplp_t mplp;
    // load contig names and BAM index.
    data = calloc(MR->numBams, sizeof(void*)); // data[i] for the i-th input
    int beg = 0, end = 1<<30;  // set the default region

    uint32_t * num_found = calloc(MR->numBams, sizeof(uint32_t));     // the number of good links found for each bam
    uint32_t *** inserts = calloc(MR->numBams, sizeof(uint32_t**));   // observed insert sizes per bam per OT
    int num_done = 0;                                                 // the number of bams that enough reads
    int ** orient_counts = calloc(MR->numBams, sizeof(int*));         // hold counts for each of the orientation types we'll see

    for (i = 0; i < MR->numBams; ++i) {
        num_found[i] = 0;
        orient_counts[i] = calloc(3, sizeof(uint32_t*));
        inserts[i] = calloc(3, sizeof(int*));
        for (j = 0; j < 3; ++j)
        {
            inserts[i][j] = calloc(BM_PAIRS_FOR_TYPE, sizeof(uint32_t));
            orient_counts[i][j] = 0;
        }

        data[i] = calloc(1, sizeof(aux_t));
        data[i]->fp = bgzf_open(bam_files[i], "r"); // open BAM
        data[i]->min_mapQ = 0;                     // set the mapQ filter
        data[i]->min_len  = 0;                     // set the qlen filter
        bam_hdr_t *htmp;
        htmp = bam_hdr_read(data[i]->fp);          // read the BAM header ( I think this must be done for each file for legitness!)
        if (i == 0) {
            h = htmp; // keep the header of the 1st BAM
        } else { bam_hdr_destroy(htmp); } // if not the 1st BAM, trash the header
    }

    // the core multi-pileup loop
    mplp = bam_mplp_init(MR->numBams, read_bam, (void**)data); // initialization
    n_plp = calloc(MR->numBams, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp = calloc(MR->numBams, sizeof(void*)); // plp[i] points to the array of covering reads (internal in mplp)

    // initialise
    int pos = 0; // current position in the contig ( 1 indexed )

    j = 0;
    // go through each of the contigs in the file, from tid == 0 --> end
    while ((num_done < MR->numBams) && (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0)) { // come to the next covered position
        if (pos < beg || pos >= end) continue; // out of range; skip
        for (i = 0; i < MR->numBams; ++i) {
            if(num_found[i] <= BM_PAIRS_FOR_TYPE)
            {
                // for each read in the pileup
                for (j = 0; j < n_plp[i]; ++j) {
                    const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
                    if (p->is_del || p->is_refskip) {} // having dels or refskips at tid:pos
                    else {
                        // now we do links if we've been asked to
                        bam1_core_t core = p->b[0].core;
                        // check to see if this is a proper linking paired read
                        if ((core.pos == pos) &&                        // make sure this is the first time we've ween this particular pair
                            (core.flag & BAM_FPAIRED) &&                // read is a paired read
                            (core.flag & BAM_FREAD1) &&                 // read is first in pair (avoid dupe links)
                            ((core.flag & BM_BAM_FMAPPED) == 0) &&      // both ends are mapped
                            ((core.flag & supp_check) == 0) &&          // is primary mapping (optional)
                            core.tid == core.mtid) {                    // hits same contig

                            // get orientation type
                            int ot = OT_NONE;
                            if((core.flag&BAM_FREVERSE) != 0) {             // <-1--
                                if((core.flag&BAM_FMREVERSE) != 0)
                                    ot = OT_SAME;                              // <-1-- <-2--
                                else
                                    ot = OT_OUT;                               // <-1-- --2->
                            }
                            else {                                          // --1->
                                if((core.flag&BAM_FMREVERSE) != 0)
                                    ot = OT_IN;                                // --1-> <-2--
                                else
                                    ot = OT_SAME;                              // --1-> --2->
                            }
                            //int isize = abs(core.isize);
                            // check to see if it's about face
                            if(core.isize < 0)
                            {
                                if(ot == OT_IN)
                                    ot = OT_OUT;
                                else if(ot == OT_OUT)
                                    ot = OT_IN;
                            }

                            if(core.isize > 0) {
                                inserts[i][ot][orient_counts[i][ot]] += core.isize;
                            } else {
                                inserts[i][ot][orient_counts[i][ot]] -= core.isize;
                            }
                            // increment thse guys
                            ++orient_counts[i][ot];
                            ++num_found[i];
                            /*
                            printf("TYPE: I: %d, CID: %d, P: %d, MP: %d, R: %d, MR: %d, BAM: %d, OT: %d, FOUND: %d\n",
                                    isize,
                                    core.mtid,
                                    core.pos,
                                    core.mpos,
                                    ((core.flag&BAM_FREVERSE) != 0),
                                    ((core.flag&BAM_FMREVERSE) != 0),
                                    i,
                                    ot,
                                    num_found[i]);
                            */
                            if(num_found[i] ==  BM_PAIRS_FOR_TYPE)
                            {
                                // this bam is done
                                ++num_done;
                            }
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i < MR->numBams; ++i) {
        int num_to_find = MR->bamFiles[i]->numTypes;
        while(num_to_find > 0) {
            int type = 0;
            int max = orient_counts[i][type];
            for(j = 1; j < 3; ++j) {
                if(orient_counts[i][j] > max) {
                    type = j;
                    max = orient_counts[i][j];
                }
            }
            int type_index = MR->bamFiles[i]->numTypes - num_to_find;
            ((MR->bamFiles[i])->types[type_index])->orientationType = type;
            ((MR->bamFiles[i])->types[type_index])->insertSize = BM_mean(inserts[i][type], max);
            ((MR->bamFiles[i])->types[type_index])->insertStdev = BM_stdDev(inserts[i][type], max, ((MR->bamFiles[i])->types[type_index])->insertSize);
            ((MR->bamFiles[i])->types[type_index])->supporting = max;
            // reset this so that we find the next best one
            orient_counts[i][type] = 0;
            --num_to_find;
        }
    }

    free(num_found);
    free(n_plp); free(plp);
    bam_mplp_destroy(mplp);
    bam_hdr_destroy(h);

    for (i = 0; i < MR->numBams; ++i) {
        bgzf_close(data[i]->fp);
        if (data[i]->iter) bam_itr_destroy(data[i]->iter);
        free(data[i]);
        free(orient_counts[i]);
        for (j = 0; j < 3; ++j) {
            free(inserts[i][j]);
        }
        free(inserts[i]);

    }
    free(inserts);
    free(data);
    free(orient_counts);
    free(bam_files);
}

void adjustPlpBp(BM_mappingResults * MR,
                 uint32_t ** positionHolder,
                 int tid
) {
    uint32_t num_bams = MR->numBams;
    uint32_t * plp_sum = calloc(num_bams, sizeof(uint32_t));
    int pos = 0, i = 0;
    if(strcmp(MR->coverage_mode, "outlier") == 0) {
        uint32_t * drops = calloc(num_bams, sizeof(uint32_t));
        for(i = 0; i < num_bams; ++i) {
            // set the cut off at a stdev either side of the mean
            float m = BM_mean(positionHolder[i], MR->contigLengths[tid]);
            float std = BM_stdDev(positionHolder[i], MR->contigLengths[tid], m);
            float lower_cut = ((m-std) < 0) ? 0 : (m-std);
            float upper_cut = m+std;
            for(pos = 0; pos < MR->contigLengths[tid]; ++pos) {
                if((positionHolder[i][pos] <= upper_cut) && (positionHolder[i][pos] >= lower_cut)) {
                    // OK
                    plp_sum[i] += positionHolder[i][pos];
                } else {
                    // DROP
                    ++drops[i];
                }
            }
            MR->plpBp[tid][i] = plp_sum[i];
            MR->contigLengthCorrectors[tid][i] = drops[i];
        }
        free(drops);
    } else {
        for(i = 0; i < num_bams; ++i) {
            for(pos = 0; pos < MR->contigLengths[tid]; ++pos) {
                plp_sum[i] += positionHolder[i][pos];
            }
            MR->plpBp[tid][i] = plp_sum[i];
        }
    }
    free(plp_sum);
}

float ** calculateCoverages(BM_mappingResults * MR) {
    int i = 0, j = 0;
    if(MR->numContigs != 0 && MR->numBams != 0) {
        if(MR->plpBp != NULL) {
            float ** ret_matrix = calloc(MR->numContigs, sizeof(float*));
            for(i = 0; i < MR->numContigs; ++i) {
                ret_matrix[i] = calloc(MR->numBams, sizeof(float));
                for(j = 0; j < MR->numBams; ++j) {
                    // print average coverages
                    if(strcmp(MR->coverage_mode, "outlier") == 0) {
                        // the counts are reduced so we should reduce the contig length accordingly
                        ret_matrix[i][j] = (float)MR->plpBp[i][j]/(float)(MR->contigLengths[i]-MR->contigLengthCorrectors[i][j]);
                    } else {
                        // otherwise it's just a straight up average
                        ret_matrix[i][j] = (float)MR->plpBp[i][j]/(float)MR->contigLengths[i];
                    }
                }
            }
            return ret_matrix;
        }
    }
    return NULL;
}

void destroyCoverages(float ** covs, int numContigs) {
    int i = 0;
    for(i = 0; i < numContigs; ++i) {
        free(covs[i]);
    }
    free(covs);
}

    /***********************
    ***      LINKS      ***
    ***********************/

int initLW(BM_LinkWalker * walker, BM_mappingResults * MR) {
    return initLinkWalker(walker, MR->links);
}

int stepLW(BM_LinkWalker * walker) {
    return stepLinkWalker(walker);
}

void destroyLW(BM_LinkWalker * walker) {
    destroyLinkWalker(walker);
}

    /***********************
    ***     BAMFILES     ***
    ***********************/

void destroyBamFiles(BM_bamFile ** BFs, int numBams) {
    int i = 0;
    int j = 0;
    if(BFs != 0) {
        for(i = 0; i < numBams; ++i) {
            free(BFs[i]->fileName);
            for(j = 0; j < BFs[i]->numTypes; ++j) {
                free(BFs[i]->types[j]);
            }
        }
        free(BFs);
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

    /***********************
    *** PRINTING AND I/O ***
    ***********************/

void printError(char* errorMessage, int line) {
    //-----
    // Consistent way to print error messages
    //
    printf("ERROR: At line: %d\n\t%s\n\n", line, errorMessage);
}

void print_MR(BM_mappingResults * MR) {
    int i = 0, j = 0;
    if(MR->numContigs != 0 && MR->numBams != 0) {
        if(MR->plpBp != NULL) {
            float ** covs = calculateCoverages(MR); // get the coverages we want
            if(covs != NULL) {
                // print away!
                printf("#contig\tlength");
                for(j = 0; j < MR->numBams; ++j) {
                    printf("\t%s",(MR->bamFiles[j])->fileName);
                }
                printf("\n");
                for(i = 0; i < MR->numContigs; ++i) {
                    printf("%s\t%d", MR->contigNames[i], MR->contigLengths[i]);
                    for(j = 0; j < MR->numBams; ++j) {
                        printf("\t%0.4f", covs[i][j]);
                    }
                    printf("\n");
                }
                printf("---\n");
                // we're responsible for cleaning up the covs structure
                destroyCoverages(covs, MR->numContigs);
                if(MR->isLinks) {
                    char** bfns = calloc(MR->numBams, sizeof(char*));
                    printf("#file\tinsert\tstdev\torientation\tsupporting\n");
                    for(j = 0; j < MR->numBams; ++j) {
                        bfns[j] = (MR->bamFiles[j])->fileName;
                        printBamFileInfo(MR->bamFiles[j]);
                    }
                    printLinks(MR->links, bfns, MR->contigNames);
                    free(bfns);
                }
            }
        }
    }
}

void printBamFileInfo(BM_bamFile * BF)
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
