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
#define BMM_BAM_FSUPP (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)
#define BMM_BAM_FMAPPED (BAM_FMUNMAP | BAM_FUNMAP)

BMM_mapping_results * create_MR(void)
{
    //-----
    // allocate space for a new MR struct
    //
    return calloc(1, sizeof(BMM_mapping_results));
}

void init_MR(BMM_mapping_results * MR,
             bam_hdr_t * BAM_header,
             int numBams,
             char * bamFiles[],
             int doLinks,
             char * coverageMode,
             int ignoreSuppAlignments
)
{
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
    MR->num_contigs = BAM_header->n_targets;
    MR->num_bams = numBams;
    MR->is_links_included = doLinks;
    MR->coverage_mode = strdup(coverageMode);
    MR->is_ignore_supps = ignoreSuppAlignments;

    if(MR->num_contigs != 0 && MR->num_bams != 0) {
        // make room to store read counts
        MR->plp_bp = calloc(MR->num_contigs, sizeof(uint32_t*));
        for(i = 0; i < MR->num_contigs; ++i) {
            MR->plp_bp[i] = calloc(MR->num_bams, sizeof(uint32_t));
        }

        // store BAM file names
        MR->bam_file_names = calloc(MR->num_bams, sizeof(char*));
        MR->bam_file_name_lengths = calloc(MR->num_bams, sizeof(uint16_t*));
        for (i = 0; i < numBams; ++i) {
            MR->bam_file_names[i] = strdup(bamFiles[i]);
            MR->bam_file_name_lengths[i] = strlen(bamFiles[i]);
        }

        // place for contig storage
        MR->contig_names = calloc(MR->num_contigs, sizeof(char*));
        MR->contig_name_lengths = calloc(MR->num_contigs, sizeof(uint16_t*));
        MR->contig_lengths = calloc(MR->num_contigs, sizeof(uint32_t));
        for(i =0; i < MR->num_contigs; ++i) {
            MR->contig_names[i] = strdup(BAM_header->target_name[i]);
            MR->contig_name_lengths[i] = strlen(BAM_header->target_name[i]);
            MR->contig_lengths[i] = (uint32_t)BAM_header->target_len[i];
        }

        //----------------------------
        // only allocate if we NEED to
        //----------------------------
        if (strcmp(MR->coverage_mode, "outlier") == 0) {
            MR->contig_length_correctors = calloc(MR->num_contigs, sizeof(uint32_t*));
            for(i = 0; i < MR->num_contigs; ++i) {
                MR->contig_length_correctors[i] = calloc(MR->num_bams, sizeof(uint32_t));
            }
        } else {
            MR->contig_length_correctors = NULL;
        }

        if(doLinks)
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

void merge_MRs(BMM_mapping_results * MR_A, BMM_mapping_results * MR_B)
{
    //----
    // Merge the contents of MR_B into MR_A
    //
    // Check to see that the number of contigs are the same
    int i = 0, j = 0, k = 0;
    if(MR_A->num_contigs != MR_B->num_contigs)
    {
        char str[80];
        sprintf(str, "Unequal number of contigs in MR structs to be merged (Num contigs: %d, %d)", MR_A->num_contigs, MR_B->num_contigs);
        printError(str, __LINE__);
        return;
    }

    // now check that the lengths of all the contigs are the same
    // just check every 100th contig
    for(i=0; i < MR_A->num_contigs; ++i)
    {
        if(i < MR_A->num_contigs)
        {
            if(MR_A->contig_lengths[i] != MR_B->contig_lengths[i])
            {
                char str[80];
                sprintf(str, "Unequal contig lengths in MR structs to be merged (Contig#: %d, Lengths: %d, %d)", i, MR_A->num_contigs, MR_B->num_contigs);
                printError(str, __LINE__);
                return;
            }
        }
        else
        {
            break;
        }
//        i += 100;
    }

    // we can assume that the headers are the same. So now time to merge the data

    // keep a backup of these guys
    uint32_t old_num_bams = MR_A->num_bams;
    char ** old_bam_file_names = MR_A->bam_file_names;
    uint32_t ** old_contig_length_correctors = MR_A->contig_length_correctors;
    uint32_t ** old_plp_bp = MR_A->plp_bp;

    //-----
    // Fix the num bams and bam file names
    MR_A->num_bams += MR_B->num_bams;
    // realloc the memory for MR_A
    MR_A->bam_file_names = calloc(MR_A->num_bams, sizeof(char*));
    for (i = 0; i < old_num_bams; ++i) {
        MR_A->bam_file_names[i] = strdup(old_bam_file_names[i]);
    }
    for (j = 0; j < MR_B->num_bams; ++j) {
        MR_A->bam_file_names[i] = strdup(MR_B->bam_file_names[j]);
        ++i;
    }
    // free these original data
    if(old_bam_file_names != 0) {
        for(i = 0; i < old_num_bams; ++i) {
            if(old_bam_file_names[i] != 0)
                free(old_bam_file_names[i]);
        }
        free(old_bam_file_names);
    }

    //-----
    // Pileups
    MR_A->plp_bp = calloc(MR_A->num_contigs, sizeof(uint32_t*));
    for(i = 0; i < MR_A->num_contigs; ++i) {
        MR_A->plp_bp[i] = calloc(MR_A->num_bams, sizeof(uint32_t));
    }
    for(k = 0; k < MR_A->num_contigs; ++k)
    {
        for (i = 0; i < old_num_bams; ++i) {
            //printf("A: %d, %d, %d\n", k, i, old_plp_bp[k][i]);
            MR_A->plp_bp[k][i] = old_plp_bp[k][i];
        }
        for (j = 0; j < MR_B->num_bams; ++j) {
            //printf("B: %d, %d, %d\n", k, i, MR_B->plp_bp[k][j]);
            MR_A->plp_bp[k][i] = MR_B->plp_bp[k][j];
            ++i;
        }
    }
    if(old_plp_bp != 0) {
        for(i = 0; i < MR_A->num_contigs; ++i) {
            if(old_plp_bp[i] != 0)
                free(old_plp_bp[i]);
        }
        free(old_plp_bp);
    }

    //-----
    // Contig length correctors
    //
    if (strcmp(MR_A->coverage_mode, "outlier") == 0) {
        MR_A->contig_length_correctors = calloc(MR_A->num_contigs, sizeof(uint32_t*));
        for(i = 0; i < MR_A->num_contigs; ++i) {
            MR_A->contig_length_correctors[i] = calloc(MR_A->num_bams, sizeof(uint32_t));
        }
        for(k = 0; k < MR_A->num_contigs; ++k)
        {
            for (i = 0; i < old_num_bams; ++i) {
                MR_A->contig_length_correctors[k][i] = old_contig_length_correctors[k][i];
            }
            for (j = 0; j < MR_B->num_bams; ++j) {
                MR_A->contig_length_correctors[k][i] = MR_B->contig_length_correctors[k][j];
                ++i;
            }
        }
        if(old_contig_length_correctors != 0) {
            for(i = 0; i < MR_A->num_contigs; ++i) {
                if(old_contig_length_correctors[i] != 0)
                    free(old_contig_length_correctors[i]);
            }
            free(old_contig_length_correctors);
        }
    }


    //-----
    // Links
    if(MR_A->is_links_included)
    {
        // break the pattern of overwriting. We can just add the links here
        char **keys = NULL;
        size_t *key_sizes = NULL;
        size_t key_count = 0;
        keys = (char **)cfuhash_keys_data(MR_B->links, &key_count, &key_sizes, 0);

        for (i = 0; i < (int)key_count; i++) {
            BMM_link_pair * LP = cfuhash_get(MR_B->links, keys[i]);
            free(keys[i]);
            BMM_link_info* LI = LP->LI;
            do {
                addLink(MR_A->links,
                        LP->cid_1,
                        LP->cid_2,
                        LI->pos_1,
                        LI->pos_2,
                        LI->orient_1,
                        LI->orient_2,
                        LI->bam_ID+old_num_bams);
            } while(getNextLinkInfo(&LI));

        }
        free(keys);
        free(key_sizes);
    }
}

void destroy_MR(BMM_mapping_results * MR)
{
    //----
    // free the MR object
    //
    int i = 0;
    if(MR->num_contigs != 0 && MR->num_bams != 0) {
        if(MR->plp_bp != 0) {
            for(i = 0; i < MR->num_contigs; ++i) {
                if(MR->plp_bp[i] != 0)
                    free(MR->plp_bp[i]);
            }
            free(MR->plp_bp);
        }

        if(MR->bam_file_names != 0) {
            for(i = 0; i < MR->num_bams; ++i) {
                if(MR->bam_file_names[i] != 0)
                    free(MR->bam_file_names[i]);
            }
            free(MR->bam_file_names);
            free(MR->bam_file_name_lengths);
        }

        if(MR->contig_names != 0) {
            for(i = 0; i < MR->num_contigs; ++i) {
                if(MR->contig_names[i] != 0)
                    free(MR->contig_names[i]);
            }
            free(MR->contig_names);
            free(MR->contig_name_lengths);
        }

        if(MR->contig_lengths != 0)
            free(MR->contig_lengths);

        if(MR->contig_length_correctors != 0) {
            for(i = 0; i < MR->num_contigs; ++i) {
                if(MR->contig_length_correctors[i] != 0)
                    free(MR->contig_length_correctors[i]);
            }
            free(MR->contig_length_correctors);
        }
    }

    free(MR->coverage_mode);

    // destroy paired links
    if(MR->is_links_included)
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
                          int doLinks,
                          int ignoreSuppAlignments,
                          char* coverageMode,
                          char* bamFiles[],
                          BMM_mapping_results * MR
) {
    //-----
    // work out coverage depths and also pairwise linkages if asked to do so
    //
    int supp_check = 0x0; // include supp mappings
    if (ignoreSuppAlignments) {
        supp_check = BMM_BAM_FSUPP;
    }

    // initialize the auxiliary data structures
    const bam_pileup1_t **plp;
    bam_hdr_t *h = 0; // BAM header of the 1st input
    aux_t **data;
    bam_mplp_t mplp;
    int i = 0, tid = 0, *n_plp;
    // load contig names and BAM index.
    data = calloc(numBams, sizeof(void*)); // data[i] for the i-th input
    int beg = 0, end = 1<<30;  // set the default region

    for (i = 0; i < numBams; ++i) {
        data[i] = calloc(1, sizeof(aux_t));
        data[i]->fp = bgzf_open(bamFiles[i], "r"); // open BAM
        data[i]->min_mapQ = mapQ;                    // set the mapQ filter
        data[i]->min_len  = minLen;                 // set the qlen filter
        bam_hdr_t *htmp;
        htmp = bam_hdr_read(data[i]->fp);         // read the BAM header ( I think this must be done for each file for legitness!)
        if (i == 0) {
            h = htmp; // keep the header of the 1st BAM
        } else { bam_hdr_destroy(htmp); } // if not the 1st BAM, trash the header
    }

    // initialise the mapping results struct
    init_MR(MR,
            h,
            numBams,
            bamFiles,
            doLinks,
            coverageMode,
            ignoreSuppAlignments
           );

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
                position_holder[i] = calloc(MR->contig_lengths[tid], sizeof(uint32_t));
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
                else if(MR->is_links_included) {
                    // now we do links if we've been asked to
                    bam1_core_t core = p->b[0].core;
                    //printf("-->%d %d %d %d %d %d\n", (core.flag & BAM_FPAIRED), (core.flag & BAM_FREAD1), ((core.flag & BMM_BAM_FMAPPED) == 0), ((core.flag & supp_check) == 0), core.tid, core.mtid);
                    // check to see if this is a proper linking paired read
                    if ((core.flag & BAM_FPAIRED) &&                // read is a paired read
                        (core.flag & BAM_FREAD1) &&                 // read is first in pair (avoid dupe links)
                        ((core.flag & BMM_BAM_FMAPPED) == 0) &&      // both ends are mapped
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
                                i);                                 // bam file ID
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

void adjustPlpBp(BMM_mapping_results * MR,
                 uint32_t ** positionHolder,
                 int tid
) {
    uint32_t num_bams = MR->num_bams;
    uint32_t * plp_sum = calloc(num_bams, sizeof(uint32_t));
    int pos = 0, i = 0;
    if(strcmp(MR->coverage_mode, "outlier") == 0) {
        uint32_t * drops = calloc(num_bams, sizeof(uint32_t));
        for(i = 0; i < num_bams; ++i) {
            // set the cut off at a stdev either side of the mean
            float m = BMM_mean(positionHolder[i], MR->contig_lengths[tid]);
            float std = BMM_stdDev(positionHolder[i], MR->contig_lengths[tid], m);
            float lower_cut = ((m-std) < 0) ? 0 : (m-std);
            float upper_cut = m+std;
            for(pos = 0; pos < MR->contig_lengths[tid]; ++pos) {
                if((positionHolder[i][pos] <= upper_cut) && (positionHolder[i][pos] >= lower_cut)) {
                    // OK
                    plp_sum[i] += positionHolder[i][pos];
                } else {
                    // DROP
                    ++drops[i];
                }
            }
            MR->plp_bp[tid][i] = plp_sum[i];
            MR->contig_length_correctors[tid][i] = drops[i];
        }
        free(drops);
    } else {
        for(i = 0; i < num_bams; ++i) {
            for(pos = 0; pos < MR->contig_lengths[tid]; ++pos) {
                plp_sum[i] += positionHolder[i][pos];
            }
            MR->plp_bp[tid][i] = plp_sum[i];
        }
    }
    free(plp_sum);
}

float ** calculateCoverages(BMM_mapping_results * MR) {
    int i = 0, j = 0;
    if(MR->num_contigs != 0 && MR->num_bams != 0) {
        if(MR->plp_bp != NULL) {
            float ** ret_matrix = calloc(MR->num_contigs, sizeof(float*));
            for(i = 0; i < MR->num_contigs; ++i) {
                ret_matrix[i] = calloc(MR->num_bams, sizeof(float));
                for(j = 0; j < MR->num_bams; ++j) {
                    // print average coverages
                    if(strcmp(MR->coverage_mode, "outlier") == 0) {
                        // the counts are reduced so we should reduce the contig length accordingly
                        ret_matrix[i][j] = (float)MR->plp_bp[i][j]/(float)(MR->contig_lengths[i]-MR->contig_length_correctors[i][j]);
                    } else {
                        // otherwise it's just a straight up average
                        ret_matrix[i][j] = (float)MR->plp_bp[i][j]/(float)MR->contig_lengths[i];
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
    *** LINKS ***
    ***********************/

int initLW(BMM_LinkWalker * walker, BMM_mapping_results * MR)
{
    return initLinkWalker(walker, MR->links);
}

int stepLW(BMM_LinkWalker * walker)
{
    return stepLinkWalker(walker);
}

void destroyLW(BMM_LinkWalker * walker)
{
    destroyLinkWalker(walker);
}

    /***********************
    *** PRINTING AND I/O ***
    ***********************/

void printError(char* errorMessage, int line)
{
    //-----
    // Consistent way to print error messages
    //
    printf("ERROR: At line: %d\n\t%s\n\n", line, errorMessage);
}

void print_MR(BMM_mapping_results * MR) {
    int i = 0, j = 0;
    if(MR->num_contigs != 0 && MR->num_bams != 0) {
        if(MR->plp_bp != NULL) {
            float ** covs = calculateCoverages(MR); // get the coverages we want
            if(covs != NULL) {
                // print away!
                printf("Contig\tLength");
                for(j = 0; j < MR->num_bams; ++j) {
                    printf("\t%s",MR->bam_file_names[j]);
                }
                printf("\n");
                for(i = 0; i < MR->num_contigs; ++i) {
                    printf("%s\t%d", MR->contig_names[i], MR->contig_lengths[i]);
                    for(j = 0; j < MR->num_bams; ++j) {
                        printf("\t%0.4f", covs[i][j]);
                    }
                    printf("\n");
                }
                // we're responsible for cleaning up the covs structure
                destroyCoverages(covs, MR->num_contigs);
                if(MR->is_links_included) {
                    printLinks(MR->links, MR->bam_file_names, MR->contig_names);
                }
            }
        }
    }
}
