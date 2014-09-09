//#############################################################################
//
//   bamParser.c
//
//   Engine for parsing BAM files - non-parallel
//   Functions include:
//   Determine "type" of bam file (insert size, orientation of reads etc...)
//   Determine average coverage values
//   Find linking read pairs
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
#include <libgen.h>

// htslib
//#include "htslib/bgzf.h"
#include "htslib/sam.h"

// cfuhash
#include "cfuhash.h"

// local includes
#include "bamParser.h"
#include "pairedLink.h"

BM_fileInfo * createBFI(void) {
    return (BM_fileInfo*) calloc(1, sizeof(BM_fileInfo));
}

void initBFI(BM_fileInfo * BFI,
             bam_hdr_t * BAM_header,
             int numBams,
             char * bamFiles[],
             int * links,
             char * coverageMode,
             int ignoreSuppAlignments
            ) {
    // check for a valid coverage mode
    int valid_coverage_mode = strcmp(coverageMode, "vanilla");
    if(valid_coverage_mode != 0)
    { valid_coverage_mode = strcmp(coverageMode, "outlier"); }
    if(valid_coverage_mode != 0)
    { valid_coverage_mode = strcmp(coverageMode, "none"); }
    if(valid_coverage_mode != 0) {
        char str[80];
        sprintf(str, "Invalid coverage mode '%s'", coverageMode);
        printError(str, __LINE__);
        return;
    }

    int i = 0; int j = 0;
    BFI->numContigs = BAM_header->n_targets;
    BFI->numBams = numBams;
    BFI->isLinks = (links != 0);
    BFI->coverage_mode = strdup(coverageMode);
    BFI->isIgnoreSupps = ignoreSuppAlignments;
    BFI->contigLengthCorrectors = 0;
    BFI->links = 0;

    if(BFI->numContigs != 0 && BFI->numBams != 0) {
        if(strcmp(coverageMode, "none") != 0) {
            // unless this is an extraction run we will
            // need to make room to store read counts
            BFI->plpBp = calloc(BFI->numContigs, sizeof(uint32_t*));
            for(i = 0; i < BFI->numContigs; ++i) {
                BFI->plpBp[i] = calloc(BFI->numBams, sizeof(uint32_t));
            }
        }

        // make BAM file structs
        BFI->bamFiles = (BM_bamFile**) calloc(BFI->numBams, sizeof(BM_bamFile*));
        for (i = 0; i < numBams; ++i) {
            BM_bamFile * BF = (BM_bamFile*) calloc(1, sizeof(BM_bamFile));
            BF->fileName = strdup(bamFiles[i]);
            BF->fileNameLength = strlen(bamFiles[i]);

            // support for mixed insert sizes in one bam (i.e. shadow library)
            BF->numTypes = 0;
            BF->types = 0;
            if((strcmp(coverageMode, "none") != 0) || BFI->isLinks) {
                if(links != 0) {
                    BF->numTypes = links[i];
                    BF->types = (BM_bamType**) calloc(BF->numTypes, sizeof(BM_bamType*));
                    for(j = 0; j < BF->numTypes; ++j) {
                        BM_bamType * BT = (BM_bamType*) calloc(1, sizeof(BM_bamType));
                        BT->orientationType = OT_NONE;
                        BT->insertSize = 0;
                        BT->insertStdev = 0;
                        BT->supporting = 0;
                        BF->types[j] = BT;
                    }
                }
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

        //----------------------------
        // only allocate if we NEED to
        //----------------------------
        if (strcmp(BFI->coverage_mode, "outlier") == 0) {
            BFI->contigLengthCorrectors = calloc(BFI->numContigs, sizeof(uint32_t*));
            for(i = 0; i < BFI->numContigs; ++i) {
                BFI->contigLengthCorrectors[i] = calloc(BFI->numBams, sizeof(uint32_t));
            }
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
        sprintf(str, "Unequal number of contigs in BFI structs to be merged (Num contigs: %d, %d)", BFI_A->numContigs, BFI_B->numContigs);
        printError(str, __LINE__);
        return;
    }

    // now check that the lengths of all the contigs are the same
    // just check every 100th contig
    for(i=0; i < BFI_A->numContigs; ++i) {
        if(i < BFI_A->numContigs) {
            if(BFI_A->contigLengths[i] != BFI_B->contigLengths[i]) {
                char str[80];
                sprintf(str, "Unequal contig lengths in BFI structs to be merged (Contig#: %d, Lengths: %d, %d)", i, BFI_A->numContigs, BFI_B->numContigs);
                printError(str, __LINE__);
                return;
            }
        }
        else {
            break;
        }
    }

    // we can assume that the headers are the same. So now time to merge the data
    // keep a backup of these guys
    uint32_t old_numBams = BFI_A->numBams;
    BM_bamFile ** old_bamFiles = BFI_A->bamFiles;
    uint32_t ** old_contigLengthCorrectors = BFI_A->contigLengthCorrectors;
    uint32_t ** old_plpBp = BFI_A->plpBp;

    //-----
    // BAM files
    BFI_A->numBams += BFI_B->numBams;
    // realloc the memory for BFI_A
    BFI_A->bamFiles = (BM_bamFile**) calloc(BFI_A->numBams, sizeof(BM_bamFile*));
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
    }

    //-----
    // Pileups
    BFI_A->plpBp = calloc(BFI_A->numContigs, sizeof(uint32_t*));
    for(i = 0; i < BFI_A->numContigs; ++i) {
        BFI_A->plpBp[i] = calloc(BFI_A->numBams, sizeof(uint32_t));
    }
    for(k = 0; k < BFI_A->numContigs; ++k)
    {
        for (i = 0; i < old_numBams; ++i) {
            //printf("A: %d, %d, %d\n", k, i, old_plpBp[k][i]);
            BFI_A->plpBp[k][i] = old_plpBp[k][i];
        }
        for (j = 0; j < BFI_B->numBams; ++j) {
            //printf("B: %d, %d, %d\n", k, i, BFI_B->plpBp[k][j]);
            BFI_A->plpBp[k][i] = BFI_B->plpBp[k][j];
            ++i;
        }
    }
    if(old_plpBp != 0) {
        for(i = 0; i < BFI_A->numContigs; ++i) {
            if(old_plpBp[i] != 0)
                free(old_plpBp[i]);
        }
        free(old_plpBp);
    }

    //-----
    // Contig length correctors
    //
    if (strcmp(BFI_A->coverage_mode, "outlier") == 0) {
        BFI_A->contigLengthCorrectors = calloc(BFI_A->numContigs, sizeof(uint32_t*));
        for(i = 0; i < BFI_A->numContigs; ++i) {
            BFI_A->contigLengthCorrectors[i] = calloc(BFI_A->numBams, sizeof(uint32_t));
        }
        for(k = 0; k < BFI_A->numContigs; ++k)
        {
            for (i = 0; i < old_numBams; ++i) {
                BFI_A->contigLengthCorrectors[k][i] = old_contigLengthCorrectors[k][i];
            }
            for (j = 0; j < BFI_B->numBams; ++j) {
                BFI_A->contigLengthCorrectors[k][i] = BFI_B->contigLengthCorrectors[k][j];
                ++i;
            }
        }
        if(old_contigLengthCorrectors != 0) {
            for(i = 0; i < BFI_A->numContigs; ++i) {
                if(old_contigLengthCorrectors[i] != 0)
                    free(old_contigLengthCorrectors[i]);
            }
            free(old_contigLengthCorrectors);
        }
    }

    //-----
    // Links
    if(BFI_A->isLinks) {
        // break the pattern of overwriting. We can just add the links here
        char **keys = NULL;
        size_t *key_sizes = NULL;
        size_t key_count = 0;
        keys = (char **)cfuhash_keys_data(BFI_B->links, &key_count, &key_sizes, 0);
        for (i = 0; i < (int)key_count; i++) {
            BM_linkPair * LP = cfuhash_get(BFI_B->links, keys[i]);
            free(keys[i]);
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
        free(keys);
        free(key_sizes);
    }
    destroyBFI(BFI_B);
}

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

int parseCoverageAndLinks(int typeOnly,
                          int numBams,
                          int baseQ,
                          int mapQ,
                          int minLen,
                          int * links,
                          int ignoreSuppAlignments,
                          char* coverageMode,
                          char* bamFiles[],
                          BM_fileInfo * BFI
) {
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
    data = (aux_t**) calloc(numBams, sizeof(aux_t*)); // data[i] for the i-th input
    int beg = 0, end = 1<<30;  // set the default region

    for (i = 0; i < numBams; ++i) {
        bam_hdr_t *htmp;
        data[i] = (aux_t*) calloc(1, sizeof(aux_t));
        data[i]->fp = bgzf_open(bamFiles[i], "r"); // open BAM
        data[i]->min_mapQ = mapQ;                  // set the mapQ filter
        data[i]->min_len  = minLen;                // set the qlen filter
        htmp = bam_hdr_read(data[i]->fp);          // read the BAM header ( I think this must be done for each file for legitness!)
        if (i == 0) {
            h = htmp; // keep the header of the 1st BAM
        } else { bam_hdr_destroy(htmp); } // if not the 1st BAM, trash the header
    }

    // initialise the BFI struct
    initBFI(BFI,
            h,
            numBams,
            bamFiles,
            links,
            coverageMode,
            ignoreSuppAlignments
           );

    // type the bam files, but only if we're doing links
    if((links != 0) || typeOnly) {
        typeBamFiles(BFI);
    }

    // if we're only here to type the file then we're done
    if(typeOnly) {
        bam_hdr_destroy(h);

        for (i = 0; i < numBams; ++i) {
            bgzf_close(data[i]->fp);
            if (data[i]->iter) bam_itr_destroy(data[i]->iter);
            free(data[i]);
        }
        free(data);
        return 0;
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
                adjustPlpBp(BFI, position_holder, prev_tid);
                for (i = 0; i < numBams; ++i) {
                    free(position_holder[i]); // free this up
                }
            }
            for (i = 0; i < numBams; ++i) {
                // reset for next contig
                position_holder[i] = calloc(BFI->contigLengths[tid], sizeof(uint32_t));
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
                else if(BFI->isLinks) {
                    // now we do links if we've been asked to
                    bam1_core_t core = p->b[0].core;
                    //printf("-->%d %d %d %d %d %d\n", (core.flag & BAM_FPAIRED), (core.flag & BAM_FREAD1), ((core.flag & BM_BAM_FMAPPED) == 0), ((core.flag & supp_check) == 0), core.tid, core.mtid);
                    // check to see if this is a proper linking paired read
                    if ((core.pos == pos) &&                        // make sure it's the first time we've seen this link
                        (core.flag & BAM_FPAIRED) &&                // read is a paired read
                        (core.flag & BAM_FREAD1) &&                 // read is first in pair (avoid dupe links)
                        ((core.flag & BM_BAM_FMAPPED) == 0) &&      // both ends are mapped
                        ((core.flag & supp_check) == 0) &&          // is primary mapping (optional)
                        core.tid != core.mtid) {                    // hits different contigs

                        // looks legit
                        // get a link info
                        BM_linkInfo * LI = makeLinkInfo(core.tid,                           // contig 1
                                                        core.mtid,                          // contig 2
                                                        core.pos,                           // pos 1
                                                        core.mpos,                          // pos 2
                                                        ((core.flag&BAM_FREVERSE) != 0),    // 1 == reversed
                                                        ((core.flag&BAM_FMREVERSE) != 0),   // 0 = agrees
                                                        i                                   // bam file ID
                                                        );

                        // add the link
                        addLink(BFI->links,
                                LI,
                                core.tid,
                                core.mtid);
                    }
                }
            }
            position_holder[i][pos] = n_plp[i] - rejects; // add this position's depth
        }
    }

    if(prev_tid != -1) {
        // at the end of a contig
        adjustPlpBp(BFI, position_holder, prev_tid);
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

void extractReads(char ** bamFiles,
                  int numBams,
                  char ** contigs,
                  int numContigs) {
    //-----
    // code uses the pattern outlined in samtools view (sam_view.c)
    // thanks lh3!
    //
    int i = 0;
    // always ignoreSuppAlignments
    int supp_check = BM_BAM_FSUPP;

    samFile *in = 0, *out = 0;
    bam_hdr_t *header = NULL;

    // let's build a BFI
    // open file handlers
    BM_fileInfo * BFI = 0;
    if ((in = sam_open(bamFiles[0], "r")) == 0) {
        fprintf(stderr, "ERROR: Failed to open \"%s\" for reading.\n", bamFiles[0]);
    } else {
        if ((header = sam_hdr_read(in)) == 0) {
            fprintf(stderr, "ERROR: Failed to read the header from \"%s\".\n", (BFI->bamFiles[i])->fileName);
        } else {
            // samfile and header are legit!
            BFI = calloc(1, sizeof(BM_fileInfo));
            initBFI(BFI,
                    header,
                    numBams,
                    bamFiles,
                    0,
                    "none",
                    1);
        }
        if (in) sam_close(in);
    }

    out = sam_open("-", "w");

    for (i = 0; i < BFI->numBams; ++i) {
        // open file handlers --> use the same header throughout
        if ((in = sam_open((BFI->bamFiles[i])->fileName, "r")) == 0) {
            fprintf(stderr, "ERROR: Failed to open \"%s\" for reading.\n", (BFI->bamFiles[i])->fileName);
            if (in) sam_close(in);
            continue;
        }

        // retrieve alignments in specified regions
        bam1_t *b;
        hts_idx_t *idx = sam_index_load(in, (BFI->bamFiles[i])->fileName); // load index
        if (idx == 0) { // index is unavailable
            fprintf(stderr, "ERROR: Random alignment retrieval only works for indexed BAM or CRAM files.\n");
            if (in) sam_close(in);
            continue;
        }
        b = bam_init1();

        cfuhash_table_t *pairs = cfuhash_new_with_initial_size(3000);
        cfuhash_set_flag(pairs, CFUHASH_FROZEN_UNTIL_GROWS);

		char * pretty_name = calloc(1000, sizeof(char));
		pretty_name = strdup(basename(((BFI)->bamFiles[i])->fileName));
		fprintf(stdout, "%d\n", (int)strlen(pretty_name)); fflush(stdout);

        int result = -1;

        int hh;
        for (hh = 0; hh < numContigs; ++hh) {
            hts_itr_t *iter = sam_itr_querys(idx, header, contigs[hh]); // parse a region in the format like `chr2:100-200'
            if (iter == NULL) { // reference name is not found
                fprintf(stderr, "WARNING: Could not find contig: [%s] in BAM: [%s].\n", contigs[hh], (BFI->bamFiles[i])->fileName);
            }
            // int contig_length = header->target_len[hh];

            // fetch alignments
            while ((result = sam_itr_next(in, iter, b)) >= 0) {
                bam1_core_t core = b->core;
                printf(">b_%s;c_%s;r_%s;\n", pretty_name, contigs[hh], bam_get_qname(b));
                printf(">%s\n", bam_get_qname(b));
                uint8_t *s = bam_get_seq(b);
                for (i = 0; i < core.l_qseq; ++i) { printf("%c", "=ACMGRSVTWYHKDBN"[bam_seqi(s, i)]); }
                printf("\n");
                s = bam_get_qual(b);
                if (s[0] != 0xff) {
                    for (i = 0; i < core.l_qseq; ++i) { printf("%c", s[i] + 33); }
                    printf("\n");
                }
                continue;

                // check to see if this is a proper linking paired read
                if ((core.flag & BAM_FPAIRED) &&                // read is a paired read
                    (core.flag & BAM_FREAD1) &&                 // read is first in pair (avoid dupe links)
                    ((core.flag & BM_BAM_FMAPPED) == 0) &&      // both ends are mapped
                    ((core.flag & supp_check) == 0) &&          // is primary mapping (optional)
                    core.tid == core.mtid) {                    // hits same contig

                    sam_write1(out, header, b);
                    /*
                    printf("FOUND: CID: %d, P: %d, MP: %d, R: %d, BFI: %d, Contig: %s, BAM: %s\n",
                            core.mtid,
                            core.pos,
                            core.mpos,
                            ((core.flag&BAM_FREVERSE) != 0),
                            ((core.flag&BAM_FMREVERSE) != 0),
                            contigs[hh],
                            (BFI->bamFiles[i])->fileName);
                    */
                }
            }

            hts_itr_destroy(iter);
            if (result < -1) {
                fprintf(stderr, "ERROR: retrieval of reads from contig: \"%s\" failed due to truncated file or corrupt BAM index file\n", header->target_name[hh]);
                break;
            }
        }

        free(pretty_name);
        bam_destroy1(b);
        hts_idx_destroy(idx); // destroy the BAM index

        cfuhash_clear(pairs);

    }
    if ( header ) bam_hdr_destroy(header);
    if ( BFI ) destroyBFI(BFI);
}

void typeBamFiles(BM_fileInfo * BFI) {
    //-----
    // code uses the pattern outlined in samtools view (sam_view.c)
    // thanks lh3!
    //
    int i = 0, j = 0;

    // always ignoreSuppAlignments
    int supp_check = BM_BAM_FSUPP;

    uint32_t * num_found = calloc(BFI->numBams, sizeof(uint32_t));     // the number of good links found for each bam
    uint32_t *** inserts = calloc(BFI->numBams, sizeof(uint32_t**));   // observed insert sizes per bam per OT
    int ** orient_counts = calloc(BFI->numBams, sizeof(int*));         // hold counts for each of the orientation types we'll see

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
        hts_idx_t *idx = sam_index_load(in, (BFI->bamFiles[i])->fileName); // load index
        if (idx == 0) { // index is unavailable
            fprintf(stderr, "ERROR: Random alignment retrieval only works for indexed BAM or CRAM files.\n");
            if (in) sam_close(in);
            if ( header ) bam_hdr_destroy(header);
            continue;
        }
        b = bam_init1();
        int result = -1;

        int hh;
        for (hh = 0; hh < header->n_targets; ++hh) {
            hts_itr_t *iter = sam_itr_querys(idx, header, header->target_name[hh]); // parse a region in the format like `chr2:100-200'
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
            while ((num_found[i] <= BM_PAIRS_FOR_TYPE) && ((result = sam_itr_next(in, iter, b)) >= 0)) {
                bam1_core_t core = b->core;
                // check to see if this is a proper linking paired read
                if ((core.flag & BAM_FPAIRED) &&                // read is a paired read
                    (core.flag & BAM_FREAD1) &&                 // read is first in pair (avoid dupe links)
                    ((core.flag & BM_BAM_FMAPPED) == 0) &&      // both ends are mapped
                    ((core.flag & supp_check) == 0) &&          // is primary mapping (optional)
                    core.tid == core.mtid) {                    // hits same contig

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
                fprintf(stderr, "ERROR: retrieval of reads from contig: \"%s\" failed due to truncated file or corrupt BAM index file\n", header->target_name[hh]);
                break;
            }

            // have we found enough?
            if(num_found[i] >= BM_PAIRS_FOR_TYPE)
                break;
            // else next contig!
        }
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
            ((BFI->bamFiles[i])->types[type_index])->insertSize = BM_mean(inserts[i][type], max);
            ((BFI->bamFiles[i])->types[type_index])->insertStdev = BM_stdDev(inserts[i][type], max, ((BFI->bamFiles[i])->types[type_index])->insertSize);
            ((BFI->bamFiles[i])->types[type_index])->supporting = max;
            // reset this so that we find the next best one
            orient_counts[i][type] = 0;
            --num_to_find;
        }
    }

    // clean up
    free(num_found);
    for (i = 0; i < BFI->numBams; ++i) {
        free(orient_counts[i]);
        for (j = 0; j < 3; ++j) {
            free(inserts[i][j]);
        }
        free(inserts[i]);
    }
    free(inserts);
    free(orient_counts);
}

void adjustPlpBp(BM_fileInfo * BFI,
                 uint32_t ** positionHolder,
                 int tid
) {
    uint32_t num_bams = BFI->numBams;
    uint32_t * plp_sum = calloc(num_bams, sizeof(uint32_t));
    int pos = 0, i = 0;
    if(strcmp(BFI->coverage_mode, "outlier") == 0) {
        uint32_t * drops = calloc(num_bams, sizeof(uint32_t));
        for(i = 0; i < num_bams; ++i) {
            // set the cut off at a stdev either side of the mean
            float m = BM_mean(positionHolder[i], BFI->contigLengths[tid]);
            float std = BM_stdDev(positionHolder[i], BFI->contigLengths[tid], m);
            float lower_cut = ((m-std) < 0) ? 0 : (m-std);
            float upper_cut = m+std;
            for(pos = 0; pos < BFI->contigLengths[tid]; ++pos) {
                if((positionHolder[i][pos] <= upper_cut) && (positionHolder[i][pos] >= lower_cut)) {
                    // OK
                    plp_sum[i] += positionHolder[i][pos];
                } else {
                    // DROP
                    ++drops[i];
                }
            }
            BFI->plpBp[tid][i] = plp_sum[i];
            BFI->contigLengthCorrectors[tid][i] = drops[i];
        }
        free(drops);
    } else {
        for(i = 0; i < num_bams; ++i) {
            for(pos = 0; pos < BFI->contigLengths[tid]; ++pos) {
                plp_sum[i] += positionHolder[i][pos];
            }
            BFI->plpBp[tid][i] = plp_sum[i];
        }
    }
    free(plp_sum);
}

float ** calculateCoverages(BM_fileInfo * BFI) {
    int i = 0, j = 0;
    if(BFI->numContigs != 0 && BFI->numBams != 0) {
        if(BFI->plpBp != NULL) {
            float ** ret_matrix = calloc(BFI->numContigs, sizeof(float*));
            for(i = 0; i < BFI->numContigs; ++i) {
                ret_matrix[i] = calloc(BFI->numBams, sizeof(float));
                for(j = 0; j < BFI->numBams; ++j) {
                    // print average coverages
                    if(strcmp(BFI->coverage_mode, "outlier") == 0) {
                        // the counts are reduced so we should reduce the contig length accordingly
                        ret_matrix[i][j] = (float)BFI->plpBp[i][j]/(float)(BFI->contigLengths[i]-BFI->contigLengthCorrectors[i][j]);
                    } else {
                        // otherwise it's just a straight up average
                        ret_matrix[i][j] = (float)BFI->plpBp[i][j]/(float)BFI->contigLengths[i];
                    }
                }
            }
            return ret_matrix;
        }
    }
    return NULL;
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
            if(BFI->plpBp != 0) {
                for(i = 0; i < BFI->numContigs; ++i) {
                    if(BFI->plpBp[i] != 0)
                        free(BFI->plpBp[i]);
                }
                free(BFI->plpBp);
            }

            if(BFI->bamFiles != 0) {
                destroyBamFiles(BFI->bamFiles, BFI->numBams);
            }

            if(BFI->contigNames != 0) {
                for(i = 0; i < BFI->numContigs; ++i) {
                    if(BFI->contigNames[i] != 0)
                        free(BFI->contigNames[i]);
                }
                free(BFI->contigNames);
                free(BFI->contigNameLengths);
            }

            if(BFI->contigLengths != 0)
                free(BFI->contigLengths);

            if(BFI->contigLengthCorrectors != 0) {
                for(i = 0; i < BFI->numContigs; ++i) {
                    if(BFI->contigLengthCorrectors[i] != 0)
                        free(BFI->contigLengthCorrectors[i]);
                }
                free(BFI->contigLengthCorrectors);
            }
        }

        free(BFI->coverage_mode);

        // destroy paired links
        if(BFI->isLinks) {
            destroyLinks(BFI->links);
            cfuhash_clear(BFI->links);
            cfuhash_destroy(BFI->links);
        }
    }
}

void destroyCoverages(float ** covs, int numContigs) {
    int i = 0;
    for(i = 0; i < numContigs; ++i) {
        free(covs[i]);
    }
    free(covs);
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
                free(BFs[i]->fileName);
                for(j = 0; j < BFs[i]->numTypes; ++j) {
                    free(BFs[i]->types[j]);
                }
            }
        }
        free(BFs);
    }
}

void destroyPairedRead(BM_pairedRead * PR) {
    if(PR) {
        destroyRead(PR->forward);
        destroyRead(PR->reverse);
    }
}

void destroyRead(BM_seqRead * SR) {
    if(SR) {
        if(SR->seqId) free(SR->seqId);
        if(SR->seq) free(SR->seq);
        if(SR->qual) free(SR->qual);
        free(SR);
    }
}

BM_seqRead * makeSeqRead(char * seqId,
                         char * seq,
                         char * qual) {
    BM_seqRead * SR = calloc(1, sizeof(BM_seqRead));
    SR->seqId = strdup(seqId);
    SR->seq = strdup(seq);
    SR->len = strlen(seq);
    if(qual)
        SR->qual = strdup(qual);
    else
        SR->qual = 0;
    return SR;
}

BM_pairedRead * makePairedRead(BM_seqRead * fSR, BM_seqRead * rSR) {
    BM_pairedRead * PR = calloc(1, sizeof(BM_pairedRead));
    PR->forward = fSR;
    PR->reverse = rSR;
    return PR;
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
        if(BFI->plpBp != NULL) {
            float ** covs = calculateCoverages(BFI); // get the coverages we want
            if(covs != NULL) {
                // print away!
                printf("#contig\tlength");
                for(j = 0; j < BFI->numBams; ++j) {
                    printf("\t%s",(BFI->bamFiles[j])->fileName);
                }
                printf("\n");
                for(i = 0; i < BFI->numContigs; ++i) {
                    printf("%s\t%d", BFI->contigNames[i], BFI->contigLengths[i]);
                    for(j = 0; j < BFI->numBams; ++j) {
                        printf("\t%0.4f", covs[i][j]);
                    }
                    printf("\n");
                }
                printf("---\n");
                // we're responsible for cleaning up the covs structure
                destroyCoverages(covs, BFI->numContigs);
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

float BM_mean(uint32_t * values, uint32_t size) {
    uint32_t sum = 0;
    int i = 0;
    for( i = 0; i < size; ++i) {
        sum += *(values + i);
    }
    return (float)sum/(float)size;
}

float BM_stdDev(uint32_t * values, uint32_t size, float m) {
    float sum = 0;
    int i = 0;
    if(m == -1)
        m = BM_mean(values, size);
    for(i = 0; i < size; ++i) {
        sum += pow((float)*(values + i) - m, 2);
    }
    return sqrt(sum/(float)size);
}

float BM_fakeStdDev(uint32_t * values, uint32_t size) {
    // everything is 3 stdevs from the mean right?
    uint32_t max = 0, min = 1<<30;
    int i = 0;
    for( i = 0; i < size; ++i) {
        if (*(values + i) > max)
            max = *(values + i);
        else if (*(values + i) < min)
            min = *(values + i);
    }
    return (float)(max-min)/6;
}
