//#############################################################################
//
//   bamParser.h
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

#ifndef BMM_BAM_PARSER_H
  #define BMM_BAM_PARSER_H

// system includes
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>

// htslib
#include "htslib/bgzf.h"
#include "htslib/sam.h"

// cfuhash
#include "cfuhash.h"

// local includes
#include "pairedLink.h"

typedef BGZF bamFile;

#ifdef __cplusplus
extern "C" {
#endif


/*! @typedef
 @abstract Auxiliary data structure used in read_bam
 @field fp the file handler
 @field iter NULL if a region not specified
 @field min_mapQ mapQ filter
 @field min_len length filter
 */
typedef struct {                    //
    bamFile *fp;                    // the file handler
    hts_itr_t *iter;                // NULL if a region not specified
    int min_mapQ, min_len;          // mapQ filter; length filter
} aux_t;

/*! @typedef
 @abstract Structure for returning mapping results
 @field plp_bp                      number of bases piled up on each contig
 @field contig_lengths              lengths of the referene sequences
 @field contig_length_correctors    corrections to contig lengths used when doing outlier coverage
 @field num_bams                    number of BAM files parsed
 @field num_contigs                 number of reference sequences
 @field contig_names                names of the reference sequences
 @field bam_file_names              names of the bam files used in this mapping result
 @field contig_name_lengths         lengths of the names of the reference sequences
 @field bam_file_name_lengths       lengths of the names of the bam files used in this mapping result
 @field is_links_included           are links being calculated
 @field coverage_mode               type of coverage to be calculate ['vanilla', 'outlier']
 @field is_ignore_supps             are supplementary alignments being ignored
 @field links                       linking pairs
 */
typedef struct {
    uint32_t ** plp_bp;
    uint32_t * contig_lengths;
    uint32_t ** contig_length_correctors;
    uint32_t num_bams;
    uint32_t num_contigs;
    char ** contig_names;
    char ** bam_file_names;
    uint16_t * contig_name_lengths;
    uint16_t * bam_file_name_lengths;
    int is_links_included;
    char * coverage_mode;
    int is_ignore_supps;
    cfuhash_table_t * links;
} BMM_mapping_results;

//#############################################################################
// BAM parsing functions
//#############################################################################


int read_bam(void *data,
             bam1_t *b);


/*!
 * @abstract Allocate space for a new MR struct
 *
 * @return BMM_mapping_results *
 *
 * @discussion Allocates the memore which should be initialised at some point in time.
 * You call destroy_MR to free this memory
 */
BMM_mapping_results * create_MR(void);

/*!
 * @abstract Initialise the mapping results struct
 *
 * @param MR                    mapping results struct to initialise
 * @param BAM_header            htslib BAM header
 * @param numBams               number of BAM files to parse
 * @param bamFiles              filenames of BAMs parsed
 * @param doLinks               set to 1 if should initialise contig_length_correctors
 * @param coverageMode          type of coverage to be calculated
 * @param ignoreSuppAlignments  only use primary alignments
 * @return void
 *
 * @discussion If you call this function then you MUST call destroy_MR
 * when you're done.
 */
void init_MR(BMM_mapping_results * MR,
             bam_hdr_t * BAM_header,
             int numBams,
             char * bamFiles[],
             int doLinks,
             char * coverageMode,
             int ignoreSuppAlignments
);

/*!
 * @abstract Merge the contents of MR_B into MR_A
 *
 * @param  MR_A  mapping results struct to copy to
 * @param  MR_B  mapping results struct to copy from
 * @return void
 *
 * @discussion MR_B remains unchanged.
 * MR_A is updated to include all the info contained in MR_B
 *
 * NOTE: We assume that all the haders of all the files are in sync. If they
 * are not, if contigs have been removed etc, then doom will swiftly follow.
 * Also, we assume that flags like do_links, do_outlier match... ...or DOOM!
 */
void merge_MRs(BMM_mapping_results * MR_A, BMM_mapping_results * MR_B);

/*!
 * @abstract Free all the memory calloced in init_MR
 *
 * @param  MR  mapping results struct to destroy
 * @return void
 */
void destroy_MR(BMM_mapping_results * MR);

/*!
 * @abstract Initialise the mapping results struct <- read in the BAM files
 *
 * @param numBams               number of BAM files to parse
 * @param baseQ                 base quality threshold
 * @param mapQ                  mapping quality threshold
 * @param minLen                min query length
 * @param doLinks               1 if links should be calculated
 * @param ignoreSuppAlignments  only use primary alignments
 * @param coverageMode          type of coverage to be calculated
 * @param bamFiles              filenames of BAM files to parse
 * @param MR                    mapping results struct to write to
 * @return 0 for success
 *
 * @discussion This function expects MR to be a null pointer. It calls
 * init_MR and stores info accordingly. TL;DR If you call this function
 * then you MUST call destroy_MR when you're done.
 *
 */
int parseCoverageAndLinks(int numBams,
                          int baseQ,
                          int mapQ,
                          int minLen,
                          int doLinks,
                          int ignoreSuppAlignments,
                          char* coverageMode,
                          char* bamFiles[],
                          BMM_mapping_results * MR);

/*!
 * @abstract Adjust (reduce) the number of piled-up bases along a contig
 *
 * @param  MR  mapping results struct to write to
 * @param  position_holder  array of pileup depths
 * @param  tid  contig currently being processed
 * @return void
 *
 * @discussion This function expects MR to be initialised.
 * it can change the values of contig_length_correctors and plp_bp
 */
void adjustPlpBp(BMM_mapping_results * MR,
                 uint32_t ** position_holder,
                 int tid);

/*!
 * @abstract Calculate the coverage for each contig for each BAM
 *
 * @param  MR  mapping results struct with mapping info
 * @return matrix of floats (rows = contigs, cols = BAMs)
 *
 * @discussion This function expects MR to be initialised.
 * NOTE: YOU are responsible for freeing the return value
 * recommended method is to use destroyCoverages
 */
float ** calculateCoverages(BMM_mapping_results * MR);

/*!
 * @abstract Destroy the coverages structure made in  calculateCoverages
 *
 * @param covs array to destroy
 * @param numContigs number of rows in array
 * @return void
 */
void destroyCoverages(float ** covs, int numContigs);


    /***********************
    *** LINKS ***
    ***********************/
/*!
 * @abstract Start moving through all of the links
 *
 * @param  MR  mapping results struct containing links
 * @return 1 if links OK or 0 otherwise
 */
int initLW(BMM_LinkWalker * walker, BMM_mapping_results * MR);

/*!
 * @abstract Move to the next LinkInfo or LinkPair
 *
 * @param  walker   pointer to LinkHolder.
 * @return 1 for step within current contig pair, 2 for new pair, 0 for end walk
 */
int stepLW(BMM_LinkWalker * walker);

/*!
 * @abstract Start moving through all of the links
 *
 * @param  walker   pointer to LinkHolder.
 * @return void
 */
void destroyLW(BMM_LinkWalker * walker);

    /***********************
    *** PRINTING AND I/O ***
    ***********************/
/*!
 * @abstract Print an error to stdout
 *
 * @param  errorMessage to print
 * @param  line Line number where called from
 * @return void
 */
 void printError(char* errorMessage, int line);

/*!
 * @abstract Print the Mapping results to stdout
 *
 * @param  MR  mapping results struct to print
 * @return void
 */
void print_MR(BMM_mapping_results * MR);

#ifdef __cplusplus
}
#endif

#endif // BMM_BAM_PARSER_H

