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

#ifndef BM_BAM_PARSER_H
  #define BM_BAM_PARSER_H

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

// Orientation types
typedef enum {OT_OUT, OT_SAME, OT_IN, OT_NONE, OT_ERROR} OT;

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
 @abstract Structure for storing information about the orientation type and insert sizes for a given bam file
 @field orientationType     orientation of the reads in this mapping
 @field insertSize          insert size of the reads (outer)
 @field insertStdev         standard deviation measured for the insert size
 @field numTested           number of reads tested to get these stats
 */
 typedef struct {
    int orientationType;                // actually going to use the ENUM, but worried about how to get around this with python ctypes
    float insertSize;
    float insertStdev;
    int supporting;
 } BM_bamType;


/*! @typedef
 @abstract Structure for storing information about a single BAM file
 @field fileName            filename of the BAM
 @field fileNameLength      length of the filename
 @field types               the orientation types and insert sizes accociated with the bam file
 @field numTypes            the number of orientation types for the bam file
 */
 typedef struct {
    char * fileName;
    uint16_t fileNameLength;
    BM_bamType ** types;
    int numTypes;
 } BM_bamFile;

/*! @typedef
 @abstract Structure for returning  mapping results
 @field plpBp                       number of bases piled up on each contig
 @field contigLengths               lengths of the referene sequences
 @field contigLengthCorrectors      corrections to contig lengths used when doing outlier coverage
 @field numBams                     number of BAM files parsed
 @field numContigs                  number of reference sequences
 @field bamFiles                    array of pointers to BM_bamFile's used in this mapping result
 @field contigNames                 names of the reference sequences
 @field contig_name_lengths         lengths of the names of the reference sequences
 @field isLinks                     are links being calculated
 @field coverage_mode               type of coverage to be calculate ['vanilla', 'outlier']
 @field isIgnoreSupps               are supplementary alignments being ignored
 @field links                       linking pairs
 */
typedef struct {
    uint32_t ** plpBp;
    uint32_t * contigLengths;
    uint32_t ** contigLengthCorrectors;
    uint32_t numBams;
    uint32_t numContigs;
    BM_bamFile ** bamFiles;
    char ** contigNames;
    uint16_t * contig_name_lengths;
    int isLinks;
    char * coverage_mode;
    int isIgnoreSupps;
    cfuhash_table_t * links;
} BM_mappingResults;

//#############################################################################
// BAM parsing functions
//#############################################################################


int read_bam(void *data,
             bam1_t *b);

/*!
 * @abstract Allocate space for a new MR struct
 *
 * @return BM_mappingResults *
 *
 * @discussion Allocates the memore which should be initialised at some point in time.
 * You call destroy_MR to free this memory
 */
BM_mappingResults * create_MR(void);

/*!
 * @abstract Initialise the mapping results struct
 *
 * @param MR                    mapping results struct to initialise
 * @param BAM_header            htslib BAM header
 * @param numBams               number of BAM files to parse
 * @param bamFiles              filenames of BAMs parsed
 * @param links                 array of OT counts per bam or 0
 * @param coverageMode          type of coverage to be calculated
 * @param ignoreSuppAlignments  only use primary alignments
 * @return void
 *
 * @discussion If you call this function then you MUST call destroy_MR
 * when you're done.
 */
void init_MR(BM_mappingResults * MR,
             bam_hdr_t * BAM_header,
             int numBams,
             char * bamFiles[],
             int * links,
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
void merge_MRs(BM_mappingResults * MR_A, BM_mappingResults * MR_B);

/*!
 * @abstract Free all the memory calloced in init_MR
 *
 * @param  MR  mapping results struct to destroy
 * @return void
 */
void destroy_MR(BM_mappingResults * MR);

/*!
 * @abstract Initialise the mapping results struct <- read in the BAM files
 *
 * @param numBams               number of BAM files to parse
 * @param baseQ                 base quality threshold
 * @param mapQ                  mapping quality threshold
 * @param minLen                min query length
 * @param links                 0 if no links, otherwise should point to an array of ints.
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
 * Each item in the links array is the number of orienation types for the corresponding bam
 *
 */
int parseCoverageAndLinks(int numBams,
                          int baseQ,
                          int mapQ,
                          int minLen,
                          int * links,
                          int ignoreSuppAlignments,
                          char* coverageMode,
                          char* bamFiles[],
                          BM_mappingResults * MR);



/*!
 * @abstract work out the orientation type and insert size for given bam files
 *
 * @param  MR  mapping results struct containing BM_bamFile structs
 * @return void
 *
 * @discussion This function expects MR to be initialised.
 */
#define BM_PAIRS_FOR_TYPE 10000    // only need to parse this many pairs to infer the type of the bam file
void typeBamFiles(BM_mappingResults * MR);

/*!
 * @abstract Adjust (reduce) the number of piled-up bases along a contig
 *
 * @param  MR  mapping results struct to write to
 * @param  position_holder  array of pileup depths
 * @param  tid  contig currently being processed
 * @return void
 *
 * @discussion This function expects MR to be initialised.
 * it can change the values of contigLengthCorrectors and plpBp
 */
void adjustPlpBp(BM_mappingResults * MR,
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
float ** calculateCoverages(BM_mappingResults * MR);

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
int initLW(BM_LinkWalker * walker, BM_mappingResults * MR);

/*!
 * @abstract Move to the next LinkInfo or LinkPair
 *
 * @param  walker   pointer to LinkHolder.
 * @return 1 for step within current contig pair, 2 for new pair, 0 for end walk
 */
int stepLW(BM_LinkWalker * walker);

/*!
 * @abstract Start moving through all of the links
 *
 * @param  walker   pointer to LinkHolder.
 * @return void
 */
void destroyLW(BM_LinkWalker * walker);


    /***********************
    ***     BAMFILES     ***
    ***********************/

/*!
 * @abstract Free all the memory used to store bamFile structs
 *
 * @param  BFs      Malloc'd array of BF structs
 * @return void
 */
 void destroyBamFiles(BM_bamFile ** BFs, int numBams);


/*!
 * @abstract Human readable orientation type
 *
 * @param  type         Orientation type
 * @return char *       Human readable string
 */
char * OT2Str(OT type);

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
void print_MR(BM_mappingResults * MR);


/*!
 * @abstract Print the bamFile-ish info
 *
 * @param  BF  bam file struct to print
 * @return void
 */
void printBamFileInfo(BM_bamFile * BF);

#ifdef __cplusplus
}
#endif

#endif // BM_BAM_PARSER_H

