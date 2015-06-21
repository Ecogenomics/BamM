//#############################################################################
//
//   bamParser.h
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

#ifndef BM_BAM_PARSER_H
  #define BM_BAM_PARSER_H

// system includes
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>

// htslib
#include "bgzf.h"
#include "sam.h"

// cfuhash
#include "cfuhash.h"

// local includes
#include "pairedLink.h"
#include "coverageEstimators.h"

typedef BGZF bamFile;

#ifdef __cplusplus
extern "C" {
#endif

// proper linking read is a properly paired
#define BM_BAM_FMAPPED (BAM_FMUNMAP | BAM_FUNMAP)

// relative orientation of paired reads
typedef enum {OT_OUT, OT_SAME, OT_IN, OT_NONE, OT_ERROR} OT;

/*! @typedef
 * @abstract Auxiliary data structure used in read_bam
 *
 * @field fp the file handler
 * @field iter NULL if a region not specified
 * @field min_mapQ mapQ filter
 * @field min_len length filter
 */
typedef struct {            //
    bamFile *fp;            // the file handler
    hts_itr_t *iter;        // NULL if a region not specified
    int min_mapQ, min_len;  // mapQ filter; length filter
} aux_t;

/*! @typedef
 * @abstract Structure for storing information about the orientation type
 *           and insert sizes for a given bam file
 *
 * @field orientationType     orientation of the reads in this mapping
 * @field insertSize          insert size of the reads (outer)
 * @field insertStdev         standard deviation measured for the insert size
 * @field numTested           number of reads tested to get these stats
 */
 typedef struct BM_bamType {
    int orientationType;                // actually going to use the ENUM,
                                        // but worried about how to get around
                                        // this with python ctypes
    float insertSize;
    float insertStdev;
    int supporting;
 } BM_bamType;


/*! @typedef
 * @abstract Structure for storing information about a single BAM file
 *
 * @field fileName            filename of the BAM
 * @field fileNameLength      length of the filename
 * @field types               the orientation types and insert sizes
 *                            accociated with the bam file
 * @field numTypes            the number of orientation types for the bam file
 */
 typedef struct BM_bamFile {
    char * fileName;
    uint16_t fileNameLength;
    BM_bamType ** types;
    int numTypes;
 } BM_bamFile;

/*! @typedef
 * @abstract Structure for storing BAM file information
 *
 * @field coverages              calculated coverage of each contig per bam
 * @field contigLengths          lengths of the referene sequences
 * @field numBams                number of BAM files parsed
 * @field numContigs             number of reference sequences
 * @field bamFiles               array of pointers to BM_bamFile's
 * @field contigNames            names of the reference sequences
 * @field contigNameLengths      lengths reference sequence names
 * @field isLinks                are links being calculated
 * @field coverageType           coverage type (see coverageEstimators.h)
 * @field isIgnoreSupps          are supplementary alignments being ignored
 * @field links                  linking pairs
 */
typedef struct BM_fileInfo {
    float ** coverages;
    uint32_t * contigLengths;
    uint32_t numBams;
    uint32_t numContigs;
    BM_bamFile ** bamFiles;
    char ** contigNames;
    uint16_t * contigNameLengths;
    int isLinks;
    BM_coverageType * coverageType;
    int isIgnoreSupps;
    cfuhash_table_t * links;
} BM_fileInfo;

    /**********************
    ***  BAM FILE DATA  ***
    **********************/
/*!
 * @abstract Allocate space for a new BFI struct
 *
 * @return BM_fileInfo *
 *
 * @discussion Allocates the memory which should be initialised
 * You call destroyBFI to free this memory
 */
BM_fileInfo * createBFI(void);

/*!
 * @abstract Initialise the mapping results struct
 *
 * @param BFI                   mapping results struct to initialise
 * @param BAM_header            htslib BAM header
 * @param numBams               number of BAM files to parse
 * @param bamFiles              filenames of BAMs parsed
 * @param types                 array of num orientation types per bam
 * @param isLinks               == 1 -> links will be calculated (alloc storage)
 * @param coverageType          type of coverage to be calculated
 *                              ("none" if we're just extracting reads)
 * @param ignoreSuppAlignments  only use primary alignments
 * @return void
 *
 * @discussion If you call this function then you MUST call destroyBFI
 * when you're done.
 */
void initBFI(BM_fileInfo * BFI,
             bam_hdr_t * BAM_header,
             int numBams,
             char * bamFiles[],
             int * types,
             int isLinks,
             BM_coverageType * coverageType,
             int ignoreSuppAlignments
);

/*!
 * @abstract Merge the contents of BFI_B into BFI_A
 *
 * @param  BFI_A  BM_bamFileInfo struct to copy to
 * @param  BFI_B  BM_bamFileInfo struct to copy from
 * @return void
 *
 * @discussion BFI_B remains unchanged.
 * BFI_A is updated to include all the info contained in BFI_B
 *
 * NOTE: We assume that all the headers of all the files are in sync. If they
 * are not, if contigs have been removed etc, then doom will swiftly follow.
 * Also, we assume that flags like do_links, do_outlier match... ...or DOOM!
 *
 * BFI_B is deleted in the course of the merge
 */
void mergeBFIs(BM_fileInfo * BFI_A, BM_fileInfo * BFI_B);

    /***********************
    *** PARSING BAMFILES ***
    ***********************/
/*!
 * @abstract Auxillary function needed as helper for libhts
*/
int read_bam(void *data,
             bam1_t *b);

/*!
 * @abstract Initialise the mapping results struct <- read in the BAM files
 *
 * @param  doLinks              == 1 -> calculate linking pairs
 * @param  doCovs               == 1 -> calculate coverage profiles
 * @param  numBams              number of BAM files to parse
 * @param  baseQ                base quality threshold
 * @param  mapQ                 mapping quality threshold
 * @param  minLen               min query length
 * @param  maxMisMatches        maximum number of mismatches to accept (NM flag)
 * @param  types                array of num orientation types per bam
 * @param  ignoreSuppAlignments == 1 -> ignore supplmentary alignments
 * @param  ignoreSecondaryAlignments  == 1 -> ignore secondary alignments
 * @param  coverageType         BM_CoverageType (cov type to be calculated)
 * @param  bamFiles             filenames of BAM files to parse
 * @param  BFI                  BM_bamFileInfo struct to write to
 * @return 0 for success
 *
 * @discussion This function expects BFI to be a null pointer. It calls
 * initBFI and stores info accordingly. TL;DR If you call this function
 * then you MUST call destroyBFI when you're done.
 *
 * Each item in the links array is the number of orienation types for the bam
 *
 */
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
                          BM_fileInfo * BFI);

/*!
 * @abstract work out the orientation type and insert size for given bam files
 *
 * @param  BFI         BM_bamFileInfo struct containing BM_bamFile structs
 * @return void
 *
 * @discussion This function expects BFI to be initialised.
 */
// only need to parse this many pairs to infer the type of the bam file
#define BM_PAIRS_FOR_TYPE 10000
// ignore reads with centers mapping this far from end of contig
#define BM_IGNORE_FROM_END 3000
void typeBamFiles(BM_fileInfo * BFI);

    /***********************
    ***      LINKS       ***
    ***********************/
/*!
 * @abstract Start moving through all of the links
 *
 * @param  BFI         BM_bamFileInfo struct containing links
 * @return 1 if links OK or 0 otherwise
 */
int initLW(BM_LinkWalker * walker, BM_fileInfo * BFI);

/*!
 * @abstract Move to the next LinkInfo or LinkPair
 *
 * @param  walker   pointer to LinkHolder.
 * @return 1 for step within current contig pair, 2 for new pair, 0 for end walk
 */
int stepLW(BM_LinkWalker * walker);

    /***********************
    *** HATIN' SEGFAULTS ***
    ***********************/
/*!
 * @abstract Free all the memory calloced in initBFI
 *
 * @param  BFI         BM_bamFileInfo struct to destroy
 * @return void
 */
void destroyBFI(BM_fileInfo * BFI);

/*!
 * @abstract Destroy the link walker instance
 *
 * @param  walker   pointer to LinkHolder.
 * @return void
 */
void destroyLW(BM_LinkWalker * walker);

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
 * @param  BFI         BM_bamFileInfo struct to print
 * @return void
 */
void printBFI(BM_fileInfo * BFI);

/*!
 * @abstract Print the bamFile-ish info
 *
 * @param  BF  bam file struct to print
 * @return void
 */
void printBamFileType(BM_bamFile * BF);

#ifdef __cplusplus
}
#endif

#endif // BM_BAM_PARSER_H

