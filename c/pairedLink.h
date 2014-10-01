//#############################################################################
//
//   pairedLink.h
//
//   Implements struct and methods for storing paired read links
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

#ifndef BM_PAIRED_LINK_H
  #define BM_PAIRED_LINK_H

// system includes
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>

// cfuhash
#include "cfuhash.h"


// Link types
// Do contigs point start to start, end to end etc...
typedef enum {LT_NONE, LT_SS, LT_SE, LT_ES, LT_EE, LT_ERROR} LT;

/*! @typedef
 * @abstract Structure for storing information about specific links
 *
 * @field reversed1      orientation of read in contig 1 (1 == reversed)
 * @field reversed2      orientation of read in contig 2 (1 == reversed)
 * @field readLength1    length of the first read
 * @field readLength2    length of the second read
 * @field pos1           position of read in contig 1
 * @field pos2           position of read in contig 2
 * @field bid            id of the BAM file link originates from
 * @field nextLink       link to next link info struct (linked list)
 */
typedef struct {
    uint16_t reversed1;
    uint16_t reversed2;
    uint16_t readLength1;
    uint16_t readLength2;
    uint32_t pos1;
    uint32_t pos2;
    uint32_t bid;
    struct BM_linkInfo * nextLink;
} BM_linkInfo;

/*! @typedef
 * @abstract Structure for storing information about specific links
 *
 * @field cid1     tid of contig 1 (from BAM header)
 * @field cid2     tid of contig 2 (cid1 < cid2)
 * @field numLinks  number of link info structures in linked list
 * @field LI        first link info struct for this contig pair
 */
typedef struct {
    uint32_t cid1;
    uint32_t cid2;
    uint32_t numLinks;
    BM_linkInfo * LI;
} BM_linkPair;


/*! @typedef
 * @abstract Structure for moving through the hash of links
 *
 * @field keys      set of keys in the hash
 * @field keyCount  index into the set of keys
 * @field numKeys   number of keys in the hash
 * @param linkHash  pointer to the hash we're working with
 * @field pair      pointer to the current BM_linkPair
 * @field LI        link info pair struct ready for processing
 */
typedef struct {
    char ** keys;
    size_t keyCount;
    size_t numKeys;
    cfuhash_table_t * linkHash;
    BM_linkPair * pair;
    BM_linkInfo * LI;
} BM_LinkWalker;


#ifdef __cplusplus
extern "C" {
#endif

        /***********************
        ***   HOUSEKEEPING   ***
        ***********************/
/*!
 * @abstract Make consistent keys for contg pairs
 *
 * @param  cid1       tid of contig 1 ( from BAM header )
 * @param  cid2       tid of contig 2
 * @param  keyStore   string to write result into
 * @return void
*/
void makeContigKey(char* keyStore, int cid1, int cid2);

BM_linkInfo * makeLinkInfo(int cid1,
                           int cid2,
                           int pos1,
                           int pos2,
                           int reversed1,
                           int reversed2,
                           int bid
                           );

BM_linkInfo * cloneLinkInfo(BM_linkInfo * LI);


/*!
 * @abstract Add a link info struct to the main link table.
 *
 * @param  linkHash   hash of all the links
 * @param  LI         link to add to the link hash
 * @param  cid1       tid of contig 1 ( from BAM header )
 * @param  cid2       tid of contig 2
 * @return void
 *
 * @discussion The contigs can be added in any order hovever the pos and
 * reversed1 variables must match this ordering.
 * The code will sort cids accordingly.
 */
void addLink(cfuhash_table_t * linkHash,
             BM_linkInfo * LI,
             int cid1,
             int cid2
             );

        /***********************
        *** NAVIGATING LINKS ***
        ***********************/
/*!
 * @abstract Start moving through all of the links
 *
 * @param linkHash  pointer to hash holding all the links
 * @return 1 if links OK or 0 otherwise
 */
int initLinkWalker(BM_LinkWalker * walker, cfuhash_table_t * linkHash);

/*!
 * @abstract Move to the next LinkInfo or LinkPair
 *
 * @param  walker   pointer to LinkHolder.
 * @return 1 for step within current contig pair, 2 for new pair, 0 for end walk
 */
int stepLinkWalker(BM_LinkWalker * walker);

/*!
 * @abstract Walk along the link info linked list
 *
 * @param  LI_ptr pointer to the current node
 * @return 1 if should keep walking, 0 is at end of list
 */
int getNextLinkInfo(BM_linkInfo** LI_ptr);

        /***********************
        *** HATIN' SEGFAULTS ***
        ***********************/
/*!
 * @abstract Start moving through all of the links
 *
 * @param  walker   pointer to LinkHolder.
 * @return void
 */
void destroyLinkWalker(BM_LinkWalker * walker);

/*!
 * @abstract Walk along the link info linked list destoying the current node
 *
 * @param  LI_ptr pointer to the current node
 * @return 1 if should keep walking, 0 is at end of list
 */
int destroyLinkInfo_andNext(BM_linkInfo** LI_ptr);

/*!
 * @abstract Destroy all links information
 *
 * @param  linkHash pointer to the hash to be desroyed
 * @return void
 */
void destroyLinks(cfuhash_table_t * linkHash);

        /***********************
        *** PRINTING AND I/O ***
        ***********************/
/*!
 * @abstract Convert a link type into human readable string
 *
 * @param  type     type of the link in question
 * @return char*    human readable version of the link tyoe
 */
char * LT2Str(LT type);

/*!
 * @abstract Prints information about contigs and number / source of links
 *
 * @param  contigNames      array of char * contig names
 * @param  bamNames         array of char * bam names
 * @param  linkHash         hash containing all the links
 * @return void
 */
void printLinks(cfuhash_table_t * linkHash,
                char ** bamNames,
                char ** contigNames);

/*!
 * @abstract Print contigs in pair and number of links joining them
 *
 * @param  LP               pair of contigs to print links for
 * @param  contigNames      array of char * contig names
 * @return void
 */
void printLinkPair(BM_linkPair* LP, char ** contigNames);

/*!
 * @abstract Prints information about contigs and number / source of links
 *
 * @param  LI               LinkInfo struct to print
 * @param  bamNames         array of char * bam names
 * @return void
 */
void printLinkInfo(BM_linkInfo* LI, char ** bamNames);

#ifdef __cplusplus
}
#endif

#endif // BM_PAIRED_LINK_H