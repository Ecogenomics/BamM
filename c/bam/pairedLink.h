//#############################################################################
//
//   pairedLink.h
//
//   Implements struct and methods for storing paired read links
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
typedef enum {LT_NONE, LT_SS, LT_SE, LT_ES, LT_EE, LT_ERROR} LT;   // Do contigs point start to start, end to end etc...

/*! @typedef
 @abstract Structure for storing information about specific links
 @field reversed1      orientation of read in contig 1 (1 == reversed)
 @field reversed2      orientation of read in contig 2 (1 == reversed)
 @field readLength1    length of the first read
 @field readLength2    length of the second read
 @field pos1           position of read in contig 1
 @field pos2           position of read in contig 2
 @field bid            id of the BAM file link originates from
 @field type           the type of this link
 @field nextLink       link to next link info struct (linked list)
 */
typedef struct {
    uint16_t reversed1;
    uint16_t reversed2;
    uint16_t readLength1;
    uint16_t readLength2;
    uint32_t pos1;
    uint32_t pos2;
    uint32_t bid;
    LT type;
    struct BM_linkInfo * nextLink;
} BM_linkInfo;

/*! @typedef
 * @abstract Structure for storing information about specific links
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

/*!
 * @abstract Add a link info struct to the main link table.
 *
 * @param  cid1       tid of contig 1 ( from BAM header )
 * @param  cid2       tid of contig 2
 * @param  pos1       position of read in contig 1
 * @param  pos2       position of read in contig 2
 * @param  reversed1  orientation of read in contig 1
 * @param  reversed2  orientation of read in contig 2
 * @param  bid        id of the BAM file link originates from
 * @param  type       type of this link (SE, ES, ...)
 * @return void
 *
 * @discussion The contigs can be added in any order hovever the pos and reversed1
 * variables must match this ordering. The code will sort cids accordingly.
 */
void addLink(cfuhash_table_t * linkTable,
             int cid1,
             int cid2,
             int pos1,
             int pos2,
             int reversed1,
             int reversed2,
             int bid,
             LT type);

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

/*!
 * @abstract Walk along the link info linked list
 *
 * @param  LI_ptr pointer to the current node
 * @return 1 if should keep walking, 0 is at end of list
 */
int getNextLinkInfo(BM_linkInfo** LI_ptr);

        /***********************
        *** PRINTING AND I/O ***
        ***********************/
char * LT2Str(LT type);
void printLinks(cfuhash_table_t * linkHash, char ** bamNames, char ** contigNames);
void printLinkPair(BM_linkPair* LP, char ** contigNames);
void printLinkInfo(BM_linkInfo* LI, char ** bamNames);

#ifdef __cplusplus
}
#endif

#endif // BM_PAIRED_LINK_H