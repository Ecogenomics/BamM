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

#ifndef BMM_PAIRED_LINK_H
  #define BMM_PAIRED_LINK_H

// system includes
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>

// cfuhash
#include "cfuhash.h"

/*! @typedef
 @abstract Structure for storing information about specific links
 @field orient_1    orientation of read in contig 1 (1 == reversed)
 @field orient_2    orientation of read in contig 2 (1 == reversed)
 @field pos_1       position of read in contig 1
 @field pos_2       position of read in contig 2
 @field bam_ID      id of the BAM file link originates from
 @field next_link   link to next link info struct (linked list)
 */
typedef struct {
    uint16_t orient_1;
    uint16_t orient_2;
    uint32_t pos_1;
    uint32_t pos_2;
    uint32_t bam_ID;
    struct BMM_link_info * next_link;
} BMM_link_info;

/*! @typedef
 * @abstract Structure for storing information about specific links
 * @field cid_1     tid of contig 1 (from BAM header)
 * @field cid_2     tid of contig 2 (cid_1 < cid_2)
 * @field numLinks  number of link info structures in linked list
 * @field LI        first link info struct for this contig pair
 */
typedef struct {
    uint32_t cid_1;
    uint32_t cid_2;
    uint32_t numLinks;
    BMM_link_info * LI;
} BMM_link_pair;


/*! @typedef
 * @abstract Structure for moving through the hash of links
 * @field keys      set of keys in the hash
 * @field keyCount  index into the set of keys
 * @field numKeys   number of keys in the hash
 * @param linkHash  pointer to the hash we're working with
 * @field pair      pointer to the current BMM_link_pair
 * @field LI        link info pair struct ready for processing
 */
typedef struct {
    char ** keys;
    size_t keyCount;
    size_t numKeys;
    cfuhash_table_t * linkHash;
    BMM_link_pair * pair;
    BMM_link_info * LI;
} BMM_LinkWalker;


#ifdef __cplusplus
extern "C" {
#endif

/*!
 * @abstract Add a link info struct to the main link table.
 *
 * @param  cid_1     tid of contig 1 ( from BAM header )
 * @param  cid_2     tid of contig 2
 * @param  pos_1     position of read in contig 1
 * @param  pos_2     position of read in contig 2
 * @param  orient_1  orientation of read in contig 1
 * @param  orient_2  orientation of read in contig 2
 * @param  bam_ID    id of the BAM file link originates from
 * @return void
 *
 * @discussion The contigs can be added in any order hovever the pos and orient_1
 * variables must match this ordering. The code will sort cids accordingly.
 */
void addLink(cfuhash_table_t * linkTable,
             int cid_1,
             int cid_2,
             int pos_1,
             int pos_2,
             int orient_1,
             int orient_2,
             int bam_ID);

/*!
 * @abstract Start moving through all of the links
 *
 * @param linkHash  pointer to hash holding all the links
 * @return 1 if links OK or 0 otherwise
 */
int initLinkWalker(BMM_LinkWalker * walker, cfuhash_table_t * linkHash);

/*!
 * @abstract Move to the next LinkInfo or LinkPair
 *
 * @param  walker   pointer to LinkHolder.
 * @return 1 for step within current contig pair, 2 for new pair, 0 for end walk
 */
int stepLinkWalker(BMM_LinkWalker * walker);

/*!
 * @abstract Start moving through all of the links
 *
 * @param  walker   pointer to LinkHolder.
 * @return void
 */
void destroyLinkWalker(BMM_LinkWalker * walker);

/*!
 * @abstract Walk along the link info linked list destoying the current node
 *
 * @param  LI_ptr pointer to the current node
 * @return 1 if should keep walking, 0 is at end of list
 */
int destroyLinkInfo_andNext(BMM_link_info** LI_ptr);

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
int getNextLinkInfo(BMM_link_info** LI_ptr);

        /***********************
        *** PRINTING AND I/O ***
        ***********************/

void printLinks(cfuhash_table_t * linkHash, char ** bamNames, char ** contigNames);
void printLinkPair(BMM_link_pair* LP, char ** contigNames);
void printLinkInfo(BMM_link_info* LI, char ** bamNames);

#ifdef __cplusplus
}
#endif

#endif // BMM_PAIRED_LINK_H