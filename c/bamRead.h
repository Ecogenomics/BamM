//#############################################################################
//
//   bamRead.h
//
//   Data structures for representing mapped reads in a BAM file
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

#ifndef BM_BAM_READ_H
#define BM_BAM_READ_H

// system includes
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>

// cfuhash
#include "cfuhash.h"

// local includes

#ifdef __cplusplus
extern "C" {
#endif

// style for printing headers, fasta and fastq sequences
#define FASTA_FORMAT ">g_%s;%s%s\n%s\n"
#define FASTQ_FORMAT "@g_%s;%s%s\n%s\n+\n%s\n"
#define HEADER_FORMAT "g_%s;%s%s\n"

/*
 *   Read pairing info
 *
 *   ERROR       just for fun
 *   FIR         means first in properly paired mapping
 *   SEC         means second "
 *   SNGL_FIR    means first paired but unpaired in mapping
 *   SNGL_SEC    means second paired but unpaired in mapping
 *   SNGL        means cannonical unpaired read
 */
typedef enum {RPI_ERROR,
              RPI_FIR,
              RPI_SEC,
              RPI_SNGL_FIR,
              RPI_SNGL_SEC,
              RPI_SNGL} RPI;
/*
 *   Extraction info
 *
 *   MI_ER_EM_EG    Unset
 *   MI_PR_PM_PG    Paired read, paired mapping, paired grouping
 *   MI_PR_PM_UG    Paired read, paired mapping, unpaired grouping
 *   MI_PR_UM_NG    Paired read, unpaired mapping
 *   MI_UR_NM_NG    Unpaired read
 */
typedef enum {MI_ER_EM_EG,
              MI_PR_PM_PG,
              MI_PR_PM_UG,
              MI_PR_UM_NG,
              MI_UR_NM_NG} MI;
/*
 * Convert mapping information to a readable string
*/
static const char MITEXT[5][12];

/*! @typedef
 * @abstract Structure for storing a mapped read (linked list style)
 *
 * @field seqId             id of the read as stored in the BAM
 * @field seq               nucleic acid sequence
 * @field qual              quality sequence
 * @field idLen             length of the sequence ID
 * @field seqLen            length of the sequence
 * @field qualLen           length of the quality sequence
 * @field rpi               identifies which read in a pair or unpaired
 * @field rpi               encodes the pairing infor for the read header
 * @field group             target group this read belongs to
 * @field nextRead          linked list structure (next in chain)
 * @field partnerRead       the pair of this read (if paired and if read 1)
 * @field nextPrintingRead  the next read that will be printed
 */
typedef struct BM_mappedRead {
    char * seqId;
    char * seq;
    char * qual;
    uint16_t idLen;
    uint16_t seqLen;
    uint16_t qualLen;
    uint8_t rpi;
    uint8_t mi;
    uint16_t group;
    struct BM_mappedRead * nextRead;
    struct BM_mappedRead * partnerRead;
    struct BM_mappedRead * nextPrintingRead;
} BM_mappedRead;

/*!
 * @abstract Initialise a BM_mappedRead struct
 *
 * @param   seqId       Id of the read as given in the BAM
 * @param   seq         read sequence
 * @param   qual        quality score of the read
 * @param   idLen       length of the sequence ID
 * @param   seqLen      length of the sequence
 * @param   qualLen     length of the quality sequence
 * @param   rpi         read pairing info
 * @param   group       target group ID to assign to this read
 * @param   prevRead    pointer to the previously made mapped read
 * @return  BM_mappedRead *
 */
 BM_mappedRead * makeMappedRead(char * seqId,
                                char * seq,
                                char * qual,
                                uint16_t idLen,
                                uint16_t seqLen,
                                uint16_t qualLen,
                                uint8_t rpi,
                                uint16_t group,
                                BM_mappedRead * prev_MR);

/*!
 * @abstract Set the MI code for a BM_mappedRead struct
 *
 * @param  BM_mappedRead *  the read to update
 * @param  code             the code to set
 * @return void
 */
 void setMICode(BM_mappedRead * BMM, MI code);

/*!
 * @abstract Get the next BM_mappedRead struct in the base chain
 *
 * @param  MR               current BM_mappedRead struct in chain
 * @return BM_mappedRead*   next BM_mappedRead
 */
BM_mappedRead * getNextMappedRead(BM_mappedRead * MR);

/*!
 * @abstract Get the next BM_mappedRead struct in the print chain
 *
 * @param  MR    current BM_mappedRead struct in print chain
 * @return BM_mappedRead*   next BM_mappedRead
 */
BM_mappedRead * getNextPrintRead(BM_mappedRead * MR);

/*!
 * @abstract Set the next BM_mappedRead struct in the print chain
 *
 * @param  baseMR    current BM_mappedRead struct in chain
 * @param  nextMR    current BM_mappedRead struct in chain
 * @return void
 */
void setNextPrintRead(BM_mappedRead * baseMR, BM_mappedRead * nextMR);

/*!
 * @abstract Get the paired partner BM_mappedRead
 *
 * @param  MR               query BM_mappedRead
 * @return BM_mappedRead*   the paired partner
 */
BM_mappedRead * getPartner(BM_mappedRead * MR);

/*!
 * @abstract Is the partner of this read in the same target group?
 *
 * @param  MR    current BM_mappedRead struct in print chain
 * @return int   0 -> different group, 1 -> same group
 */
int partnerInSameGroup(BM_mappedRead * MR);

/*!
 * @abstract Delete a chain of BM_mappedRead structs
 *
 * @param  root_MR    first BM_mappedRead struct in chain to delete
 * @return void
 */
void destroyMappedReads(BM_mappedRead * root_MR);

/*!
 * @abstract Delete a print chain of BM_mappedRead structs
 *
 * @param  root_MR    first BM_mappedRead struct in chain to delete
 * @return void
 */
void destroyPrintChain(BM_mappedRead * root_MR);

/*!
 * @abstract Pretty print a BM_mappedRead struct
 *
 * @param  MR          BM_mappedRead struct to print
 * @param  f           file handle to print to or 0 for stdout
 * @param  groupName    string representing the mapping target's name
 * @param  headerOnly  1 if should print only read header
 * @param  pairedOutput 1 if printing to a paired output file
 * @return void
 */
void printMappedRead(BM_mappedRead * MR,
                     FILE * f,
                     char * groupName,
                     int headerOnly,
                     int pairedOutput);

/*!
 * @abstract Pretty print a BM_mappedRead struct to a string
 *
 * @param  MR           BM_mappedRead struct to print
 * @param  buffer       character buffer to write to
 * @param  count        number of characters written into buffer
 * @param  groupName    string representing the mapping target's name
 * @param  headerOnly   1 if should print only read headers
 * @param  pairedOutput 1 if printing to a paired output file
 * @return void
 */
void sprintMappedRead(BM_mappedRead * MR,
                      char * buffer,
                      int * count,
                      char * groupName,
                      int headerOnly,
                      int pairedOutput);

/*!
 * @abstract Pretty print BM_mappedRead structs
 *
 * @param  root_MR      first BM_mappedRead struct to print
 * @param  f            file handle to print to or 0 for stdout
 * @param  groupNames   [ string ] representing the mapping targets' name
 * @param  headersOnly  1 if should print only read headers
 * @param  pairedOutput 1 if printing to a paired output file
 * @return void
 */
void printMappedReads(BM_mappedRead * root_MR,
                      FILE * f,
                      char ** groupNames,
                      int headersOnly,
                      int pairedOutput);

#ifdef __cplusplus
}
#endif

#endif // BM_BAM_READ_H
