//#############################################################################
//
//   bamRead.c
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

// system includes
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <libgen.h>

// cfuhash
#include "cfuhash.h"

// local includes
#include "bamRead.h"

BM_mappedRead * makeMappedRead(char * seqId,
                                char * seq,
                                char * qual,
                                uint16_t idLen,
                                uint16_t seqLen,
                                uint16_t qualLen,
                                uint8_t rpi,
                                uint16_t group,
                                BM_mappedRead * prev_MR) {
    BM_mappedRead * MR = calloc(1, sizeof(BM_mappedRead));
    MR->seqId = strdup(seqId);
    if( 0 != seq ) {
        MR->seq = strdup(seq);
     }
     else {
         MR->seq = 0;
     }
    MR->rpi = rpi;
    MR->idLen = idLen;
    MR->seqLen = seqLen;
    MR->qualLen = qualLen;
    MR->group = group; // fprintf(stdout, "%d\n", group); fflush(stdout);
    if( 0 != qual )
        MR->qual = strdup(qual);
    else
        MR->qual = 0;

    if(prev_MR)
        prev_MR->nextRead = MR;

    MR->nextRead = 0;
    MR->partnerRead = 0;
    MR->nextPrintingRead = 0;

    return MR;
}

BM_mappedRead * getNextMappedRead(BM_mappedRead * MR) {
    return MR->nextRead;
}

BM_mappedRead * getNextPrintRead(BM_mappedRead * MR) {
    return MR->nextPrintingRead;
}

void setNextPrintRead(BM_mappedRead * baseMR, BM_mappedRead * nextMR) {
    baseMR->nextPrintingRead = nextMR;
}

BM_mappedRead * getPartner(BM_mappedRead * MR) {
    return MR->partnerRead;
}

int partnerInSameGroup(BM_mappedRead * MR) {
    if(MR->partnerRead) {
        // there is a partner
        return ((MR->partnerRead)->group == MR->group);
    }
    return 0;
}

void destroyMappedReads(BM_mappedRead * root_MR) {
    BM_mappedRead * tmp_MR = 0;
    while(root_MR) {
        if(root_MR->seqId) free(root_MR->seqId);
        if(root_MR->seq) free(root_MR->seq);
        if(root_MR->qual) free(root_MR->qual);
        tmp_MR = root_MR->nextRead;
        free(root_MR);
        root_MR = tmp_MR;
    }
}

void destroyPrintChain(BM_mappedRead * root_MR) {
    BM_mappedRead * tmp_MR = 0;
    while(root_MR) {
        if(root_MR->seqId) free(root_MR->seqId);
        if(root_MR->seq) free(root_MR->seq);
        if(root_MR->qual) free(root_MR->qual);
        tmp_MR = root_MR->nextPrintingRead;
        free(root_MR);
        root_MR = tmp_MR;
    }
}

void printMappedRead(BM_mappedRead * MR, FILE * f) {
    if(0 == f) { f = stdout; }
    if(MR->qual) {
        fprintf(f, FASTQ_FORMAT, MR->group, MR->seqId, MR->seq, MR->qual);
    }
    else {
        fprintf(f, FASTA_FORMAT, MR->group, MR->seqId, MR->seq);
    }
}

void sprintMappedRead(BM_mappedRead * MR, char * buffer, int * count) {
    if(MR->qual) {
        *count = sprintf(buffer, FASTQ_FORMAT, MR->group, MR->seqId, MR->seq, MR->qual);
    }
    else {
        *count = sprintf(buffer, FASTA_FORMAT, MR->group, MR->seqId, MR->seq);
    }
}

void printMappedReads(BM_mappedRead * root_MR, FILE * f) {

    if(0 == f) { f = stdout; }

    BM_mappedRead * MR = root_MR;
    if(MR->qual) {
        while(MR) {
            fprintf(f, FASTQ_FORMAT, MR->group, MR->seqId, MR->seq, MR->qual);
            MR = MR->nextRead;
        }
    }
    else {
        while(MR) {
            fprintf(f, FASTA_FORMAT, MR->group, MR->seqId, MR->seq);
            MR = MR->nextRead;
        }
    }
}
