//#############################################################################
//
//   pairedLink.c
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

// system includes
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

// cfuhash
#include "cfuhash.h"

// local includes
#include "pairedLink.h"

void makeContigKey(char* keyStore, int cid1, int cid2) {
    if(cid1 < cid2) {sprintf(keyStore, "%d,%d", cid1, cid2);}
    else {sprintf(keyStore, "%d,%d",cid2, cid1);}
}

BM_linkInfo * makeLinkInfo(int cid1,
                           int cid2,
                           int pos1,
                           int pos2,
                           int reversed1,
                           int reversed2,
                           int bid
                           ) {
    // store the link info, swap order of cid1 and cid2 if needed
    BM_linkInfo* LI = (BM_linkInfo*) calloc(1, sizeof(BM_linkInfo));
    if(cid1 < cid2) {
        LI->reversed1 = reversed1;
        LI->reversed2 = reversed2;
        LI->pos1 = pos1;
        LI->pos2 = pos2;
    }
    else
    {
        LI->reversed1 = reversed2;
        LI->reversed2 = reversed1;
        LI->pos1 = pos2;
        LI->pos2 = pos1;
    }

    LI->readLength1 = 0;
    LI->readLength2 = 0;

    LI->bid = bid;
    return LI;
}

BM_linkInfo * cloneLinkInfo(BM_linkInfo * LI) {
    BM_linkInfo* new_LI = (BM_linkInfo*) calloc(1, sizeof(BM_linkInfo));
    new_LI->reversed1 = LI->reversed1;
    new_LI->reversed2 = LI->reversed2;
    new_LI->readLength1 = LI->readLength1;
    new_LI->readLength2 = LI->readLength2;
    new_LI->pos1 = LI->pos1;
    new_LI->pos2 = LI->pos2;
    new_LI->bid = LI->bid;
    return new_LI;
}


void addLink(cfuhash_table_t * linkHash,
             BM_linkInfo * LI,
             int cid1,
             int cid2
             ) {
    BM_linkInfo** nextLink_ptr = (BM_linkInfo**) &LI->nextLink;
    // see if the key is in the hash already
    char * key = calloc(30, sizeof(char)); // allocate room for the key
    makeContigKey(key, cid1, cid2);
    BM_linkPair * base_LP = cfuhash_get(linkHash, key);
    if (base_LP != NULL) {
        // exists in the hash -> daisy chain it on
        *nextLink_ptr = base_LP->LI;
        base_LP->LI = LI;
        base_LP->numLinks++;
    }
    else {
        // we'll need to build a bit of infrastructure
        // store the contig ids once only
        BM_linkPair * LP = (BM_linkPair*) calloc(1, sizeof(BM_linkPair));
        if(cid1 < cid2) {
            LP->cid1 = cid1;
            LP->cid2 = cid2;
        }
        else {
            LP->cid1 = cid2;
            LP->cid2 = cid1;
        }
        LP->LI = LI;
        LP->numLinks = 1;
        *nextLink_ptr = LI; // point to self means end of list

        // finally, add the lot to the hash
        cfuhash_put(linkHash, key, LP);
    }
    if (key != 0) {
        free(key);
        key = 0;
    }
}

int initLinkWalker(BM_LinkWalker * walker, cfuhash_table_t * linkHash) {
    // get some memory
    walker->linkHash = linkHash;
    walker->keyCount = 0;
    walker->numKeys = 0;
    size_t *key_sizes = NULL;
    // get the keys from the hash
    walker->keys = (char **)cfuhash_keys_data(linkHash,
                                              &(walker->numKeys),
                                              &key_sizes,
                                              0);
    if((walker->keys)[walker->keyCount] != 0) {
        // we don't need this
        if (key_sizes != 0) {
            free(key_sizes);
            key_sizes = 0;
        }
        // load the first linkPair
        walker->pair = cfuhash_get(linkHash, (walker->keys)[walker->keyCount]);
        if ((walker->keys)[walker->keyCount] != 0) {
            free((walker->keys)[walker->keyCount]);
            (walker->keys)[walker->keyCount] = 0;
        }
        // tee-up the first linkInfo
        walker->LI = (walker->pair)->LI;
        return 2;
    }
    else {
        // no links
        if (walker->keys != 0) {
            free(walker->keys);
            walker->keys = 0;
        }
    }
    return 0;
}

int stepLinkWalker(BM_LinkWalker * walker) {
    if(getNextLinkInfo(&(walker->LI))) {
        // still on the same pair
        return 1;
    }
    else {
        // at the end of the links of this contig pair
        walker->keyCount += 1;
        if(walker->keyCount < walker->numKeys) {
            walker->pair = cfuhash_get(walker->linkHash,
                                       walker->keys[walker->keyCount]);
            if (walker->keys[walker->keyCount] != 0) {
                free(walker->keys[walker->keyCount]);
                walker->keys[walker->keyCount] = 0;
            }
            walker->LI = (walker->pair)->LI;
            return 2;
        }
    }
    return 0;
}

int getNextLinkInfo(BM_linkInfo** LI_ptr) {
    BM_linkInfo* nextLink = (BM_linkInfo*) (*LI_ptr)->nextLink;
    if(*LI_ptr ==  nextLink) { return 0; }
    else {
        *LI_ptr = nextLink;
        return 1;
    }
}

void destroyLinkWalker(BM_LinkWalker * walker) {
    if(walker != 0) {
        if(walker->keys != 0) {
            free(walker->keys);
            walker->keys = 0;
        }
    }
}

int destroyLinkInfo_andNext(BM_linkInfo** LI_ptr) {
    BM_linkInfo* nextLink = (BM_linkInfo*) (*LI_ptr)->nextLink;
    if(*LI_ptr ==  nextLink) {
        // at the end of the chain
        if (*LI_ptr != 0 ) {
            free(*LI_ptr);
            *LI_ptr = 0;
        }
        return 0;
    }
    else {
        // not done yet
        BM_linkInfo* tmp_link = *LI_ptr;
        if (tmp_link != 0) {
            free(tmp_link);
            tmp_link = 0;
        }
        *LI_ptr = nextLink;
        return 1;
    }
}

void destroyLinks(cfuhash_table_t * linkHash) {
    char **keys = NULL;
    size_t *key_sizes = NULL;
    size_t key_count = 0;
    int i = 0;
    keys = (char **)cfuhash_keys_data(linkHash, &key_count, &key_sizes, 0);

    for (i = 0; i < (int)key_count; i++) {
        BM_linkPair * base_LP = cfuhash_get(linkHash, keys[i]);
        if(keys[i] != 0)
        {
            free(keys[i]);
            keys[i] = 0;
        }
        BM_linkInfo* LI = base_LP->LI;
        if (LI != 0)
            while(destroyLinkInfo_andNext(&LI));
        if(base_LP !=0) {
            free(base_LP);
            base_LP = 0;
        }
    }

    if(keys != 0) {
        free(keys);
        keys = 0;
    }

    if(key_sizes != 0) {
        free(key_sizes);
        key_sizes = 0;
    }
}

char * LT2Str(LT type) {
    switch(type) {
        case LT_NONE:
            return "NONE";
            break;
        case LT_SS:
            return "SS";
            break;
        case LT_SE:
            return "SE";
            break;
        case LT_ES:
            return "ES";
            break;
        case LT_EE:
            return "EE";
            break;
        case LT_ERROR:
            return "ERROR";
            break;
    }
    return "UNKNOWN";
}

void printLinks(cfuhash_table_t * linkHash,
                char ** bamNames,
                char ** contigNames) {
    BM_LinkWalker walker;
    int ret_val = initLinkWalker(&walker, linkHash);
    if(ret_val)
        printf("#cid_1\tcid_2\tpos_1\trev_1\tpos_2\trev_2\tfile\n");
        while(ret_val != 0) {
            printf("%s\t%s\t",
                   contigNames[walker.pair->cid1],
                   contigNames[walker.pair->cid2]);
            printLinkInfo(walker.LI, bamNames);
            ret_val = stepLinkWalker(&walker);
        }
        destroyLinkWalker(&walker);
}

void printLinkPair(BM_linkPair* LP, char ** contigNames) {
    printf("#%s\t%s\t%d\n",
           contigNames[LP->cid1],
           contigNames[LP->cid2],
           LP->numLinks);
}

void printLinkInfo(BM_linkInfo* LI, char ** bamNames) {
    printf("%d\t%d\t%d\t%d\t%s\n",
           LI->pos1,
           LI->reversed1,
           LI->pos2,
           LI->reversed2,
           bamNames[LI->bid]);
}
