//#############################################################################
//
//   pairedLink.c
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

// system includes
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

// cfuhash
#include "cfuhash.h"

// local includes
#include "pairedLink.h"

void makeContigKey(char* keyStore, int cid1, int cid2)
{
    // force cid1 < cid2 for consistent keys
    if(cid1 < cid2) {sprintf(keyStore, "%d,%d", cid1, cid2);}
    else {sprintf(keyStore, "%d,%d",cid2, cid1);}
}

void addLink(cfuhash_table_t * linkHash,
             int cid1,
             int cid2,
             int pos1,
             int pos2,
             int reversed1,
             int reversed2,
             int bid,
             LT type
            )
{
    // store the link info, swap order of cid1 and cid2 if needed
    BM_linkInfo* LI = (BM_linkInfo*) calloc(1, sizeof(BM_linkInfo));
    if(cid1 < cid2){
        LI->reversed1 = reversed1;
        LI->reversed2 = reversed2;
        LI->pos1 = pos1;
        LI->pos2 = pos2;
        LI->type = type;
    }
    else
    {
        LI->reversed1 = reversed2;
        LI->reversed2 = reversed1;
        LI->pos1 = pos2;
        LI->pos2 = pos1;
        // take care of the types if we're swapping...
        if(type == LT_SE) { LI->type = LT_ES; }
        else if(type == LT_ES) { LI->type = LT_SE; }
        else { LI->type = type; }
    }
    LI->bid = bid;

    BM_linkInfo** nextLink_ptr = (BM_linkInfo**) &LI->nextLink;

    // see if the key is in the hash already
    char * key = calloc(30, sizeof(char)); // allocate room for the key
    makeContigKey(key, cid1, cid2);
    BM_linkPair * base_LP = cfuhash_get(linkHash, key);
    if (base_LP != NULL)
    {
        // exists in the hash -> daisy chain it on
        *nextLink_ptr = base_LP->LI;
        base_LP->LI = LI;
        base_LP->numLinks++;
    }
    else
    {
        // we'll need to build a bit of infrastructure
        // store the contig ids once only
        BM_linkPair * LP = (BM_linkPair*) calloc(1, sizeof(BM_linkPair));
        if(cid1 < cid2)
        {
            LP->cid1 = cid1;
            LP->cid2 = cid2;
        }
        else
        {
            LP->cid1 = cid2;
            LP->cid2 = cid1;
        }
        LP->LI = LI;
        LP->numLinks = 1;
        *nextLink_ptr = LI; // point to self means end of list

        // finally, add the lot to the hash
        cfuhash_put(linkHash, key, LP);
    }
    free(key);
}

int destroyLinkInfo_andNext(BM_linkInfo** LI_ptr)
{
    BM_linkInfo* nextLink = (BM_linkInfo*) (*LI_ptr)->nextLink;
    if(*LI_ptr ==  nextLink) // at the end of the chain
    {
        free(*LI_ptr);
        return 0;
    }
    else // not done yet
    {
        BM_linkInfo* tmp_link = *LI_ptr;
        free(tmp_link);
        *LI_ptr = nextLink;
        return 1;
    }
}

void destroyLinks(cfuhash_table_t * linkHash)
{
    char **keys = NULL;
    size_t *key_sizes = NULL;
    size_t key_count = 0;
    int i = 0;
    keys = (char **)cfuhash_keys_data(linkHash, &key_count, &key_sizes, 0);

    for (i = 0; i < (int)key_count; i++) {
        BM_linkPair * base_LP = cfuhash_get(linkHash, keys[i]);
        if(keys[i] != 0)
            free(keys[i]);
        BM_linkInfo* LI = base_LP->LI;
        while(destroyLinkInfo_andNext(&LI));
        if(base_LP !=0)
            free(base_LP);
    }
    if(keys != 0)
        free(keys);

    if(key_sizes != 0)
        free(key_sizes);
}

int initLinkWalker(BM_LinkWalker * walker, cfuhash_table_t * linkHash)
{
    // get some memory
    //BM_LinkWalker * walker = calloc(1, sizeof(BM_LinkWalker));
    walker->linkHash = linkHash;
    walker->keyCount = 0;
    walker->numKeys = 0;
    size_t *key_sizes = NULL;
    // get the keys from the hash
    walker->keys = (char **)cfuhash_keys_data(linkHash, &(walker->numKeys), &key_sizes, 0);
    if((walker->keys)[walker->keyCount] != 0)
    {
        // we don't need this
        free(key_sizes);
        // load the first linkPair
        walker->pair = cfuhash_get(linkHash, (walker->keys)[walker->keyCount]);
        free((walker->keys)[walker->keyCount]);
        // tee-up the first linkInfo
        walker->LI = (walker->pair)->LI;
        return 1;
    }
    else
    {
        // no links
        free(walker->keys);
    }
    return 0;
}

int stepLinkWalker(BM_LinkWalker * walker)
{
    if(getNextLinkInfo(&(walker->LI)))
    {
        // still on the same pair
        return 1;
    }
    else
    {
        // at the end of the links of this contig pair
        walker->keyCount += 1;
        if(walker->keyCount < walker->numKeys)
        {
            walker->pair = cfuhash_get(walker->linkHash, walker->keys[walker->keyCount]);
            free(walker->keys[walker->keyCount]);
            walker->LI = (walker->pair)->LI;
            return 2;
        }
    }
    return 0;
}

void destroyLinkWalker(BM_LinkWalker * walker)
{
    if(walker != 0)
    {
        if(walker->keys != 0)
        {
            free(walker->keys);
        }
        //free(walker);
    }
}


int getNextLinkInfo(BM_linkInfo** LI_ptr)
{
    BM_linkInfo* nextLink = (BM_linkInfo*) (*LI_ptr)->nextLink;
    if(*LI_ptr ==  nextLink) {return 0;}
    else
    {
        *LI_ptr = nextLink;
        return 1;
    }
}

char * LT2Str(LT type)
{
    switch(type)
    {
        case LT_NONE:
            return "UNSET";
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


void printLinks(cfuhash_table_t * linkHash, char ** bamNames, char ** contigNames)
{
    BM_LinkWalker walker;
    initLinkWalker(&walker, linkHash);
    int ret_val = 2;
    printf("#cid_1\tcid_2\tpos_1\trev_1\tpos_2\trev_2\ttype\tfile\n");
    while(ret_val != 0)
    {
        /*
        if(ret_val == 2)
        {
            // new contig pair
            printLinkPair(walker.pair, contigNames);
        }
        */
        printf("%s\t%s\t", contigNames[walker.pair->cid1], contigNames[walker.pair->cid2]);
        printLinkInfo(walker.LI, bamNames);
        ret_val = stepLinkWalker(&walker);
    }
    destroyLinkWalker(&walker);
}

void printLinkPair(BM_linkPair* LP, char ** contigNames)
{
    printf("#%s\t%s\t%d\n",  contigNames[LP->cid1], contigNames[LP->cid2], LP->numLinks);
}

void printLinkInfo(BM_linkInfo* LI, char ** bamNames)
{
    printf("%d\t%d\t%d\t%d\t%s\t%s\n",  LI->pos1, LI->reversed1, LI->pos2, LI->reversed2, LT2Str(LI->type), bamNames[LI->bid]);
}
