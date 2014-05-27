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

void makeContigKey(char* keyStore, int cid_1, int cid_2)
{
    // force cid_1 < cid_2 for consistent keys
    if(cid_1 < cid_2) {sprintf(keyStore, "%d,%d", cid_1, cid_2);}
    else {sprintf(keyStore, "%d,%d",cid_2, cid_1);}
}

void addLink(cfuhash_table_t * linkHash,
             int cid_1,
             int cid_2,
             int pos_1,
             int pos_2,
             int orient_1,
             int orient_2,
             int bam_ID
            )
{
    // store the link info, swap order of cid_1 and cid_2 if needed
    BMM_link_info* LI = (BMM_link_info*) calloc(1, sizeof(BMM_link_info));
    if(cid_1 < cid_2){
        LI->orient_1 = orient_1;
        LI->orient_2 = orient_2;
        LI->pos_1 = pos_1;
        LI->pos_2 = pos_2;
    }
    else
    {
        LI->orient_1 = orient_2;
        LI->orient_2 = orient_1;
        LI->pos_1 = pos_2;
        LI->pos_2 = pos_1;
    }
    LI->bam_ID = bam_ID;

    BMM_link_info** next_link_ptr = (BMM_link_info**) &LI->next_link;

    // see if the key is in the hash already
    char * key = calloc(30, sizeof(char)); // allocate room for the key
    makeContigKey(key, cid_1, cid_2);
    BMM_link_pair * base_LP = cfuhash_get(linkHash, key);
    if (base_LP != NULL)
    {
        // exists in the hash -> daisy chain it on
        *next_link_ptr = base_LP->LI;
        base_LP->LI = LI;
        base_LP->numLinks++;
    }
    else
    {
        // we'll need to build a bit of infrastructure
        // store the contig ids once only
        BMM_link_pair * LP = (BMM_link_pair*) calloc(1, sizeof(BMM_link_pair));
        if(cid_1 < cid_2)
        {
            LP->cid_1 = cid_1;
            LP->cid_2 = cid_2;
        }
        else
        {
            LP->cid_1 = cid_2;
            LP->cid_2 = cid_1;
        }
        LP->LI = LI;
        LP->numLinks = 1;
        *next_link_ptr = LI; // point to self means end of list

        // finally, add the lot to the hash
        cfuhash_put(linkHash, key, LP);
    }
    free(key);
}

int destroyLinkInfo_andNext(BMM_link_info** LI_ptr)
{
    BMM_link_info* next_link = (BMM_link_info*) (*LI_ptr)->next_link;
    if(*LI_ptr ==  next_link) // at the end of the chain
    {
        free(*LI_ptr);
        return 0;
    }
    else // not done yet
    {
        BMM_link_info* tmp_link = *LI_ptr;
        free(tmp_link);
        *LI_ptr = next_link;
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
        BMM_link_pair * base_LP = cfuhash_get(linkHash, keys[i]);
        if(keys[i] != 0)
            free(keys[i]);
        BMM_link_info* LI = base_LP->LI;
        while(destroyLinkInfo_andNext(&LI));
        if(base_LP !=0)
            free(base_LP);
    }
    if(keys != 0)
        free(keys);

    if(key_sizes != 0)
        free(key_sizes);
}

int initLinkWalker(BMM_LinkWalker * walker, cfuhash_table_t * linkHash)
{
    // get some memory
    //BMM_LinkWalker * walker = calloc(1, sizeof(BMM_LinkWalker));
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

int stepLinkWalker(BMM_LinkWalker * walker)
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

void destroyLinkWalker(BMM_LinkWalker * walker)
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


int getNextLinkInfo(BMM_link_info** LI_ptr)
{
    BMM_link_info* next_link = (BMM_link_info*) (*LI_ptr)->next_link;
    if(*LI_ptr ==  next_link) {return 0;}
    else
    {
        *LI_ptr = next_link;
        return 1;
    }
}

void printLinks(cfuhash_table_t * linkHash, char ** bamNames, char ** contigNames)
{
    BMM_LinkWalker walker;
    initLinkWalker(&walker, linkHash);
    int ret_val = 2;
    while(ret_val != 0)
    {
        if(ret_val == 2)
        {
            // new contig pair
            printLinkPair(walker.pair, contigNames);
        }
        printLinkInfo(walker.LI, bamNames);
        ret_val = stepLinkWalker(&walker);
    }
    destroyLinkWalker(&walker);
}

void printLinkPair(BMM_link_pair* LP, char ** contigNames)
{
    printf("  (%s, %s, %d links)\n",  contigNames[LP->cid_1], contigNames[LP->cid_2], LP->numLinks);
}

void printLinkInfo(BMM_link_info* LI, char ** bamNames)
{
    printf("    (%d,%d -> %d,%d, %s)\n",  LI->pos_1, LI->orient_1, LI->pos_2, LI->orient_2, bamNames[LI->bam_ID]);
}
