//#############################################################################
//
//   parser.c
//
//   An example showing how to use the c components of BamM (for parsing).
//   Also works as a handy little program too!
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

// local includes
#include "bamParser.h"
#include "pairedLink.h"

int main(int argc, char *argv[])
{
    // parse the command line
    float bp_version = 0.1;
    int n = 0, baseQ = 0, mapQ = 0, min_len = 0;
    char * coverage_mode = 0;
    char * linkstr = 0;
    int * links = 0;
    while ((n = getopt(argc, argv, "q:Q:l:m:L:")) >= 0) {
        switch (n) {
            case 'l': min_len = atoi(optarg); break; // minimum query length
            case 'q': baseQ = atoi(optarg); break;   // base quality threshold
            case 'Q': mapQ = atoi(optarg); break;    // mapping quality threshold
            case 'L': linkstr = strdup(optarg); break;
            case 'm': coverage_mode = strdup(optarg); break;
        }
    }
    if (optind == argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "[bamParser] version: %0.1f\n", bp_version);
        fprintf(stderr, "Usage: bamParser [options] in1.bam [in2.bam [...]]\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "   -L <int>,[<int>,[...]]  find pairing links*\n");
        fprintf(stderr, "   -l <int>                minQLen\n");
        fprintf(stderr, "   -q <int>                base quality threshold\n");
        fprintf(stderr, "   -Q <int>                mapping quality threshold\n");
        fprintf(stderr, "   -m <string>             coverage mode [vanilla, outlier]\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "*For links: Specify the number of orientation types needed for each BAM in the\n");
        fprintf(stderr, "same order as the BAMS are specified. EX: for PE library use 1, for MP use 2\n");
        fprintf(stderr, "(this will allow 2 insert sizes to incorporate shadow library shenanigans)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "NOTE: This is a comma separated parameter.\n");
        fprintf(stderr, "\n");
        return 1;
    }

    int num_bams = argc - optind; // the number of BAMs on the command line

    // work out links!
    int i = 0;
    if(linkstr != 0) {
        // user specified links to parse. Split the commas and
        // convert to ints
        links = calloc(num_bams, sizeof(int));
        char *ch;
        ch = strtok(linkstr, ",");
        while (ch != NULL) {
            links[i] = atoi(ch);
            ++i;
            ch = strtok(NULL, " ,");
        }
        free(linkstr);
    }

    // set the default coverage_mode
    if(coverage_mode == 0)
        coverage_mode = strdup("vanilla");

    // an BFI can handle multiple BAMs but we make one BAM == one BFI
    // because we want to demonstrate the "merge" functionality here.
    BM_fileInfo ** BFIs = calloc(num_bams, sizeof(BM_fileInfo*));

    for (i = 0; i < num_bams; ++i) {
        int * l = 0;
        if(linkstr != 0) {
            l = calloc(1, sizeof(int));
            l[0] = links[i];
        }
        char **bam_files = calloc(1, sizeof(char*));             // bam file names
        bam_files[0] = strdup(argv[optind+i]);
        BM_fileInfo * bfi = calloc(1, sizeof(BM_fileInfo));
        BFIs[i] = bfi;
        int ignore_supps = 1;

        // this is where the magick happens
        parseCoverageAndLinks(0,
                              1,
                              baseQ,
                              mapQ,
                              min_len,
                              l,
                              ignore_supps,
                              coverage_mode,
                              bam_files,
                              bfi);
        free(bam_files[0]);
        free(bam_files);
    }
    if(links != 0) {
        free(links);
    }
    free(coverage_mode);

    // merge all the BAMs together
    for(i = 1; i < num_bams; ++i) {
        mergeBFIs(BFIs[0], BFIs[i]);
        BFIs[i] = 0;
    }

    // print some stuff to the screen
    printBFI(BFIs[0]);

    // clean up memory
    for(i = 0; i < num_bams; ++i) {
        destroyBFI(BFIs[i]);
    }

    return 0;
}
