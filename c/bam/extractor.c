//#############################################################################
//
//   extrator.c
//
//   An example showing how to use the c components of BamM (for extraction).
//   Also works as a handy little program too!
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
#include <string.h>
#include <stdio.h>
#include <unistd.h>

// local includes
#include "bamParser.h"
#include "pairedLink.h"

int main(int argc, char *argv[])
{
    // parse the command line
    float ex_version = 0.1;
    char * contig_list_file_name = 0;   // file containing a list of contig names
    char * out_file_prefix = 0;         // prefix for all output files
    int n, i;
    while ((n = getopt(argc, argv, "l:p:")) >= 0) {
        switch (n) {
            case 'l': contig_list_file_name = strdup(optarg); break;
            case 'p': out_file_prefix = strdup(optarg); break;
        }
    }
    if (optind == argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "[bamExtractor] version: %0.1f\n", ex_version);
        fprintf(stderr, "Usage: bamExtractor [options] in1.bam [in2.bam [...]]\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "   -l <string>             file containing list of contigs to extract reads for\n");
        fprintf(stderr, "   -p <string>             prefix for all output files\n");
        fprintf(stderr, "\n");
        return 1;
    }

    if(out_file_prefix == 0) {
        fprintf(stderr, "ERROR: you must supply a out file prefix with -p");
        return 1;
    }
    if(contig_list_file_name == 0) {
        fprintf(stderr, "ERROR: you must supply a contig list with -l");
        return 1;
    }

    // populate the list of contigNames
    int num_contigs = 0;
    // lets start with 10,000,000 contigs. If this isn't enouh we'll have
    // to expand. This is 40MB, which should be trivial to waste :(
    char ** contig_list = calloc(10000000, sizeof(char*));
    FILE *file = fopen ( contig_list_file_name, "r" );
    if (file != NULL) {
        char line [1000];
        while(fgets(line,sizeof line,file)!= NULL) {
            // for each line
            contig_list[num_contigs] = strdup(line);
            // remove newline
            size_t len = strlen(contig_list[num_contigs]);
            if( contig_list[num_contigs][len-1] == '\n' )
                contig_list[num_contigs][len-1] = 0;
            ++num_contigs;
        }

        fclose(file);
    }
    else {
        perror(contig_list_file_name); //print the error message on stderr.
    }

    int num_bams = argc - optind; // the number of BAMs on the command line

    // bam file names
    char **bam_files = calloc(num_bams, sizeof(char*));

    for (i = 0; i < num_bams; ++i) {
        bam_files[i] = strdup(argv[optind+i]);
    }

    // this is where the magick happens
    extractReads(bam_files, num_bams, contig_list, num_contigs);

    // clean up
    if(contig_list != 0) {
        for (i = 0; i < num_contigs; ++i) {
            if(contig_list[i] != 0)
                free(contig_list[i]);
        }
        free(contig_list);
    }

    free(contig_list_file_name);
    free(out_file_prefix);

    if(bam_files != 0) {
        for (i = 0; i < num_bams; ++i) {
            if(bam_files[i] != 0)
                free(bam_files[i]);
        }
        free(bam_files);
    }
    return 0;
}
