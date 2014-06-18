#BamM

## Overview

BamM is a c library, wrapped in python, that parses BAM files.
The code is intended to provide a faster, more stable interface to parsing BAM files than PySam, but doesn't implement all/any of PySam's features.
Do you want all the links that join two contigs in a BAM? Do you need to get coverage? Would you like to just work out the insert size and orientation of some mapped reads?

Then BamM is for you!

## Installation

Dependencies:

The BAM parsing is done using c and a few external libraries. This slightly complicates the makes BamM installation process but fortunately not too much.

If you're installing system-wide then you can use your favourite package manager to install htslib and libcfu. For local installs, or installs that will work with the linux "modules" system, you need to be a bit trickier. This is the type of install I'll be documenting here. Read on.

The notes here are for installing on Ubuntu, but should be transferrable to any other Linux system. I'm not sure if these notes are transferabble to fashionable overpriced systems with rounded rectangles and retina displays and I'm almost cetain you'll need some sysadmin-fu to get things going on a Windows system. If you do ge it all set up then please let me know how and I'll buy you a 5 shot venti, 2/5th decaf, ristretto shot, 1pump Vanilla, 1pump Hazelnut, breve,1 sugar in the raw, with whip, carmel drizzle on top, free poured, 4 pump mocha.

First, you need git, zlib and a C-compiler. On Ubuntu this looks like:

    sudo apt-get -y install git build-essential zlib1g-dev

Next you'll need htslib (Samtools guts) and libcfu (hash objects) (if you haven't already installed them system wide)

###Notes on installing htslib:
Get the latest htslib from github:

    git clone https://github.com/samtools/htslib.git

For various resons we need to install a statically linked version of htslib. When making use this command instead of just 'make':

    make CFLAGS='-g -Wall -O2 -fPIC -static-libgcc -shared'

###Notes on installing libcfu:
I have built this librray around libcfu 0.03. It is available here:

    http://downloads.sourceforge.net/project/libcfu/libcfu/libcfu-0.03/libcfu-0.03.tar.gz

On my system I have trouble installing it becuase of inconsistencies with print statements. I fix this with sed like so, and then install locally:

    tar -xvf libcfu-0.03.tar.gz
    cd libcfu-0.03
    sed -i -e "s/%d/%zd/g" examples/*.c
    sed -i -e "s/%u/%zu/g" examples/*.c
    # Remove the '--prefix=' part to install system wide
    ./configure --prefix=`pwd`/build
    # Code in the examples folder breaks compilation so nuke this from the Makefile
    sed -i -e "s/src examples doc/src doc/" Makefile
    # make and install
    make CFLAGS='-g -Wno-unused-variable -O2 -fPIC -static-libgcc -shared'
    make install

If you install these libraries to local folders (e.g. in your home folder) then you need to take note of where you installed them. If you installed them system-wide then it *should* be no hassle.

###Install BamM
Get the latest version from github (pip hates this code for some reason...):

    git clone https://github.com/minillinim/BamM.git

If you installed htslib and libcfu system-wide then installation is very straight forward:

    python setup.py install

If you installed one or more of these libraries locally then you need to tell setup.py where they are:

    python setup.py install --with-libhts-lib /path/to/htslib --with-libhts-inc /path/to/htslib --with-libcfu-inc /path/to/libcfu/include/ --with-libcfu-lib path/to/libcfu/lib/

Relative paths are OK. You can add the --prefix flag to setup.py to install BamM locally. Once done, don't forget to add BamM to your PYTHONPATH. Also, if htslib and libcfu are in non-standard places you'll need to mess with your LD_LIBRARY_PATH.

## Example usage

    #-----------------------
    # first import it
    from bamm.BamParser import BamParser

    #-----------------------
    # use a BAM parser to get the orientation type of a BAM file only

    BP = BamParser(coverageMode="none")             # note you must set coverage to None

    BP.typeBams(['file1.bam', 'file2.bam'],         # list of file names
                types=[2,1],                        # number of insert types per BAM
                threads=2)                          # number of threads

    BP.printBamTypes()                              # tell the world!

    #-----------------------
    # use a BAM parser to get coverage profiles and links

    BP = BamParser()

    BP.parseBams(['file1.bam', 'file2.bam'],        # list of file names
                 doLinks=True,                      # work out links during parse
                 types=[2,1],                       # number of insert types per BAM
                 threads=3)                         # number of threads

    BP.printBamTypes()
    BP.printCoverages()
    BP.printLinks()

    #-----------------------
    # the full range of options at initialisation are:

     baseQuality                    # quality score threshold of reads to accept during pileup coverage calculation ( default = 0 )
     minLength                      # minimum length threshold of a mapped read ( default = 0 )
     mappingQuality                 # BWA/BAM mapping quality threshold
     coverageMode                   # 'vanilla' (standard pileup), 'outlier' (truncated mean) or 'none' (use for typing only)
     doLinks                        # calculate links during parsing ( default = False )
     ignoreSuppAlignments           # ignore supplementary alignments ( default = True )


## Using the C directly

Here is an example C program which uses this library to parse bam files (fast but not multithreaded).
see "parser.c" in the src directory

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

## Help

If you experience any problems using BamM, open an [issue](https://github.com/minillinim/BamM/issues) on GitHub and tell us about it.

## Licence and referencing

Project home page, info on the source tree, documentation, issues and how to contribute, see http://github.com/minillinim/BamM

This software is currently unpublished

## Copyright

Copyright (c) Michael Imelfort, See LICENSE.txt for further details.
