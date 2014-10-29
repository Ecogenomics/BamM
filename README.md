#BamM

## Overview

BamM is a c library, wrapped in python, that parses BAM files.
The code is intended to provide a faster, more stable interface to parsing BAM files than PySam, but doesn't implement all/any of PySam's features.
Do you want all the links that join two contigs in a BAM? Do you need to get coverage? Would you like to just work out the insert size and orientation of some mapped reads?

Then BamM is for you!

## Installation

Dependencies:

The BAM parsing is done using c and a few external libraries. This slightly complicates the BamM installation process but fortunately not too much.

If you're running 'bamm make' you'll need to have bwa and samtools installed. Installation of these tools is really straightforward. You can find the code and instructions at:

Samtools:   https://github.com/samtools/samtools
BWA:        https://github.com/lh3/bwa

If you're installing system-wide then you can use your favourite package manager to install htslib and libcfu. For local installs, or installs that will work with the linux "modules" system, you need to be a bit trickier. This is the type of install I'll be documenting here. Read on.

The notes here are for installing on Ubuntu, but should be transferrable to any other Linux system. I'm not sure if these notes are transferabble to fashionable overpriced systems with rounded rectangles and retina displays and I'm almost cetain you'll need some sysadmin-fu to get things going on a Windows system. If you do ge it all set up then please let me know how and I'll buy you a 5 shot venti, 2/5th decaf, ristretto shot, 1 pump Vanilla, 1 pump Hazelnut, breve, 1 sugar in the raw, with whip, carmel drizzle on top, free poured, 4 pump mocha.

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

Relative paths are OK. You can add the --prefix flag to setup.py to install BamM locally. Once done, don't forget to add BamM to your PYTHONPATH.
Also, if htslib and libcfu are in non-standard places and you plan to access the C code, you'll need to mess with your LD_LIBRARY_PATH.

## Example usage

I've wrapped the python in a script / library called bamm.

## Using BamM on the command line

BamM has 3 modes; 'make', 'parse' and 'extract'. The first option allows you to make BAM files. The second option lets you derive coverage profiles or linking information. The final option lets you extract reads that map to a set(s) of contigs.

#-----------------------
# EX: Map several reads files to a common reference

    $ bamm make -d my_assembly.fa -i interleaved_1.fastq.gz interleaved_2.fastq.gz -c paired_R1.fastq.gz paired_R2.fastq.gz -s unpaired.fastq.gz [-t 20] [-v]

This command will make 4 BAM files by mapping the two interleaved readsets, the one paired readset and the singleton read set onto the reference m_assembly.fa
The code calls BWA and samtools to produce a set of sorted and indexed BAM files. If you specify -t <threads> then bamm will pass this onto BWA and samtools.
Use -v to get more verbose output.

NOTE: To save space, the final BAM files contain only mapped reads.
NOTE: Output files are automatically named based on the names of the read files, however you can specify a prefix to append to the beginning of all output files.

#-----------------------
# EX: Find the insert size and read orientation associated with a BAM file

    $ bamm parse file.bam

This produces output like this:

    #file   insert  stdev   orientation supporting
    file.bam    899.7514    14.7167 IN  10000

The 'IN' orientation indicates that this is a paired-end (PE) library with an insert of ~900 bp and a standard deviation of ~15 bp.
Mate pair (MP) libraries will typically have orientation 'OUT'.

Many MP libraries also have a shadow library which looks like someone added some PE reads to the mix. You can tell BamM to look for more than one insert type
by specifying the -n option

    $ bamm parse -n 2 mate_pair_file.bam

    #file   insert  stdev   orientation supporting
    mate_pair_file.bam    2524.4540    729.1291 OUT  10000
    mate_pair_file.bam    251.6253    44.7241 IN  10000

Multiple BAM files are separated using spaces. The -n argument is space separated too.

    $ bamm parse pe_file.bam mp_file.bam -n 1 2

#-----------------------
# EX: Create a coverage profile or read-linking information from several BAM files

    $ bamm parse -c coverage.tsv -m <COV_MODE> f1.bam f2.bam f3.bam [-t 3]

Produces this output in the file 'coverage.tsv'

    #contig         Length  f1.bam          f2.bam          f3.bam
    contig_1        946     103.0000        327.0000        369.0000
    contig_3        1147    130.0000        492.0000        778.0000
    contig_5        1465    228.0000        643.0000        970.0000
    contig_7        168     34.0000         82.0000         102.0000
    contig_9        4045    899.0000        1756.0000       2649.0000

The -t option indicates the maximum number of threads bamm will use. This option speeds up the process but you should only use as may threads as you have bam files. If you have 6 bam files then you'll see no improvement when using -t 6, -t 7 or -t 700.
BamM implements several coverage calculation methods:

    Mode        Description

    counts      Total number of reads that START mapping on each contig (see above output).
    cmean       Number of reads that START mapping on each contig divided by contig length.
    pmean       Average pileup depth along the length of the contig.
    pmedian     Median pileup depth along the length of the contig.
    tpmean      Average pileup depth along the length of the contig after trimming upper and lower (-r) percent of pileup depths.
    opmean      Average pileup depth along the length of the contig after trimming upper and lower (-r) stdevs of pileup depths.

BamM will calculate links between contigs if passed the '-l <filename>' argument.

    $ bamm parse -l linke.tsv f1.bam f2.bam f3.bam

Produces this output in the file links.tsv:

    #cid_1          cid_2           len_1   pos_1   rev_1   len_2   pos_2   rev_2   file
    contig_2203     contig_3479     1664    334     0       3873    2866    0       f2.bam
    contig_2203     contig_3479     1664    384     0       3873    2818    0       f2.bam
    contig_2203     contig_3479     1664    383     0       3873    2831    0       f2.bam
    contig_2203     contig_3479     1664    349     0       3873    2864    0       f2.bam
    contig_2203     contig_3479     1664    338     0       3873    2862    0       f2.bam

The first (non-header) line is interpreted like this:

    contig_2203 is linked to contig_3479.
    The first read is towards the start of contig_2203 (len_1 == 1664, pos_1 == 334) and is in the same orientation as the contig (rev_1 == 0)
    The second (paired) read is towards the end of contig_3479 (len_2 == 3873, pos_2 == 2866) and is also in the same orientation as the contig (rev_2 == 0)
    The linking information was extracted from file: f2.bam.

The lines following this one describe other links between the two contigs.
Finally, note that the -c and -l options are not mutually exclusive and can be run at the same time.

#-----------------------
# EX: Extract all reads mapping to a particular set of contigs

    $ bamm extract -g group1.file group2.file -b f1.bam f2.bam f3.bam

Will extract all reads from each of the three BAM files that map to the contigs in group1 or group2.
The 'group' files can be multiple (gzipped) FASTA (like the fna files you can extract from GroopM) or lists of contig headers (one sequence per line).

Unless specified otherwise, BamM differentiates between paired and unpaired reads (from a mapping and group perspective), reads from different BAM files and reads mapping to contigs in different groups.


static const char MITEXT[6][2][12] = {{"p_PR_PM_UG;\0",   // FIR, U
                                       "p_PR_PM_PG;\0"},  // FIR, P
                                      {"p_PR_PM_UG;\0",  // SEC, U
                                       "p_PR_PM_PG;\0"},  // SEC, P
                                      {"p_PR_UM_NG;\0",  // SNGL_FIR, U
                                       "p_PR_EM_NG;\0"},  // SNGL_FIR, P
                                      {"p_PR_UM_NG;\0",  // SNGL_SEC, U
                                       "p_PR_EM_NG;\0"},  // SNGL_SEC, P
                                      {"p_UR_NM_NG;\0",  // SNGL, U
                                       "p_UR_EM_NG;\0"},  // SNGL, P
                                      {"p_ER_NM_NG;\0",  // ERR, U
                                       "p_ER_NM_NG;\0"}}; // ERR, P


NOTE: this command can produce A LOT of output files.

## Using the modules in your own Python code

    #-----------------------
    # First import it
    from bamm.BamParser import BamParser


    def parseBams(self,
                  bamFiles,
                  doLinks=False,
                  doCovs=False,
                  types=None,
                  threads=1,
                  verbose=False):

    #-----------------------
    # EX: use a BAM parser to get the orientation type of a BAM file only

    BP = BamParser(coverageMode="none")             # note you must set coverage to None

    BP.parseBams(['file1.bam', 'file2.bam'],        # list of file names
                 types=[2,1],                       # number of insert types per BAM
                 threads=2,                         # number of threads
                 verbose=True)                      # be more verbose

    BP.printBamTypes()                              # tell the world!

    #-----------------------
    # EX: use a BAM parser to get coverage profiles and links

    BP = BamParser()

    BP.parseBams(['file1.bam', 'file2.bam'],        # list of file names
                 doCovs=True,                       # indicate coverages should be calculated
                 doLinks=True,                      # work out links during parse
                 types=[2,1],                       # number of insert types per BAM
                 threads=3,                         # number of threads
                 verbose=True)                      # be more verbose

	# NOTE: doLinks implies doTypes!

    BP.printBamTypes()
    BP.printCoverages()
    BP.printLinks("links.tsv")						# can specify a filename for all the print functions

    #-----------------------
    # the full range of options at initialisation are:

     baseQuality                    # quality score threshold of reads to accept during pileup coverage calculation ( default = 0 )
     minLength                      # minimum length threshold of a mapped read ( default = 0 )
     mappingQuality                 # BWA/BAM mapping quality threshold
     coverageMode                   # 'vanilla' (standard pileup), 'outlier' (truncated mean)
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
