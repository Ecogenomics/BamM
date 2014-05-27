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
    int n = 0, do_links = 0, baseQ = 0, mapQ = 0, min_len = 0, do_outlier_coverage = 0;
    while ((n = getopt(argc, argv, "q:Q:l:Lo")) >= 0) {
        switch (n) {
            case 'l': min_len = atoi(optarg); break; // minimum query length
            case 'q': baseQ = atoi(optarg); break;   // base quality threshold
            case 'Q': mapQ = atoi(optarg); break;    // mapping quality threshold
            case 'L': do_links = 1; break;
            case 'o': do_outlier_coverage = 1; break;
        }
    }
    if (optind == argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: samtools depth [options] in1.bam [in2.bam [...]]\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "   -L                  find pairing links\n");
        fprintf(stderr, "   -l <int>            minQLen\n");
        fprintf(stderr, "   -q <int>            base quality threshold\n");
        fprintf(stderr, "   -Q <int>            mapping quality threshold\n");
        fprintf(stderr, "   -o                  do outlier coverage corrections\n");
        fprintf(stderr, "\n");
        return 1;
    }

    int num_bams = argc - optind; // the number of BAMs on the command line
    int i = 0;
    char **bam_files = calloc(num_bams, sizeof(char*));             // bam file names
    for (i = 0; i < num_bams; ++i) {
        bam_files[i] = strdup(argv[optind+i]);
    }

    BMM_mapping_results * mr = calloc(1, sizeof(BMM_mapping_results));
    int ignore_supps = 1;
    int ret_val = parseCoverageAndLinks(num_bams,
                                        baseQ,
                                        mapQ,
                                        min_len,
                                        do_links,
                                        ignore_supps,
                                        do_outlier_coverage,
                                        bam_files,
                                        mr);
    print_MR(mr);
    destroy_MR(mr);

    for (i = 0; i < num_bams; ++i) {
        free(bam_files[i]);
    }
    free(bam_files);
    return ret_val;
}
