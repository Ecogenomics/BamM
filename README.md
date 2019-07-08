

**_DEPRECATION NOTICE_**

BamM is no longer being maintained. Instead try [CoverM](https://github.com/wwood/CoverM) which is easier to both install and use, and is faster.



# BamM (deprecated)

## Overview

BamM is a c library, wrapped in python, that parses BAM files.
The code is intended to provide a faster, more stable interface to parsing BAM files than PySam, but doesn't implement all/any of PySam's features.
Do you want all the links that join two contigs in a BAM? Do you need to get coverage? Would you like to just work out the insert size and orientation of some mapped reads?

Then BamM is for you!

## Installation

There are 3 ways to install BamM:

1. [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/bamm/README.html)
2. Through GNU Guix with `guix package --install bamm`
3. From source

### Instructions for installation from source

Dependencies:

The BAM parsing is done using c and a few external libraries. This slightly complicates the BamM installation process but fortunately not too much.

BamM is not available via pip, so you'll need to install it manually. Complete instructions are available in the manual at:

http://ecogenomics.github.io/BamM/manual/BamM_manual.pdf

If you're running 'bamm make' you'll need to have bwa and samtools installed. Installation of these tools is really straightforward. You can find the code and instructions at:

* Samtools:   http://sourceforge.net/projects/samtools  (tested with version: 1.2)
* BWA:        https://github.com/lh3/bwa            (tested with version: 0.7.12)

BamM depends on two libraries: htslib and libcfu. Several users have reported difficulty installing these dependencies so they are now packaged with the BamM source. It it still possible to use your local version of these libraries, see the manual for more details.

We have tested BamM on CentOS and Ubuntu and recommend these or similar systems. Several people have reported success with installing BamM on OSX varieties however a greater number have not been able to get it installed on their MACs. This is something we're aware of and are working on however **it is probably safe to assume that BamM is not supported on OSX**.

## Example usage

Please refer to the manual or http://ecogenomics.github.io/BamM for more information about example usage.

## Help

If you experience any problems using BamM, open an [issue](https://github.com/ecogenomics/BamM/issues) on GitHub and tell us about it.

## Licence and referencing

Project home page, info on the source tree, documentation, issues and how to contribute, see http://github.com/ecogenomics/BamM

This software is currently unpublished

## Copyright

Copyright (c) Ben Woodcroft, Tim Lamberton, Michael Imelfort, See LICENSE.txt for further details.
