#BamM

## Overview

BamM is a c library, wrapped in python, that parses BAM files.
The code is intended to provide a faster, more stable interface to parsing BAM files than PySam, but doesn't implement all/any of PySam's features.
Do you want all the links that join two contigs in a BAM? Do you need to get coverage? Would you like to just work out the insert size and orientation of some mapped reads?

Then BamM is for you!

## Installation

Dependencies:

The BAM parsing is done using c and a few external libraries. This slightly complicates the BamM installation process but fortunately not too much.

BamM is not available via pip, so you'll need to install it manually. Complete instructions are available in the manual at:

http://minillinim.github.io/BamM/manual/BamM_manual.pdf

If you're running 'bamm make' you'll need to have bwa and samtools installed. Installation of these tools is really straightforward. You can find the code and instructions at:

Samtools:   https://github.com/samtools/samtools
BWA:        https://github.com/lh3/bwa

If you're installing system-wide then you can use your favourite package manager to install htslib and libcfu. For local installs, or installs that will work with the linux "modules" system, you need to be a bit trickier. This is the type of install that's documented in the manual.

## Example usage

Please refer to the manual or http://minillinim.github.io/BamM for more information about example usage.

## Help

If you experience any problems using BamM, open an [issue](https://github.com/minillinim/BamM/issues) on GitHub and tell us about it.

## Licence and referencing

Project home page, info on the source tree, documentation, issues and how to contribute, see http://github.com/minillinim/BamM

This software is currently unpublished

## Copyright

Copyright (c) Michael Imelfort, See LICENSE.txt for further details.
