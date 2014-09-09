#!/usr/bin/env python
###############################################################################
#                                                                             #
#    BamExtractor.py                                                          #
#                                                                             #
#    Class for extracting reads from BAM files                                #
#                                                                             #
#    Copyright (C) Michael Imelfort                                           #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2014"
__credits__ = ["Michael Imelfort"]
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Beta"

###############################################################################

# system imports
import os
import ctypes as c
from multiprocessing import Pool, Manager
import numpy as np
import sys
import gzip
import mimetypes

# local imports
from cWrapper import *
from bamFile import *
from bammExceptions import *

###############################################################################
###############################################################################
###############################################################################
# Multiprocessing requires that all passed items be pickleable. That is they
# must be vanilla variables or functions defined in the file itself, ie. not
# within a class. We get around this by writing an external function which calls
# a class function. Hacky, but it works.
###############################################################################
###############################################################################
###############################################################################

def externalParseWrapper(bAMeXTRACTOR, bid, gBFI, func, doContigNames):
    """ctypes pointers are unpickleable -- what we need is a hack!

    See BamParser._parseOneBam for what this function should be doing
    """
    pass
    # go back into the class to do the work
    """
    coverages = []
    contig_lengths = None
    contig_names = None
    links = {}

    BFI = bAMpARSER._parseOneBam(bid)

    # only do this if we are doing covs or links (or both)
    if bAMpARSER.doCovs or bAMpARSER.doLinks:
        contig_lengths = np.array([int(i) for i in c.cast(BFI.contigLengths, c.POINTER(c.c_uint32*BFI.numContigs)).contents])
        plpBp = np.array([[int(j) for j in c.cast(i, c.POINTER(c.c_uint32*BFI.numBams)).contents] for i in c.cast(BFI.plpBp,c.POINTER(c.POINTER(c.c_uint32*BFI.numBams)*BFI.numContigs)).contents])

        # transfer the coverages over
        coverages = np.zeros((BFI.numContigs, BFI.numBams))
        if bAMpARSER.coverageMode == 'outlier':
            contig_length_correctors = np.array([[int(j) for j in c.cast(i, c.POINTER(c.c_uint32*BFI.numBams)).contents] for i in c.cast(BFI.contigLengthCorrectors,c.POINTER(c.POINTER(c.c_uint32*BFI.numBams)*BFI.numContigs)).contents])
            for c_idx in range(int(BFI.numContigs)):
                for b_idx in range(int(BFI.numBams)):
                    coverages[c_idx,b_idx] = float(plpBp[c_idx,b_idx])/float(contig_lengths[c_idx] - contig_length_correctors[c_idx])
        else:
            for c_idx in range(BFI.numContigs):
                for b_idx in range(BFI.numBams):
                    coverages[c_idx,b_idx] = float(plpBp[c_idx,b_idx])/float(contig_lengths[c_idx])

        # we only need to do the contig names for one of the threads
        if doContigNames:
            contig_names = []
            contig_name_lengths = np.array([int(i) for i in c.cast(BFI.contigNameLengths, c.POINTER(c.c_uint16*BFI.numContigs)).contents])
            contig_name_array = c.cast(BFI.contigNames, c.POINTER(c.POINTER(c.c_char)*BFI.numContigs)).contents
            for i in range(BFI.numContigs):
                contig_names.append("".join([j for j in c.cast(contig_name_array[i], c.POINTER(c.c_char*contig_name_lengths[i])).contents]))

    # we always populate the bam file type information classes
    bam_file_name = bAMpARSER.bamFiles[bid]
    BF = BM_bamFile(bid, bam_file_name)
    BF_C = (c.cast(BFI.bamFiles, c.POINTER(c.POINTER(BM_bamFile_C)*1)).contents)[0].contents
    num_types = BF_C.numTypes
    BTs_C = c.cast(BF_C.types, c.POINTER(c.POINTER(BM_bamType_C)*num_types)).contents
    for bt_c in BTs_C:
        BT = BM_bamType((bt_c.contents).orientationType,
                        (bt_c.contents).insertSize,
                        (bt_c.contents).insertStdev,
                        (bt_c.contents).supporting)
        BF.types.append(BT)

    if bAMpARSER.doLinks:
        links = pythonizeLinks(BFI, BF, contig_lengths)
    else:
        links = {}

    # make the python object
    BBFI = BM_fileInfo(coverages,
                       contig_lengths,
                       BFI.numBams,
                       BFI.numContigs,
                       contig_names,
                       [BF],
                       links)

    # append onto the global list
    gBFI.append(BBFI)

    # destroy the C-allocateed memory
    pBFI = c.POINTER(BM_fileInfo_C)
    pBFI = c.pointer(BFI)
    CW = CWrapper()
    CW._destroyBFI(pBFI)
    """

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamExtractor:
    """Main class for reading in and parsing contigs"""
    def __init__(self,
                targets,                    # list of contig IDs or fasta file (used as a filter)
                bamfiles,                   # list of bamfiles to extract reads from
                prefix="",                  # append thiss to all output files
                outFolder=".",              # wriate output to this folder
                shuffle=False,              # use shuffled format for paired reads
                mixBams=False,              # use one file for all bams
                combineReadsFalse=False,    # combine paired and unpaired into one file
                ignoreUnpaired=False,       # ignore all upaired reads
                bigFile=False,              # do NOT gzip outputs
                headersOnly=False           # write read headers only
                ):

        # make sure the output folder exists
        self.outFolder = outFolder
        self.makeSurePathExists(self.outFolder)

        # work out how we'll write files
        if bigFile:
            self.writeOpen = open
        else:
            self.writeOpen = gzip.open

        self.prefix = prefix
        self.shuffle = shuffle
        self.mixBams = mixBams
        self.ignoreUnpaired = ignoreUnpaired
        self.headersOnly = headersOnly

        # munge the targets
        self.targets = []
        try:
            read_open = open
            # handle gzipped files
            mime = mimetypes.guess_type(targets)
            if mime[1] == 'gzip':
                read_open = gzip.open
        except:
            raise InvalidParameterSetException('Error when guessing targets file mimetype')
        with read_open(targets, "r") as t_fh:
            self.makeTargetList(t_fh)
        if len(self.targets) == 0:
            raise InvalidParameterSetException('No targets supplied')

    def makeTargetList(self, t_fh):
        """Get the list of targets to hit"""
        # work out if the targets are lists of contig IDs or just contigs
        # assume that if the file is fasta then the first character will be ">"
        # otherwise it must be a list
        first_line = t_fh.readline()
        try:
            if first_line[0] == ">":
                t = first_line.rstrip()[1:]
                if t != "":
                    self.targets.append(t)
                for line in t_fh:
                    if line[0] == ">":
                        t = line.rstrip()[1:]
                        if t != "":
                            self.targets.append(t)
            else:
                t = first_line.rstrip()
                if t != "":
                    self.targets.append(t)
                for line in t_fh:
                    t = line.rstrip()
                    if t != "":
                        self.targets.append(t)
        except:
            raise InvalidParameterSetException('Something is wrong with the supplied targets file')

    def makeOutputFiles(self):
        """Open all the output file handles we'll need"""
        pass

    def extract(self):
        """Extract all the reads"""
        self.makeOutputFiles()

    def makeSurePathExists(self, path):
        try:
            os.makedirs(path)
        except OSError as exception:
            import errno
            if exception.errno != errno.EEXIST:
                raise
