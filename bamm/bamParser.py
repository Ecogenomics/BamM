#!/usr/bin/env python
###############################################################################
#                                                                             #
#    BamParser.py                                                             #
#                                                                             #
#    Class for parsing BAM files                                              #
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

# local imports
from cWrapper import *
from bamLink import *
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

def pythonizeLinks(BFI, bamFile, contigLengths):
    """Unwrap the links-associated C structs and return a python-ized dict"""
    links = {}
    CW = CWrapper()
    pBFI = c.POINTER(BM_fileInfo_C)
    pBFI = c.pointer(BFI)

    LW = BM_LinkWalker_C()
    pLW = c.POINTER(BM_LinkWalker_C)
    pLW = c.pointer(LW)
    success = CW._initLW(pLW, pBFI)
    if(success == 2):
        ret_val = 2
        LP = None
        while(ret_val != 0):
            if ret_val == 2:
                # need a new contig pair
                LP = BM_linkPair(((LW.pair).contents).cid1, ((LW.pair).contents).cid2)
                key = "%d,%d" % (((LW.pair).contents).cid1, ((LW.pair).contents).cid2)
                links[key] = LP
            # add a link
            LI = (LW.LI).contents
            LP.addLink(LI.reversed1,
                       LI.reversed2,
                       LI.pos1,
                       LI.pos2,
                       bamFile)
            ret_val = CW._stepLW(pLW)
        CW._destroyLW(pLW)

    return links

def externalParseWrapper(bAMpARSER, bid, gBFI, func, doContigNames):
    """ctypes pointers are unpickleable -- what we need is a hack!

    See BamParser._parseOneBam for what this function should be doing
    """
    # go back into the class to do the work
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

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self,
                 baseQuality=0,
                 minLength=0,
                 mappingQuality=0,
                 coverageMode='vanilla',
                 ignoreSuppAlignments=True
                 ):
        #---------------------------------
        # information about how the parser will be used
        #---------------------------------
        self.baseQuality = baseQuality
        self.mappingQuality = mappingQuality
        self.minLength = minLength

        if ignoreSuppAlignments:
            self.ignoreSuppAlignments = 1
        else:
            self.ignoreSuppAlignments = 0

        if coverageMode not in ['vanilla', 'outlier', 'none']:
             raise InvalidCoverageModeException("Unknown coverage mode '%s' supplied" % coverageMode)
        self.coverageMode = coverageMode

        #---------------------------------
        # internal variables
        #---------------------------------
        self.BFI = None          # internal mapping results object

        # these are set when we make the call to parse
        self.bamFiles = []
        self.types = []
        self.doLinks = False
        self.doTypes = False
        self.doCovs = False

#------------------------------------------------------------------------------
# Bam parseratering

    def typeBams(self, bamFiles, types=None, threads=1):
        # set these now
        self.bamFiles = bamFiles
        if types is None:
            self.types = [1]*len(self.bamFiles)
        else:
            self.types = types

        # check that the bam files and their indexes exist
        for bam in bamFiles:
            if not os.path.isfile(bam):
                raise BAMFileNotFoundException("BAM file %s could not be found" % bam)
            elif not os.path.isfile("%s.bai" % bam):
                raise BAMIndexNotFoundException("Index file %s could not be found" % ("%s.bai" % bam))

        global gBFI
        gBFI = Manager().list()
        pool = Pool(processes=threads)
        for bid in range(len(bamFiles)):
            pool.apply_async(func=externalParseWrapper, args=(self, bid, gBFI, "type", False))
        pool.close()
        pool.join()

        # merge all the separate mapping results
        self.BFI = gBFI[0]
        for i in range(1, len(gBFI)):
            self.BFI.consume(gBFI[i])

    def parseBams(self,
                  bamFiles,
                  doLinks=False,
                  doTypes=False,
                  doCovs=False,
                  types=None,
                  threads=1):
        """Parse bam files to get coverage and linking reads

        stores results in internal mapping results list
        """
        # set these now
        self.bamFiles = bamFiles

        # how may insert types for each bam file?
        if types is None:
            self.types = [1]*len(self.bamFiles)
        else:
            self.types = types

        if len(self.types) != len(self.bamFiles):
            raise InvalidNumberOfTypesException("%d types for %d BAM files" % (len(self.types), len(self.bamFiles)))

        # make sure (again) that we're doing something
        self.doLinks = doLinks
        self.doCovs = doCovs
        if not (self.doCovs or self.doLinks):
            self.doTypes = True
        else:
            self.doTypes = doTypes

        # check that the bam files and their indexes exist
        for bam in bamFiles:
            if not os.path.isfile(bam):
                raise BAMFileNotFoundException("BAM file %s could not be found" % bam)
            elif not os.path.isfile("%s.bai" % bam):
                raise BAMIndexNotFoundException("Index file %s could not be found" % ("%s.bai" % bam))

        global gBFI
        gBFI = Manager().list()
        pool = Pool(processes=threads)
        do_contig_names = True
        for bid in range(len(bamFiles)):
            pool.apply_async(func=externalParseWrapper, args=(self, bid, gBFI, "parse", do_contig_names))
            if do_contig_names:
                # we only need to parse the contig names once
                do_contig_names = False
        pool.close()
        pool.join()

        baseBFI_index = 0
        if self.doCovs or self.doLinks:
            # all the BFIs are made. Only one has the contig IDs. find it's index
            for i in range(len(gBFI)):
                if len(gBFI[i].contigNames) > 0:
                    baseBFI_index = i
                    break

        # merge all the separate mapping results
        self.BFI = gBFI[baseBFI_index]
        for i in range(len(gBFI)):
            if i != baseBFI_index:
                self.BFI.consume(gBFI[i])

    def _typeOneBam(self, bid):
        """Parse a single BAM to get it's type"""
        BFI = BM_fileInfo_C()        # destroy needs to be called on this -> it should be called by the calling function
        pBFI = c.POINTER(BM_fileInfo_C)
        pBFI = c.pointer(BFI)

        bamfiles_c_array = (c.c_char_p * 1)()
        bamfiles_c_array[:] = [self.bamFiles[bid]]

        types_c_arr = (c.c_int * 1)()
        types_c_arr[:] = [self.types[bid]]

        CW = CWrapper()
        CW._parseCoverageAndLinks(1,        # set typeOnly flag
                                  1,
                                  0,
                                  0,
                                  0,
                                  types_c_arr,
                                  1,
                                  c.create_string_buffer("none"),
                                  bamfiles_c_array,
                                  pBFI)
        return BFI

    def _parseOneBam(self, bid):
        """Parse a single BAM file and append the result to the internal mapping results list"""
        BFI = BM_fileInfo_C()        # destroy needs to be called on this -> it should be called by the calling function
        pBFI = c.POINTER(BM_fileInfo_C)
        pBFI = c.pointer(BFI)

        bamfiles_c_array = (c.c_char_p * 1)()
        bamfiles_c_array[:] = [self.bamFiles[bid]]

        types_c_arr = (c.c_int * 1)()
        types_c_arr[:] = [self.types[bid]]

        CW = CWrapper()
        if self.doLinks or self.doCovs:
            CW._parseCoverageAndLinks(0,        # unset typeOnly flag
                                      1,
                                      self.baseQuality,
                                      self.mappingQuality,
                                      self.minLength,
                                      types_c_arr,
                                      self.ignoreSuppAlignments,
                                      c.create_string_buffer(self.coverageMode),
                                      bamfiles_c_array,
                                      pBFI)
        else:
            # types only
            CW._parseCoverageAndLinks(1,        # set typeOnly flag
                                      1,
                                      0,
                                      0,
                                      0,
                                      types_c_arr,
                                      1,
                                      c.create_string_buffer("none"),
                                      bamfiles_c_array,
                                      pBFI)

        return BFI

#------------------------------------------------------------------------------
# Printing and IO

    def printBamTypes(self, fileName=""):
        if self.BFI is None:
            raise NoBAMSFoundException
        else:
            if fileName == "":
                self.BFI.printBamTypes(sys.stdout)
            else:
                with open(fileName, "w") as fh:
                    self.BFI.printBamTypes(fh)

    def printCoverages(self, fileName=""):
        if self.BFI is None:
            raise NoBAMSFoundException
        else:
            if fileName == "":
                self.BFI.printCoverages(sys.stdout)
            else:
                with open(fileName, "w") as fh:
                    self.BFI.printCoverages(fh)

    def printLinks(self, fileName=""):
        if self.BFI is None:
            raise NoBAMSFoundException
        else:
            if fileName == "":
                self.BFI.printLinks(sys.stdout)
            else:
                with open(fileName, "w") as fh:
                    self.BFI.printLinks(fh)


###############################################################################
###############################################################################
###############################################################################
###############################################################################

def makeSurePathExists(path):
    try:
        os_makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

###############################################################################
###############################################################################
###############################################################################
###############################################################################
