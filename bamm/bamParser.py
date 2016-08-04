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
#    This library is free software; you can redistribute it and/or            #
#    modify it under the terms of the GNU Lesser General Public               #
#    License as published by the Free Software Foundation; either             #
#    version 3.0 of the License, or (at your option) any later version.       #
#                                                                             #
#    This library is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        #
#    Lesser General Public License for more details.                          #
#                                                                             #
#    You should have received a copy of the GNU Lesser General Public         #
#    License along with this library.                                         #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2014,2015"
__credits__ = ["Michael Imelfort"]
__license__ = "LGPLv3"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"

###############################################################################

# system imports
import os
import ctypes as c
from multiprocessing import Pool, Manager
import multiprocessing as mp
import numpy as np
import sys

# local imports
from cWrapper import (CWrapper,
                      BM_bamFile_C,
                      BM_fileInfo_C,
                      BM_LinkWalker_C,
                      BM_coverageType_C,
                      BM_bamType_C,
                      CT)
from bamLink import BM_linkPair
from bamFile import BM_bamFile, BM_bamType, BM_fileInfo, BM_coverageType
from bamMaker import BamValidator
from bammExceptions import (InvalidNumberOfTypesException,
                            BAMFileNotFoundException,
                            BAMIndexNotFoundException,
                            NoBAMSFoundException)

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

def externalParseWrapper(bAMpARSER,
                         parseQueue,
                         BFI_list,
                         verbose,
                         doContigNames):
    '''Single-process BAMfile parsing

    cTypes pointers are unpickleable unless they are top level, so this function
    lives outside the class. In this case we reduce the number of member
    variables passed to it by passing the class instead. Any implicit copy
    operations do not affect the workflow as it stands now. If you modify this
    function you need to be aware of the limitations of python multiprocessing,
    Queues, pickling and shared memory.

    Extra logic is also contained in BamParser._parseOneBam

    Inputs:
     bAMpARSER - BamParser instance, a valid BamParser instance
     parseQueue - Manager.Queue, bids (BAMs) yet to be parsed
     BFI_list - Manager.List, place all processed BFIs on this list
     verbose - == True -> be verbose
     doContigNames - == True -> load contigs names from the C-land BFI struct
    '''
    CW = CWrapper()
    while True:
        # get the next one off the list
        bid = parseQueue.get(block=True, timeout=None)
        if bid is None: # poison pill
            break

        if verbose:
            print "Parsing file: %s" % bAMpARSER.bamFiles[bid]

        # go back into the class to do the work
        coverages = []
        contig_lengths = []
        contig_names = []
        links = {}

        BFI = bAMpARSER._parseOneBam(bid)

        # only do this if we are doing covs or links (or both)
        if bAMpARSER.doCovs or bAMpARSER.doLinks:
            contig_lengths = \
                np.array([int(i) for i in
                          c.cast(BFI.contigLengths,
                                 c.POINTER(c.c_uint32*BFI.numContigs)).contents
                          ])

            coverages = np.array([[float(j) for j in c.cast(i, c.POINTER(c.c_float*BFI.numBams)).contents] for i in
                                  c.cast(BFI.coverages, c.POINTER(c.POINTER(c.c_float*BFI.numBams)*BFI.numContigs)).contents])

            # we only need to do the contig names for one of the threads
            if doContigNames:
                contig_names = []
                contig_name_lengths = \
                    np.array([int(i) for i in
                              c.cast(BFI.contigNameLengths,
                                     c.POINTER(c.c_uint16*BFI.numContigs)
                                     ).contents
                              ])

                contig_name_array = \
                    c.cast(BFI.contigNames,
                           c.POINTER(c.POINTER(c.c_char)*BFI.numContigs)
                           ).contents

                for i in range(BFI.numContigs):
                    contig_names.append((c.cast(contig_name_array[i],
                                                c.POINTER(c.c_char * \
                                                          contig_name_lengths[i]
                                                          )
                                                ).contents
                                         ).value
                                        )

        # we always populate the bam file type information classes
        bam_file_name = bAMpARSER.bamFiles[bid]
        BF = BM_bamFile(bid, bam_file_name)
        BF_C = \
            (c.cast(BFI.bamFiles,
                    c.POINTER(c.POINTER(BM_bamFile_C)*1)).contents)[0].contents

        num_types = BF_C.numTypes
        BTs_C = c.cast(BF_C.types,
                       c.POINTER(c.POINTER(BM_bamType_C)*num_types)).contents

        for bt_c in BTs_C:
            BT = BM_bamType((bt_c.contents).orientationType,
                            (bt_c.contents).insertSize,
                            (bt_c.contents).insertStdev,
                            (bt_c.contents).supporting)
            BF.types.append(BT)

        if bAMpARSER.doLinks:
            links = pythonizeLinks(BFI, BF)
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
        BFI_list.append(BBFI)

        # destroy the C-allocateed memory
        pBFI = c.POINTER(BM_fileInfo_C)
        pBFI = c.pointer(BFI)
        CW._destroyBFI(pBFI)

        if doContigNames:
            # we only need to parse the contig names once
            doContigNames = False

def pythonizeLinks(BFI, bamFile):
    '''Unpeel the links-associated C structs and return a python dictionary
    of LinkPair instances

    Inputs:
     BFI - BM_fileInfo_C, C-land BamFileInfo struct
     bamFile - uid of the bamFile associated with the BFI

    Outputs:
     A python dictionary of LinkPair instances
    '''
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
                LP = BM_linkPair(((LW.pair).contents).cid1,
                                 ((LW.pair).contents).cid2)
                # makeKey should return unique ID
                links[LP.makeKey()] = LP
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

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamParser:
    ''' class for generating coverage profiles and linking read information
    from several BAM files

    Uses python multiprocessing to parallelize processing
    '''

    def __init__(self,
                 coverageType,
                 minLength=0,
                 baseQuality=0,
                 mappingQuality=0,
                 useSuppAlignments=False,
                 useSecondaryAlignments=False,
                 maxMisMatches=1000
                 ):
        '''
        Default constructor.

        Set quality thresholds used in later parsing

        Inputs:
         coverageType - BM_coverageType, stores type of coverage to calculate
         minLength - int, ignore contigs shorter than this length
         mappingQuality - int, ignore positions with a lower base quality score
         mappingQuality - int, skip all reads with a lower mapping quality score
         useSuppAlignments - == True -> DON'T skip supplementary alignments
         useSecondaryAlignments - == True -> DON'T skip secondary alignments
         maxMisMatches - int, skip all reads with more mismatches (NM aux files)

        Outputs:
         None
        '''
        #---------------------------------
        # information about how the parser will be used
        #---------------------------------
        self.baseQuality = baseQuality
        self.mappingQuality = mappingQuality
        self.minLength = minLength
        self.maxMisMatches = maxMisMatches

        if useSuppAlignments:
            self.ignoreSuppAlignments = 0
        else:
            self.ignoreSuppAlignments = 1

        if useSuppAlignments:
            self.ignoreSecondaryAlignments = 0
        else:
            self.ignoreSecondaryAlignments = 1

        self.coverageType = coverageType

        #---------------------------------
        # internal variables
        #---------------------------------
        self.BFI = None          # internal mapping results object

        # these are set when we make the call to parse
        self.bamFiles = []
        self.types = []
        self.doLinks = False
        self.doCovs = False

#------------------------------------------------------------------------------
# Bam parseratering

    def parseBams(self,
                  bamFiles,
                  doLinks=False,
                  doCovs=False,
                  types=None,
                  threads=1,
                  verbose=False):
        '''Parse bam files to get coverage and linking reads.

        Manages threading
        Stores results in internal mapping results list

        Inputs:
         bamFiles - [string], full paths to BAMs
         doLinks - == True -> find linking pairs
         doCovs - == True -> calculate coverage profiles
         types - [int] or None, number of insert types per bamfile
         threads - int, max number of threads to use
         verbose - == True -> be verbose

        Outputs:
         0 if the parsing worked, 1 otherwise
        '''
        # set these now
        self.bamFiles = bamFiles

        # how may insert types for each bam file?
        if types is None:
            self.types = [1]*len(self.bamFiles)
        else:
            self.types = types

        if len(self.types) != len(self.bamFiles):
            raise InvalidNumberOfTypesException("%d types for %d BAM files" % \
                                                (len(self.types),
                                                 len(self.bamFiles)))

        # make sure (again) that we're doing something
        self.doLinks = doLinks
        self.doCovs = doCovs

        BV = BamValidator(silent=not verbose)
        # check that the bam files and their indexes exist
        for bam in bamFiles:
            if not os.path.isfile(bam):
                raise BAMFileNotFoundException("BAM file %s not found" % bam)
            elif not os.path.isfile("%s.bai" % bam):
                raise BAMIndexNotFoundException("Index file %s not found" % \
                                                ("%s.bai" % bam))
            BV.validate_bam(bam)

        # start running the parser in multithreaded mode
        parse_queue = mp.Queue()
        # each thread can place their new BFIs on a single global list
        BFI_list = mp.Manager().list()

        # place the bids on the queue
        for bid in range(len(bamFiles)):
            parse_queue.put(bid)

        # place one None on the queue for each thread we have access to
        for _ in range(threads):
            parse_queue.put(None)

        try:
            # only the first thread and the first job should parse contig names
            parse_proc = [mp.Process(target=externalParseWrapper,
                                     args = (self,
                                             parse_queue,
                                             BFI_list,
                                             verbose,
                                             True)
                                     )
                          ]
            # all the other threads will not parse contig names
            parse_proc += [mp.Process(target=externalParseWrapper,
                                      args = (self,
                                              parse_queue,
                                              BFI_list,
                                              verbose,
                                              True)
                                      ) for _ in range(threads-1)
                           ]

            for p in parse_proc:
                p.start()

            for p in parse_proc:
                p.join()

            # all processes are finished, collapse the BFI_list
            self.collapseBFIs(BFI_list)

            # success
            return 0

        except  (KeyboardInterrupt, SystemExit):
            # ctrl-c! Make sure all processes are terminated
            for p in parse_proc:
                p.terminate()

            # dismal failure
            return 1
            
    def _parseOneBam(self, bid):
        '''Parse a single BAM file and append the result
        to the internal mapping results list

        Called from the ExternalParseWrapper

        Inputs:
         bid - unique identifier of the BAM to parse

        Outputs:
         A populated BM_FileInfo_C  struct containing the parsing results
        '''
        # destroy needs to be called on this
        # -> it should be called by the calling function
        BFI = BM_fileInfo_C()
        pBFI = c.POINTER(BM_fileInfo_C)
        pBFI = c.pointer(BFI)

        BCT = BM_coverageType_C()
        BCT.type = self.coverageType.cType
        BCT.upperCut = float(self.coverageType.cUpper)
        BCT.lowerCut = float(self.coverageType.cLower)
        pBCT = c.POINTER(BM_coverageType_C)
        pBCT = c.pointer(BCT)

        bamfiles_c_array = (c.c_char_p * 1)()
        bamfiles_c_array[:] = [self.bamFiles[bid]]

        types_c_array = (c.c_int * 1)()
        types_c_array[:] = [self.types[bid]]

        CW = CWrapper()
        if self.doLinks or self.doCovs:
            CW._parseCoverageAndLinks(self.doLinks,
                                      self.doCovs,
                                      1,        # numBams always one here
                                      self.baseQuality,
                                      self.mappingQuality,
                                      self.minLength,
                                      self.maxMisMatches,
                                      types_c_array,
                                      self.ignoreSuppAlignments,
                                      self.ignoreSecondaryAlignments,
                                      pBCT,
                                      bamfiles_c_array,
                                      pBFI)
        else:
            # types only
            BCT.type = CT.NONE # just to be sure!
            CW._parseCoverageAndLinks(self.doLinks,
                                      self.doCovs,
                                      1,        # numBams always one here
                                      0,
                                      0,
                                      0,
                                      0,
                                      types_c_array,
                                      1,
                                      1,
                                      pBCT,
                                      bamfiles_c_array,
                                      pBFI)

        return BFI

    def collapseBFIs(self, BFI_list):
        '''Collapse multiple BFI objects into one and
        set it as member variable BFI

        Only one of the threads will bother to parse contig names. Find it's
        BFI and then merge all the other BFIs (from other threads) into it

        Inputs:
         BFI_list - Manager.List, list of all BFIs created during parsing

        Outputs:
         None
        '''
        baseBFI_index = 0
        if self.doCovs or self.doLinks:
            # all the BFIs are made. Only one has the contig IDs. find it.
            for i in range(len(BFI_list)):
                if len(BFI_list[i].contigNames) > 0:
                    baseBFI_index = i
                    break

        # merge all the separate mapping results
        self.BFI = BFI_list[baseBFI_index]
        for i in range(len(BFI_list)):
            if i != baseBFI_index:
                self.BFI.consume(BFI_list[i])

#------------------------------------------------------------------------------
# Printing and IO

    def printBamTypes(self, fileName=""):
        '''Print template size and orientation information
        to a file or to stdout.

        Inputs:
         fileName - string, full path to output file or ""

        Outputs:
         None
        '''
        if self.BFI is None:
            raise NoBAMSFoundException
        else:
            if fileName == "":
                self.BFI.printBamTypes(sys.stdout)
            else:
                with open(fileName, "w") as fh:
                    self.BFI.printBamTypes(fh)

    def printCoverages(self, fileName=""):
        '''Print coverage profile information to a file or to stdout

        Inputs:
         fileName - string, full path to output file or ""

        Outputs:
         None
        '''
        if self.BFI is None:
            raise NoBAMSFoundException
        else:
            if fileName == "":
                self.BFI.printCoverages(sys.stdout)
            else:
                with open(fileName, "w") as fh:
                    self.BFI.printCoverages(fh)

    def printLinks(self, fileName=""):
        '''Print paired links information to a file or to stdout

        Inputs:
         fileName - string, full path to output file or ""

        Outputs:
         None
        '''
        if self.BFI is None:
            raise NoBAMSFoundException
        else:
            if fileName == "":
                self.BFI.printLinks(sys.stdout)
            else:
                with open(fileName, "w") as fh:
                    self.BFI.printLinks(dict(zip(range(len(self.bamFiles)),self.bamFiles)), fh)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
