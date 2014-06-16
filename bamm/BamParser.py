#!/usr/bin/env python
###############################################################################
#                                                                             #
#    BamParser.py                                                             #
#                                                                             #
#    Class for parsing BAM files                                              #
#                                                                             #
#    Copyright (C) Michael Imelfort, Donovan Parks                            #
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
__credits__ = ["Michael Imelfort, Donovan Parks"]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Dev"

###############################################################################

import os
import ctypes as c
import pkg_resources
from multiprocessing import Pool, Manager
import numpy as np
from string import maketrans as s_maketrans
from re import sub as re_sub

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamMException(BaseException): pass
class InvalidCoverageModeException(BamMException): pass

###############################################################################
###############################################################################
###############################################################################
###############################################################################

#------------------------------------------------------------------------------
# Managing orientation and linking types

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

# Read orientations
# type 0 <--- --->
# type 1 ---> --->
# type 2 ---> <---
global OT
OT = enum('OUT', 'SAME', 'IN', 'NONE', 'ERROR')

# Links are stored as triples (contig1, contig2, linktype)
# There are 4 linktypes:
# SS  <--1--- ---2-->
# SE  <--1--- <--2---
# ES  ---1--> ---2-->
# EE  ---1--> <--2---
global LT
LT = enum('SS','SE','ES','EE','ERROR')

def OT2Str(OT):
    """For the humans!"""
    if OT == OT.OUT:
        return 'OUT'
    if OT == OT.SAME:
        return 'SAME'
    if OT == OT.IN:
        return 'IN'
    if OT == OT.NONE:
        return 'NONE'

    return 'ERROR'

def LT2Str(cid1, cid2, gap, linkType, terse=False):
    """For the humans!"""
    if terse:
        if linkType == LT.SS:
            return "SS"
        if linkType == LT.SE:
            return "SE"
        if linkType == LT.ES:
            return "ES"
        if linkType == LT.EE:
            return "EE"
        return "??"
    if linkType == LT.SS:
        return str(cid1) + " lies before " + str(cid2) + " in the opposite direction with gap "+str(gap)
    if linkType == LT.SE:
        return str(cid1) + " lies after " + str(cid2) + " in the same direction with gap "+str(gap)
    if linkType == LT.ES:
        return str(cid1) + " lies before " + str(cid2) + " in the same direction with gap "+str(gap)
    if linkType == LT.EE:
        return str(cid1) + " lies after " + str(cid2) + " in the opposite direction with gap "+str(gap)

    return 'Who knows?'

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamFile(object):
    """Store information about a BAM file"""
    def __init__(self,
                 bid,
                 fileName,
                 orientationType = OT.NONE,
                 insertSize = 0.,
                 insertStdev = 0.):
         self.bid = bid,
         self.fileName = fileName
         self.orientationType = orientationType
         self.insertSize = insertSize
         self.insertStdev = insertStdev

    def __str__(self):
        return "FileName: %s Orientation: %s Insert mean: %0.4f Stdev: %0.4f" % (self.fileName, OT2Str(OT), self.insertSize, self.insertStdev)


###############################################################################
###############################################################################
###############################################################################
###############################################################################


# fields defined in cfuhash.c but not accessed at this level
class cfuhash_table_t(c.Structure):
    pass

# links-associated structures "C land"
"""
typedef struct {
    uint16_t reversed1;
    uint16_t reversed2;
    uint16_t readLength1;
    uint16_t readLength2;
    uint32_t pos1;
    uint32_t pos2;
    uint32_t bam_ID;
    struct BMM_linkInfo * nextLink;
} BMM_linkInfo;

typedef struct {
    uint32_t cid1;
    uint32_t cid2;
    uint32_t numLinks;
    BMM_linkInfo * LI;
} BMM_linkPair;

typedef struct {
    char ** keys;
    size_t keyCount;
    size_t numKeys;
    cfuhash_table_t * linkHash;
    BMM_linkPair * pair;
    BMM_linkInfo * LI;
} BMM_LinkWalker;

"""
class BMM_linkInfo_C(c.Structure):
    pass

class BMM_linkInfo_C(c.Structure):
    _fields_ = [("reversed1", c.c_uint16),
                ("reversed2", c.c_uint16),
                ("readLength1", c.c_uint16),
                ("readLength2", c.c_uint16),
                ("pos1", c.c_uint32),
                ("pos2", c.c_uint32),
                ("bam_ID", c.c_uint32),
                ("nextLink",c.POINTER(BMM_linkInfo_C))
                ]

class BMM_linkPair_C(c.Structure):
    _fields_ = [("cid1", c.c_uint32),
                ("cid2", c.c_uint32),
                ("numLinks", c.c_uint32),
                ("LI",c.POINTER(BMM_linkInfo_C))
                ]

class BMM_LinkWalker_C(c.Structure):
    _fields_= [("keys", c.POINTER(c.POINTER(c.c_char))),
               ("keyCount", c.c_size_t),
               ("numKeys", c.c_size_t),
               ("links",c.POINTER(cfuhash_table_t)),
               ("pair",c.POINTER(BMM_linkPair_C)),
               ("LI",c.POINTER(BMM_linkInfo_C))
               ]

# links-associated structures "Python land"
class BMM_linkInfo(object):
    def __init__(self,
                 r1,
                 r2,
                 p1,
                 p2,
                 ot,
                 gap = 0,
                 bid = None):
        self.reversed1 = r1
        self.reversed2 = r2
        self.pos1 = p1
        self.pos2 = p2
        self.orientationType = ot
        self.gap = gap
        self.bamID = bid

    def __str__(self):
        return "    (%d,%d -> %d,%d, (%s, %d), %s)\n" % (self.pos1, self.reversed1, self.pos2, self.reversed2, LT2Str(self.orientationType), self.gap, self.bamID)

    def printMore(self, bamFileNames):
        return "    (%d,%d -> %d,%d, (%s, %d), %s)\n" % (self.pos1, self.reversed1, self.pos2, self.reversed2, LT2Str(self.orientationType), self.gap, bamFileNames[self.bamID])

class BMM_linkPair(object):
    def __init__(self,
                 cid1,
                 cid2):
        self.cid1 = cid1
        self.cid2 = cid2
        self.numLinks = 0
        self.links = []

    def addLink(self,
                r1,
                r2,
                p1,
                p2,
                rl,
                contigLengths,
                bamFile):

        LI = BMM_linkInfo(r1,
                          r2,
                          p1,
                          p2,
                          bamFile.bid)
        # work out the gap and the orientation type of this link
        (gap, ot) = self.determineOTDifferentRefs(LI, rl, bamFile, contigLengths)

        if ot != LT.ERROR:
            LI.gap = gap
            LI.orientationType = ot
            self.links.append(LI)
            self.numLinks += 1

    def determineOTDifferentRefs(self, LI, rl, bamFile, contigLengths):
        """Determine the orientation type and insert size of two reads

        We assume that:
           1. both reads are on different references
           2. ar1 is the first read in the pair

        I swear this is the last time I write this code!
        """
        # first check to see if the start read lies in the right position
        # max_ins should* be: mean + 3 * stdev but assemblers like SaSSY and Velvet
        # can produce nodes which overlap at the ends. So we need to account for negative gaps
        # thus we choose max ins = mean + 3 * stdev + rl
        max_ins = bamFile.insertSize + (3 * bamFile.insertStdev) + rl
        gap = bamFile.insertSize

        #print max_ins, cid1, LI.pos1, contigLengths[cid1], LI.reversed1, cid2, LI.pos2, contigLengths[cid2], LI.reversed2
        if LI.pos1 <= ( max_ins - 2 * rl ):
            # read 1 lies at the start of its contig
            r1_at_start = True
            gap -= (LI.pos1 + rl)
            max_ins -= (LI.pos1 + rl)
        elif LI.pos1 >= (contigLengths[cid1] - max_ins + rl ):
            # read 1 lies at the end of its contig
            r1_at_start = False
            gap -= (contigLengths[cid1] - LI.pos1)
            max_ins -= (contigLengths[cid1] - LI.pos1)
        else:
            # read 1 lies in the middle of its contig
            return (0, LT.ERROR)
        #print gap, ( max_ins - rl ), (contigLengths[cid2] - max_ins)

        # now check read 2
        if LI.pos2 <= ( max_ins - rl ):
            # read 2 lies at the start of its contig
            r2_at_start = True
            gap -= (LI.pos2 + rl)
        elif LI.pos2 >= (contigLengths[cid2] - max_ins ):
            # read 2 lies at the end of its contig
            r2_at_start = False
            gap -= (contigLengths[cid2] - LI.pos2)
        else:
            # read 2 lies in the middle of its contig
            return (0, LT.ERROR)
        #print gap, max_ins

        # now put it all together!
        # print r1_at_start, LI.pos1, LI.reversed1, "|", r2_at_start, LI.pos2, LI.reversed2, "|",
        if r1_at_start:
            # -x-1-->
            if LI.reversed1:
                # -<-1-->
                if r2_at_start:
                    # <--1->- -x-2-->
                    if LI.reversed2:
                        # <--1->- -<-2--> ==> IN + SS
                        if bamFile.orientationType == OT.IN:
                            #print "0 IN + SS"
                            return (gap, LT.SS)
                    else:
                        # <--1->- ->-2--> ==> SAME + SS
                        if bamFile.orientationType == OT.SAME:
                            #print "1 SAME + SS"
                            return (gap, LT.SS)
                else:
                    # <--1->- <x-2---
                    if LI.reversed2:
                        # <--1->- <>-2--- ==> SAME + SE
                        if bamFile.orientationType == OT.SAME:
                            #print "2 SAME + SE"
                            return (gap, LT.SE)
                    else:
                        # <--1->- <<-2--- ==> IN + SE
                        if bamFile.orientationType == OT.IN:
                            #print "3 IN + SE"
                            return (gap, LT.SE)
            else: # r1 agrees
                # ->-1-->
                if r2_at_start:
                    # <--2-x- ->-1-->
                    if LI.reversed2:
                        # <--2->- ->-1--> ==> SAME + SS
                        if bamFile.orientationType == OT.SAME:
                            #print "4 SAME + SS"
                            return (gap, LT.SS)
                    else:
                        # <--2-<- ->-1--> ==> OUT + SS
                        if bamFile.orientationType == OT.OUT:
                            #print "5 OUT + SS"
                            return (gap, LT.SS)
                else:
                    # ---2-x> ->-1-->
                    if LI.reversed2:
                        # ---2-<> ->-1--> ==> OUT + SE
                        if bamFile.orientationType == OT.OUT:
                            #print "6 OUT + SE"
                            return (gap, LT.SE)
                    else:
                        # ---2->> ->-1--> ==> SAME + SE
                        if bamFile.orientationType == OT.SAME:
                            #print "7 SAME + SE"
                            return (gap, LT.SE)
        else: # r1 at end
            # ---1-x>
            if LI.reversed1:
                # ---1-<>
                if r2_at_start:
                    # ---1-<>- -x-2-->
                    if LI.reversed2:
                        # ---1-<> -<-2--> ==> SAME + ES
                        if bamFile.orientationType == OT.SAME:
                            #print "8 SAME + ES"
                            return (gap, LT.ES)
                    else:
                        # ---1-<> ->-2--> ==> OUT + ES
                        if bamFile.orientationType == OT.OUT:
                            #print "9 OUT + ES"
                            return (gap, LT.ES)
                else:
                    # ---1-<> <x-2---
                    if LI.reversed2:
                        # ---1-<> <>-2--- ==> OUT + EE
                        if bamFile.orientationType == OT.OUT:
                            #print "a OUT + EE"
                            return (gap, LT.EE)
                    else:
                        # ---1-<> <<-2--- ==> SAME + EE
                        if bamFile.orientationType == OT.SAME:
                            #print "b SAME + EE"
                            return (gap, LT.EE)
            else: # r1 agrees
                # ---1->>
                if r2_at_start:
                    # ---1->> -x-2-->
                    if LI.reversed2:
                        # ---1->> -<-2--> ==> IN + ES
                        if bamFile.orientationType == OT.IN:
                            #print "c IN + SS"
                            return (gap, LT.ES)
                    else:
                        # ---1->> ->-2--> ==> SAME + ES
                        if bamFile.orientationType == OT.SAME:
                            #print "d SAME + ES"
                            return (gap, LT.ES)
                else:
                    # ---1->> <x-2---
                    if LI.reversed2:
                        # ---1->> <>-2--- ==> SAME + EE
                        if bamFile.orientationType == OT.SAME:
                            #print "e SAME + EE"
                            return (gap, LT.EE)
                    else:
                        # ---1->> <<-2--- ==> IN + EE
                        if bamFile.orientationType == OT.IN:
                            #print "f IN + EE"
                            return (gap, LT.EE)
        #print
        return (0, LT.ERROR)

    def __str__(self):
        str = "  (%d, %d, %d links)\n" % (self.cid1, self.cid2, len(self.links))
        for link in self.links:
            str += "%s" % link
        return str

    def printMore(self, contigNames, bamFileNames):
        str = "  (%s, %s, %d links)\n" % (contigNames[self.cid1], contigNames[self.cid2], len(self.links))
        for link in self.links:
            str += link.printMore(bamFileNames)
        return str

# mapping results structure "C land"
"""
typedef struct {
    uint32_t ** plpBp;
    uint32_t * contigLengths;
    uint32_t ** contigLengthCorrectors;
    uint32_t numBams;
    uint32_t numContigs;
    char ** contigNames;
    char ** bamFileNames;
    int isLinks;
    char * coverage_mode;
    int isIgnoreSupps;
    cfuhash_table_t * links;
} BM_mappingResults;
"""
class BM_mappingResults_C(c.Structure):
    _fields_ = [("plpBp", c.POINTER(c.POINTER(c.c_uint32))),
                ("contigLengths",c.POINTER(c.c_uint32)),
                ("contigLengthCorrectors",c.POINTER(c.POINTER(c.c_uint32))),
                ("numBams",c.c_uint32),
                ("numContigs",c.c_uint32),
                ("contigNames",c.POINTER(c.POINTER(c.c_char))),
                ("bamFileNames",c.POINTER(c.POINTER(c.c_char))),
                ("contig_name_lengths",c.POINTER(c.c_uint16)),
                ("bam_file_name_lengths",c.POINTER(c.c_uint16)),
                ("isLinks",c.c_int),
                ("coverage_mode",c.POINTER(c.c_char)),
                ("isIgnoreSupps",c.c_int),
                ("links",c.POINTER(cfuhash_table_t))
                ]
# mapping results structure "Python land"
class BMM_mappingResults(object):
    def __init__(self,
                 coverages,
                 contigLengths,
                 numBams,
                 numContigs,
                 contigNames,
                 bamFileName,
                 links):
        self.coverages = np.array(coverages)
        self.contigLengths = np.array(contigLengths)
        self.numBams = numBams
        self.numContigs = numContigs
        self.contigNames = contigNames
        self.bamFileNames = [bamFileName]
        self.links = links

    def consume(self, MRb):
        """Merge MR's internals with this one"""
        tmp = np.zeros((self.numContigs, (self.numBams+1)))
        tmp[:,:-1] = self.coverages
        tmp[:,-1] = np.reshape(MRb.coverages, (1, self.numContigs))
        self.coverages = tmp
        self.numBams += 1
        self.bamFileNames.append(MRb.bamFileNames[0])
        for key in MRb.links.keys():
            try:
                (self.links[key]).links += (MRb.links[key]).links
            except KeyError:
                self.links[key] = MRb.links[key]

    def __str__(self):
        str = "Contig\tLength"
        for i in range(self.numBams):
            str += "\t%s" % self.bamFileNames[i]
        str += "\n"
        for i in range(self.numContigs):
            str += "%s\t%d" % (self.contigNames[i], self.contigLengths[i])
            for j in range(self.numBams):
                str += "\t%0.2f" % (self.coverages[i,j])
            str += "\n"
        #return str
        if len(self.links) != 0:
            for key in self.links.keys():
                str += self.links[key].printMore(self.contigNames, self.bamFileNames)

        return str

###############################################################################
###############################################################################
# Multiprocessing requires that all passed items be pickleable. That is they
# must be vanilla variables or functions defined in the file itself, ie. not
# within a class. We get around this by writing an external function which calls
# a class function. Hacky, but it works.
###############################################################################
###############################################################################

def pythonizeLinks(MR, bamFile, contigLengths):
    """Unwrap the links-associated C structs and return a python-ized dict"""
    links = {}
    CW = CWrapper()
    pMR = c.POINTER(BM_mappingResults_C)
    pMR = c.pointer(MR)

    LW = BMM_LinkWalker_C()
    pLW = c.POINTER(BMM_LinkWalker_C)
    pLW = c.pointer(LW)
    success = CW._initLW(pLW, pMR)
    if(success == 1):
        ret_val = 2
        LP = None
        while(ret_val != 0):
            if ret_val == 2:
                # need a new contig pair
                LP = BMM_linkPair(((LW.pair).contents).cid1, ((LW.pair).contents).cid2)
                key = "%d,%d" % (((LW.pair).contents).cid1, ((LW.pair).contents).cid2)
                links[key] = LP
            # add a link
            LI = (LW.LI).contents
            LP.addLink(LI.reversed1,
                       LI.reversed2,
                       LI.pos1,
                       LI.pos2,
                       readLength,
                       contigLengths,
                       bamFile)
            ret_val = CW._stepLW(pLW)
        CW._destroyLW(pLW)

    return links

def externalParseWrapper(bAMpARSER, bamFileName, bid, MR, doContigNames):
    """ctypes pointers are unpickleable -- what we need is a hack!

    See BamParser._parseOneBam for what this function should be doing
    """
    print "HHH"
    MR = bAMpARSER._parseOneBam(bamFileName)
    print "HHsH"


    contig_names = []
    contig_lengths = np.array([int(i) for i in c.cast(MR.contigLengths, c.POINTER(c.c_uint32*MR.numContigs)).contents])
    plpBp = np.array([[int(j) for j in c.cast(i, c.POINTER(c.c_uint32*MR.numBams)).contents] for i in c.cast(MR.plpBp,c.POINTER(c.POINTER(c.c_uint32*MR.numBams)*MR.numContigs)).contents])

    coverages = np.zeros((MR.numContigs, MR.numBams))
    if bAMpARSER.coverageMode == 'outlier':
        contig_length_correctors = np.array([[int(j) for j in c.cast(i, c.POINTER(c.c_uint32*MR.numBams)).contents] for i in c.cast(MR.contigLengthCorrectors,c.POINTER(c.POINTER(c.c_uint32*MR.numBams)*MR.numContigs)).contents])
        for c_idx in range(int(MR.numContigs)):
            for b_idx in range(int(MR.numBams)):
                coverages[c_idx,b_idx] = float(plpBp[c_idx,b_idx])/float(contig_lengths[c_idx] - contig_length_correctors[c_idx])
    else:
        for c_idx in range(MR.numContigs):
            for b_idx in range(MR.numBams):
                coverages[c_idx,b_idx] = float(plpBp[c_idx,b_idx])/float(contig_lengths[c_idx])

    if doContigNames:
        contig_name_lengths = np.array([int(i) for i in c.cast(MR.contig_name_lengths, c.POINTER(c.c_uint16*MR.numContigs)).contents])
        contig_name_array = c.cast(MR.contigNames, c.POINTER(c.POINTER(c.c_char)*MR.numContigs)).contents
        for i in range(MR.numContigs):
            contig_names.append("".join([j for j in c.cast(contig_name_array[i], c.POINTER(c.c_char*contig_name_lengths[i])).contents]))

    BF = BamFile(bid, bamFileName)
    if bAMpARSER.doLinks:
        links = pythonizeLinks(MR, BF, contig_lengths)
    else:
        links = {}

    MRR = BMM_mappingResults(coverages,
                             contig_lengths,
                             MR.numBams,
                             MR.numContigs,
                             contig_names,
                             BF,
                             links)
    MR.append(MRR)

    # we need to call some C on this guy
    pMR = c.POINTER(BM_mappingResults_C)
    pMR = c.pointer(MR)
    CW = CWrapper()
    CW._destroyMR(pMR)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class CWrapper:
    """Can't pickle cTypes pointers and functions. Use this CWrap C-Wrapper as a hack"""
    def __init__(self, numBams=1, numContigs=0):
        #---------------------------------
        # load the c library
        #---------------------------------
        package_dir, filename = os.path.split(__file__)
        package_dir = os.path.abspath(package_dir)
        c_lib = os.path.join(package_dir, 'c', 'bam', 'libBamM.a')
        self.libPMBam = c.cdll.LoadLibrary(c_lib)

        #---------------------------------
        # import C functions
        #---------------------------------
        self._mergeMR = self.libPMBam.mergeMRs
        """
        @abstract Merge the contents of MR_B into MR_A

        @param  MR_A  mapping results struct to copy to
        @param  MR_B  mapping results struct to copy from
        @return void
        @discussion MR_B remains unchanged.
        MR_A is updated to include all the info contained in MR_B

        void mergeMRs(BM_mappingResults_C * MR_A, BM_mappingResults_C * MR_B);
        """

        self._destroyMR = self.libPMBam.destroyMR
        """
        @abstract Free all the memory calloced in initMR

        @param  MR  mapping results struct to destroy
        @return void

        void destroyMR(BM_mappingResults_C * MR)
        """

        self._parseCoverageAndLinks = self.libPMBam.parseCoverageAndLinks
        """
        @abstract Initialise the mapping results struct <- read in the BAM files

        @param numBams               number of BAM files to parse
        @param baseQ                 base quality threshold
        @param mapQ                  mapping quality threshold
        @param minLen                min query length
        @param doLinks               1 if links should be calculated
        @param ignoreSuppAlignments  only use primary alignments
        @param coverageMode          type of coverage to be calculated
        @param bamFiles              filenames of BAM files to parse
        @param MR                    mapping results struct to write to
        @return 0 for success

        @discussion This function expects MR to be a null pointer. It calls
        initMR and stores info accordingly. TL;DR If you call this function
        then you MUST call destroyMR when you're done.

        int parseCoverageAndLinks(int numBams,
                                  int baseQ,
                                  int mapQ,
                                  int minLen,
                                  int doLinks,
                                  int ignoreSuppAlignments,
                                  char* coverageMode,
                                  char* bamFiles[],
                                  BM_mappingResults_C * MR
                                 )
        """

        self._adjustPlpBp = self.libPMBam.adjustPlpBp
        """
        @abstract Adjust (reduce) the number of piled-up bases along a contig

        @param  MR  mapping results struct to write to
        @param  position_holder  array of pileup depths
        @param  tid  contig currently being processed
        @return void

        @discussion This function expects MR to be initialised.
        it can change the values of contigLengthCorrectors and plpBp

        void adjustPlpBp(BM_mappingResults_C * MR,
                         uint32_t ** position_holder,
                         int tid)
        """

        self._calculateCoverages = self.libPMBam.calculateCoverages
        """
        @abstract Calculate the coverage for each contig for each BAM

        @param  MR  mapping results struct with mapping info
        @return matrix of floats (rows = contigs, cols = BAMs)

        @discussion This function expects MR to be initialised.
        NOTE: YOU are responsible for freeing the return value
        recommended method is to use destroyCoverages

        float ** calculateCoverages(BM_mappingResults_C * MR);
        """

        self._destroyCoverages = self.libPMBam.destroyCoverages
        """
        @abstract Destroy the coverages structure

        @param covs array to destroy
        @param numContigs number of rows in array
        @return void

        void destroyCoverages(float ** covs, int numContigs)
        """

        self._initLW = self.libPMBam.initLW
        """
        @abstract Start moving through all of the links

        @param  MR  mapping results struct containing links
        @return pointer to LinkHolder if links exist or NULL
        BMM_LinkWalker * initLW(BM_mappingResults * MR);
        """

        self._stepLW = self.libPMBam.stepLW
        """
        @abstract Move to the next LinkInfo or LinkPair

        @param  walker   pointer to LinkHolder.
        @return 1 for step within current contig pair, 2 for new pair, 0 for end walk

        int stepLW(BMM_LinkWalker * walker);
        """

        self._destroyLW = self.libPMBam.destroyLW
        """
        @abstract Start moving through all of the links

        @param  walker   pointer to LinkHolder.
        @return void

        void destroyLW(BMM_LinkWalker * walker);
        """

        self._printMR = self.libPMBam.printMR
        """
        @abstract Print the contents of the MR struct

        @param  MR   mapping results struct with mapping info

        void printMR(BM_mappingResults_C * MR)
        """

class BamParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self,
                 baseQuality,
                 minLength,
                 mappingQuality=0,
                 coverageMode='vanilla',
                 doLinks=False,
                 ignoreSuppAlignments=True
                 ):
        #---------------------------------
        # information about how the parser will be used
        #---------------------------------
        self.baseQuality = baseQuality
        self.mappingQuality = mappingQuality
        self.minLength = minLength

        if doLinks:
            self.doLinks = 1
        else:
            self.doLinks = 0

        if ignoreSuppAlignments:
            self.ignoreSuppAlignments = 1
        else:
            self.ignoreSuppAlignments = 0

        if coverageMode not in ['vanilla', 'outlier']:
             raise InvalidCoverageModeException("Unknown coverage mode '%s' supplied" % coverageMode)
        self.coverageMode = coverageMode

        #---------------------------------
        # internal variables
        #---------------------------------
        self.MR = None          # internal mapping results object

        LTInverter = {LT.SS:LT.SS,
                           LT.SE:LT.ES,
                           LT.ES:LT.SE,
                           LT.EE:LT.EE,
                           LT.ERROR:LT.ERROR}

        self.compl = s_maketrans('ACGT', 'TGCA')

#------------------------------------------------------------------------------
# Orientation stuffz

    def determineOTSameRef(self, ar1, ar2):
        """Determine the orientation type and insert size of two reads

        We assume that:
           1. both reads are on the same reference
           2. ar1 comes before ar2
        """
        isize = LI.pos2 - LI.pos1 + ar1.rlen
        if LI.reversed1:
            # <-1--
            if LI.reversed2:
                # <-1-- <-2--
                return (OT.SAME, isize)
            else:
                # <-1-- -2-->
                return (OT.OUT, isize)
        else:
            # -->
            if LI.reversed2:
                # --1-> <-2--
                return (OT.IN, isize)
            else:
                # --1-> --2->
                return (OT.SAME, isize)

    def revComp(self, seq):
        """Return the reverse complement of a sequence"""
        # build a dictionary to know what letter to switch to
        return seq.translate(self.compl)[::-1]

#------------------------------------------------------------------------------
# Bam parseratering

    def parseBams(self, bamFiles, numThreads=1):
        """Parse bam files to get coverage and linking reads

        stores results in internal mapping results list
        """
        global MR
        MR = Manager().list()
        pool = Pool(processes=numThreads)
        do_contig_names = True
        bid = 0
        for bamFile in bamFiles:
            pool.apply_async(func=externalParseWrapper, args=(self, bamFile, bid, MR, do_contig_names))
            bid += 1
            if do_contig_names:
                # we only need to parse the contig names once
                do_contig_names = False
        pool.close()
        pool.join()

        # all the MRs are made. Only one has the contig IDs. find it's index
        baseMR_index = 0
        for i in range(len(MR)):
            if len(MR[i].contigNames) > 0:
                baseMR_index = i
                break
        # merge all the separate mapping results
        self.MR = MR[baseMR_index]
        for i in range(len(MR)):
            if i != baseMR_index:
                self.MR.consume(MR[i])

    def _parseOneBam(self, bamFile):
        """Parse a single BAM file and append the result to the internal mapping results list"""
        MR = BM_mappingResults_C()        # destroy needs to be called on this -> it should be called by the calling function
        pMR = c.POINTER(BM_mappingResults_C)
        pMR = c.pointer(MR)
        bamfiles_c_array = (c.c_char_p * 1)()
        bamfiles_c_array[:] = [bamFile]
        CW = CWrapper()
        CW._parseCoverageAndLinks(1,
                                  self.baseQuality,
                                  self.mappingQuality,
                                  self.minLength,
                                  self.doLinks,
                                  self.ignoreSuppAlignments,
                                  c.create_string_buffer(self.coverageMode),
                                  bamfiles_c_array,
                                  pMR)
        return MR


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
