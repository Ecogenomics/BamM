#!/usr/bin/env python
###############################################################################
#                                                                             #
#    BamFiles.py                                                              #
#                                                                             #
#    Class for storing information about BAM files                            #
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
import ctypes as c
import numpy as np
import sys

# local imports
from bamLink import *

###############################################################################
###############################################################################
###############################################################################
###############################################################################
# BAM files "C land"
"""
typedef struct BM_bamType {
   int orientationType;
   float insertSize;
   float insertStdev;
   int supporting;
} BM_bamType;

typedef struct BM_bamFile {
   char * fileName;
   uint16_t fileNameLength;
   BM_bamType ** types;
   int numTypes;
} BM_bamFile;

"""
class BM_bamType_C(c.Structure):
    _fields_ = [("orientationType", c.c_int),
                ("insertSize", c.c_float),
                ("insertStdev", c.c_float),
                ("supporting", c.c_int)
                ]

class BM_bamFile_C(c.Structure):
    _fields_ = [("fileName", c.POINTER(c.c_char)),
                ("fileNameLength", c.c_uint16),
                ("types", c.POINTER(c.POINTER(BM_bamType_C))),
                ("numTypes", c.c_int)
                ]

# BAM files "Python land"
class BM_bamType(object):
    """Store information about the inserts and orientation types
    of a BAM file"""
    def __init__(self,
                 orientationType = OT.NONE,
                 insertSize = 0.,
                 insertStdev = 0.,
                 supporting = 0):

         self.orientationType = orientationType
         self.insertSize = insertSize
         self.insertStdev = insertStdev
         self.supporting  = supporting

    def __str__(self):
        return "%0.4f\t%0.4f\t%s\t%d" % (self.insertSize,
                                         self.insertStdev,
                                         OT2Str(self.orientationType),
                                         self.supporting)

class BM_bamFile(object):
    """Store information about a BAM file"""
    def __init__(self,
                 bid,
                 fileName):
         self.bid = bid
         self.fileName = fileName
         self.types = []

    def __str__(self):
        return "\n".join(["%s\t%s" % (self.fileName, type) for type in self.types])

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Bam info structure "C land"
"""
typedef struct BM_fileInfo {
    uint32_t ** plpBp;
    uint32_t * contigLengths;
    uint32_t ** contigLengthCorrectors;
    uint32_t numBams;
    uint32_t numContigs;
    BM_bamFile ** bamFiles;
    char ** contigNames;
    uint16_t * contigNameLengths;
    int isLinks;
    char * coverage_mode;
    int isIgnoreSupps;
    cfuhash_table_t * links;
} BM_fileInfo;

"""
class BM_fileInfo_C(c.Structure):
    _fields_ = [("plpBp", c.POINTER(c.POINTER(c.c_uint32))),
                ("contigLengths",c.POINTER(c.c_uint32)),
                ("contigLengthCorrectors",c.POINTER(c.POINTER(c.c_uint32))),
                ("numBams",c.c_uint32),
                ("numContigs",c.c_uint32),
                ("bamFiles",c.POINTER(c.POINTER(BM_bamFile_C))),
                ("contigNames",c.POINTER(c.POINTER(c.c_char))),
                ("contigNameLengths",c.POINTER(c.c_uint16)),
                ("isLinks",c.c_int),
                ("coverage_mode",c.POINTER(c.c_char)),
                ("isIgnoreSupps",c.c_int),
                ("links",c.POINTER(cfuhash_table_t))
                ]

# mapping results structure "Python land"
class BM_fileInfo(object):
    def __init__(self,
                 coverages,
                 contigLengths,
                 numBams,
                 numContigs,
                 contigNames,
                 bamFiles,
                 links):
        self.coverages = np.array(coverages)
        self.contigLengths = np.array(contigLengths)
        self.numBams = numBams
        self.numContigs = numContigs
        self.contigNames = contigNames
        self.bamFiles = bamFiles
        self.links = links

    def consume(self, BFIb):
        """Merge BFI's internals with this one"""
        if len(self.coverages) > 0:
            tmp = np.zeros((self.numContigs, (self.numBams+1)))
            tmp[:,:-1] = self.coverages
            tmp[:,-1] = np.reshape(BFIb.coverages, (1, self.numContigs))
            self.coverages = tmp

        self.numBams += BFIb.numBams
        for BF in BFIb.bamFiles:
            self.bamFiles.append(BF)

        for key in BFIb.links.keys():
            try:
                (self.links[key]).links += (BFIb.links[key]).links
            except KeyError:
                self.links[key] = BFIb.links[key]

    def printBamTypes(self, fileHandle=None):
        if fileHandle is None:
            fileHandle = sys.stdout
        # header
        fileHandle.write("#file\tinsert\tstdev\torientation\tsupporting\n")
        # values
        fileHandle.write("\n".join([str(BF) for BF in self.bamFiles]) + "\n")

    def printCoverages(self, fileHandle=None):
        """print the coverage profiles"""
        if fileHandle is None:
            fileHandle = sys.stdout
        # header
        fileHandle.write("#contig\tLength")
        for i in range(self.numBams):
            fileHandle.write("\t%s" % self.bamFiles[i].fileName)
        fileHandle.write("\n")
        # values
        fileHandle.write("\n".join(["%s\t%d\t" % (self.contigNames[i], self.contigLengths[i]) +
                                    "\t".join(["%0.4f" % self.coverages[i,j] for j in range(self.numBams)])
                                    for i in range(self.numContigs)]) + "\n")

    def printLinks(self, fileHandle=None):
        """print all the links"""
        if fileHandle is None:
            fileHandle = sys.stdout
        if len(self.links) != 0:
            # header
            fileHandle.write("%s\n" % "\t".join(["#cid_1","cid_2","len_1","pos_1","rev_1","len_2","pos_2","rev_2","file"]))
            # values
            fileHandle.write("\n".join([self.links[key].printMore(self.contigNames, self.contigLengths, [self.bamFiles[i].fileName for i in range(self.numBams)]) for key in self.links.keys()]) + "\n")

    def __str__(self):
        """print everything"""
        retstr = self.printBamTypes(returnStr = True) + "\n"
        retstr += self.printCoverages(returnStr = True)
        if len(self.links) != 0:
            retstr += self.printLinks(returnStr = True)
        return retstr

###############################################################################
###############################################################################
###############################################################################
###############################################################################
