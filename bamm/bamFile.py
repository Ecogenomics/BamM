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
import ctypes as c
import numpy as np
import sys

# local imports
from bamLink import BM_linkPair, BM_linkInfo
from cWrapper import OT2Str, OT

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Coverage types "Python Land"
class BM_coverageType(object):
    '''Container class for storing the type of coverage to calculate'''

    def __init__(self,
                 cType,
                 cUpper,
                 cLower):
        '''
        Default constructor.

        Inputs:
         cType -  enum, from the CT enum (see cWrapper.py)
         cUpper - float, percent of coverage values or number of stdevs. used
                  to determine upper cutoff when calculating coverage values
         cLower - float, percent of coverage values or number of stdevs. used
                  to determine lower cutoff when calculating coverage values

        Outputs:
         None

        '''
        self.cType = cType
        self.cUpper = cUpper
        self.cLower = cLower

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# BAM files "Python land"
class BM_bamType(object):
    ''' Container class for Storing information about the inserts
    and orientation types of a BAM file'''

    def __init__(self,
                 orientationType = OT.NONE,
                 insertSize = 0.,
                 insertStdev = 0.,
                 supporting = 0):
        '''
        Default constructor.

        Inputs:
         orientationType - enum, from the OT enum (see cWrapper.py)
         insertSize - float, estimated mean template size for the library
         insertStdev - float, estimated stdev of the template size
         supporting - int, number or reads used to define stats

        Outputs:
         None

        '''
        self.orientationType = orientationType
        self.insertSize = insertSize
        self.insertStdev = insertStdev
        self.supporting  = supporting

    def __str__(self):
        '''override basic string function'''
        return "%0.4f\t%0.4f\t%s\t%d" % (self.insertSize,
                                         self.insertStdev,
                                         OT2Str(self.orientationType),
                                         self.supporting)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BM_bamFile(object):
    '''Container class for storing information about a BAM file'''
    def __init__(self,
                 bid,
                 fileName):
        '''
        Default constructor.

        Inputs:
         bid - unique identifier for the bameFile
         fileName - long fileName

        Outputs:
         None
        '''
        self.bid = bid
        self.fileName = fileName
        self.types = []

    def __str__(self):
        '''override basic string function'''
        return "\n".join(["%s\t%s" % (self.fileName, type)
                          for type in self.types])

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BM_fileInfo(object):
    '''High level management of BamFiles, contigs and associated
    mapping statistics'''

    def __init__(self,
                 coverages,
                 contigLengths,
                 numBams,
                 numContigs,
                 contigNames,
                 bamFiles,
                 links):
        '''
        Default constructor.

        Inputs:
         coverages - (C x B) array of floats; C = num contigs and B = num bams
         contigLengths -(C) array of ints storing contig lengths
         numBams - The number of BAMs the instance is managing ( == B )
         numContigs - The number of unique contigs across all BAMs ( == C )
         contigNames - (C) array of strings storing long contig names
         bamFiles - (B) array of BM_BamFile instances
         links - { uid, BM_linkPair }, all links joining the contigs, uid made
                 by BM_BM_linkPair.makeKey()

        Outputs:
         None
        '''
        self.coverages = np.array(coverages)
        self.contigLengths = np.array(contigLengths)
        self.numBams = numBams
        self.numContigs = numContigs
        self.contigNames = contigNames
        self.bamFiles = bamFiles
        self.links = links

    def consume(self, BFIb):
        '''Merge another BFI's internals with this one

        You should destory the merged BFI after merging to avoid
        duplication of results

        Inputs:
         BFIb, BM_fileInfo, that will be merged

        Outputs:
         None
        '''
        if len(self.coverages) > 0:
            tmp = np.zeros((self.numContigs, (self.numBams+1)))
            tmp[:,:-1] = self.coverages
            try:
                tmp[:,-1] = np.reshape(BFIb.coverages, (1, self.numContigs))
            except ValueError:
                print "Error combining results from different BAMs. Are you " \
                      "sure they're from the same mapping?"
                raise
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
        ''' Write bam type information to the given filehandle or stdout

        used for printing output to types file

        Inputs:
         fileHandle - open FILE, == None implies writing to stdout

        Outputs:
         None
        '''
        if fileHandle is None:
            fileHandle = sys.stdout
        # header
        fileHandle.write("%s\n" % "\t".join(["#file",
                                             "insert",
                                             "stdev",
                                             "orientation",
                                             "supporting"])
                         )
        # values
        fileHandle.write("%s\n" % "\n".join([str(BF) for BF in self.bamFiles]))

    def printCoverages(self, fileHandle=None):
        '''Write coverage profile information to the given filehandle or stdout

        used for printing output to types file

        Inputs:
         fileHandle - open FILE, == None implies writing to stdout

        Outputs:
         None
        '''
        if fileHandle is None:
            fileHandle = sys.stdout
        # header
        fileHandle.write("#contig\tLength\t%s\n" % \
            "\t".join([self.bamFiles[i].fileName for i in range(self.numBams)]))
        # values
        fileHandle.write("%s\n" % "\n".join(["%s\t%d\t" % \
                                             (self.contigNames[i],
                                              self.contigLengths[i]) +
                                             "\t".join(["%0.4f" % self.coverages[i,j]
                                                        for j in range(self.numBams)])
                                             for i in range(self.numContigs)]))

    def printLinks(self, bamFileNames, fileHandle=None):
        '''Write contig linking information to the given filehandle or stdout

        used for printing output to types file

        Inputs:
         bamFileNames - { bamId : string }, storage for long bam file names
         fileHandle - open FILE, == None implies writing to stdout

        Outputs:
         None
        '''
        if fileHandle is None:
            fileHandle = sys.stdout
        if len(self.links) != 0:
            # header
            fileHandle.write("%s\n" % "\t".join(["#cid_1",
                                                 "cid_2",
                                                 "len_1",
                                                 "pos_1",
                                                 "rev_1",
                                                 "len_2",
                                                 "pos_2",
                                                 "rev_2",
                                                 "file"]))
            # values
            fileHandle.write("%s\n" % \
                "\n".join([self.links[key].printMore(self.contigNames,
                                                     self.contigLengths,
                                                     bamFileNames)
                           for key in self.links.keys()]))

    def __str__(self):
        '''Override basic string function

        calls printBamTypes, printCoverages and printLinks
        '''
        retstr = self.printBamTypes(returnStr = True) + "\n"
        retstr += self.printCoverages(returnStr = True)
        if len(self.links) != 0:
            retstr += self.printLinks(returnStr = True)
        return retstr

###############################################################################
###############################################################################
###############################################################################
###############################################################################
