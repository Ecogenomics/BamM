#!/usr/bin/env python
###############################################################################
#                                                                             #
#    bamRead.py                                                               #
#                                                                             #
#    Class for storing information about mapped reads                         #
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
__copyright__ = "Copyright 2014"
__credits__ = ["Michael Imelfort"]
__license__ = "LGPLv3"
__version__ = "0.1.0"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Beta"

###############################################################################

# system imports
import ctypes as c
import sys
import multiprocessing as mp
import re
import gzip
import Queue
import random

# local import
from bamm.bammExceptions import *

###############################################################################
###############################################################################
###############################################################################
###############################################################################

#
# This re is used to strip pair information from a read
#
pairStripper = re.compile( '(_1$|_2$|/1$|/2$)' )
metaStripper = re.compile( '(.*;r_)' )

#------------------------------------------------------------------------------
# Managing paired and unparied reads (and relative ordering)

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

global RPI
RPI = enum('FIR', 'SEC', 'SNGL', 'ERROR')

def RPI2Str(rpi):
    """For the humans!"""
    if rpi == RPI.FIR:
        return 'First'
    if rpi == RPI.SEC:
        return 'Second'
    if rpi == RPI.SNGL:
        return 'Single'
    return 'ERROR'

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# mapped read structure "C land"
"""
typedef struct BM_mappedRead {
    char * seqId,
    char * seq,
    char * qual,
    uint16_t idLen,
    uint16_t seqLen,
    uint16_t qualLen,
    uint8_t rpi,
    uint16_t bin,
    BM_mappedRead * prev_MR
} BM_mappedRead;

"""
class BM_mappedRead_C(c.Structure):
    pass

class BM_mappedRead_C(c.Structure):
    _fields_ = [("seqId", c.POINTER(c.c_char)),
                ("seq", c.POINTER(c.c_char)),
                ("qual", c.POINTER(c.c_char)),
                ("idLen", c.c_uint16),
                ("seqLen", c.c_uint16),
                ("qualLen", c.c_uint16),
                ("rpi", c.c_uint8),
                ("bin", c.c_uint16),
                ("nextRead",c.POINTER(BM_mappedRead_C))
                ]

# mapped read structure "Python land"
class BM_mappedRead(object):
    def __init__(self,
                 seqId,
                 seq=None,
                 qual=None,
                 rpi=RPI.SNGL,
                 bin=0
                 ):
        self.seqId = seqId
        self.seq = seq
        self.qual=qual
        self.rpi = rpi
        self.bin = bin

    def getUniversalId(self):
        """Strips off any pesky _1 _2 /1 /2"""
        return metaStripper.sub('', pairStripper.sub('', self.seqId))

    def print_MR(self, targetNames, fileHandle=None):
        """Print this read"""
        if fileHandle is None:
            fileHandle = sys.stdout

        if self.seq is None:
            # headers only
            fileHandle.write("b_%s;%s\n" %(targetNames[self.bin], self.seqId))
        else:
            if self.qual is not None:
                fileHandle.write("@b_%s;%s\n%s\n+\n%s\n" %(targetNames[self.bin], self.seqId, self.seq, self.qual))
            else:
                fileHandle.write(">b_%s;%s\n%s\n" % (targetNames[self.bin], self.seqId, self.seq))

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ReadSet(object):
    """Container class to encapsulate the concept of a(n un)paired read set"""
    def __init__(self,
                 fileName1,                     # always need at least one file name
                 fileName2 = None,              # second filename for paired non-shuffled
                 paired=False,                  # paired flag indicates how we should deplete the queue
                 zipped=True                    # gzip output
                 ):
        self.unopened = True                    # is this file unopened?
        self.fileName1 = fileName1              # the file to write to
        self.fileName2 = fileName2              # the other file to write to... ...perhaps
        self.buffer = mp.Manager().dict()       # thread safe dictionary to act as a buffer
        self.reads = mp.Manager().Queue()       # the reads that we can write to the file
        self.isPaired = paired                  # pop off reads two at a time?

        # work out how we'll write files
        self.zipOutput = zipped                 # totes
        if zipped:
            self.writeOpen = gzip.open
        else:
            self.writeOpen = open

        # we need this to throw an attribute error
        # on the first time it's read -> ie. leave undefined
        #self.isFasta = #$@!:

    def add(self, BMM):
        """Add a mapped read to the queue, respect read pairings"""

        # first check that we are adding fasta to fasta and fastq to fastq
        try:
            if self.isFasta ^ (BMM.qual is None):
                raise MixedFileTypesException("You cannot mix Fasta and Fastq reads together in an output file")
        except AttributeError:
            # we have not defined self.isFasta in the __init__
            # Now we can set the type of the file and we should only get
            # here on the first read
            self.isFasta = BMM.qual is None
            if self.isFasta:
                ext = ".fna"
            else:
                ext = ".fq"

            if self.zipOutput:
                ext += ".gz"
            self.fileName1 += ext
            if self.fileName2 is not None:
                 self.fileName2 += ext

            print "Here", self.fileName1, self.fileName2

        return
        # if the rpi is unpaired then we don't need to worry about the
        # order of the reads. Otherwise we should be sure there is a paired read coming
        # at some stage (or is here now) and these must be placed sequentially
        # onto the queue
        if BMM.rpi == RPI.SNGL:
            self.reads.put(BMM)
        else:
            uid = BMM.getUniversalId()
            try:
                stored_BMM = self.buffer[uid]
                # put read 1 on first
                if BMM.rpi == RPI.FIR:
                    self.reads.put(BMM)
                    self.reads.put(stored_BMM)
                else:
                    self.reads.put(stored_BMM)
                    self.reads.put(BMM)

                # free this memory
                del self.buffer[uid]

            except KeyError:
                # first on the scene, store it here
                self.buffer[uid] = BMM

    def getSize(self):
        """return the size of the read queue and buffer"""
        return(len(self.buffer), self.reads.qsize())

    def write(self, targetNames, maxDump=None):
        """Write a chunk of reads to file"""
        qs = self.reads.qsize()
        if maxDump is None:
            maxDump = qs
        elif maxDump < qs:  # don't want to read off the end of the queue
            maxDump = qs

        if self.isPaired:
            # enforce always writing even numbers of reads
            if maxDump % 2 != 0:
                maxDump -= 1

            # open files
            fh1 = self.writeOpen(self.fileName1, 'a')
            if self.fileName2 is None:
                print "Writing shuffled file:", self.fileName1
                fh2 = fh1
            else:
                print "Writing coupled files: ", self.fileName1, self.fileName2
                fh2 = self.writeOpen(self.fileName2, 'a')
            # write
            while maxDump > 0:
                print maxDump
                print "going to get A"
                BMM = self.reads.get(block=True, timeout=None)
                print "got"
                BMM.print_MR(targetNames, fh1)
                print "going to get B"
                BMM = self.reads.get(block=True, timeout=None)
                print "got"
                BMM.print_MR(targetNames, fh2)
                maxDump -= 2

            # and close
            fh1.close()
            if self.fileName2 is not None:
                fh2.close()
        else:
            with self.writeOpen(self.fileName1, "a") as fh:
                while maxDump > 0:
                    print maxDump
                    print "going to get C"
                    BMM = self.reads.get(block=True, timeout=None)
                    print "got"
                    BMM.print_MR(targetNames, fh)
                    maxDump -= 1

###############################################################################
###############################################################################
###############################################################################
###############################################################################


