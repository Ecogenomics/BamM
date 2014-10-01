#!/usr/bin/env python
###############################################################################
#                                                                             #
#    CWrapper.py                                                              #
#                                                                             #
#    Class for exposing the basic c functions                                 #
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

import os
import ctypes as c

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# fields defined in cfuhash.c but not accessed at this level
class cfuhash_table_t(c.Structure):
    pass

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BM_mappedRead_C(c.Structure):
    pass

class BM_mappedRead_C(c.Structure):
    '''
    typedef struct BM_mappedRead {
        char * seqId,
        char * seq,
        char * qual,
        uint16_t idLen,
        uint16_t seqLen,
        uint16_t qualLen,
        uint8_t rpi,
        uint16_t group,
        BM_mappedRead * nextRead,
        BM_mappedRead * partnerRead,
        BM_mappedRead * nextPrintingRead
    } BM_mappedRead;
    '''
    _fields_ = [("seqId", c.POINTER(c.c_char)),
                ("seq", c.POINTER(c.c_char)),
                ("qual", c.POINTER(c.c_char)),
                ("idLen", c.c_uint16),
                ("seqLen", c.c_uint16),
                ("qualLen", c.c_uint16),
                ("rpi", c.c_uint8),
                ("group", c.c_uint16),
                ("nextRead",c.POINTER(BM_mappedRead_C)),
                ("partnerRead",c.POINTER(BM_mappedRead_C)),
                ("nextPrintingRead",c.POINTER(BM_mappedRead_C))
                ]

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BM_linkInfo_C(c.Structure):
    pass

class BM_linkInfo_C(c.Structure):
    '''
    typedef struct {
        uint16_t reversed1;
        uint16_t reversed2;
        uint16_t readLength1;
        uint16_t readLength2;
        uint32_t pos1;
        uint32_t pos2;
        uint32_t bam_ID;
        struct BM_linkInfo * nextLink;
    } BM_linkInfo;
    '''
    _fields_ = [("reversed1", c.c_uint16),
                ("reversed2", c.c_uint16),
                ("readLength1", c.c_uint16),
                ("readLength2", c.c_uint16),
                ("pos1", c.c_uint32),
                ("pos2", c.c_uint32),
                ("bam_ID", c.c_uint32),
                ("nextLink",c.POINTER(BM_linkInfo_C))
                ]

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BM_linkPair_C(c.Structure):
    '''
    typedef struct {
        uint32_t cid1;
        uint32_t cid2;
        uint32_t numLinks;
        BM_linkInfo * LI;
    } BM_linkPair;
    '''
    _fields_ = [("cid1", c.c_uint32),
                ("cid2", c.c_uint32),
                ("numLinks", c.c_uint32),
                ("LI",c.POINTER(BM_linkInfo_C))
                ]

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BM_LinkWalker_C(c.Structure):
    '''
    typedef struct {
        char ** keys;
        size_t keyCount;
        size_t numKeys;
        cfuhash_table_t * linkHash;
        BM_linkPair * pair;
        BM_linkInfo * LI;
    } BM_LinkWalker;
    '''
    _fields_= [("keys", c.POINTER(c.POINTER(c.c_char))),
               ("keyCount", c.c_size_t),
               ("numKeys", c.c_size_t),
               ("links",c.POINTER(cfuhash_table_t)),
               ("pair",c.POINTER(BM_linkPair_C)),
               ("LI",c.POINTER(BM_linkInfo_C))
               ]

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BM_bamType_C(c.Structure):
    '''
    typedef struct BM_bamType {
       int orientationType;
       float insertSize;
       float insertStdev;
       int supporting;
    } BM_bamType;
    '''
    _fields_ = [("orientationType", c.c_int),
                ("insertSize", c.c_float),
                ("insertStdev", c.c_float),
                ("supporting", c.c_int)
                ]

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BM_bamFile_C(c.Structure):
    '''
    typedef struct BM_bamFile {
       char * fileName;
       uint16_t fileNameLength;
       BM_bamType ** types;
       int numTypes;
    } BM_bamFile;
    '''
    _fields_ = [("fileName", c.POINTER(c.c_char)),
                ("fileNameLength", c.c_uint16),
                ("types", c.POINTER(c.POINTER(BM_bamType_C))),
                ("numTypes", c.c_int)
                ]

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BM_fileInfo_C(c.Structure):
    '''
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
    '''
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

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class CWrapper:
    ''' Multiprocessing can't pickle cTypes pointers and functions.
    This class is a hack which quarantines the cTypes functions.
    '''
    def __init__(self):
        '''Default constructor.

        Loads libBamM.a and instantiates wrappers to the functions we wish
        to export

        Inputs:
         None

        Outputs;
         None
        '''
        #---------------------------------
        # load the c library
        #---------------------------------
        package_dir, filename = os.path.split(__file__)
        package_dir = os.path.abspath(package_dir)
        c_lib = os.path.join(package_dir, 'c', 'libBamM.a')
        self.libPMBam = c.cdll.LoadLibrary(c_lib)

        #---------------------------------
        # import C functions
        #---------------------------------

        self._mergeBFI = self.libPMBam.mergeBFIs

        #-----------------
        self._destroyBFI = self.libPMBam.destroyBFI

        #-----------------
        self._extractReads = self.libPMBam.extractReads
        self._extractReads.argtypes = [c.POINTER(c.c_char),
                                       c.POINTER(c.c_char_p),
                                       c.c_int,
                                       c.POINTER(c.c_uint16),
                                       c.POINTER(c.c_char),
                                       c.c_int,
                                       c.c_int,
                                       c.c_int,
                                       c.c_int,
                                       c.c_int]
        self._extractReads.restype = c.POINTER(BM_mappedRead_C)

        #-----------------
        self._getNextMappedRead = self.libPMBam.getNextMappedRead
        self._getNextMappedRead.argtypes = [c.POINTER(BM_mappedRead_C)]
        self._getNextMappedRead.restype = c.POINTER(BM_mappedRead_C)

        #-----------------
        self._setNextPrintRead = self.libPMBam.setNextPrintRead
        self._setNextPrintRead.argtypes = [c.POINTER(BM_mappedRead_C),
                                           c.POINTER(BM_mappedRead_C)]
        self._setNextPrintRead.restype = None

        #-----------------
        self._getNextPrintRead = self.libPMBam.getNextPrintRead
        self._getNextPrintRead.argtypes = [c.POINTER(BM_mappedRead_C)]
        self._getNextPrintRead.restype = c.POINTER(BM_mappedRead_C)

        #-----------------
        self._getPartner = self.libPMBam.getPartner
        self._getPartner.argtypes = [c.POINTER(BM_mappedRead_C)]
        self._getPartner.restype = c.POINTER(BM_mappedRead_C)

        #-----------------
        self._partnerInSameGroup = self.libPMBam.partnerInSameGroup

        #-----------------
        self._destroyMappedReads = self.libPMBam.destroyMappedReads
        self._destroyMappedReads.argtypes = [c.POINTER(BM_mappedRead_C)]
        self._destroyMappedReads.restype = None

        #-----------------
        self._destroyPrintChain = self.libPMBam.destroyPrintChain
        self._destroyPrintChain.argtypes = [c.POINTER(BM_mappedRead_C)]
        self._destroyPrintChain.restype = None

        #-----------------
        self._printMappedRead = self.libPMBam.printMappedRead
        self._sprintMappedRead = self.libPMBam.sprintMappedRead
        self._printMappedReads = self.libPMBam.printMappedReads

        #-----------------
        self._parseCoverageAndLinks = self.libPMBam.parseCoverageAndLinks

        #-----------------
        self._adjustPlpBp = self.libPMBam.adjustPlpBp

        #-----------------
        self._calculateCoverages = self.libPMBam.calculateCoverages

        #-----------------
        self._destroyCoverages = self.libPMBam.destroyCoverages

        #-----------------
        self._initLW = self.libPMBam.initLW

        #-----------------
        self._stepLW = self.libPMBam.stepLW

        #-----------------
        self._destroyLW = self.libPMBam.destroyLW

        #-----------------
        self._printBFI = self.libPMBam.printBFI

###############################################################################
###############################################################################
###############################################################################
###############################################################################
