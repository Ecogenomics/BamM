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
__copyright__ = "Copyright 2014,2015"
__credits__ = ["Michael Imelfort"]
__license__ = "LGPLv3"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"

###############################################################################

import os
import ctypes as c
from pkg_resources import working_set, resource_filename

from bamm.bammExceptions import printError

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# C-style enums FTW!
def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

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
# Managing paired and unparied reads (and relative ordering)
#
# NOTE: This RPI definition corresponds to the definition in bamRead.h
#
# FIR   means first in properly paired mapping
# SEC   means second "
# SNGLP means paired in BAM but unpaired in mapping
# SNGL  means unpaired in BAM
# ERROR just for fun
#
global RPI
RPI = enum('ERROR', 'FIR', 'SEC', 'SNGL_FIR', 'SNGL_SEC', 'SNGL')

def RPI2Str(rpi):
    '''Convert an RPI into a human readable string

    Inputs:
     rpi - RPI to convert

    Outputs:
     Human readable string
    '''
    if rpi == RPI.FIR:
        return 'First'
    elif rpi == RPI.SEC:
        return 'Second'
    elif rpi == RPI.SNGL_FIR:
        return 'First_single'
    elif rpi == RPI.SNGL_SEC:
        return 'Second_single'
    elif rpi == RPI.SNGL:
        return 'Single'
    return 'ERROR'

# RPI.SNGL_FIR needs to be written to the singles file etc...
RPIConv = {RPI.ERROR:RPI.ERROR,
           RPI.FIR:RPI.FIR,
           RPI.SEC:RPI.SEC,
           RPI.SNGL_FIR:RPI.SNGL,
           RPI.SNGL_SEC:RPI.SNGL,
           RPI.SNGL:RPI.SNGL}

###############################################################################
###############################################################################
###############################################################################
###############################################################################

global MI
MI = enum('ER_EM_EG', 'PR_PM_PG', 'PR_PM_UG', 'PR_UM_NG', 'UR_NM_NG')

def MI2Str(mi):
    '''Convert an MI into a human readable string

    Inputs:
     mi - MI to convert

    Outputs:
     Human readable string
    '''
    if mi == MI.PR_PG_PM:
        return 'PR_PG_PM'
    elif mi == MI.PR_PM_UG:
        return 'PR_PM_UG'
    elif mi == MI.PR_UM_NG:
        return 'PR_UM_NG'
    elif mi == MI.UR_NM_NG:
        return 'UR_NM_NG'
    return 'ER_EM_EG'

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
        uint8_t mi,
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
                ("mi", c.c_uint8),
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
#------------------------------------------------------------------------------
# Managing orientation and linking types
#
# NOTE: This OT definition corresponds to the definition in bamParser.h
#
# Read orientations
# type 0 <--- --->
# type 1 ---> --->
# type 2 ---> <---
global OT
OT = enum('OUT', 'SAME', 'IN', 'NONE', 'ERROR')

def OT2Str(ot):

    '''Convert an orientation type into a human readable string

    Inputs:
     ot - OT to convert

    Outputs:
     Human readable string
    '''
    if ot == OT.OUT:
        return 'OUT'
    if ot == OT.SAME:
        return 'SAME'
    if ot == OT.IN:
        return 'IN'
    if ot == OT.NONE:
        return 'NONE'
    return 'ERROR'

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
# Managing different ways to calculate coverage
#
# NOTE: This CT definition corresponds to the definition in coverageEstimators.h
#
# CT_NONE             do not calculate coverage
# CT_COUNT            read counts, unaffected by contig length
# CT_C_MEAN           read counts, divided by contig length
# CT_P_MEAN           mean pileup depth
# CT_P_MEDIAN         median pileup depth
# CT_P_MEAN_TRIMMED   pileup mean trancated based on upper lower %
# CT_P_MEAN_OUTLIER   pileup mean trancated based on distributions
# P_VARIANCE          Variance of pileup depth
# MAPPED_COUNT        reference positions mapped
# MAPPED_MEAN         percentage reference positions mapped
# MAPPED_MEAN_TRIMMED percentage reference positions mapped excluding upper lower %

global CT
CT = enum('NONE',
          'COUNT',
          'C_MEAN',
          'P_MEAN',
          'P_MEDIAN',
          'P_MEAN_TRIMMED',
          'P_MEAN_OUTLIER',
          'P_VARIANCE',
          'MAPPED_COUNT',
          'MAPPED_MEAN',
          'MAPPED_MEAN_TRIMMED')

def CT2Str(ct):
    '''Convert a CT into a human readable string

    Inputs:
     ct - CT to convert

    Outputs:
     Human readable string
    '''
    if ct == CT.NONE:
        return 'None'
    elif ct == CT.COUNT:
        return 'Read_counts'
    elif ct == CT.C_MEAN:
        return 'Mean_read_counts'
    elif ct == CT.P_MEAN:
        return 'Mean_pileup_depth'
    elif ct == CT.P_MEDIAN:
        return 'Median_pileup_depth'
    elif ct == CT.P_MEAN_TRIMMED:
        return 'Mean_trimmed_pileup_depth'
    elif ct == CT.P_MEAN_OUTLIER:
        return 'Mean_outlier_pileup_depth'
    elif ct == CT.P_VARIANCE:
        return 'Variance_pileup_depth'
    elif ct == CT.MAPPED_COUNT:
        return 'Length_mapped'
    elif ct == CT.MAPPED_COUNT_MEAN:
        return 'Percent_mapped'
    elif ct == CT.MAPPED_COUNT_MEAN_TRIMMED:
        return 'Trimmed_percent_mapped'

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BM_coverageType_C(c.Structure):
    '''
    typedef struct BM_coverageType {
        CT type;
        float upperCut;
        float lowerCut;
    } BM_coverageType;
    '''
    _fields_ = [("type", c.c_int),
                ("upperCut", c.c_float),
                ("lowerCut", c.c_float)
                ]

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BM_fileInfo_C(c.Structure):
    '''
    typedef struct BM_fileInfo {
        float ** coverages;
        uint32_t * contigLengths;
        uint32_t numBams;
        uint32_t numContigs;
        BM_bamFile ** bamFiles;
        char ** contigNames;
        uint16_t * contigNameLengths;
        int isLinks;
        BM_coverageType * coverageType;
        int isIgnoreSupps;
        cfuhash_table_t * links;
    } BM_fileInfo;
    '''
    _fields_ = [("coverages", c.POINTER(c.POINTER(c.c_float))),
                ("contigLengths",c.POINTER(c.c_uint32)),
                ("numBams",c.c_uint32),
                ("numContigs",c.c_uint32),
                ("bamFiles",c.POINTER(c.POINTER(BM_bamFile_C))),
                ("contigNames",c.POINTER(c.POINTER(c.c_char))),
                ("contigNameLengths",c.POINTER(c.c_uint16)),
                ("isLinks",c.c_int),
                ("coverageType",c.POINTER(BM_coverageType_C)),
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
    def __init__(self, UT=False):
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
        c_lib = os.path.abspath(resource_filename('bamm', 'libBamM.a'))
        if UT:
            # unit tests are run from within the install dir which confuses
            # pkg_resources as there is a folder there called bamm
            # this is a hack to get around that -> prefer resource_filename
            for dist_path, dist in working_set.entry_keys.items():
                for d in dist:
                    if d == 'bamm':
                        c_lib = os.path.join(dist_path, d, 'libBamM.a')

        try:
            self.libPMBam = c.cdll.LoadLibrary(c_lib)
        except OSError:
            printError("Problem importing the BamM c library. This typically " \
                  "means that BamM is not installed correctly.\nPlease check " \
                  "the installation logs for more details.\nIf you don't have "\
                  "the installation logs then please try to reinstall BamM " \
                  "and look at the output.\nLooking for the c library at: %s" %\
                  c_lib)
            raise

        #---------------------------------
        # import C functions
        #---------------------------------

        self._filterReads = self.libPMBam.filterReads
        self._filterReads.argtypes = [c.POINTER(c.c_char),
                                      c.POINTER(c.c_char),
                                      c.c_int,
                                      c.c_int,
                                      c.c_int,
                                      c.c_float,
                                      c.c_float,
                                      c.c_int,
                                      c.c_int,
                                      c.c_int]
        self._filterReads.restype = None

        #----------------

        self._profileReads = self.libPMBam.profileReads
        self._profileReads.argtypes = [c.POINTER(c.c_char),
                                       c.c_int,
                                       c.c_int]
        self._profileReads.restype = None

        #----------------

        
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
        self._setMICode = self.libPMBam.setMICode
        self._setMICode.argTypes = [c.POINTER(BM_mappedRead_C), c.c_int]
        self._setMICode.restype = None

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
        self._initLW = self.libPMBam.initLW

        #-----------------
        self._stepLW = self.libPMBam.stepLW

        #-----------------
        self._destroyLW = self.libPMBam.destroyLW

        #-----------------
        self._printBFI = self.libPMBam.printBFI

        #----------------- [Note: importing for unit tests]
        self._estimate_COUNT_Coverage = self.libPMBam.estimate_COUNT_Coverage
        self._estimate_COUNT_Coverage.argtypes = [c.POINTER(c.c_uint32), c.POINTER(BM_coverageType_C), c.c_uint32]
        self._estimate_COUNT_Coverage.restype = c.c_float

        self._estimate_C_MEAN_Coverage = self.libPMBam.estimate_C_MEAN_Coverage
        self._estimate_C_MEAN_Coverage.argtypes = [c.POINTER(c.c_uint32), c.POINTER(BM_coverageType_C), c.c_uint32]
        self._estimate_C_MEAN_Coverage.restype = c.c_float

        self._estimate_P_MEAN_Coverage = self.libPMBam.estimate_P_MEAN_Coverage
        self._estimate_P_MEAN_Coverage.argtypes = [c.POINTER(c.c_uint32), c.POINTER(BM_coverageType_C), c.c_uint32]
        self._estimate_P_MEAN_Coverage.restype = c.c_float

        self._estimate_P_MEDIAN_Coverage = self.libPMBam.estimate_P_MEDIAN_Coverage
        self._estimate_P_MEDIAN_Coverage.argtypes = [c.POINTER(c.c_uint32), c.POINTER(BM_coverageType_C), c.c_uint32]
        self._estimate_P_MEDIAN_Coverage.restype = c.c_float

        self._estimate_P_MEAN_TRIMMED_Coverage = self.libPMBam.estimate_P_MEAN_TRIMMED_Coverage
        self._estimate_P_MEAN_TRIMMED_Coverage.argtypes = [c.POINTER(c.c_uint32), c.POINTER(BM_coverageType_C), c.c_uint32]
        self._estimate_P_MEAN_TRIMMED_Coverage.restype = c.c_float

        self._estimate_P_MEAN_OUTLIER_Coverage = self.libPMBam.estimate_P_MEAN_OUTLIER_Coverage
        self._estimate_P_MEAN_OUTLIER_Coverage.argtypes = [c.POINTER(c.c_uint32), c.POINTER(BM_coverageType_C), c.c_uint32]
        self._estimate_P_MEAN_OUTLIER_Coverage.restype = c.c_float
        
        self._estimate_P_VARIANCE_Coverage = self.libPMBam.estimate_P_VARIANCE_Coverage
        self._estimate_P_VARIANCE_Coverage.argtypes = [c.POINTER(c.c_uint32), c.POINTER(BM_coverageType_C), c.c_uint32]
        self._estimate_P_VARIANCE_Coverage.restype = c.c_float
        """        
        self._estimate_MAPPED_COUNT_Coverage = self.libPMBam.estimate_MAPPED_COUNT_Coverage
        self._estimate_MAPPED_COUNT_Coverage.argtypes = [c.POINTER(c.c_uint32), c.POINTER(BM_coverageType_C), c.c_uint32]
        self._estimate_MAPPED_COUNT_Coverage.restype = c.c_float
        
        self._estimate_MAPPED_MEAN_Coverage = self.libPMBam.estimate_MAPPED_MEAN_Coverage
        self._estimate_MAPPED_MEAN_Coverage.argtypes = [c.POINTER(c.c_uint32), c.POINTER(BM_coverageType_C), c.c_uint32]
        self._estimate_MAPPED_MEAN_Coverage.restype = c.c_float
        
        self._estimate_MAPPED_MEAN_TRIMMED_Coverage = self.libPMBam.estimate_MAPPED_MEAN_TRIMMED_Coverage
        self._estimate_MAPPED_MEAN_TRIMMED_Coverage.argtypes = [c.POINTER(c.c_uint32), c.POINTER(BM_coverageType_C), c.c_uint32]
        self._estimate_MAPPED_MEAN_TRIMMED_Coverage.restype = c.c_float
        """
        self._BM_median = self.libPMBam.BM_median
        self._BM_median.argtypes = [c.POINTER(c.c_uint32), c.c_uint32]
        self._BM_median.restype = c.c_float

        self._BM_mean = self.libPMBam.BM_mean
        self._BM_mean.argtypes = [c.POINTER(c.c_uint32), c.c_uint32]
        self._BM_mean.restype = c.c_float

        self._BM_stdDev = self.libPMBam.BM_stdDev
        self._BM_stdDev.argtypes = [c.POINTER(c.c_uint32), c.c_uint32, c.c_float]
        self._BM_stdDev.restype = c.c_float



###############################################################################
###############################################################################
###############################################################################
###############################################################################
