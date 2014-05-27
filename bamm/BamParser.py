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
    uint16_t orient_1;
    uint16_t orient_2;
    uint32_t pos_1;
    uint32_t pos_2;
    uint32_t bam_ID;
    struct BMM_link_info * next_link;
} BMM_link_info;

typedef struct {
    uint32_t cid_1;
    uint32_t cid_2;
    uint32_t numLinks;
    BMM_link_info * LI;
} BMM_link_pair;

typedef struct {
    char ** keys;
    size_t keyCount;
    size_t numKeys;
    cfuhash_table_t * linkHash;
    BMM_link_pair * pair;
    BMM_link_info * LI;
} BMM_LinkWalker;

"""
class BMM_link_info_C(c.Structure):
    pass

class BMM_link_info_C(c.Structure):
    _fields_ = [("orient_1", c.c_uint16),
                ("orient_2", c.c_uint16),
                ("pos_1", c.c_uint32),
                ("pos_2", c.c_uint32),
                ("bam_ID", c.c_uint32),
                ("next_link",c.POINTER(BMM_link_info_C))
                ]

class BMM_link_pair_C(c.Structure):
    _fields_ = [("cid_1", c.c_uint32),
                ("cid_2", c.c_uint32),
                ("numLinks", c.c_uint32),
                ("LI",c.POINTER(BMM_link_info_C))
                ]

class BMM_LinkWalker_C(c.Structure):
    _fields_= [("keys", c.POINTER(c.POINTER(c.c_char))),
               ("keyCount", c.c_size_t),
               ("numKeys", c.c_size_t),
               ("links",c.POINTER(cfuhash_table_t)),
               ("pair",c.POINTER(BMM_link_pair_C)),
               ("LI",c.POINTER(BMM_link_info_C))
               ]

# links-associated structures "Python land"
class BMM_linkInfo(object):
    def __init__(self,
                 o1,
                 o2,
                 p1,
                 p2,
                 bid):
        self.orient1 = o1
        self.orient2 = o2
        self.pos1 = p1
        self.pos2 = p2
        self.bamID = bid

    def __str__(self):
        return "    (%d,%d -> %d,%d, %s)\n" % (self.pos1, self.orient1, self.pos2, self.orient2, self.bamID)

    def printMore(self, bamFileNames):
        return "    (%d,%d -> %d,%d, %s)\n" % (self.pos1, self.orient1, self.pos2, self.orient2, bamFileNames[self.bamID])

class BMM_linkPair(object):
    def __init__(self,
                 cid1,
                 cid2):
        self.cid1 = cid1
        self.cid2 = cid2
        self.numLinks = 0
        self.links = []

    def addLink(self,
                o1,
                o2,
                p1,
                p2,
                bid):
        LI = BMM_linkInfo(o1,
                         o2,
                         p1,
                         p2,
                         bid)
        self.links.append(LI)
        self.numLinks += 1

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
    uint32_t ** plp_bp;
    uint32_t * contig_lengths;
    uint32_t ** contig_length_correctors;
    uint32_t num_bams;
    uint32_t num_contigs;
    char ** contig_names;
    char ** bam_file_names;
    int is_links_included;
    int is_outlier_coverage;
    int is_ignore_supps;
    cfuhash_table_t * links;
} BMM_mapping_results;
"""
class BMM_mapping_results_C(c.Structure):
    _fields_ = [("plp_bp", c.POINTER(c.POINTER(c.c_uint32))),
                ("contig_lengths",c.POINTER(c.c_uint32)),
                ("contig_length_correctors",c.POINTER(c.POINTER(c.c_uint32))),
                ("num_bams",c.c_uint32),
                ("num_contigs",c.c_uint32),
                ("contig_names",c.POINTER(c.POINTER(c.c_char))),
                ("bam_file_names",c.POINTER(c.POINTER(c.c_char))),
                ("contig_name_lengths",c.POINTER(c.c_uint16)),
                ("bam_file_name_lengths",c.POINTER(c.c_uint16)),
                ("is_links_included",c.c_int),
                ("is_outlier_coverage",c.c_int),
                ("is_ignore_supps",c.c_int),
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

def pythonizeLinks(MR):
    """Unwrap the links-associated C structs and return a python-ized dict"""
    links = {}
    CW = CWrapper()
    pMR = c.POINTER(BMM_mapping_results_C)
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
                LP = BMM_linkPair(((LW.pair).contents).cid_1, ((LW.pair).contents).cid_2)
                key = "%d,%d" % (((LW.pair).contents).cid_1, ((LW.pair).contents).cid_2)
                links[key] = LP
            # add a link
            LI = (LW.LI).contents
            LP.addLink(LI.orient_1, LI.orient_2, LI.pos_1, LI.pos_2, LI.bam_ID)
            ret_val = CW._stepLW(pLW)
        CW._destroyLW(pLW)

    return links

def externalParseWrapper(bAMpARSER, bamFile, _MR, doContigNames):
    """ctypes pointers are unpickleable -- what we need is a hack!

    See BamParser._parseOneBam for what this function should be doing
    """
    MR = bAMpARSER._parseOneBam(bamFile)

    contig_names = []
    contig_lengths = np.array([int(i) for i in c.cast(MR.contig_lengths, c.POINTER(c.c_uint32*MR.num_contigs)).contents])
    plp_bp = np.array([[int(j) for j in c.cast(i, c.POINTER(c.c_uint32*MR.num_bams)).contents] for i in c.cast(MR.plp_bp,c.POINTER(c.POINTER(c.c_uint32*MR.num_bams)*MR.num_contigs)).contents])

    coverages = np.zeros((MR.num_contigs, MR.num_bams))
    if MR.is_outlier_coverage:
        contig_length_correctors = np.array([[int(j) for j in c.cast(i, c.POINTER(c.c_uint32*MR.num_bams)).contents] for i in c.cast(MR.contig_length_correctors,c.POINTER(c.POINTER(c.c_uint32*MR.num_bams)*MR.num_contigs)).contents])
        for c_idx in range(int(MR.num_contigs)):
            for b_idx in range(int(MR.num_bams)):
                coverages[c_idx,b_idx] = float(plp_bp[c_idx,b_idx])/float(contig_lengths[c_idx] - contig_length_correctors[c_idx])
    else:
        for c_idx in range(MR.num_contigs):
            for b_idx in range(MR.num_bams):
                coverages[c_idx,b_idx] = float(plp_bp[c_idx,b_idx])/float(contig_lengths[c_idx])

    if doContigNames:
        contig_name_lengths = np.array([int(i) for i in c.cast(MR.contig_name_lengths, c.POINTER(c.c_uint16*MR.num_contigs)).contents])
        contig_name_array = c.cast(MR.contig_names, c.POINTER(c.POINTER(c.c_char)*MR.num_contigs)).contents
        for i in range(MR.num_contigs):
            contig_names.append("".join([j for j in c.cast(contig_name_array[i], c.POINTER(c.c_char*contig_name_lengths[i])).contents]))

    if bAMpARSER.doLinks:
        links = pythonizeLinks(MR)
    else:
        links = {}

    MRR = BMM_mappingResults(coverages,
                            contig_lengths,
                            MR.num_bams,
                            MR.num_contigs,
                            contig_names,
                            bamFile,
                            links)
    _MR.append(MRR)

    # we need to call some C on this guy
    pMR = c.POINTER(BMM_mapping_results_C)
    pMR = c.pointer(MR)
    CW = CWrapper()
    CW._destroy_MR(pMR)

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
        package_dir = package_dir.replace("bamm","" )
        c_lib = os.path.join(package_dir, 'c', 'bam', 'libPMBam.a')
        self.libPMBam = c.cdll.LoadLibrary(c_lib)

        #---------------------------------
        # import C functions
        #---------------------------------
        self._merge_MR = self.libPMBam.merge_MRs
        """
        @abstract Merge the contents of MR_B into MR_A

        @param  MR_A  mapping results struct to copy to
        @param  MR_B  mapping results struct to copy from
        @return void
        @discussion MR_B remains unchanged.
        MR_A is updated to include all the info contained in MR_B

        void merge_MRs(BMM_mapping_results_C * MR_A, BMM_mapping_results_C * MR_B);
        """

        self._destroy_MR = self.libPMBam.destroy_MR
        """
        @abstract Free all the memory calloced in init_MR

        @param  MR  mapping results struct to destroy
        @return void

        void destroy_MR(BMM_mapping_results_C * MR)
        """

        self._parseCoverageAndLinks = self.libPMBam.parseCoverageAndLinks
        """
        @abstract Initialise the mapping results struct <- read in the BAM files

        @param numBams  number of BAM files to parse
        @param baseQ  base quality threshold
        @param mapQ  mapping quality threshold
        @param minLen  min query length
        @param doLinks  1 if links should be calculated
        @param ignoreSuppAlignments  only use primary alignments
        @param doOutlierCoverage  set to 1 if should initialise contig_length_correctors
        @param bamFiles  filenames of BAM files to parse
        @param MR  mapping results struct to write to
        @return 0 for success

        @discussion This function expects MR to be a null pointer. It calls
        init_MR and stores info accordingly. TL;DR If you call this function
        then you MUST call destroy_MR when you're done.

        int parseCoverageAndLinks(int numBams,
                                  int baseQ,
                                  int mapQ,
                                  int minLen,
                                  int doLinks,
                                  int ignoreSuppAlignments,
                                  int doOutlierCoverage,
                                  char* bamFiles[],
                                  BMM_mapping_results_C * MR
                                 )
        """

        self._adjustPlpBp = self.libPMBam.adjustPlpBp
        """
        @abstract Adjust (reduce) the number of piled-up bases along a contig

        @param  MR  mapping results struct to write to
        @param  position_holder  array of pileup depths
        @param  tid  contig currently being processed
        @param  doOutlierCoverage  remove effects fo very high or very low regions
        @return void

        @discussion This function expects MR to be initialised.
        it can change the values of contig_length_correctors and plp_bp

        void adjustPlpBp(BMM_mapping_results_C * MR,
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

        float ** calculateCoverages(BMM_mapping_results_C * MR);
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
        BMM_LinkWalker * initLW(BMM_mapping_results * MR);
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

        self._print_MR = self.libPMBam.print_MR
        """
        @abstract Print the contents of the MR struct

        @param  MR   mapping results struct with mapping info

        void print_MR(BMM_mapping_results_C * MR)
        """

class BamParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self,
                 baseQuality,
                 mappingQuality,
                 minLength,
                 doLinks,
                 ignoreSuppAlignments,
                 doOutliers
                 ):
        #---------------------------------
        # information about how the parser will be used
        #---------------------------------
        self.baseQuality = baseQuality
        self.mappingQuality = mappingQuality
        self.minLength = minLength
        self.doLinks = doLinks
        self.ignoreSuppAlignments = ignoreSuppAlignments
        self.doOutliers = doOutliers

        #---------------------------------
        # internal variables
        #---------------------------------
        self.MR = None          # internal mapping results object

    def parseBams(self, bamFiles, numThreads=1):
        """Parse bam files to get coverage and linking reads

        stores results in internal mapping results list
        """

        global _MR
        _MR = Manager().list()
        pool = Pool(processes=numThreads)
        do_contig_names = True
        for bamFile in bamFiles:
            pool.apply_async(func=externalParseWrapper, args=(self, bamFile, _MR, do_contig_names))
            if do_contig_names:
                # we only need to parse the contig names once
                do_contig_names = False
        pool.close()
        pool.join()

        # all the MRs are made. Only one has the contig IDs. find it's index
        base_MR_index = 0
        for i in range(len(_MR)):
            if len(_MR[i].contigNames) > 0:
                base_MR_index = i
                break
        # merge all the separate mapping results
        self.MR = _MR[base_MR_index]
        for i in range(len(_MR)):
            if i != base_MR_index:
                self.MR.consume(_MR[i])

    def _parseOneBam(self, bamFile):
        """Parse a single BAM file and append the result to the internal mapping results list"""
        MR = BMM_mapping_results_C()
        pMR = c.POINTER(BMM_mapping_results_C)
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
                                  self.doOutliers,
                                  bamfiles_c_array,
                                  pMR)
        return MR

    def destroy(self):
        """Clean up c-malloc'd memory"""
        if self.MR is not None:
            CW = CWrapper()
            pMR = c.POINTER(BMM_mapping_results_C)
            pMR = c.pointer(self.MR)
            CW._destroy_MR(pMR)

    def printMappings(self):
        """print a mapping results structure"""
        if self.MR is not None:
            CW = CWrapper()
            pMR = c.POINTER(BMM_mapping_results_C)
            pMR = c.pointer(self.MR)
            CW._print_MR(pMR)



###############################################################################
###############################################################################
###############################################################################
###############################################################################
