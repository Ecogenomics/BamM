###############################################################################
#                                                                             #
#    BamFilter.py                                                             #
#                                                                             #
#    Class for filtering BAM file reads                                       #
#                                                                             #
#    Copyright (C) Tim Lamberton                                              #
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

__author__ = "Tim Lamberton"
__copyright__ = "Copyright 2015"
__credits__ = ["Tim Lamberton"]
__license__ = "LGPLv3"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"

###############################################################################

# system imports
import os
import ctypes as c

# local imports
from cWrapper import CWrapper


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamFilter:
    '''Class used to manage filtering reads from multiple BAM files'''
    def __init__(self,
                 bamFile,
                 outFolder=".",
                 minMapQual=0,
                 minLength=0,
                 maxMisMatches=1000,
                 minPcId=0.,
                 minPcAln=0.,
                 invertMatch=False,
                 useSuppAlignments=False,
                 useSecondaryAlignments=False,
                 ):
        '''
        Default constructor.

        Set all the instance variables, make ReadSets, organise output files

        Inputs:
         bamFiles - string, list of bamfiles to extract reads from
         outFolder - path, write output to this folder
         minMapQual - int, skip all reads with a lower mapping quality score
         minLength - int
         maxMisMatches - int, skip all reads with more mismatches (NM aux files)
         minPcId - float
         minPcAln - float
         useSuppAlignments - == True -> DON'T skip supplementary alignments
         useSecondaryAlignments - == True -> DON'T skip secondary alignments

        Outputs:
         None
        '''
        # make sure the output folder exists
        self.outFolder = outFolder
        # it's a waste if we abort but I like to check if write permissions
        # are intact before I do lots of work.
        self.makeSurePathExists(self.outFolder)

        self.bamFile = bamFile
        self.prettyBamFileName = os.path.basename(self.bamFile).replace(".bam", "")

        self.minMapQual = minMapQual
        self.minLength = minLength
        self.maxMisMatches = maxMisMatches
        self.minPcId = minPcId
        self.minPcAln = minPcAln

        if useSuppAlignments:
            self.ignoreSuppAlignments = 0
        else:
            self.ignoreSuppAlignments = 1

        if useSuppAlignments:
            self.ignoreSecondaryAlignments = 0
        else:
            self.ignoreSecondaryAlignments = 1
            
        if invertMatch:
            self.invertMatch = 1
        else:
            self.invertMatch = 0


    def filter(self):
        '''Start extracting reads from the BAM files

        This function is responsible for starting and stopping all threads and
        processes used in bamm extract. Due to python multiprocessing's need to
        pickle everything the actual work of extraction is carried out in the
        first level function called externalExtractWrapper. See there for actual
        extraction details. This function is primarily concerned with thread
        and process management.

        Inputs:
         threads - int, the number of threads / processes to use
         verbose - bool, True if lot's of stuff should be printed to screen

        Outputs:
         None
        '''
        # first we need to C-ify variables
        bamfile_c = c.c_char_p()
        bamfile_c = self.bamFile

        outfile_c = c.c_char_p()
        outfile_c = os.path.join(os.path.abspath(self.outFolder), "%s_%s.bam" % (self.prettyBamFileName, 'filtered'))

        min_mapping_quality_c = c.c_uint32()
        min_mapping_quality_c = self.minMapQual

        min_query_length_c = c.c_uint32()
        min_query_length_c = self.minLength

        max_mismatches_c = c.c_uint32()
        max_mismatches_c = self.maxMisMatches

        min_percentage_id_c = c.c_float()
        min_percentage_id_c = self.minPcId

        min_percentage_aln_c = c.c_float()
        min_percentage_aln_c = self.minPcAln

        # call the C function to filter the reads
        CW = CWrapper()
        CW._filterReads(bamfile_c,
                        outfile_c,
                        min_mapping_quality_c,
                        min_query_length_c,
                        max_mismatches_c,
                        min_percentage_id_c,
                        min_percentage_aln_c,
                        self.invertMatch,
                        self.ignoreSuppAlignments,
                        self.ignoreSecondaryAlignments)

    def makeSurePathExists(self, path):
        '''Make sure that a path exists, make it if necessary

        Inputs:
         path - string, full or relative path to create

        Outputs:
         None
        '''
        try:
            os.makedirs(path)
        except OSError as exception:
            import errno
            if exception.errno != errno.EEXIST:
                raise

###############################################################################
###############################################################################
###############################################################################
###############################################################################

