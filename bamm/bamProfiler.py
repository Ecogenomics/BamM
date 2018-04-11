###############################################################################
#                                                                             #
#    BamProfiler.py                                                             #
#                                                                             #
#    Class for profiling BAM file reads                                       #
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
import ctypes as c

# local imports
from cWrapper import CWrapper


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamProfiler:
    '''Class used to manage profiling reads from a BAM file'''
    def __init__(self,
                 bamFile,
                 useSuppAlignments=False,
                 useSecondaryAlignments=False
                 ):
        '''
        Default constructor.

        Set all the instance variables, make ReadSets, organise output files

        Inputs:
         useSuppAlignments - == True -> DON'T skip supplementary alignments
         useSecondaryAlignments - == True -> DON'T skip secondary alignments

        Outputs:
         None
        '''

        self.bamFile = bamFile

        if useSuppAlignments:
            self.ignoreSuppAlignments = 0
        else:
            self.ignoreSuppAlignments = 1

        if useSuppAlignments:
            self.ignoreSecondaryAlignments = 0
        else:
            self.ignoreSecondaryAlignments = 1

		
    def profile(self):
        '''Start filtering reads from BAM files according to
		configured quality metrics.
        '''
        # first we need to C-ify variables
        bamfile_c = c.c_char_p()
        bamfile_c = self.bamFile

        # call the C function to filter the reads
        CW = CWrapper()
        CW._profileReads(bamfile_c,
                         self.ignoreSuppAlignments,
                         self.ignoreSecondaryAlignments)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

