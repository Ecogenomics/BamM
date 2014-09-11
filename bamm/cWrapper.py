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

class CWrapper:
    """Can't pickle cTypes pointers and functions. Use this CWrap C-Wrapper as a hack"""
    def __init__(self, numBams=1, numContigs=0):
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
        self._destroyBFI = self.libPMBam.destroyBFI

        self._extractReads = self.libPMBam.extractReads

        self._parseCoverageAndLinks = self.libPMBam.parseCoverageAndLinks
        self._adjustPlpBp = self.libPMBam.adjustPlpBp
        self._calculateCoverages = self.libPMBam.calculateCoverages
        self._destroyCoverages = self.libPMBam.destroyCoverages

        self._initLW = self.libPMBam.initLW
        self._stepLW = self.libPMBam.stepLW
        self._destroyLW = self.libPMBam.destroyLW

        self._printBFI = self.libPMBam.printBFI

###############################################################################
###############################################################################
###############################################################################
###############################################################################
