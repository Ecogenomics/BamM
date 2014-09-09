#!/usr/bin/env python
###############################################################################
#                                                                             #
#    bammExceptions.py                                                        #
#                                                                             #
#    Handy dandy place for storing exceptions                                 #
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
__version__ = "1.0"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Alpha"

###############################################################################
###############################################################################
###############################################################################
###############################################################################
import sys

def printError(error):
    """Make a string kinda errory"""
    sys.stderr.write("*******************************************************************************\n\nERROR:\n\n" + str(error) + "\n\n*******************************************************************************\n\n")

class BamMException(BaseException): pass
class InvalidCoverageModeException(BamMException): pass
class InvalidNumberOfTypesException(BamMException): pass
class BAMFileNotFoundException(BamMException): pass
class BAMIndexNotFoundException(BamMException): pass
class NoBAMSFoundException(BamMException): pass

class InvalidParameterSetException(BaseException): pass

###############################################################################
###############################################################################
###############################################################################
###############################################################################
