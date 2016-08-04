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
###############################################################################
###############################################################################
###############################################################################
import sys

def printError(error):
    """Make a string kinda errory"""
    sys.stderr.write("%s\n\nERROR:\n\n%s\n\n%s\n\n" % (''.join(["*"]*80),
                                                       str(error),
                                                       ''.join(["*"]*80)))
def printShortUsage(mode=None):
    """Point the user to a fuller help with a short message on STDERR
    mode: bamm mode, or None for general help
    """
    
    if mode:
        mode_and_space = ' %s' % mode
    else:
        mode_and_space = ''
        
    sys.stderr.write("For further usage details:\n\nbamm%s -h\n\n" % mode_and_space)

class BamMException(BaseException): pass

class InvalidInstallationException(BamMException): pass


class InvalidCoverageModeException(BamMException): pass
class InvalidNumberOfTypesException(BamMException): pass
class BAMFileNotFoundException(BamMException): pass
class BAMIndexNotFoundException(BamMException): pass
class NoBAMSFoundException(BamMException): pass

class InvalidParameterSetException(BaseException): pass

class MixedFileTypesException(BamMException): pass

class DuplicateSequenceNameException(BamMException): pass
###############################################################################
###############################################################################
###############################################################################
###############################################################################
