#!/usr/bin/env python
###############################################################################
#
# sysInfoBamM.py - get some infor about the system we're running on
#
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

if __name__ == "__main__":

    try:
        import numpy
        print "Numpy version: %s" % numpy.__version__
        print "Location: %s" % numpy.__file__
    except ImportError:
        print "Error importing numpy"

    import platform
    import sys
    print platform.system(), platform.release()
    print(sys.version)
    print(sys.path)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
