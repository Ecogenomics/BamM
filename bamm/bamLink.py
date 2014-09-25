#!/usr/bin/env python
###############################################################################
#                                                                             #
#    bamLink.py                                                              #
#                                                                             #
#    Class for storing information about links between contigs                #
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

# local imports
from cWrapper import *

###############################################################################
###############################################################################
###############################################################################
###############################################################################

#------------------------------------------------------------------------------
# Managing orientation and linking types

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

# Read orientations
# type 0 <--- --->
# type 1 ---> --->
# type 2 ---> <---
global OT
OT = enum('OUT', 'SAME', 'IN', 'NONE', 'ERROR')

def OT2Str(ot):
    """For the humans!"""
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

# links-associated structures "Python land"
class BM_linkInfo(object):
    def __init__(self,
                 r1,
                 r2,
                 p1,
                 p2,
                 bid = None):
        self.reversed1 = r1
        self.reversed2 = r2
        self.pos1 = p1
        self.pos2 = p2
        self.bamID = bid

    def __str__(self):
        return "%d\t%d\t%d\t%d\t%s" % (self.pos1, self.reversed1, self.pos2, self.reversed2, self.bamID)

    def printMore(self, bamFileNames, len1, len2):
        return "%d\t%d\t%d\t%d\t%d\t%d\t%s" % (len1, self.pos1, self.reversed1, len2, self.pos2, self.reversed2, bamFileNames[self.bamID])

class BM_linkPair(object):
    def __init__(self,
                 cid1,
                 cid2):
        self.cid1 = cid1
        self.cid2 = cid2
        self.numLinks = 0
        self.links = []

    def addLink(self,
                r1,
                r2,
                p1,
                p2,
                bamFile):

        LI = BM_linkInfo(r1,
                         r2,
                         p1,
                         p2,
                         bamFile.bid)
        self.links.append(LI)
        self.numLinks += 1

    def __str__(self):
        return "\n".join(["%s\t%s\t%s" % (self.cid1, self.cid2, link) for link in self.links])

    def printMore(self,
                  contigNames,
                  contigLengths,
                  bamFileNames):
        """Print function used to export a tabular format"""
        return "\n".join(["%s\t%s\t%s" % (contigNames[self.cid1],
                                          contigNames[self.cid2],
                                          link.printMore(bamFileNames, contigLengths[self.cid1], contigLengths[self.cid2])
                                          )
                          for link in self.links])

###############################################################################
###############################################################################
###############################################################################
###############################################################################
