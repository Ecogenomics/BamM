#!/usr/bin/env python
###############################################################################
#                                                                             #
#    BamLinks.py                                                              #
#                                                                             #
#    Class for storing information about links between contigs                #
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
__version__ = "0.1.0"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Beta"

###############################################################################

# system imports
import ctypes as c

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

# fields defined in cfuhash.c but not accessed at this level
class cfuhash_table_t(c.Structure):
    pass

# links-associated structures "C land"
"""
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

typedef struct {
    uint32_t cid1;
    uint32_t cid2;
    uint32_t numLinks;
    BM_linkInfo * LI;
} BM_linkPair;

typedef struct {
    char ** keys;
    size_t keyCount;
    size_t numKeys;
    cfuhash_table_t * linkHash;
    BM_linkPair * pair;
    BM_linkInfo * LI;
} BM_LinkWalker;

"""
class BM_linkInfo_C(c.Structure):
    pass

class BM_linkInfo_C(c.Structure):
    _fields_ = [("reversed1", c.c_uint16),
                ("reversed2", c.c_uint16),
                ("readLength1", c.c_uint16),
                ("readLength2", c.c_uint16),
                ("pos1", c.c_uint32),
                ("pos2", c.c_uint32),
                ("bam_ID", c.c_uint32),
                ("nextLink",c.POINTER(BM_linkInfo_C))
                ]

class BM_linkPair_C(c.Structure):
    _fields_ = [("cid1", c.c_uint32),
                ("cid2", c.c_uint32),
                ("numLinks", c.c_uint32),
                ("LI",c.POINTER(BM_linkInfo_C))
                ]

class BM_LinkWalker_C(c.Structure):
    _fields_= [("keys", c.POINTER(c.POINTER(c.c_char))),
               ("keyCount", c.c_size_t),
               ("numKeys", c.c_size_t),
               ("links",c.POINTER(cfuhash_table_t)),
               ("pair",c.POINTER(BM_linkPair_C)),
               ("LI",c.POINTER(BM_linkInfo_C))
               ]

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
