#!/usr/bin/env python
###############################################################################
#                                                                             #
#    bamLink.py                                                               #
#                                                                             #
#    Class for storing information about paired links between contigs         #
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

# system imports
import ctypes as c

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# links-associated structures "Python land"
class BM_linkInfo(object):
    '''A single link joining two contigs'''

    def __init__(self,
                 r1,
                 r2,
                 p1,
                 p2,
                 bid = None):
        '''Default constructor.

        Initializes a BM_LinkInfo instance with the provided set of properties.

        Inputs:
         r1 - int, == 1 if the read maps to contig 1 in reverse orientation
         r2 - int, == 1 if the read maps to contig 2 in reverse orientation
         pos1 - int, leftmost position of read on contig 1
         pos2 - int, leftmost position of read on contig 2
         bid - unique identifier for the bam file describing this link

        Outputs:
         None
        '''
        self.reversed1 = r1
        self.reversed2 = r2
        self.pos1 = p1
        self.pos2 = p2
        self.bamID = bid

    def printMore(self,
                  bamFileNames,
                  len1,
                  len2):
        '''Advanced string function

        used for creating output to links file

        Inputs:
         bamFileNames - { bamId : string }, storage for long bam file names
         len1 - int, length of contig 1
         len1 - int, length of contig 2

        Outputs:
         String describing the link
        '''
        return "%d\t%d\t%d\t%d\t%d\t%d\t%s" % (len1,
                                               self.pos1,
                                               self.reversed1,
                                               len2,
                                               self.pos2,
                                               self.reversed2,
                                               bamFileNames[self.bamID])

    def __str__(self):
        '''override basic string function'''
        return "%d\t%d\t%d\t%d\t%s" % (self.pos1,
                                       self.reversed1,
                                       self.pos2,
                                       self.reversed2,
                                       self.bamID)



###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BM_linkPair(object):
    '''Container class for storing all links joining two contigs'''

    def __init__(self,
                 cid1,
                 cid2):
        '''Default constructor.

        Initializes a BM_LinkPair instance with the provided set of properties.

        Inputs:
         cid1 - unique identifier for contig 1
         cid2 - unique identifier for contig 2

        Outputs:
         None
        '''
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
        '''Add a link between the two contigs

        Inputs:
         r1 - int, == 1 if the read maps to contig 1 in reverse orientation
         r2 - int, == 1 if the read maps to contig 2 in reverse orientation
         pos1 - int, leftmost position of read on contig 1
         pos2 - int, leftmost position of read on contig 2
         bamFile - BM_bamFile instance that describes this link

        Outputs:
         None
        '''
        LI = BM_linkInfo(r1,
                         r2,
                         p1,
                         p2,
                         bamFile.bid)
        self.links.append(LI)
        self.numLinks += 1

    def makeKey(self):
        '''Return a unique key for this link pair'''
        return "%d,%d" % (self.cid1, self.cid2)

    def printMore(self,
                  contigNames,
                  contigLengths,
                  bamFileNames):
        '''Advanced string function

        used for creating output to links file
        calls link.printmore()

        Inputs:
         contigNames - { cid : string }, storage for long contig names
         contigLengths - { cid : int }, storage for contig lengths
         bamFileNames - { bamId : string }, storage for long bam file names

        Outputs:
         Multi-line string descrtibing all the links between the two contigs
        '''

        """Print function used to export a tabular format"""
        return "\n".join(["%s\t%s\t%s"%(contigNames[self.cid1],
                                        contigNames[self.cid2],
                                        link.printMore(bamFileNames,
                                                       contigLengths[self.cid1],
                                                       contigLengths[self.cid2]
                                                       )
                                          )
                          for link in self.links])

    def __str__(self):
        '''override basic string function'''
        return "\n".join(["%s\t%s\t%s" % (self.cid1,
                                          self.cid2,
                                          link) for link in self.links]
                         )

###############################################################################
###############################################################################
###############################################################################
###############################################################################
