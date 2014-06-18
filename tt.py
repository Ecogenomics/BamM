#!/usr/bin/env python
###############################################################################
#
# StoreM.py - Metagenomic data storage
#
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
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

import argparse
import sys

from bamm.BamParser import BamParser

# Links are stored as triples (contig1, contig2, linktype)
# There are 4 linktypes:
# SS  <--1--- ---2-->
# SE  <--1--- <--2---
# ES  ---1--> ---2-->
# EE  ---1--> <--2---
global LT
LT = enum('SS','SE','ES','EE','ERROR')

def LT2Str(cid1, cid2, gap, linkType, terse=False):
    """For the humans!"""
    if terse:
        if linkType == LT.SS:
            return "SS"
        if linkType == LT.SE:
            return "SE"
        if linkType == LT.ES:
            return "ES"
        if linkType == LT.EE:
            return "EE"
        return "??"
    if linkType == LT.SS:
        return str(cid1) + " lies before " + str(cid2) + " in the opposite direction with gap "+str(gap)
    if linkType == LT.SE:
        return str(cid1) + " lies after " + str(cid2) + " in the same direction with gap "+str(gap)
    if linkType == LT.ES:
        return str(cid1) + " lies before " + str(cid2) + " in the same direction with gap "+str(gap)
    if linkType == LT.EE:
        return str(cid1) + " lies after " + str(cid2) + " in the opposite direction with gap "+str(gap)

    return 'Who knows?'

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def doWork( args ):
    """ Main wrapper"""
    BP = BamParser(args.base_quality,
                   args.length,
                   mappingQuality=args.mapping_quality,
                   coverageMode=args.coverage_mode,
                   doLinks=args.links)

    BP.parseBams(args.bams, numThreads=3)
    print BP.MR

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('bams', nargs='+', help="BAM files to parse")
    parser.add_argument('--links', '-L', help="find pairing links", action='store_true', default=False)
    parser.add_argument('--length', '-l', help="minimum Q length", type=int, default=50)
    parser.add_argument('--base_quality', '-q', help="base quality threshold (Qscore)", type=int, default=20)
    parser.add_argument('--mapping_quality', '-Q', help="mapping quality threshold", type=int, default=0)
    parser.add_argument('--coverage_mode' ,'-m', help="how to calculate coverage [vanilla, outlier]", default='vanilla')

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
