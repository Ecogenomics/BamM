#!/usr/bin/env python
###############################################################################
#                                                                             #
#    BamExtractor.py                                                          #
#                                                                             #
#    Class for extracting reads from BAM files                                #
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
import os
import ctypes as c
from multiprocessing import Lock, Manager, Process
import numpy as np
import sys
import gzip
import Queue

# local imports
from cWrapper import *
from bamFile import *
from bammExceptions import *
from bamRead import *

###############################################################################
###############################################################################
###############################################################################
# Multiprocessing requires that all passed items be pickleable. That is they
# must be vanilla variables or functions defined in the file itself, ie. not
# within a class. We get around this by writing an external function which calls
# a class function. Hacky, but it works.
###############################################################################
###############################################################################
###############################################################################

def externalExtractWrapper(bAMeXTRACTOR, extractQueue, readSetQueue, verbose):
    """ctypes pointers are unpickleable -- what we need is a hack!

    See BamExtractor._parseOneBam for what this function should be doing
    """
    CW = CWrapper()
    files_to_parse = True
    timeout = 0

    # how many reads should accrue before we write them out?
    cache_size = 500

    while True:
        # see if there is anything to write to disk
        try:
            read_set = readSetQueue.get(block=True, timeout=timeout)
            replace = True

            # put it back on the of the queue
            with bAMeXTRACTOR.Lock:
                if bAMeXTRACTOR.numReadingThreads == 0:
                    # all threads have finished adding to the lists
                    # we don't need to put this one back on the end of the queue
                    replace = False

            # write some stuff
            if not replace:
                # we are going to write the rest of this guy to file regardless
                read_set.write(bAMeXTRACTOR.targetNames)
            else:
                (buffer_size, num_reads_queued) = read_set.getSize()
                if num_reads_queued >= cache_size:
                    print "going to write"
                    read_set.write(bAMeXTRACTOR.targetNames, maxDump=cache_size)

            # put it back on the end of the queue if necessary
            if replace:
                readSetQueue.put(read_set)

        except Queue.Empty:
            # empty queue only matters when the extract Q is depleted
            if not files_to_parse:
                break

        timeout=2

        # get the next one off the list if we've not seen our None already
        if files_to_parse:
            bid = extractQueue.get(block=True, timeout=None)
            if bid is None: # poison pill
                files_to_parse = False
                with bAMeXTRACTOR.Lock:
                    bAMeXTRACTOR.numReadingThreads -= 1

            else:
                if verbose:
                    print "Extracting reads from file: %s" % bAMeXTRACTOR.bamFiles[bid]
                bAMeXTRACTOR._extractFromOneBam(bid)
                if verbose:
                    print "Done extracting reads from file: %s" % bAMeXTRACTOR.bamFiles[bid]


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamExtractor:
    """Main class for reading in and parsing contigs"""
    def __init__(self,
                targets,                    # list of contig IDs or fasta file (used as a filter)
                bamFiles,                   # list of bamfiles to extract reads from
                prefix="",                  # append this string to the beginning of all output files
                targetNames=[],             # list of names of the groups in the targets list
                outFolder=".",              # wriate output to this folder
                mixBams=False,              # use one file for all bams
                mixTargets=False,           # use one file for all groups
                mixReads=False,             # use one file for paired / unpaired reads
                shuffle=False,              # use shuffled format for paired reads
                ignoreUnpaired=False,       # ignore all upaired reads
                bigFile=False,              # do NOT gzip outputs
                headersOnly=False           # write read headers only
                ):

        # make sure the output folder exists
        self.outFolder = outFolder
        self.makeSurePathExists(self.outFolder) # it's a waste if we abort but I like to check if write permissions are intact before I do lots of work.

        # the number of threads that are still reading BAM files
        self.numReadingThreads = 0
        # a global lock to protect numReadingThreads
        self.Lock = Lock()

        self.bamFiles = bamFiles
        self.prettyBamNames = []
        for bam in self.bamFiles:
            self.prettyBamNames.append(os.path.basename(bam).replace(".bam", ""))
        self.prefix = prefix

        self.mixBams = mixBams
        self.mixTargets = mixTargets
        self.mixReads = mixReads

        self.ignoreUnpaired = ignoreUnpaired
        self.shuffle = shuffle
        if headersOnly:
            self.headersOnly = 1
        else:
            self.headersOnly = 0

        # are we going to zip the output?
        if bigFile:
            self.zipped = False
        else:
            self.zipped = True

        # munge the targets
        if targetNames == []:
            # no names specified, just use "List_1", "List_2" etc...
            targetNames = ["list_%d" % i for i in range(1, len(targets)+1)]
        self.targetNames = targetNames

        # initialise to the first set of targets
        self.targets = targets[0]
        self.bins = [0]*len(self.targets)

        for targ in range(1, len(targets)):
            self.targets += targets[targ]
            self.bins += [targ] * len(targets[targ])

        self.manager = Manager()

        self.outFiles = {}          # target, bam, rpi -> (filename, readQueue)
        self.readSetsQueue = self.organiseOutFiles()

    def organiseOutFiles(self):
        """Open all the output file handles we'll need

        full path determined
        extension and zippability not determined
        """
        # we need to make a filename for every eventuality
        prefix = os.path.join(os.path.abspath(self.outFolder), self.prefix)

        # place all the file queues onto one global queue. This is where we'll do the writing from
        read_set_queue = self.manager.Queue()

        for bid in range(len(self.bamFiles)):
            if self.mixBams:
                bam_str = "allMapd"
            else:
                bam_str = self.prettyBamNames[bid]

            if bam_str not in self.outFiles:
                self.outFiles[bid] = {}

            for tid in range(len(self.targetNames)):
                if self.mixTargets:
                    tar_str = "allTargets"
                else:
                    tar_str = self.targetNames[tid]

                if tid not in self.outFiles[bid]:
                    self.outFiles[bid][tid] = {}
                if self.prefix == "":
                    fn = "%s%s.%s" % (prefix, bam_str, tar_str)
                else:
                    fn = "%s.%s.%s" % (prefix, bam_str, tar_str)
                if self.mixReads:
                    # all reads should go into the one file
                    read_set_P = ReadSet(fn + ".allReads", zipped=self.zipped)
                    read_set_S = read_set_P

                    read_set_queue.put(read_set_P)

                elif self.shuffle:
                    # one file for pairs and one for singles
                    read_set_P = ReadSet(fn + ".pairedReads", paired=True, zipped=self.zipped)
                    read_set_S = ReadSet(fn + ".unpairedReads", zipped=self.zipped)

                    read_set_queue.put(read_set_P)
                    read_set_queue.put(read_set_S)

                else:
                    # each in their own file
                    paired_fn1 = fn + ".1"
                    paired_fn2 = fn + ".2"
                    read_set_P = ReadSet(paired_fn1, fileName2=paired_fn2, paired=True, zipped=self.zipped)
                    read_set_S = ReadSet(fn + ".unpairedReads", zipped=self.zipped)

                    read_set_queue.put(read_set_P)
                    read_set_queue.put(read_set_S)

                # we use the filenames to link everything up below
                self.outFiles[bid][tid][RPI.FIR] = read_set_P
                self.outFiles[bid][tid][RPI.SEC] = read_set_P
                if self.ignoreUnpaired:
                    self.outFiles[bid][tid][RPI.SNGL] = None
                else:
                    self.outFiles[bid][tid][RPI.SNGL] = read_set_S

        return read_set_queue

    def extract(self, threads=1, verbose=False):
        """Extract all the reads"""
        # start running the parser in multithreaded mode
        extract_queue = self.manager.Queue()

        # place the bids on the queue
        for bid in range(len(self.bamFiles)):
            extract_queue.put(bid)

        # place one None on the queue for each thread we have access to
        for _ in range(threads):
            extract_queue.put(None)

        # eventually we need to stop putting file objects back onto the queue
        # as each thread gets to the end of it's parsing it decrements this value
        # when it is 0 then we know that there is nothing more to read
        with self.Lock:
            self.numReadingThreads = threads

        try:
            # make the processes
            extract_proc = [Process(target=externalExtractWrapper,
                                    args=(self, extract_queue, self.readSetsQueue, verbose)
                                    )
                            for _ in range(threads)]
            for p in extract_proc:
                p.start()

            for p in extract_proc:
                p.join()

            # success
            return 0

        except:
            # ctrl-c! Make sure all processes are terminated
            for p in extract_proc:
                p.terminate()

            # dismal failure
            return 1

    def _extractFromOneBam(self, bid):
        """Extract reads mapping to contigs from a single BAM"""
        bamfile_c = c.c_char_p()
        bamfile_c = self.bamFiles[bid]

        pretty_name_c = c.c_char_p()
        pretty_name_c = self.prettyBamNames[bid]

        num_contigs = len(self.targets)
        contigs_c_array = (c.c_char_p * num_contigs)()
        contigs_c_array[:] = self.targets

        bins_c_array = (c.c_uint16 * num_contigs)()
        bins_c_array[:] = self.bins

        headers_only_c = c.c_uint32()
        headers_only_c = self.headersOnly

        pBMM = c.POINTER(BM_mappedRead_C)
        BMM = BM_mappedRead_C

        # call the C function
        CW = CWrapper()
        pBMM = CW._extractReads(bamfile_c,
                                contigs_c_array,
                                num_contigs,
                                bins_c_array,
                                pretty_name_c,
                                headers_only_c)

        print "fin C part for:",  self.prettyBamNames[bid]

        for_destruction = pBMM      # need to remember the root of the linked list so we can destroy it

        num_made = 0
        coutt = 0
        # now munge the c linked list into something more pythonic
        while pBMM != 0:
            # get hold of the next item in the linked list
            root = c.cast(pBMM, c.POINTER(BM_mappedRead_C))
            # is, bin and rpi are always there
            id = (c.cast(root.contents.seqId, c.POINTER(c.c_char*root.contents.idLen)).contents).value
            bin = root.contents.bin
            rpi = root.contents.rpi
            qual = None
            seq = None
            # we only stored these guys if asked to
            if 0 == self.headersOnly:
                seq = (c.cast(root.contents.seq, c.POINTER(c.c_char*root.contents.seqLen)).contents).value
                qual_len = root.contents.qualLen
                if qual_len > 0:
                     qual (c.cast(root.contents.qual, c.POINTER(c.c_char*qual_len)).contents).value

            # make a mapped read and place on the correct queue
            new_BMM = BM_mappedRead(id,
                                    seq=seq,
                                    qual=qual,
                                    rpi=rpi,
                                    bin=bin
                                    )

            # place the read onto the management queue
            self.outFiles[bid][bin][rpi].add(new_BMM)

            coutt += 1
            num_made += 1
            if coutt >= 100000:
                coutt = 0
                print "Added %d for: %s" %(num_made, self.prettyBamNames[bid])


            # next!
            pBMM = CW._nextMappedRead(pBMM)

        # clean up all the C allocated memory now.
        CW._destroyMappedReads(for_destruction)

    def makeSurePathExists(self, path):
        try:
            os.makedirs(path)
        except OSError as exception:
            import errno
            if exception.errno != errno.EEXIST:
                raise
