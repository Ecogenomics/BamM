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
__version__ = "1.0.0-b.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Beta"

###############################################################################

# system imports
import os
import ctypes as c
from multiprocessing import Manager, Process, Value
import numpy as np
import sys
import gzip
import Queue
import time
from threading import Thread

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

def externalExtractWrapper(threadId,
                           outFilePrefixes,
                           bamPaths,
                           prettyBamFileNames,
                           numGroups,
                           perContigGroups,
                           contigs,
                           printQueue,
                           extractQueue,
                           requestQueue,
                           freeQueue,
                           responseQueue,
                           headersOnly,
                           mixGroups,
                           verbose=False
                           ):
    """ctypes pointers are unpickleable -- what we need is a hack!

    See BamExtractor._parseOneBam for what this function should be doing
    """
    while True:
        if verbose:
            printQueue.put("%s Begin" % (threadId) )
        p_bid = extractQueue.get(block=True, timeout=None)
        if p_bid is None: # poison pill
            break
        else:
            if verbose:
                printQueue.put("%s Extracting reads from file: %s" % (threadId, prettyBamFileNames[p_bid] ) )
            extractFromOneBam(p_bid,
                              outFilePrefixes,
                              bamPaths,
                              prettyBamFileNames,
                              numGroups,
                              perContigGroups,
                              contigs,
                              requestQueue,
                              freeQueue,
                              responseQueue,
                              printQueue,
                              headersOnly,
                              mixGroups,
                              threadId=threadId,
                              verbose=verbose)
            if verbose:
                printQueue.put("%s Done extracting reads from file: %s" % (threadId, prettyBamFileNames[p_bid] ) )

def extractFromOneBam(bid,
                      outFilePrefixes,
                      bamPaths,
                      prettyBamFileNames,
                      numGroups,
                      perContigGroups,
                      contigs,
                      requestQueue,
                      freeQueue,
                      responseQueue,
                      printQueue,
                      headersOnly,
                      mixGroups,
                      threadId="",
                      verbose=False):
    """Extract reads mapping to contigs from a single BAM"""
    bamfile_c = c.c_char_p()
    bamfile_c = bamPaths[bid]

    pretty_name_c = c.c_char_p()
    pretty_name_c = prettyBamFileNames[bid]

    num_contigs = len(contigs)
    contigs_c_array = (c.c_char_p * num_contigs)()
    contigs_c_array[:] = contigs

    groups_c_array = (c.c_uint16 * num_contigs)()
    groups_c_array[:] = perContigGroups

    headers_only_c = c.c_uint32()
    headers_only_c = headersOnly

    pBMM = c.POINTER(BM_mappedRead_C)

    # call the C function to extract the reads
    CW = CWrapper()
    pBMM = CW._extractReads(bamfile_c,
                            contigs_c_array,
                            num_contigs,
                            groups_c_array,
                            pretty_name_c,
                            headers_only_c)

    if verbose:
        printQueue.put("%s Finished C-based extraction for: %s" % (threadId, prettyBamFileNames[bid]))

    overlapper = {} # helper, lets us know which readSets have the same filename
    chain_info = {} # store info about the start / end and count of a printing read chain
    for gid in range(numGroups):
        chain_info[gid] = {}
        for rpi in [RPI.FIR, RPI.SNGL]: # ReadSets exist for only FIR and SNGL
            file_name = outFilePrefixes[bid][gid][rpi]
            try:
                storage = overlapper[file_name]
            except KeyError:
                storage = [None, None, 0] # [start of chain, end of chain, chain length]
                overlapper[file_name] = storage
            chain_info[gid][rpi] = {'storage' : storage}

    # now munge the c linked list into something more pythonic
    while pBMM:
        # get hold of the next item in the linked list

        rpi = pBMM.contents.rpi
        c_rpi = RPIConv[rpi]
        addable = []

        if c_rpi != RPI.SEC:  # RPI.FIR or RPI.SNGL
            # append RPI.FIR and RPI.SNGL, SEC is handled below
            addable.append([c.addressof(pBMM.contents), c_rpi])

        # use raw rpi here!
        if rpi == RPI.FIR:
            # We know this guys has a partner however
            # we may need to treat this as a single read
            # or we may have to step up the order of it's partner
            r2_rpi = RPI.ERROR
            if (1 == CW._partnerInSameGroup(pBMM)):
                # partner is in same group.
                # RPI.FIR and RPI.SEC ALWAYS point to the same ReadSet
                r2_rpi = RPI.FIR
            else:
                # partner is in a different group
                # we should treat as a single, unless we don't care (mixGroups == True)
                if mixGroups:
                    # we don't care, print it now as a pair
                    # RPI.FIR and RPI.SEC ALWAYS point to the same ReadSet
                    r2_rpi = RPI.FIR
                else:
                    # we'll treat both paired reads as singles
                    r2_rpi = RPI.SNGL
                    addable[0][1] = RPI.SNGL # update this guy
            addable.append([CW._getPartner(pBMM), r2_rpi])

        # update the printing chain
        for mappedRead in addable:

            working_rpi = mappedRead[1]
            tmp_pBMM = c.cast(mappedRead[0], c.POINTER(BM_mappedRead_C))
            has_qual = (tmp_pBMM.contents.qualLen != 0)
            group = tmp_pBMM.contents.group

            # set and check the quality info
            try:
                if chain_info[group][working_rpi]['isFastq'] ^ has_qual:
                    raise MixedFileTypesException("You cannot mix Fasta and Fastq reads together in an output file")
            except KeyError:
                # Now we can set the type of the file.
                # We will only get here on the first read for each group, rpi
                chain_info[group][working_rpi]['isFastq'] = has_qual

            # build or maintain the chain
            if chain_info[group][working_rpi]['storage'][1] is None:
                # this is the first time we've seen this print chain
                chain_info[group][working_rpi]['storage'][0] = mappedRead[0]
                chain_info[group][working_rpi]['storage'][1] = mappedRead[0]
                chain_info[group][working_rpi]['storage'][2] = 1
            else:
                # join this pBMM onto the end of the existing chain
                CW._setNextPrintRead(c.cast(chain_info[group][working_rpi]['storage'][1],
                                            c.POINTER(BM_mappedRead_C)),
                                     c.cast(mappedRead[0],
                                            c.POINTER(BM_mappedRead_C)))
                chain_info[group][working_rpi]['storage'][1] = mappedRead[0]
                chain_info[group][working_rpi]['storage'][2] += 1

        # next!
        pBMM = CW._getNextMappedRead(pBMM)

    # place the start of each read onto the appropriate management queue
    printQueue.put("%s Preparing to write" % (threadId))
    for gid in range(numGroups):
        for rpi in [RPI.FIR, RPI.SNGL]:
            if chain_info[gid][rpi]['storage'][1] is not None:
                pBMM_chain = c.cast(chain_info[gid][rpi]['storage'][0], c.POINTER(BM_mappedRead_C))

                # we need to print here, so what we will do is make a request to the RSM
                # for a fileName etc. that we can write to. We block on this call
                # so we may have to wait for a bit BUT... it's either this, or go single
                # threaded. So this is what we'll do.
                requestQueue.put((threadId,
                                  bid,
                                  gid,
                                  rpi,
                                  chain_info[gid][rpi]['isFastq']))
                # wait for the RSM to return us a copy of a ReadSet
                RS = responseQueue.get(block=True, timeout=None)
                if RS is None:
                    # free the memory, it is useless to me!
                    CW._destroyPrintChain(pBMM_chain)
                else:
                    # we can print stuff
                    pBMM_destroy = c.POINTER(BM_mappedRead_C)
                    pBMM_destroy = pBMM_chain
                    RS.writeChain(pBMM_chain, chain_info[gid][rpi]['isFastq'])
                    CW._destroyPrintChain(pBMM_destroy)

                    # free the RS now
                    freeQueue.put((threadId,
                                   bid,
                                   gid,
                                   rpi))

                # set this to None so it's not added twice
                chain_info[gid][rpi]['storage'][1] = None

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamExtractor:
    """Main class for reading in and parsing contigs"""
    def __init__(self): pass

    def initialise(self,
                   contigs,                # list of contig IDs or fasta file (used as a filter)
                   bamFiles,              # list of bamfiles to extract reads from
                   prefix="",             # append this string to the beginning of all output files
                   groupNames=[],         # list of names of the groups in the groups list
                   outFolder=".",         # wriate output to this folder
                   mixBams=False,         # use one file for all bams
                   mixGroups=False,       # use one file for all groups
                   mixReads=False,        # use one file for paired / unpaired reads
                   interleaved=False,         # use interleaved format for paired reads
                   ignoreUnpaired=False,  # ignore all upaired reads
                   bigFile=False,         # do NOT gzip outputs
                   headersOnly=False      # write read headers only
                   ):
        """Set all instance variables, make the ReadSets"""
        # make sure the output folder exists
        self.outFolder = outFolder
        self.makeSurePathExists(self.outFolder) # it's a waste if we abort but I like to check if write permissions are intact before I do lots of work.

        self.bamFiles = bamFiles
        self.prettyBamFileNames = []
        for bam in self.bamFiles:
            self.prettyBamFileNames.append(os.path.basename(bam).replace(".bam", ""))
        self.prefix = prefix

        self.mixBams = mixBams
        self.mixGroups = mixGroups
        self.mixReads = mixReads

        self.ignoreUnpaired = ignoreUnpaired
        self.interleaved = interleaved
        if headersOnly:
            self.headersOnly = 1
        else:
            self.headersOnly = 0

        # are we going to zip the output?
        if bigFile:
            self.zipped = False
        else:
            self.zipped = True

        # munge the groups
        if groupNames == []:
            # no names specified, just use "group_1", "group_2" etc...
            groupNames = ["group_%d" % i for i in range(1, len(contigs)+1)]
        self.groupNames = groupNames

        # initialise to the first set of groups
        self.contigs = contigs[0]
        self.perContigGroups = [0]*len(self.contigs)

        for i in range(1, len(contigs)):
            self.contigs += contigs[i]
            self.perContigGroups += [i] * len(contigs[i])

        self.manager = Manager()

        # make sure printing to stdout is handled in a threadsafe manner
        self.outputStream = sys.stderr
        self.printQueue = self.manager.Queue()
        self.printDelay = 0.5   # delay between checking for new print statements

        self.RSM = ReadSetManager(self.manager)

        # make sure the RSM can talk to us
        self.RSM.setPrintQueue(self.printQueue)

        self.outFilePrefixes = self.RSM.organiseOutFiles(self.prettyBamFileNames,
                                                         self.groupNames,
                                                         self.zipped,
                                                         self.interleaved,
                                                         self.ignoreUnpaired,
                                                         self.mixBams,
                                                         self.mixGroups,
                                                         self.mixReads,
                                                         self.outFolder,
                                                         self.prefix)

    def extract(self, threads=1, verbose=False):
        '''Start extracting reads from the BAM files

        This function is responsible for starting and stopping all threads and
        processes used in bamm extract. Due to python multiprocessing's need to
        pickle everything the actual work of extraction is carried out in the
        first level function called externalExtractWrapper. See there for actual
        extraction details. This function is primarily concerned with thread
        and process management.

        Inputs:
         threads - int, the number of threads / processes to use
         verbose - bool, True if lot's of stuff should be printed to screen

        Outputs:
         None
        '''
        # make a queue containing all the bids to extract reads from
        extract_queue = self.manager.Queue()
        for bid in range(len(self.bamFiles)):
            extract_queue.put(bid)

        # place one None on the extract queue for each thread we have access to
        # AKA poison pill
        for _ in range(threads):
            extract_queue.put(None)

        # each thread gets a unique identifier
        thread_ids = ["Thread_%s" % str(tt) for tt in range(threads)]

        # start the Queue management processes and threads

        # printing process
        print_process = Process(target=self.managePrintQueue)
        print_process.start()

        # several threads for writing to disk
        request_management_threads = [Thread(target=self.RSM.manageRequests)]
        #                                   args=(self.groupNames, verbose)) for _ in range(threads)]
        for w in request_management_threads:
            w.start()

        # each thread gets its own queue for recieving ReadSet instances on
        response_queues = dict(zip(thread_ids,
                                   [self.manager.Queue() for _ in range(threads)]
                                   )
                               )
        # The RSM is waiting wor this queue too
        self.RSM.setResponseQueues(response_queues)

        # start the machine
        try:
            # make the extraction processes
            extract_proc = [Process(target=externalExtractWrapper,
                                    args=(thread_ids[tt],
                                          self.outFilePrefixes,
                                          self.bamFiles,
                                          self.prettyBamFileNames,
                                          len(self.groupNames),
                                          self.perContigGroups,
                                          self.contigs,
                                          self.printQueue,
                                          extract_queue,
                                          self.RSM.requestQueue,
                                          self.RSM.freeQueue,
                                          response_queues[thread_ids[tt]],
                                          self.headersOnly,
                                          self.mixGroups,
                                          verbose
                                          )
                                    )
                            for tt in range(threads)]

            # start the extraction processes
            for p in extract_proc:
                p.start()

            for p in extract_proc:
                p.join()

            # stop any rogue file writing action
            self.RSM.invalidateThreads()
            for w in request_management_threads:
                self.RSM.requestQueue.put(None)
            for w in request_management_threads:
                w.join()

            # stop the printer
            self.printQueue.put(None)
            print_process.join()

            # success
            return 0

        except:
            # ctrl-c! Make sure all processes are terminated
            sys.stderr.write("\nEXITING...\n")

            for p in extract_proc:
                p.terminate()

            # stop any rogue file writing action
            self.RSM.invalidateThreads()
            for w in request_management_threads:
                self.RSM.requestQueue.put(None)
            for w in request_management_threads:
                w.join()

            # stop the printer
            print_process.terminate()

            # dismal failure
            return 1

    def managePrintQueue(self):
        '''Write all the print requests to stdout / stderr

        This function is run as a process and so can be terminated.
        Place a None on the printQueue to terminate the process.

        Change self.outputStream to determine where text will be written to.

        Inputs:
         None

        Outputs:
         None
        '''
        while True:
            stuff = self.printQueue.get(timeout=None, block=True)
            if stuff is None:
                break
            else:
                self.outputStream.write("%s\n" % stuff)

    def makeSurePathExists(self, path):
        '''Make sure that a path exists, make it if necessary

        Inputs:
         path - string, full or relative path to create

        Outputs:
         None
        '''
        try:
            os.makedirs(path)
        except OSError as exception:
            import errno
            if exception.errno != errno.EEXIST:
                raise

