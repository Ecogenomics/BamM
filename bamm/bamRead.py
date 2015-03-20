###############################################################################
#                                                                             #
#    bamRead.py                                                               #
#                                                                             #
#    Class for storing information about mapped reads                         #
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
import sys
import multiprocessing as mp
import gzip
import Queue
import time
from copy import deepcopy
import os

# local import
from bammExceptions import InvalidParameterSetException
from cWrapper import RPI, CWrapper, BM_mappedRead_C

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ReadSetManager(object):
    '''The principle manager of a collection of ReadSet objects

    Determines the exact names of the output files when extracting reads.
    Ensures that only one thread can write to one file at a time.
    '''

    def __init__(self, manager):
        '''Default constructor.

        Initializes a ReadSet instance with the provided set of properties.

        Inputs:
         manager - multiprocessing.Manager() instance owned by the
                   BamExtractor use this to make all Queues.
        Outputs:
         None
        '''
        # We use two data structures to manage ReadSets
        # self.outFiles is a nested hash that links bamFile Id
        # group Id and read-pair information (rpi) to a ReadSet
        # self.outFiles = { bamId : { groupId : { rpi : ReadSet } } }
        self.outFiles = {}
        # self.fnPrefix2ReadSet is a simple hash that links the
        # file name prefix to a ReadSet
        # self.fnPrefix2ReadSet = { outPrefix : ReadSet }
        self.fnPrefix2ReadSet = {}

        # The two variables ensure that only one thread has access to a
        # ReadSet at a time
        self.lock = manager.Lock()
        self.readSetInUse = {}

        # The program uses threads which cannot be terminated without
        # using a trick like altering a global variable.
        # NOTE: Use self.invalidateThreads to modify this variable
        self._threadsAreValid = True

        # The RSM communicates with the parsing threads using
        # a collection of queues.
        # The request queue is used by threads to make requests
        # for access to ReadSets. The freeQueue is used to indicate
        # that access to a resource is no loneger needed.
        # See self.manageRequests for details.
        self.requestQueue = manager.Queue()
        self.freeQueue = manager.Queue()

        # These Queues must be passed in by the BamExtractor
        # The RSM gives threads access to the ReadSets by placing
        # copies of them on the appropriate responseQueue.
        # self.responseQueues is a hash of Queues:
        # { outPrefix : mp.Manager().Queue() }
        self.responseQueues = None

        # All strings for printing to std out are placed on this queue
        # to ensure print statements are not garbled
        self.printQueue = None

    def invalidateThreads(self):
        '''Stop all the threads from running

        Also stops any looping processes in the ReadSets

        Inputs:
         None
        Outputs:
         None
        '''
        self._threadsAreValid = False
        for file_name in self.fnPrefix2ReadSet.keys():
            self.fnPrefix2ReadSet[file_name].threadsAreValid = False

    def setResponseQueues(self, queues):
        '''Set response queues used to pass ReadSets to the extract threads

        Inputs:
         queues - dict { outPrefix : mp.Manager().Queue() }, one for each
                  thread that's parsing BAM files for the BAmExtractor
        Outputs:
         None
        '''
        self.responseQueues = queues

    def setPrintQueue(self, queue):
        '''Set the print queue for communcating upwards

        Inputs:
         queue - mp.Manager().Queue(), print queue managed by the BamExtractor
        Outputs:
         None
        '''
        self.printQueue = queue

    def getReadSet(self, bid, gid, rpi, threadId):
        '''Ensure that only one thread at a time has access to a read set

        Inputs:
         bid - int, the unique identifier for a BAM file
         gid - int, the unique identifier for a target group
         rpi - enum, describes the type of read (RPI.FIR etc.)
         threadId - string, identifies the thread that is requesting
                    the resource.
        Outputs:
         The requested read set or None if it is not available
        '''
        read_set = self.outFiles[bid][gid][rpi]
        file_prefix = read_set.getConstFP()
        with self.lock:
            if file_prefix in self.readSetInUse:
                read_set = None
            else:
                self.readSetInUse[file_prefix] = threadId
        return read_set

    def freeReadSet(self, bid, gid, rpi, threadId):
        '''Indicate that a thread is finished with a read set

        Inputs:
         bid - int, the unique identifier for a BAM file
         gid - int, the unique identifier for a target group
         rpi - enum, describes the type of read (RPI.FIR etc.)
         threadId - string, identifies the thread that is freeing
                    the resource.
        Outputs:
         None

        Raises:
         InvalidParameterSetException - if the supplied information
                                        doesn't make sense
        '''
        read_set = self.outFiles[bid][gid][rpi]
        file_prefix = read_set.getConstFP()
        with self.lock:
            try:
                if self.readSetInUse[file_prefix] == threadId:
                    del self.readSetInUse[file_prefix]
                else:
                    raise InvalidParameterSetException( \
                         "%s owned by %s, not %s" % \
                             (file_prefix,
                              self.readSetInUse[file_prefix],
                              threadId))
            except KeyError:
                raise InvalidParameterSetException("%s not owned by anyone" % \
                                                   file_prefix)
                pass

    def manageRequests(self):
        '''Manage requests by parsing threads to access ReadSets

        This process runs on it's own thread and continues until it discovers
        a None on the requestQueue or self._threadsAreValid is set to False

        Inputs:
         None

        Outputs:
         None
        '''
        # loop on a global, so that way we can kill threads as needed
        while self._threadsAreValid:
            item = self.requestQueue.get(timeout=None, block=True)
            if item is None:
                break
            else:
                # is this a request for a readset
                (thread_id, bid, gid, rpi, is_fastq) = item
                return_queue = self.responseQueues[thread_id]
                read_set = None
                while read_set is None and self._threadsAreValid:
                    read_set = self.getReadSet(bid, gid, rpi, thread_id)

                    if read_set is None:
                        # the read_set is in use by another thread. Pop an entry
                        # off the top of the freeQueue and see if that helps
                        try:
                            (f_thread_id,
                             f_bid,
                             f_gid,
                             f_rpi) = self.freeQueue.get(block=True,
                                                         timeout=2)
                            try:
                                self.freeReadSet(f_bid,
                                                 f_gid,
                                                 f_rpi,
                                                 f_thread_id)
                            except InvalidParameterSetException: pass
                        except Queue.Empty:
                            # avoid wheel spinning
                            time.sleep(2)
                    # else, we have the RS and we're good to go

                # read_set should NOT be None here
                if read_set is None:
                    # free the thread
                    return_queue.put(None)
                    break

                # let's start trying to get ready for writing
                # filename should be set and fasta/fastq should be checked
                return_queue.put(deepcopy(read_set))

                # this read set is now "opened". This change will
                # be used on the next usage of the ReadSet
                if is_fastq:
                    read_set._fastqWritten = True
                else:
                    read_set._fastaWritten = True

    def organiseOutFiles(self,
                         prettyBamFileNames,
                         groupNames,
                         zipped,
                         interleaved,
                         mixBams,
                         mixGroups,
                         mixReads,
                         headersOnly,
                         outFolder,
                         prefix,
                         ):
        '''Determine all the outFile prefixes needed for extraction.

        The RSM manages a collection of output file objects called ReadSets
        These are made here and placed into hashes for retrieval later. This
        function also populates instance variables self.fnPrefix2ReadSet and
        self.outFiles; the two main ways that ReadSets can be accessed.

        File names vary wildly depending on the values of flags such as:
        mixGroups, mixBams etc. Take care when modifying this code.

        Inputs:
         prettyBamFileNames - [ string ] simple versions of the BAM filenames
         groupNames - [ string ] identify contig groups. EX: [bin1, bin2]
         zipped - bool, True if the output should be zipped
         interleaved - bool, True if the output should be interleaved
         mixBams - bool, True if BAM file origin should be ignored
         mixGroups - bool, True if group origin should be ignored
         mixReads - bool, True if (un)paired distinction should be ignored
         headersOnly - bool, True if only headers should be printed
         outFolder - string, folder to write output to
         prefix - string, prefix to apply to all output files

        Outputs:
         of_prefixes - { bamId : { groupId : { rpi : outFile prefix } } },
                       this hash can be used to get the file name prefix for a
                       read set based on BAM file origin, group origin and
                       pairing information.
        '''
        of_prefixes = {}
        base_path = os.path.join(os.path.abspath(outFolder), prefix)

        # we need to make a filename for every eventuality
        for bid in range(len(prettyBamFileNames)):
            if mixBams:
                bam_str = "allMapped"
            else:
                bam_str = prettyBamFileNames[bid]

            if bam_str not in self.outFiles:
                self.outFiles[bid] = {}
                of_prefixes[bid] = {}

            for gid in range(len(groupNames)):
                if mixGroups:
                    grp_str = "allGroups"
                else:
                    grp_str = groupNames[gid]

                if gid not in self.outFiles[bid]:
                    self.outFiles[bid][gid] = {}
                    of_prefixes[bid][gid] = {}

                if prefix == "":
                    fn = "%s%s.%s" % (base_path, bam_str, grp_str)
                else:
                    fn = "%s.%s.%s" % (base_path, bam_str, grp_str)

                if mixReads:
                    # all reads should go into the one file
                    working_fn = fn + ".allReads"
                    try:
                        read_set_P = self.fnPrefix2ReadSet[working_fn]
                    except KeyError:
                        read_set_P = ReadSet(groupNames,
                                             working_fn,
                                             zipped=zipped,
                                             headersOnly=headersOnly)
                        self.fnPrefix2ReadSet[working_fn] = read_set_P
                    read_set_S = read_set_P

                elif interleaved:
                    # one file for pairs and one for singles
                    working_fn_P = fn + ".pairedReads"
                    working_fn_S = fn + ".unpairedReads"
                    try:
                        read_set_P = self.fnPrefix2ReadSet[working_fn_P]
                    except KeyError:
                        read_set_P = ReadSet(groupNames,
                                             working_fn_P,
                                             paired=True,
                                             zipped=zipped,
                                             headersOnly=headersOnly)
                        self.fnPrefix2ReadSet[working_fn_P] = read_set_P

                    try:
                        read_set_S = self.fnPrefix2ReadSet[working_fn_S]
                    except KeyError:
                        read_set_S = ReadSet(groupNames,
                                             working_fn_S,
                                             zipped=zipped,
                                             headersOnly=headersOnly)
                        self.fnPrefix2ReadSet[working_fn_S] = read_set_S

                else:
                    # each in their own file
                    working_fn_1 = fn + ".1"
                    working_fn_2 = fn + ".2"
                    working_fn_S = fn + ".unpairedReads"
                    try:
                        read_set_P = self.fnPrefix2ReadSet[working_fn_1]
                    except KeyError:
                        read_set_P = ReadSet(groupNames,
                                             working_fn_1,
                                             outPrefix2=working_fn_2,
                                             paired=True,
                                             zipped=zipped,
                                             headersOnly=headersOnly)
                        self.fnPrefix2ReadSet[working_fn_1] = read_set_P

                    try:
                        read_set_S = self.fnPrefix2ReadSet[working_fn_S]
                    except KeyError:
                        read_set_S = ReadSet(groupNames,
                                             working_fn_S,
                                             zipped=zipped,
                                             headersOnly=headersOnly)
                        self.fnPrefix2ReadSet[working_fn_S] = read_set_S

                # we use the filenames to link everything up below
                self.outFiles[bid][gid][RPI.FIR] = read_set_P
                self.outFiles[bid][gid][RPI.SEC] = read_set_P
                of_prefixes[bid][gid][RPI.FIR] = read_set_P.getConstFP()
                of_prefixes[bid][gid][RPI.SEC] = read_set_P.getConstFP()

                self.outFiles[bid][gid][RPI.SNGL] = read_set_S
                self.outFiles[bid][gid][RPI.SNGL_FIR] = read_set_S
                self.outFiles[bid][gid][RPI.SNGL_SEC] = read_set_S
                of_prefixes[bid][gid][RPI.SNGL] = read_set_S.getConstFP()
                of_prefixes[bid][gid][RPI.SNGL_FIR] = read_set_S.getConstFP()
                of_prefixes[bid][gid][RPI.SNGL_SEC] = read_set_S.getConstFP()

        return of_prefixes

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ReadSet(object):
    '''Container class to encapsulate the concept of a(n un)paired read set

    This class manages file properties and writing functionality.'''

    def __init__(self,
                 groupNames,
                 outPrefix1,
                 outPrefix2 = None,
                 paired=False,
                 zipped=True,
                 headersOnly=False
                 ):
        '''Default constructor.

        Initializes a ReadSet instance with the provided set of properties.

        Inputs:
         groupNames - [ string ], target names, one for each target group.
         outPrefix1 - string prefix of the first output file (read1).
         outPrefix2 - string prefix of the second output file or None for
                      interleaved or unpaired files.
         isPaired - bool, True if the file is a paired read file.
         zipped - bool, True if the data should be compressed for writing.
         headersOnly - bool, True if only headers should be printed.
        '''
        # output read file properties
        self.groupNames = groupNames
        self.isPaired = paired
        self.zipOutput = zipped
        self.headersOnly = headersOnly

        # prefixes of the files we'll be writing to
        self._outPrefix1 = outPrefix1
        # prefix2 may well be None
        self._outPrefix2 = outPrefix2

        # work out how we'll write files
        if zipped:
            self._writeOpen = gzip.open
        else:
            self._writeOpen = open

        # If the RSM catches a ^C during a write operation
        # we need a global variable to exit and kill the thread
        # Loop as long as threads are valid
        self._threadsAreValid = True

        # we could potentially be writing both fasta and fastq
        # these variables keep track of such things and make sure that
        # the right data goes to the right place
        self._fastqWritten = False
        self._fastaWritten = False

    def getConstFP(self, fNumber=1):
        '''return the unchanging filename prefix associated with this ReadSet

        Inputs:
         fNumber - which file prefix is needed? [1 or 2]
        Outputs:
         The corresponding prefix
        '''
        if 1 == fNumber:
            return self._outPrefix1
        else:
            return self._outPrefix2

    def determineFileSuffix(self, isFastq):
        '''Determine the suffix of the file depending on if we're writing
        fasta or fastq

        Inputs:
         isFastq - bool, is the file fastq (or fasta)
        Outputs:
         Filenames to be written to.
         (fileName1, fileName2)
         fileName2 will be None for unpaired or interleaved-paired files
        '''
        if self.headersOnly:
            ext = ".list"
        else:
            if isFastq:
                ext = ".fq"
            else:
                ext = ".fna"

        if self.zipOutput:
            ext += ".gz"

        file_name1 = self.getConstFP()+ext
        file_name2 = self.getConstFP(fNumber=2)
        if file_name2 is not None:
            file_name2 += ext

        return (file_name1, file_name2)

    def writeChain(self,
                   pBMM,
                   isFastq,
                   printQueue=None
                   ):
        '''Write a single print chain to disk

        A print chain is a linked list of mapped reads that have been
        pre-ordered and are ready to write (or print). The print chain can
        contain either Fasta or Fastq reads but never both. File names are
        determined on the fly based on the presence or absence of quality info
        of the first read in the chain (determined by the BamExtractor) and
        passed to this function as isFastq.

        NOTE: This function does NOT free any memory associated with pBMM.

        Inputs:
         pBMM - c.POINTER(BM_mappedRead_C), the start of a linked list of
                mapped reads, pre-ordered for printing by the BamExtractor
         isFastq - bool, True if reads have quality information.
         printQueue - Managed by the BamExtractor. Place all printing strings
                      here. Acts as a verbose flag.
        Outputs:
         None
        '''
        CW = CWrapper()

        # reads are written (in C land) to this string buffer
        # is 20000 bases enough for PAC-bio?
        buffer_c = c.create_string_buffer(20000)
        pbuffer_c = c.POINTER(c.c_char_p)
        pbuffer_c = c.pointer(buffer_c)

        # this variable records how much of the buffer is used for each read
        str_len_c = c.c_int(0)
        pstr_len_c = c.cast(c.addressof(str_len_c), c.POINTER(c.c_int))

        paired_c = c.c_int(1)
        unpaired_c = c.c_int(0)
        headers = c.c_int(self.headersOnly)

        # buffer to hold the group name in C format
        # it's a bit of a waste of time to pass a string to C only to have it
        # passed right back, but this approach reduces complexity and makes the
        # C code more useful, so it's preferred.
        group_name_c = c.c_char_p()

        # get the fileNames to write to
        (out_file1, out_file2) = self.determineFileSuffix(isFastq)

        # determine file write mode. This instance is likely a copy
        # of the main one managed by the RSM. so there is no need
        # to update the value of self._fastXWritten here. Just use it.
        opened = False
        if isFastq and self._fastqWritten:
            opened = True
        elif not isFastq and self._fastaWritten:
            opened = True

        if opened:
            # we will append to an existing file
            open_mode = "a"
            mode_desc = "Appending to"
        else:
            # overwrite any existing file
            open_mode = "w"
            mode_desc = "Writing"

        if self.isPaired:
            # swap writing to file 1 and file 2.
            # always start writing to fh1 first!
            isFh1 = True

            # open files
            fh1 = self._writeOpen(out_file1, open_mode)
            if out_file2 is None:
                if printQueue:
                    printQueue.put(" %s interleaved file: %s" % (mode_desc,
                                                                 out_file1))
                fh2 = fh1
            else:
                if printQueue:
                    printQueue.put(" %s coupled files: %s %s" % (mode_desc,
                                                                 out_file1,
                                                                 out_file2))
                fh2 = self._writeOpen(out_file2, open_mode)

            # write
            while pBMM and self._threadsAreValid:
                # get C to write the read into the string buffer
                group_name_c = self.groupNames[pBMM.contents.group]
                CW._sprintMappedRead(pBMM,
                                     pbuffer_c,
                                     pstr_len_c,
                                     group_name_c,
                                     headers,
                                     paired_c)
                # unwrap the buffer and transport into python land
                printable_string = \
                    (c.cast(pbuffer_c,
                            c.POINTER(c.c_char*str_len_c.value)).contents).value
                if isFh1:
                    fh1.write(printable_string)
                    isFh1 = False
                else:
                    fh2.write(printable_string)
                    isFh1 = True

                # be sure that we're going to the next PRINT read
                pBMM = CW._getNextPrintRead(pBMM)

            # and close
            fh1.close()
            if out_file2 is not None:
                fh2.close()
        else:
            fh = self._writeOpen(out_file1, open_mode)
            if printQueue:
                printQueue.put(" %s unpaired file: %s (%s)" % (mode_desc,
                                                               out_file1,
                                                               self))
            while pBMM and self._threadsAreValid:
                group_name_c = self.groupNames[pBMM.contents.group]
                CW._sprintMappedRead(pBMM,
                                     pbuffer_c,
                                     pstr_len_c,
                                     group_name_c,
                                     headers,
                                     unpaired_c)
                printable_string = \
                  (c.cast(pbuffer_c,
                  c.POINTER(c.c_char*str_len_c.value)).contents).value
                fh.write(printable_string)
                pBMM = CW._getNextPrintRead(pBMM)

            fh.close()

###############################################################################
###############################################################################
###############################################################################
###############################################################################
