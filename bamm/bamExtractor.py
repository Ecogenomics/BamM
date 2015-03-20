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
__copyright__ = "Copyright 2014,2015"
__credits__ = ["Michael Imelfort"]
__license__ = "LGPLv3"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"

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
from cWrapper import (BM_mappedRead_C,
                      CWrapper,
                      RPI,
                      RPIConv,
                      RPI2Str,
                      MI)
from bammExceptions import MixedFileTypesException
from bamRead import ReadSetManager

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
                           minMapQual,
                           maxMisMatches,
                           ignoreSuppAlignments,
                           ignoreSecondaryAlignments,
                           verbose=False
                           ):
    '''Single-process BAMfile read extraction.

    cTypes pointers are unpickleable unless they are top level, so this function
    lives outside the class and has 1,000,000 member variables passed to it.
    Life would be easier if we could pass the class but any implicit copy
    operations that follow are somewhat difficult to detect and can cause WOE.
    Lot's of WOE, believe me...

    Inputs:
     threadId - string, a unique Id for this process / thread
     outFilePrefixes - 3D dict for finding outFilePrefixes based on bamFile,
                       group and pairing information
     bamPaths - { bid : string }, full paths to the BAM files
     prettyBamFileNames - { bid : string }, short, print-friendly BAM names
     numGroups - int, the number of groups reads are split into
     perContigGroups - [int], contains groups Ids, insync with contigs array
     contigs - [string], contig ids as written in the BAM
     printQueue - Manager.Queue, thread-safe communication with users
     extractQueue - Manager.Queue, bids (BAMs) yet to be extracted from
     requestQueue - Manager.Queue, make requests for ReadSets for printing
     freeQueue - Manager.Queue, tell the RSM when finished with a ReadSet
     responseQueue - Manager.Queue, recieve copies of ReadSets from the RSM
     headersOnly - == True -> write read headers only
     mixGroups - == True -> use one file for all groups
     minMapQual - int, skip all reads with a lower mapping quality score
     maxMisMatches - int, skip all reads with more mismatches (NM aux files)
     useSuppAlignments - == True -> skip supplementary alignments
     useSecondaryAlignments - == True -> skip secondary alignments
     verbose - == True -> be verbose

    Outputs:
     None
    '''
    while True:
        p_bid = extractQueue.get(block=True, timeout=None)
        if p_bid is None: # poison pill
            break
        else:
            if verbose:
                printQueue.put("%s Preparing to extract reads from file: %s" % \
                                (threadId, prettyBamFileNames[p_bid] ) )

            # first we need to C-ify variables
            bamfile_c = c.c_char_p()
            bamfile_c = bamPaths[p_bid]

            pretty_name_c = c.c_char_p()
            pretty_name_c = prettyBamFileNames[p_bid]

            num_contigs = len(contigs)
            contigs_c_array = (c.c_char_p * num_contigs)()
            contigs_c_array[:] = contigs

            groups_c_array = (c.c_uint16 * num_contigs)()
            groups_c_array[:] = perContigGroups

            headers_only_c = c.c_uint32()
            headers_only_c = headersOnly

            min_mapping_quality_c = c.c_uint32()
            min_mapping_quality_c = minMapQual

            max_mismatches_c = c.c_uint32()
            max_mismatches_c = maxMisMatches

            pBMM = c.POINTER(BM_mappedRead_C)

            # call the C function to extract the reads
            CW = CWrapper()
            pBMM = CW._extractReads(bamfile_c,
                                    contigs_c_array,
                                    num_contigs,
                                    groups_c_array,
                                    pretty_name_c,
                                    headers_only_c,
                                    min_mapping_quality_c,
                                    max_mismatches_c,
                                    ignoreSuppAlignments,
                                    ignoreSecondaryAlignments)

            if verbose:
                printQueue.put("%s Finished C-based extraction for: %s" \
                               % (threadId, prettyBamFileNames[p_bid]))
                printQueue.put("%s Re-ordering reads before printing" % \
                               (threadId))

            # pBMM is one large linked list consisting of all mapped reads that
            # could be extracted from the BAM file. We have information about
            # the group and rpi of each read. The destination for each read is
            # encapsulated in the structure of the chain_info hash and
            # corresponding "storage" hash. We will re-order the linked list so
            # that adjacent connections indicate adjacency in the output file.
            # This is done by setting the "nextPrintRead" pointer in each BMM

            overlapper = {} # keep track of readSets with the same filename
            chain_info = {} # store start / end and count of a printing chain

            # initialise the helper data structures
            for gid in range(numGroups):
                chain_info[gid] = {}
                # ReadSets exist for only FIR and SNGL
                for rpi in [RPI.FIR, RPI.SNGL]:
                    file_name = outFilePrefixes[p_bid][gid][rpi]
                    try:
                        storage = overlapper[file_name]
                    except KeyError:
                        # [start of chain, end of chain, chain length]
                        storage = [None, None, 0]
                        overlapper[file_name] = storage
                    chain_info[gid][rpi] = {'storage' : storage}

            while pBMM:
                '''
                USE THIS CODE TO GET THE READ ID WHEN DEBUGGING
                buffer_c = c.create_string_buffer(20000)
                pbuffer_c = c.POINTER(c.c_char_p)
                pbuffer_c = c.pointer(buffer_c)

                # this variable records how much of the buffer is used for each read
                str_len_c = c.c_int(0)
                pstr_len_c = c.cast(c.addressof(str_len_c), c.POINTER(c.c_int))

                paired_c = c.c_int(1)
                headers = c.c_int(1)

                group_name_c = c.c_char_p()
                group_name_c = "THIS__"

                CW._sprintMappedRead(pBMM,
                                     pbuffer_c,
                                     pstr_len_c,
                                     group_name_c,
                                     headers,
                                     paired_c)
                # unwrap the buffer and transport into python land
                read_ID_debug = \
                    (c.cast(pbuffer_c,
                            c.POINTER(c.c_char*str_len_c.value)).contents).value

                read_ID_debug = read_ID_debug.split(";")[-1].rstrip()
                '''

                # get hold of the next item in the linked list
                rpi = pBMM.contents.rpi
                c_rpi = RPIConv[rpi]

                # we may need to add one or two reads, depending on pairing
                # always add pairs together to keep output files in sync
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
                        # we should treat as a single, unless we don't care
                        # i.e. (mixGroups == True)
                        if mixGroups:
                            # we don't care, print it now as a pair
                            # RPI.FIR and RPI.SEC ALWAYS point to same ReadSet
                            r2_rpi = RPI.FIR
                        else:
                            # we'll treat both paired reads as singles
                            r2_rpi = RPI.SNGL
                            addable[0][1] = RPI.SNGL # update this guy
                            # the storage for this rpi may remain == problems

                    addable.append([c.addressof((CW._getPartner(pBMM)).contents), r2_rpi])

                # update the printing chain
                for mappedRead in addable:
                    tmp_pBMM = c.cast(mappedRead[0], c.POINTER(BM_mappedRead_C))
                    has_qual = (tmp_pBMM.contents.qualLen != 0)
                    group = tmp_pBMM.contents.group

                    # set the MI code here
                    working_rpi = mappedRead[1]
                    stored_rpi = tmp_pBMM.contents.rpi
                    mi = MI.ER_EM_EG
                    if working_rpi == RPI.FIR:
                        mi = MI.PR_PM_PG
                    elif working_rpi == RPI.SNGL:
                        if stored_rpi == RPI.FIR or stored_rpi == RPI.SEC:
                            mi = MI.PR_PM_UG
                        if stored_rpi == RPI.SNGL_FIR or stored_rpi == RPI.SNGL_SEC:
                            mi = MI.PR_UM_NG
                        elif stored_rpi == RPI.SNGL:
                            mi = MI.UR_NM_NG

                    CW._setMICode(tmp_pBMM, mi)

                    #sys.stderr.write("%s -- %s\n" % (RPI2Str(working_rpi), RPI2Str(stored_rpi)))

                    # set and check the quality info
                    try:
                        # 'isFastq' is not set above, so on the first go
                        # this will raise a KeyError
                        if chain_info[group][working_rpi]['isFastq'] ^ has_qual:
                            # this will happen when people have merged BAMs with
                            # and without quality information
                            raise MixedFileTypesException( \
                                "You cannot mix Fasta and Fastq reads " \
                                "together in an output file")
                    except KeyError:
                        # Now we can set the type of the file.
                        # Only get here on the first read for each group, rpi
                        # Because of the way that the same storage object can be
                        # linked to multiple rpis, there's a chance that
                        # we won't set 'isFastq' for some rpis. Further down we
                        # need to be aware of this and just pass on the KeyError
                        # CODE==RPI_SKIP
                        chain_info[group][working_rpi]['isFastq'] = has_qual

                    # build or maintain the chain
                    if chain_info[group][working_rpi]['storage'][1] is None:
                        # this is the first time we've seen this print chain
                        chain_info[group][working_rpi]['storage'][0] = \
                            mappedRead[0]
                        chain_info[group][working_rpi]['storage'][1] = \
                            mappedRead[0]
                        chain_info[group][working_rpi]['storage'][2] = 1
                    else:
                        # join this pBMM onto the end of the existing chain
                        CW._setNextPrintRead( \
                            c.cast(chain_info[group][working_rpi]['storage'][1],
                                   c.POINTER(BM_mappedRead_C)),
                                             c.cast(mappedRead[0],
                                                    c.POINTER(BM_mappedRead_C)
                                                    )
                                             )
                        chain_info[group][working_rpi]['storage'][1] = \
                            mappedRead[0]
                        chain_info[group][working_rpi]['storage'][2] += 1

                # next!
                pBMM = CW._getNextMappedRead(pBMM)

            # Write the newly created chains to disk
            if verbose:
                printQueue.put("%s Re-ordering complete. Preparing to write" % \
                               (threadId))
            # search the chain_info hash for printable chains
            for gid in range(numGroups):
                for rpi in [RPI.FIR, RPI.SNGL]:
                    if chain_info[gid][rpi]['storage'][1] is not None:
                        # if we got here then there should be a chain to print
                        pBMM_chain = \
                            c.cast(chain_info[gid][rpi]['storage'][0],
                                   c.POINTER(BM_mappedRead_C)
                                   )
                        # we need to print here, so what we will do is make a
                        # request to the RSM for a fileName etc. that we can
                        # write to. We block on this call so we may have to
                        # wait for a bit BUT... it's either this, or go single
                        # threaded. So this is what we'll do.
                        try:
                            requestQueue.put((threadId,
                                              p_bid,
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
                                RS.writeChain(pBMM_chain,
                                              chain_info[gid][rpi]['isFastq'])
                                CW._destroyPrintChain(pBMM_destroy)

                                # free the RS now
                                freeQueue.put((threadId,
                                               p_bid,
                                               gid,
                                               rpi))

                            # set this to None so it's not added twice
                            chain_info[gid][rpi]['storage'][1] = None

                        except KeyError:
                            # this will happen when we have chosen to mix reads.
                            # it's no problem and I can't see that it hides any
                            # other bug. The "best" way to handle this is to set
                            # up a new variable that works out if we've set the
                            # 'isFastq' for a particular group and rpi. But this
                            # is really the same as checking chain_info[gid][rpi]
                            # for a KeyError here. So this is what we'll do...
                            # see: CODE==RPI_SKIP
                            pass

            if verbose:
                printQueue.put("%s Read extraction complete for file: %s" % \
                               (threadId, prettyBamFileNames[p_bid])
                               )

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamExtractor:
    '''Class used to manage extracting reads from multiple BAM files'''
    def __init__(self,
                 contigs,
                 bamFiles,
                 prefix="",
                 groupNames=[],
                 outFolder=".",
                 mixBams=False,
                 mixGroups=False,
                 mixReads=False,
                 interleaved=False,
                 bigFile=False,
                 headersOnly=False,
                 minMapQual=0,
                 maxMisMatches=1000,
                 useSuppAlignments=False,
                 useSecondaryAlignments=False,
                 ):
        '''
        Default constructor.

        Set all the instance variables, make ReadSets, organise output files

        Inputs:
         contigs - [[string]], list of list contig IDs (used as a filter)
         bamFiles - [string], list of bamfiles to extract reads from
         prefix - string, append this string to the start of all output files
         groupNames - [string], list of names of the groups in the contigs list
         outFolder - path, write output to this folder
         mixBams - == True -> use one file for all bams
         mixGroups - == True -> use one file for all groups
         mixReads - == True -> use one file for paired / unpaired reads
         interleaved - == True -> use interleaved format for paired reads
         bigFile - == True -> do NOT gzip outputs
         headersOnly - == True -> write read headers only
         minMapQual - int, skip all reads with a lower mapping quality score
         maxMisMatches - int, skip all reads with more mismatches (NM aux files)
         useSuppAlignments - == True -> DON'T skip supplementary alignments
         useSecondaryAlignments - == True -> DON'T skip secondary alignments

        Outputs:
         None
        '''
        # make sure the output folder exists
        self.outFolder = outFolder
        # it's a waste if we abort but I like to check if write permissions
        # are intact before I do lots of work.
        self.makeSurePathExists(self.outFolder)

        self.bamFiles = bamFiles
        self.prettyBamFileNames = [os.path.basename(bam).replace(".bam", "")
                                   for bam in self.bamFiles]

        self.prefix = prefix

        self.mixBams = mixBams
        self.mixGroups = mixGroups
        self.mixReads = mixReads

        self.interleaved = interleaved
        if headersOnly:
            self.headersOnly = 1
        else:
            self.headersOnly = 0

        self.minMapQual = minMapQual
        self.maxMisMatches = maxMisMatches

        if useSuppAlignments:
            self.ignoreSuppAlignments = 0
        else:
            self.ignoreSuppAlignments = 1

        if useSuppAlignments:
            self.ignoreSecondaryAlignments = 0
        else:
            self.ignoreSecondaryAlignments = 1

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
        self.printDelay = 0.5   # delay between checks for new print statements

        self.RSM = ReadSetManager(self.manager)

        # make sure the RSM can talk to us
        self.RSM.setPrintQueue(self.printQueue)

        self.outFilePrefixes= self.RSM.organiseOutFiles(self.prettyBamFileNames,
                                                        self.groupNames,
                                                        self.zipped,
                                                        self.interleaved,
                                                        self.mixBams,
                                                        self.mixGroups,
                                                        self.mixReads,
                                                        self.headersOnly,
                                                        self.outFolder,
                                                        self.prefix)


        '''
        for bid in range(len(self.bamFiles)):
            for gid in range(len(self.groupNames)):
                for rpi in [RPI.FIR, RPI.SEC, RPI.SNGL, RPI.SNGL_FIR, RPI.SNGL_SEC]:
                    sys.stderr.write("%s %s %s %s\n" % (self.prettyBamFileNames[bid], self.groupNames[gid], RPI2Str(rpi), str(self.outFilePrefixes[bid][gid][rpi])))
        '''

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
        request_management_threads = [Thread(target=self.RSM.manageRequests)
                                      for _ in range(threads)]
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
                                          self.minMapQual,
                                          self.maxMisMatches,
                                          self.ignoreSuppAlignments,
                                          self.ignoreSecondaryAlignments,
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

###############################################################################
###############################################################################
###############################################################################
###############################################################################