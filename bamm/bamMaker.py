#!/usr/bin/env python
###############################################################################
#                                                                             #
#    bamMaker.py                                                              #
#                                                                             #
#    Wrapper to produce a .[bt]am file with BWA for use in general chicanery  #
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
__credits__ = ["Michael Imelfort", "Ben Woodcroft", "Connor Skennerton"]
__license__ = "LGPLv3"
__version__ = "1.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Beta"

###############################################################################

# system imports
import subprocess
import os
import sys
import tempfile

# local imports
from bammExceptions import *

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def checkForDatabase(database_basename):
    return os.path.isfile(database_basename+'.amb')

class BamScheduler:
    """Schedule making multiple BAM and TAM files"""

    def __init__(self,
                 database,
                 alignmentAlgorithm,
                 indexAlgorithm,
                 paired=[],
                 interleaved=[],
                 singleEnded=[],
                 keptFiles=False,
                 keepFiles=False,
                 outputTam=False,
                 numThreads=1,
                 maxMemory=None,
                 forceOverwriting=False,
                 verbose=False
                 ):

        # the main thing is to make sure that the input paramters make sense
        self.database = database
        self.paired = paired
        self.interleaved = interleaved
        self.singleEnded = singleEnded
        if self.database is None:
            raise InvalidParameterSetException('Nothing to map reads onto, you need to supply a database')
        if self.singleEnded == [] and self.paired == [] and self.interleaved == []:
            raise InvalidParameterSetException("Nothing to map, please specify coupled, interleaved or single ended reads files")

        self.alignmentAlgorithm = alignmentAlgorithm
        self.indexAlgorithm = indexAlgorithm
        self.keptFiles = keptFiles
        self.keepFiles = keepFiles
        self.numThreads = str(numThreads)
        self.maxMemory = maxMemory
        self.outputTam = outputTam
        self.forceOverwriting = forceOverwriting

        if self.maxMemory is None:
          self.maxMemory = str(self.numThreads*2)+'G' #Default to 2GBs per number of threads

        # --kept sanity check
        if checkForDatabase(self.database):
            # dbs are there, has the user specified 'kept'
            if self.keptFiles is False and not self.forceOverwriting:
                raise InvalidParameterSetException("You didn't specify that index files have been kept but there appears to be bwa index files present.\nI'm cowardly refusing to run so as not to risk overwriting.\nForce overwriting to create new indices")
        elif self.keptFiles:
            # user specified 'kept' but there are no DBs there
            raise InvalidParameterSetException("You specified that index files have been kept but there doesn't appear to be any suitable bwa index files present")

        self.BMs = []
        self.outFiles = {}  # use this to check for repeated output files

        # the BamMaker class can check validity of it's own paramters quite well
        # we can just make a list of these guys here and then set it all going once
        # we're sure everything is going to be OK.

        # paired first
        l_paired = len(self.paired)
        if l_paired % 2 != 0:
            raise InvalidParameterSetException("Use of the -p option requires an even number of reads (ordered as pairs)")

        for p_index in range(l_paired/2):
            # make the output file name and check that it's going to be unique
            out_file = self.makeOutFileName(self.paired[2*p_index])
            if out_file in self.outFiles:
                raise InvalidParameterSetException('Output filename: %s for read set (Paired: %s %s) conflicts with previously calculated filename for read set %s',
                                                   out_file,
                                                   self.paired[2*p_index],
                                                   self.paired[2*p_index+1],
                                                   self.outFiles[out_file])
            self.outFiles[out_file] = "(Paired: %s %s)" % (self.paired[2*p_index], self.paired[2*p_index+1])
            BM = BamMaker(self.database,
                          self.alignmentAlgorithm,
                          self.indexAlgorithm,
                          out_file,
                          self.paired[2*p_index],
                          readFile2=self.paired[2*p_index+1],
                          interleaved=False,
                          singleEnded=False,
                          keptFiles=True,       # always keep these
                          keepFiles=True,
                          outputTam=self.outputTam,
                          numThreads=self.numThreads,
                          maxMemory=self.maxMemory,
                          forceOverwriting=self.forceOverwriting
                          )
            self.BMs.append(BM)

        # interleaved next
        for file in self.interleaved:
            out_file = self.makeOutFileName(file)
            if out_file in self.outFiles:
                raise InvalidParameterSetException('Output filename: %s for read set (Interleaved: %s) conflicts with previously calculated filename for read set %s',
                                                   out_file,
                                                   file,
                                                   self.outFiles[out_file])
            self.outFiles[out_file] = "(Interleaved: %s)" % (file)
            BM = BamMaker(self.database,
                          self.alignmentAlgorithm,
                          self.indexAlgorithm,
                          out_file,
                          file,
                          interleaved=True,
                          singleEnded=False,
                          keptFiles=True,
                          keepFiles=True,
                          outputTam=self.outputTam,
                          numThreads=self.numThreads,
                          maxMemory=self.maxMemory,
                          forceOverwriting=self.forceOverwriting
                          )
            self.BMs.append(BM)

        # singletons last
        for file in self.singleEnded:
            out_file = self.makeOutFileName(file)
            if out_file in self.outFiles:
                raise InvalidParameterSetException('Output filename: %s for read set (Single: %s) conflicts with previously calculated filename for read set %s',
                                                   out_file,
                                                   file,
                                                   self.outFiles[out_file])
            self.outFiles[out_file] = "(Single: %s)" % (file)
            BM = BamMaker(self.database,
                          self.alignmentAlgorithm,
                          self.indexAlgorithm,
                          out_file,
                          file,
                          interleaved=False,
                          singleEnded=True,
                          keptFiles=True,
                          keepFiles=True,
                          outputTam=self.outputTam,
                          numThreads=self.numThreads,
                          maxMemory=self.maxMemory,
                          forceOverwriting=self.forceOverwriting
                          )
            self.BMs.append(BM)

        # we've made it this far. Lets tell the user what we intend to do
        if verbose:
            for BM in self.BMs:
                print BM
            print "-------------------------------------------"

    def makeBams(self):
        """Make the bams"""
        for BM in self.BMs:
            BM.makeBam()

    def makeOutFileName(self, readsFile):
        """Consistent way to make output files"""
        # strip off the ".fa, .fa.gz etc from the end of the reads file
        out_file_name = readsFile.replace(".fasta.gz","").replace(".fa.gz","").replace(".fq.gz","").replace(".fastq.gz","").replace(".fna.gz","")
        out_file_name = out_file_name.replace(".fasta","").replace(".fa","").replace(".fq","").replace(".fastq","").replace(".fna","")
        return out_file_name

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamMaker:
    """Neat utilities for making BAM and TAM files"""

    def __init__(self,
                 database,
                 alignmentAlgorithm,
                 indexAlgorithm,
                 outFileName,
                 readFile1,
                 readFile2=None,
                 interleaved=False,
                 singleEnded=False,
                 keptFiles=False,
                 keepFiles=False,
                 outputTam=False,
                 numThreads=1,
                 maxMemory=None,
                 forceOverwriting=False
                 ):

        self.database = database
        self.readFile1 = readFile1
        self.readFile2 = readFile2
        self.outFileName = outFileName
        self.outputTam = outputTam
        self.forceOverwriting = forceOverwriting

        if self.database is None or self.readFile1 is None or self.outFileName is None:
            raise InvalidParameterSetException('You need to specify a multiple fasta file (database), at least one read file and an output file')

        self.isInterleaved = interleaved
        self.isSingleEnded = singleEnded

        self.alignmentAlgorithm = alignmentAlgorithm
        self.indexAlgorithm = indexAlgorithm

        self.keptFiles = keptFiles
        self.keepFiles = keepFiles

        self.numThreads = str(numThreads)
        self.maxMemory = maxMemory
        if self.maxMemory is None:
          self.maxMemory = str(self.numThreads*2)+'G' #Default to 2GBs per number of threads


        # intermediate files used during BAM production
        self.sai1 = None
        self.sai2 = None

        # make sure the file names are purdy
        if self.outputTam:
            if self.alignmentAlgorithm == 'mem':
                raise InvalidParameterSetException('Sorry, tam output file format for bwa-mem is not supported at this time\n(though it relatively easy to implement)')

            if not (self.outFileName.endswith('.sam') or self.outFileName.endswith('.tam')):
                self.outFileName += '.tam'
        else:
            # bam output
            if self.outFileName.endswith('.bam'):
                self.outFileName = self.outFileName[:-4] # samtools renames the output file with .bam already


        # make sure the parameters make sense

        # we know we have a database and one read file
        if self.isSingleEnded:
            # single ended!
            if self.isInterleaved:
                raise InvalidParameterSetException('The interleaved option is incompatible with single ended alignment')
            elif self.readFile2 is not None:
                raise InvalidParameterSetException('Two read files specified for a single ended alignment')
        else:
            # paired reads
            if self.readFile2 is None:
                # only OK if interleaved is set
                if self.isInterleaved is False:
                    raise InvalidParameterSetException('You must specify two read files and a database or one read file with the interleaved flag for a paired alignment or explicitly use the single ended flag')
                elif self.alignmentAlgorithm=='aln':
                    raise InvalidParameterSetException('You cannot use the "aln" algorithm with interleaved reads')
            # else we have -d, -1 and -2 --> OK

        # last sanity check
        if (self.keptFiles is False and checkForDatabase(self.database)) and not self.forceOverwriting:
            raise InvalidParameterSetException("You didn't specify that index files have been kept but there appears to be bwa index files present.\nI'm cowardly refusing to run so as not to risk overwriting.\nForce overwriting to create new indices")

        # OK we, know what we're doing

    #---------------------------------------------------------------
    # database management

    def makeDatabase(self):
        """Wrapper to BWA make index command"""
        sys.stderr.write('making database'+"\n")
        sys.stderr.flush
        if self.indexAlgorithm is None:
            subprocess.check_call('bwa index ' + self.database, shell=True)
        else:
            subprocess.check_call('bwa index -a ' + self.indexAlgorithm + ' ' + self.database, shell=True)

    def removeDatabase(self):
        """clean up any index files we no longer need"""
        sys.stderr.write('deleting indices'+"\n")
        sys.stderr.flush
        self.safeRemove(self.database+'.amb')
        self.safeRemove(self.database+'.ann')
        self.safeRemove(self.database+'.bwt')
        self.safeRemove(self.database+'.pac')
        self.safeRemove(self.database+'.rbwt')
        self.safeRemove(self.database+'.rpac')
        self.safeRemove(self.database+'.rsa')
        self.safeRemove(self.database+'.sa')

    #---------------------------------------------------------------
    # main wrapper

    def makeBam(self):
        """actually make the bam files"""
        # run the actual alignment
        if self.alignmentAlgorithm == 'bwasw':

            if self.outputTam:
                self.bwasw()
            else:
                self.bwasw_to_sorted_indexed_bam()

        elif self.alignmentAlgorithm == 'aln':

            self.sai1 = tempfile.mkstemp(suffix='.sai')[1]
            self.sai2 = tempfile.mkstemp(suffix='.sai')[1]

            # always do the first read
            self.aln(self.readFile1, self.sai1)
            if self.isSingleEnded:
                if self.outputTam:
                    self.samse()
                else:
                    self.samse_to_sorted_indexed_bam()
            else:
                self.aln(self.readFile2, self.sai2)
                if self.outputTam:
                    self.sampe()
                else:
                    self.sampe_to_sorted_indexed_bam()

            self.safeRemove(self.sai1)
            self.safeRemove(self.sai2)

        else:
            # algorithm is mem
            if self.isSingleEnded:
                self.mem_single_to_sorted_indexed_bam()
            else:
                self.mem_to_sorted_indexed_bam()

    #---------------------------------------------------------------
    # aln algorithm

    def aln(self, readFile, saiFile):
        subprocess.check_call('bwa aln -t ' + self.numThreads + ' ' + self.database + ' ' + readFile + ' > ' + saiFile, shell=True)

    def sampe(self):
        subprocess.check_call('bwa sampe ' + self.database + ' ' + self.sai1 + ' ' + self.sai2 + ' ' + self.readFile1 + ' ' + self.readFile2 + ' > ' + self.outFileName, shell=True)

    def samse(self):
        subprocess.check_call('bwa samse ' + self.database + ' ' + self.sai1 + ' ' + self.readFile1 + ' > ' + self.outFileName, shell=True)

    def sampe_to_sorted_indexed_bam(self):
        cmd =  'bwa sampe ' + self.database + ' ' + self.sai1 + ' ' + self.sai2 + ' ' + self.readFile1 + ' ' + self.readFile2
        cmd += ' | samtools view -SubhF 4 - | samtools sort -m ' + self.maxMemory + ' - ' + self.outFileName
        subprocess.check_call(cmd, shell=True)
        self.samtoolsIndex(self.outFileName)

    def samse_to_sorted_indexed_bam(self):
        cmd =  'bwa samse ' + self.database + ' ' + self.sai1 + ' ' + self.readFile1
        cmd += ' | samtools view -SubhF 4 - | samtools sort -m ' + self.maxMemory + ' - ' + self.outFileName
        subprocess.check_call(cmd, shell=True)
        self.samtoolsIndex(self.outFileName)

    #---------------------------------------------------------------
    # mem algorithm

    def mem_single_to_sorted_indexed_bam(self):
        """run bwa mem mapping with single ended reads"""
        cmd =  'bwa mem -t ' + self.numThreads + ' ' + self.database + ' ' + self.readFile1
        cmd += ' | samtools view -SubhF 4 - | samtools sort -m ' + self.maxMemory + ' - ' + self.outFileName
        subprocess.check_call(cmd, shell=True)
        self.samtoolsIndex(self.outFileName)

    def mem_to_sorted_indexed_bam(self):
        """run bwa mem. Assume -p for bwa if reads2 is None, otherwise specify reads1 and reads2"""
        bwa_cmd = 'bwa mem -t ' + self.numThreads + ' ' + self.database + ' '
        if self.isInterleaved:
          bwa_cmd += '-p ' + self.readFile1
        else:
          bwa_cmd += self.readFile1 + ' ' + self.readFile2

        cmd = bwa_cmd + ' | samtools view -SubhF 4 - | samtools sort -m ' + self.maxMemory + ' - ' + self.outFileName
        subprocess.check_call(cmd, shell=True)
        self.samtoolsIndex(self.outFileName)

    #---------------------------------------------------------------
    # bwasw algorithm

    def bwasw(self):
        if self.isSingleEnded:
            subprocess.check_call('bwa bwasw -t ' + self.numThreads + ' ' + self.database + ' ' + self.readFile1 + ' > ' + self.outFileName, shell=True)
        else:
            subprocess.check_call('bwa bwasw -t ' + self.numThreads + ' ' + self.database + ' ' + self.readFile1 + ' ' + self.readFile2 + ' > ' + self.outFileName, shell=True)

    def bwasw_to_sorted_indexed_bam(self):
        cmd =  'bwa bwasw -t ' + self.numThreads + ' ' + self.database + ' ' + self.readFile1
        if not self.isSingleEnded:
            cmd += ' ' + self.readFile2
        cmd += ' | samtools view -SubhF 4 - | samtools sort -m ' + self.maxMemory + ' - ' + self.outFileName
        subprocess.check_call(cmd, shell=True)
        self.samtoolsIndex(self.outFileName)

    #---------------------------------------------------------------
    # index a bam file

    def samtoolsIndex(self, sortedBamFile):
        # samtools index cannot be piped, so a tmpfile is required
        subprocess.check_call('samtools index ' + sortedBamFile + '.bam', shell=True)

    #---------------------------------------------------------------
    # utilities

    def safeRemove(self, fileName):
        if os.path.isfile(fileName):
            os.system('rm ' + fileName)

    def __str__(self):
        """Print out the operations and outputs that will be made"""
        str = "-------------------------------------------\n"
        if self.isSingleEnded:
            str += "  Input: Single ended\n"
            str += "    %s\n" % self.readFile1
        elif self.isInterleaved:
            str += "  Input: Interleaved\n"
            str += "    %s\n" % self.readFile1
        else:
            str += "  Input: Paired\n"
            str += "    %s, %s\n" % (self.readFile1, self.readFile2)

        if self.outputTam:
            suffix = ""
        else:
            suffix = ".bam (sorted + indexed)"

        str += "  Database: %s\n  Output: %s%s\n  Threads: %s\n" % (self.database, self.outFileName, suffix, self.numThreads)
        str += "  Alignment algorithm: %s" % self.alignmentAlgorithm
        return str
