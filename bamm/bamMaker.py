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
__version__ = "1.0"
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

class BamMaker:
    """Neat utilities for making BAM and SAM files"""

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
                    raise InvalidParameterSetException('You must specify both -1 and -d, as well as -2 or -i for a paired alignment. For single ended just use -1, -d and -s')
                elif self.alignmentAlgorithm=='aln':
                    raise InvalidParameterSetException('You cannot use --bwa-aln with interleaved reads')
            # else we have -d, -1 and -2 --> OK

        # last sanity check
        if (self.keptFiles is False and self.checkForDatabase(self.database)) and not self.forceOverwriting:
            raise InvalidParameterSetException("You didn't specify --kept but there appears to be bwa index files present.\nI'm cowardly refusing to run so as not to risk overwriting.\nUse -f to force overwriting and -k to keep new indices")

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

    def checkForDatabase(self, database_basename):
        return os.path.isfile(database_basename+'.amb')

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
