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
__version__ = "1.0.0-b.1"
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

def checkForDatabase(databaseBaseName):
    '''Check to see if bwa indexed database is present

    Inputs:
     databaseBaseName - string, full path to the prefix of the DB

    Outputs:
     True if database files exists, false otherwise
    '''
    return os.path.isfile(databaseBaseName+'.amb')

class BamScheduler:
    '''Class for scheduling making multiple BAM and TAM files

    The class implements a somewhat parallel interface for running multiple
    BWA jobs. It is a little slack in that it sequentially calls BWA with BWA's
    threading flag. On larger systems this will probably be slower than running
    gnu parallel through a bash script. but if you can do that then you probablu
    don't need this.
    '''

    def __init__(self,
                 database,
                 alignmentAlgorithm,
                 indexAlgorithm,
                 outFolder,
                 paired=[],
                 interleaved=[],
                 singleEnded=[],
                 keptFiles=False,
                 keepFiles=False,
                 outputTam=False,
                 prefix='',
                 numThreads=1,
                 maxMemory=None,
                 forceOverwriting=False,
                 verbose=False
                 ):
        '''Default constructor.

        Initializes a BamScheduler instance with the provided set of properties.

        Inputs:
         database - full path tofasta file of contigs (may be gzipped),
         alignmentAlgorithm - one of BWA's alignment algorithms,
         indexAlgorithm - one of BWA's index algorithms,
         paired - [fileName pairs], always even numbered in length, in order
                  [R_A_1, R_A_2, R_B_1, R_B_2, ...]
         interleaved - [fileNames], containing interleaved paired reads
         singleEnded - [fileNames], containing single ended reads
         keptFiles - == True -> indexes for the db already exist,
         keepFiles - == True -> don't delete indexes at the end,
         outputTam - == True -> you love text files to bits,
         numThreads - int, the maximum number of threads to use
         maxMemory - int, maximum memory program will use
         forceOverwriting - == True -> force overwriting index files,
         verbose - == True -> be verbose

        Outputs:
         None
        '''
        # the main thing is to make sure that the input parameters make sense
        self.outFolder = outFolder
        self.database = database
        self.dbBaseName = self.stripFaCrud(os.path.basename(database))
        self.prefix = prefix
        self.paired = paired
        self.interleaved = interleaved
        self.singleEnded = singleEnded
        if self.database is None:
            raise InvalidParameterSetException('Nothing to map reads onto, ' \
                                               'you need to supply a database')
        if self.singleEnded == [] and \
           self.paired == [] and \
           self.interleaved == []:
            raise InvalidParameterSetException( \
                'Nothing to map, please specify coupled, interleaved or ' \
                'single ended reads files')

        self.alignmentAlgorithm = alignmentAlgorithm
        self.indexAlgorithm = indexAlgorithm
        self.keptFiles = keptFiles
        self.keepFiles = keepFiles
        self.numThreads = str(numThreads)
        self.maxMemory = maxMemory
        self.outputTam = outputTam
        self.forceOverwriting = forceOverwriting

        if self.maxMemory is None:
            #Default to 2GBs per number of threads
            self.maxMemory = str(self.numThreads*2)+'G'

        # --kept sanity check
        if checkForDatabase(self.database):
            # dbs are there, has the user specified 'kept'
            if self.keptFiles is False and not self.forceOverwriting:
                raise InvalidParameterSetException( \
                    "You didn't specify that index files have been kept but " \
                    "there appears to be bwa index files present.\nI'm " \
                    "cowardly refusing to run so as not to risk overwriting.\n"\
                    "Force overwriting to create new indices")

        elif self.keptFiles:
            # user specified 'kept' but there are no DBs there
            raise InvalidParameterSetException( \
                "You specified that index files have been kept but there " \
                "doesn't appear to be any suitable bwa index files present")

        self.BMs = []
        self.outFiles = {}  # use this to check for repeated output files

        # the BamMaker class can check validity of it's own paramters quite well
        # we can just make a list of these guys here and then set it all going
        # when we're sure everything is going to be OK.

        # paired first
        l_paired = len(self.paired)
        if l_paired % 2 != 0:
            raise InvalidParameterSetException( \
                "Use of the -p option requires an even number of reads " \
                "(ordered as pairs)")

        for p_index in range(l_paired/2):
            # make the output file name and check that it's going to be unique
            out_file = self.makeOutFileName(self.paired[2*p_index])
            if out_file in self.outFiles:
                raise InvalidParameterSetException( \
                    'Output filename: %s for read set (Paired: %s %s) '\
                    'conflicts with previously calculated filename for ' \
                    'read set %s',
                    out_file,
                    self.paired[2*p_index],
                    self.paired[2*p_index+1],
                    self.outFiles[out_file])

            self.outFiles[out_file] = "(Paired: %s %s)" % \
                (self.paired[2*p_index], self.paired[2*p_index+1])
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
                raise InvalidParameterSetException( \
                    'Output filename: %s for read set (Interleaved: %s) ' \
                    'conflicts with previously calculated filename for ' \
                    'read set %s',
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
                raise InvalidParameterSetException( \
                    'Output filename: %s for read set (Single: %s) conflicts ' \
                    'with previously calculated filename for read set %s',
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
        '''sequentially make the bam files

        Inputs:
         None

        Outputs:
         None
        '''
        for BM in self.BMs:
            BM.makeBam()

    def stripFaCrud(self, fileName):
        # strip off the ".fa, .fa.gz etc from the end of the reads file
        out_file_name = fileName.replace(".fasta.gz","") \
                                .replace(".fa.gz","") \
                                .replace(".fq.gz","") \
                                .replace(".fastq.gz","") \
                                .replace(".fna.gz","")

        out_file_name = out_file_name.replace(".fasta","") \
                                     .replace(".fa","") \
                                     .replace(".fq","") \
                                     .replace(".fastq","") \
                                     .replace(".fna","")
        return out_file_name

    def makeOutFileName(self, readsFile):
        '''Consistent way to make output file name prefixes

        Inputs:
         readsFile - path-free read file name

        Outputs:
         pretty prefix
        '''
        return os.path.join(self.outFolder,
                            self.prefix + self.dbBaseName + '.' + \
                            os.path.basename(self.stripFaCrud(readsFile)))

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamMaker:
    '''Class that wraps bwa and samtools. Can call shell
    commands and make BAM and TAM files

    Not threaded, use the BamSceduler to start multiple BAM runs in parallel
    '''
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
        '''Default constructor.

        Initializes a BamMaker instance with the provided set of properties.

        Inputs:
         database - full path to fasta file of contigs (may be gzipped),
         alignmentAlgorithm - one of BWA's alignment algorithms,
         indexAlgorithm - one of BWA's index algorithms,
         outFileName - string, what to name the BAM
         readFile1 - string, full path to 1st read file
         readFile2 - string, full path to 2nd read file or None
         interleaved - == True -> readFile1 is interleaved pairs
         singleEnded - == True -> readFile1 is single ended reads
         keptFiles - == True -> indexes for the db already exist,
         keepFiles - == True -> don't delete indexes at the end,
         outputTam - == True -> you love text files to bits,
         numThreads - int, the maximum number of threads to use
         maxMemory - int, maximum memory program will use
         forceOverwriting - == True -> force overwriting index files,

        Outputs:
         None
        '''
        self.database = database
        self.readFile1 = readFile1
        self.readFile2 = readFile2
        self.outFileName = outFileName
        self.outputTam = outputTam
        self.forceOverwriting = forceOverwriting

        if self.database is None or \
           self.readFile1 is None or \
           self.outFileName is None:
            raise InvalidParameterSetException( \
                'You need to specify a multiple fasta file (database), ' \
                'at least one read file and an output file')

        self.isInterleaved = interleaved
        self.isSingleEnded = singleEnded

        self.alignmentAlgorithm = alignmentAlgorithm
        self.indexAlgorithm = indexAlgorithm

        self.keptFiles = keptFiles
        self.keepFiles = keepFiles

        self.numThreads = str(numThreads)
        self.maxMemory = maxMemory
        if self.maxMemory is None:
            #Default to 2GBs per number of threads
            self.maxMemory = str(self.numThreads*2)+'G'


        # intermediate files used during BAM production
        self.sai1 = None
        self.sai2 = None

        # make sure the file names are purdy
        if self.outputTam:
            if self.alignmentAlgorithm == 'mem':
                raise InvalidParameterSetException( \
                    'Sorry, tam output file format for bwa-mem is not ' \
                    'supported at this time\n(though it relatively easy ' \
                    'to implement)')

            if not (self.outFileName.endswith('.sam') or \
                    self.outFileName.endswith('.tam')):
                self.outFileName += '.tam'
        else:
            # bam output
            if self.outFileName.endswith('.bam'):
                # samtools renames the output file with .bam already
                self.outFileName = self.outFileName[:-4]


        # make sure the parameters make sense

        # we know we have a database and one read file
        if self.isSingleEnded:
            # single ended!
            if self.isInterleaved:
                raise InvalidParameterSetException( \
                    'The interleaved option is incompatible ' \
                    'with single ended alignment')

            elif self.readFile2 is not None:
                raise InvalidParameterSetException( \
                    'Two read files specified for a single ended alignment')
        else:
            # paired reads
            if self.readFile2 is None:
                # only OK if interleaved is set
                if self.isInterleaved is False:
                    raise InvalidParameterSetException( \
                        'You must specify two read files and a database ' \
                        'or one read file with the interleaved flag for a ' \
                        'paired alignment or explicitly use the single ' \
                        'ended flag')

                elif self.alignmentAlgorithm=='aln':
                    raise InvalidParameterSetException( \
                        'You cannot use the "aln" algorithm '\
                        'with interleaved reads')

            # else we have -d, -1 and -2 --> OK

        # last sanity check
        if (self.keptFiles is False and checkForDatabase(self.database)) and \
            not self.forceOverwriting:
            raise InvalidParameterSetException(\
                "You didn't specify that index files have been kept but " \
                "there appears to be bwa index files present.\nI'm cowardly " \
                "refusing to run so as not to risk overwriting.\n" \
                "Force overwriting to create new indices")

        # OK we, know what we're doing

    #---------------------------------------------------------------
    # database management

    def makeDatabase(self):
        ''' call BWA make index command

        Inputs:
         None

        Outputs:
         None
        '''
        sys.stderr.write('making database'+"\n")
        sys.stderr.flush
        if self.indexAlgorithm is None:
            subprocess.check_call(' '.join(['bwa index',
                                            self.database]),
                                  shell=True)
        else:
            subprocess.check_call(' '.join(['bwa index -a',
                                            self.indexAlgorithm,
                                            self.database]),
                                  shell=True)

    def removeDatabase(self):
        ''''Remove any index files that are no longer needed

        Inputs:
         None

        Outputs:
         None
        '''
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
        '''Use BWA ans samtools to make a BAM/TAM file

        Inputs:
         None

        Outputs:
         None
        '''
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
        '''call bwa aln

        Inputs:
         readFile - full path to reads file
         saiFile - full path to corresponding sai file

        Outputs:
         None
        '''
        subprocess.check_call(' '.join(['bwa aln -t',
                                        self.numThreads,
                                        self.database,
                                        readFile,
                                        '>',
                                        saiFile]),
                              shell=True)

    def sampe(self):
        '''call bwa sampe

        Inputs:
         None

        Outputs:
         None
        '''
        subprocess.check_call(' '.join(['bwa sampe',
                                        self.database,
                                        self.sai1,
                                        self.sai2,
                                        self.readFile1,
                                        self.readFile2,
                                        '>',
                                        self.outFileName]),
                              shell=True)

    def samse(self):
        '''call bwa samse

        Inputs:
         None

        Outputs:
         None
        '''
        subprocess.check_call(' '.join(['bwa samse',
                                        self.database,
                                        self.sai1,
                                        self.readFile1,
                                        '>',
                                        self.outFileName]),
                              shell=True)

    def sampe_to_sorted_indexed_bam(self):
        '''call bwa sampe and sort + index the result

        Inputs:
         None

        Outputs:
         None
        '''
        cmd = ' '.join(['bwa sampe',
                        self.database,
                        self.sai1,
                        self.sai2,
                        self.readFile1,
                        self.readFile2])

        cmd += ' '.join([' | samtools view -SubhF 4 - | samtools sort -m',
                         self.maxMemory,
                         '-',
                         self.outFileName])

        subprocess.check_call(cmd, shell=True)
        self.samtoolsIndex(self.outFileName)

    def samse_to_sorted_indexed_bam(self):
        '''call bwa samse and sort + index the result

        Inputs:
         None

        Outputs:
         None
        '''
        cmd = ' '.join(['bwa samse',
                        self.database,
                        self.sai1,
                        self.readFile1])

        cmd += ' '.join([' | samtools view -SubhF 4 - | samtools sort -m',
                         self.maxMemory,
                         '-',
                         self.outFileName])

        subprocess.check_call(cmd, shell=True)
        self.samtoolsIndex(self.outFileName)

    #---------------------------------------------------------------
    # mem algorithm

    def mem_single_to_sorted_indexed_bam(self):
        ''' call bwa mem mapping with single ended reads and sort + index

        Inputs:
         None

        Outputs:
         None
        '''
        cmd = ' '.join(['bwa mem -t',
                        self.numThreads,
                        self.database,
                        self.readFile1])

        cmd += ' '.join([' | samtools view -SubhF 4 - | samtools sort -m',
                         self.maxMemory,
                         '-',
                         self.outFileName])

        subprocess.check_call(cmd, shell=True)
        self.samtoolsIndex(self.outFileName)

    def mem_to_sorted_indexed_bam(self):
        ''' call bwa mem mapping with paired reads and sort + index the result

        Assume -p for bwa if reads2 is None, otherwise specify reads1 and reads2

        Inputs:
         None

        Outputs:
         None
        '''
        bwa_cmd = ' '.join(['bwa mem -t',
                            self.numThreads,
                            self.database,
                            ''])
        if self.isInterleaved:
          bwa_cmd += ' '.join(['-p',self.readFile1])
        else:
          bwa_cmd += ' '.join([self.readFile1,self.readFile2])

        cmd = bwa_cmd + ' '.join([' | samtools view -SubhF 4 -',
                                  '| samtools sort -m',
                                  self.maxMemory,
                                  '-',
                                  self.outFileName])

        subprocess.check_call(cmd, shell=True)
        self.samtoolsIndex(self.outFileName)

    #---------------------------------------------------------------
    # bwasw algorithm

    def bwasw(self):
        '''call bwasw

        Inputs:
         None

        Outputs:
         None
        '''
        if self.isSingleEnded:
            subprocess.check_call(' '.join(['bwa bwasw -t',
                                            self.numThreads,
                                            self.database,
                                            self.readFile1,
                                            '>',
                                            self.outFileName]),
                                  shell=True)
        else:
            subprocess.check_call(' '.join(['bwa bwasw -t',
                                            self.numThreads,
                                            self.database,
                                            self.readFile1,
                                            self.readFile2,
                                            '>',
                                            self.outFileName]),
                                  shell=True)

    def bwasw_to_sorted_indexed_bam(self):
        '''call bwasw and sort and index the result

        Inputs:
         None

        Outputs:
         None
        '''
        cmd = ' '.join(['bwa bwasw -t',
                        self.numThreads,
                        self.database,
                        self.readFile1])

        if not self.isSingleEnded:
            cmd += ' ' + self.readFile2

        cmd += ' '.join([' | samtools view -SubhF 4 - | samtools sort -m',
                         self.maxMemory,
                         '-',
                         self.outFileName])

        subprocess.check_call(cmd, shell=True)
        self.samtoolsIndex(self.outFileName)

    #---------------------------------------------------------------
    # index a bam file

    def samtoolsIndex(self, sortedBamFile):
        '''call samtools index on a sorted BAM file

        Inputs:
         sortedBamFile - full path to sorted bam file

        Outputs:
         None
        '''
        # samtools index cannot be piped, so a tmpfile is required
        subprocess.check_call(' '.join(['samtools index',
                                        sortedBamFile+'.bam']),
                              shell=True)

    #---------------------------------------------------------------
    # utilities

    def safeRemove(self, fileName):
        '''Delete a file without raising an exception

        Inputs:
         fileName - full path to file to be removed

        Outputs:
         None
        '''
        if os.path.isfile(fileName):
            os.system('rm ' + fileName)

    def __str__(self):
        '''Print out the operations and outputs that will be made

        Used when the scheduler is in verbose mode

        Inputs:
         None

        Outputs:
         None
        '''
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

        str += "  Database: %s\n  Output: %s%s\n  Threads: %s\n" % \
            (self.database, self.outFileName, suffix, self.numThreads)
        str += "  Alignment algorithm: %s" % self.alignmentAlgorithm
        return str

###############################################################################
###############################################################################
###############################################################################
###############################################################################
