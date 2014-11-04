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

import random
from math import isnan
from numpy import mean, median, std, array
from nose.tools import assert_equals, assert_true
import sys
import os
import json

from bamm.cWrapper import *

class TestBamParser:
    @classmethod
    def setup_class(self):
        """Setup class variables before any tests."""

        # the following files already exist
        self.fastaFile = 'MG1655.fa'
        self.readKey = 'readkey.csv'
        self.bins = ['bin1', 'bin2']

        # the following data files will be in created
        self.peFile = ['pe.1.fa', 'pe.2.fa']    # coupled
        self.mpFile = 'mp.fa'                   # interleaved
        self.upFile = 'up.fa'                   # single
        self.contigsFile = 'contigs.fa'         # contigs
        self.poFile = "predicted_outputs.json"  # predicted outputs

        # we need to make some test data based on out test model
        self.model_dir = os.path.join(os.path.split(__file__)[0], "modelling")
        self.model_data_dir = os.path.join(self.model_dir, "data")
        os.system("%s -o %s -r %s -f %s -g %s %s" % \
                  (os.path.join(self.model_dir,"makeTestData.py"),
                   self.model_dir,
                   os.path.join(self.model_data_dir, self.readKey),
                   os.path.join(self.model_data_dir, self.fastaFile),
                   os.path.join(self.model_data_dir, self.bins[0]),
                   os.path.join(self.model_data_dir, self.bins[1])
                   )
                  )

        # now we can load the predicted outputs
        with open(os.path.join(self.model_dir, self.poFile), "r") as fh:
            self.predictedOutputs = json.loads(fh.readline())

        # these files will exist after the mapping is done
        self.peBamFile = ['contigs.pe.1.bam', 'contigs.pe.1.bam.bai']
        self.mpBamFile = ['contigs.mp.bam', 'contigs.mp.bam.bai']
        self.upBamFile = ['contigs.up.bam', 'contigs.up.bam.bai']

        # predicted outputs is a hash with three keys:
        self.pOutKeys = ['coverages', 'links', 'extracts']

        # keep the names of some stuff we'll be testing
        self.pOutCovTypes = ['cmean',
                             'counts',
                             'pmean',
                             'opmean',
                             'tpmean',
                             'pmedian'
                             ]
        self.pOutContigs = ['A', 'B', 'C', 'E']

    @classmethod
    def teardown_class(self):
        self.rmTestFile(self.mpFile)
        self.rmTestFile(self.peFile[0])
        self.rmTestFile(self.peFile[1])
        self.rmTestFile(self.upFile)
        self.rmTestFile(self.contigsFile)
        self.rmTestFile(self.poFile)

        if os.path.exists(os.path.join(self.model_dir,
                                       self.peBamFile[0])):
            self.rmTestFile(self.peBamFile[0])
            self.rmTestFile(self.peBamFile[1])
            self.rmTestFile(self.mpBamFile[0])
            self.rmTestFile(self.mpBamFile[1])
            self.rmTestFile(self.upBamFile[0])
            self.rmTestFile(self.upBamFile[1])

    @classmethod
    def rmTestFile(self, file):
        os.system("rm %s" % os.path.join(self.model_dir, file))

    def test_make(self):
        cmd = 'bamm make -q -d %s -i %s -c %s %s -s %s -o %s' % (os.path.join(self.model_dir, self.contigsFile),
                                                                 os.path.join(self.model_dir, self.mpFile),
                                                                 os.path.join(self.model_dir, self.peFile[0]),
                                                                 os.path.join(self.model_dir, self.peFile[1]),
                                                                 os.path.join(self.model_dir, self.upFile),
                                                                 self.model_dir)
        os.system(cmd)
        EX = os.path.exists(os.path.join(self.model_dir,
                                         self.peBamFile[0]))
        EX &= os.path.exists(os.path.join(self.model_dir,
                                          self.peBamFile[1]))
        EX &= os.path.exists(os.path.join(self.model_dir,
                                          self.mpBamFile[0]))
        EX &= os.path.exists(os.path.join(self.model_dir,
                                          self.mpBamFile[1]))
        EX &= os.path.exists(os.path.join(self.model_dir,
                                          self.upBamFile[0]))
        EX &= os.path.exists(os.path.join(self.model_dir,
                                          self.upBamFile[1]))
        assert(EX)

#def test_parse(self):

