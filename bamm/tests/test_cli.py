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

__author__ = "Donovan Parks"
__copyright__ = "Copyright 2014"
__credits__ = ["Donovan Parks"]
__license__ = "LGPLv3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"

###############################################################################

# system imports
import os
import random

from nose.tools import assert_equals, assert_true

# local imports
from bamm.bamMaker import BamScheduler
from bamm.bamParser import BamParser
from bamm.bammExceptions import InvalidParameterSetException

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class TestCommandLineInterface:
    @classmethod
    def setup_class(self):
        """Setup class variables before any tests."""
        self.scriptDir = os.path.dirname(os.path.realpath(__file__))
        self.testDataDir = os.path.join(self.scriptDir, 'cli_test_data')
        self.database = os.path.join(self.testDataDir, 'test_database.fna')
        self.reads_1 = os.path.join(self.testDataDir, 'test_reads.1.fna')
        self.reads_2 = os.path.join(self.testDataDir, 'test_reads.2.fna')

    #########################################################################
    # make interface
    #########################################################################
    def testIncompletePairess(self):
        """Test odd number of paired files."""
        try:
            BamScheduler(self.database,
                            'mem',
                            None,
                            self.testDataDir,
                            paired=[self.reads_1],
                            interleaved=[],
                            singleEnded=[],
                            keptFiles=False,
                            keepFiles=False,
                            outputTam=False,
                            prefix='',
                            numThreads=1,
                            maxMemory=None,
                            forceOverwriting=False,
                            showCommands=False)
        except InvalidParameterSetException:
            return True

        assert_true(False, "Make interface failed to report an uneven number of paired read files.")

    def testInvalidDatabase(self):
        """Test invalid database."""
        try:
            BamScheduler('<invalid>',
                            'mem',
                            None,
                            self.testDataDir,
                            paired=[self.reads_1, self.reads_2],
                            interleaved=[],
                            singleEnded=[],
                            keptFiles=False,
                            keepFiles=False,
                            outputTam=False,
                            prefix='',
                            numThreads=1,
                            maxMemory=None,
                            forceOverwriting=False,
                            showCommands=False)
        except InvalidParameterSetException:
            return True

        assert_true(False, "Make interface failed to report invalid database file.")

    def testNoReads(self):
        """Test input with no read files."""
        try:
            BamScheduler(self.database,
                            'mem',
                            None,
                            self.testDataDir,
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
                            showCommands=False)
        except InvalidParameterSetException:
            return True

        assert_true(False, "Make interface failed to report invalid database file.")

    def testOutputTam(self):
        """Test creation of TAM file with output prefix."""
        bs = BamScheduler(self.database,
                          'bwasw',
                          None,
                          self.testDataDir,
                          paired=[self.reads_1, self.reads_2],
                          interleaved=[],
                          singleEnded=[],
                          keptFiles=True,
                          keepFiles=False,
                          outputTam=True,
                          prefix='test_results.',
                           numThreads=1,
                          maxMemory=None,
                          forceOverwriting=False,
                          showCommands=False,
                          quiet=True)

        bs.makeBams()
        outputFile = os.path.join(self.testDataDir, 'test_results.test_database.test_reads.1.tam')
        assert_true(os.path.exists(outputFile), "Make interface failed to produce expected TAM file.")
        os.remove(outputFile)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

