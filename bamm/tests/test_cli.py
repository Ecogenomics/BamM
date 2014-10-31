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

import os
import random

from nose.tools import assert_equals, assert_true

from bamm.bamMaker import BamScheduler
from bamm.bammExceptions import *

class TestCommandLineInterface:
    @classmethod
    def setup_class(self):
        """Setup class variables before any tests."""
        self.scriptDir = os.path.dirname(os.path.realpath(__file__))
        self.testDataDir = os.path.join(self.scriptDir, '..', '..', 'test', 'data')
        self.database = os.path.join(self.testDataDir, 'test_database.fna')
        self.reads_1 = self.database = os.path.join(self.testDataDir, 'test_reads.1.fna')
        self.reads_2 = self.database = os.path.join(self.testDataDir, 'test_reads.2.fna')

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
                            verbose=False)
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
                            verbose=False)
        except InvalidParameterSetException:
            return True
            
        assert_true(False, "Make interface failed to report invalid database file.")
            
    #########################################################################
    # parse interface
    #########################################################################
   
    #########################################################################
    # extract interface
    #########################################################################
