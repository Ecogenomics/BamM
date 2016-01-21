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

__author__ = "Tim Lamberton"
__copyright__ = "Copyright 2015"
__credits__ = ["Tim Lamberton"]
__license__ = "LGPLv3"
__maintainer__ = "Tim Lamberton"
__email__ = "tim.lamberton@gmail.com"

###############################################################################

# system imports
from nose.tools import assert_equals, assert_true
import sys
import os
import subprocess
import pysam

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class TestBamFilter:
    @classmethod
    def setup_class(self):
        """Setup class variables before any tests."""

        #self.bamm = os.path.join("~", "git", "BamM", "bin", "bamm")
        self.dataDir = os.path.join(os.path.split(__file__)[0], "filter_test_data")
        self.bamNames = ["1", "2"]

        # the following files already exist
        self.bamFiles = dict(zip(self.bamNames,
                                 [os.path.join(self.dataDir, "%s.bam" % name) for name in self.bamNames]))
        self.testDataDirs = dict(zip(self.bamNames,
                                     [os.path.join(self.dataDir, name) for name in self.bamNames]))


        # generated files
        self.outputBamFnames = dict(zip(self.bamNames,
                                           ["%s_filtered.bam" % name for name in self.bamNames]))


        # if True tests should fail
        if False:
            self.bamFiles = dict(zip(self.bamNames,
                                     [os.path.join(self.dataDir, "f.bam") for _ in self.bamNames]))
            self.outputBamFnames = dict(zip(self.bamNames,
                                           ["f_filtered.bam" for _ in self.bamNames]))


        # test parameters
        self.params = {
          "none": ['--use_secondary',
                   '--use_supplementary',
                   '--percentage_aln', "0",
                   '--percentage_id', "0"],
          "aln_only_90": ['--use_secondary',
                          '--use_supplementary',
                          '--percentage_aln', "0.9",
                          '--percentage_id', "0"],
          "aln_only_101": ['--use_secondary',
                           '--use_supplementary',
                           '--percentage_aln', "1.01",
                           '--percentage_id', "0"],
          "id_only_90": ['--use_secondary',
                         '--use_supplementary',
                         '--percentage_aln', "0",
                         '--percentage_id', "0.9"],
          "id_only_101": ['--use_secondary',
                          '--use_supplementary',
                          '--percentage_aln', "0",
                          '--percentage_id', "1.01"],
          "no_secondary_only": ['--use_supplementary',
                                '--percentage_aln', "0",
                                '--percentage_id', "0"],
          "no_supp_only": ['--use_secondary',
                           '--percentage_aln', "0",
                           '--percentage_id', "0"],
          "all_conds": ['--percentage_aln', "0.9",
                        '--percentage_id', "0.9"]
        }

    @classmethod
    def teardown_class(self):
        for name in self.bamNames:
            self.rmTestFile(name)

    @classmethod
    def rmTestFile(self, name):
        path = os.path.join(self.dataDir, self.outputBamFnames[name])
        if os.path.exists(path):
            os.remove(path)
        else:
            sys.stderr.write("No file: %s\n" % path)


    def generate_bam(self, name, args):
        cmd = " ".join(["bamm filter -b",
                        self.bamFiles[name],
                        "-o",
                        self.dataDir,
                        " ".join(args)])
        subprocess.check_call(cmd, shell=True)


    def assert_equal_query_sequences(self, out, expected):
        try:
            aln_expected = pysam.AlignmentFile(expected, "rb")
        except:
            raise
            raise AssertionError('File of expected reads "%s" exists and is readable.' % expected)

        try:
            aln_out = pysam.AlignmentFile(out, "rb")
        except:
            raise AssertionError('File of filtered reads "%s" exists and is readable.' % out)

        while True:
            try:
                expected_read = aln_expected.next()
            except StopIteration:
                expected_read = None

            try:
                out_read = aln_out.next()
            except StopIteration:
                out_read = None

            if expected_read is None:
                assert_true(out_read is None, 'Filtered file "%s" contains expected number of reads.' %out)
                break

            assert_true(expected_read.compare(out_read) == 0, 'Filtered file "%s" queries match expected queries.' % out)

    def assert_equal_together_query_sequences(self, (out_a, out_b), expected):
        try:
            aln_expected = pysam.AlignmentFile(expected, "rb")
        except:
            raise
            raise AssertionError('File of expected reads "%s" exists and is readable.' % expected)

        try:
            aln_out_a = pysam.AlignmentFile(out_a, "rb")
        except:
            raise AssertionError('File of filtered reads "%s" exists and is readable.' % out_a)

        try:
            aln_out_b = pysam.AlignmentFile(out_b, "rb")
        except:
            raise AssertionError('File of filtered reads "%s" exists and is readable.' % out_b)

        current_is_a = True
        out_read_b = aln_out_b.next()
        while True:
            try:
                expected_read = aln_expected.next()
            except StopIteration:
                expected_read = None

            match = False
            if current_is_a:
                try:
                    out_read_a = aln_out_a.next()
                except StopIteration:
                    out_read_a = None
                    
                match = out_read_a is not None and expected_read.compare(out_read_a)
                if not match:
                    current_is_a = False
                    match = out_read_b is not None and expected_read.compare(out_read_b)
            else:
                try:
                    out_read_b = aln_out_b.next()
                except StopIteration:
                    out_read_b = None
                    
                match = out_read_b is not None and expected_read.compare(out_read_b)
                if not match:
                    current_is_a = True
                    match = out_read_a is not None and expected_read.compare(out_read_a)
            

            if expected_read is None:
                assert(out_read_a is None and out_read_b is None, 'Filtered file "%s" contains expected number of reads.' %out)
                break

            assert_true(match, 'Filtered file "%s" queries match expected queries.' % out)


    def testFilter(self):
        for bamName in self.bamNames:

            for (testName, args) in self.params.iteritems():
                self.generate_bam(bamName, args)
                out = os.path.join(self.dataDir, self.outputBamFnames[bamName])
                test = os.path.join(self.testDataDirs[bamName], "%s_%s.bam" % (bamName, testName))
                self.assert_equal_query_sequences(out, test)
                #self.rmTestFile(bamName)
                
    def testInverseFilter(self):
        for bamName in self.bamNames:
            
            for (testName, args) in self.params.iteritems():
                args = args + ["-v"]
                self.generate_bam(bamName, args)
                out = os.path.join(self.dataDir, self.outputBamFnames[bamName])
                test = os.path.join(self.testDataDirs[bamName], "%s_%s.bam" % (bamName, testName))
                orig = os.path.join(self.dataDir, self.bamFiles[bamName])
                self.assert_equal_together_query_sequences((out, test), orig)


###############################################################################
###############################################################################
###############################################################################
###############################################################################
