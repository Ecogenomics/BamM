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

try:
    import pysam
except ImportError:
    print """ERROR: Some tests for `bamm filter` requires that pysam be installed.
    See 'http://pysam.readthedocs.io/en/latest/installation.html#installation' for 
    installation details."""
    raise

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class TestBamProfiler:
    @classmethod
    def setup_class(self):
        """Setup class variables before any tests."""

        #self.bamm = os.path.join("~", "git", "BamM", "bin", "bamm")
        self.dataDir = os.path.join(os.path.split(__file__)[0], "filter_test_data")
        self.bamNames = ["1", "2"]

        # the following files already exist
        self.bamFiles = dict(zip(self.bamNames,
                                 [os.path.join(self.dataDir, "%s.bam" % name) for name in self.bamNames]))

        # generated files
        self.outputProfiles = dict(zip(self.bamNames,
                                        ["%s_profile.txt" % name for name in self.bamNames]))

        # if True tests should fail
        if False:
            self.bamFiles = dict(zip(self.bamNames,
                                     [os.path.join(self.dataDir, "f.bam") for _ in self.bamNames]))
            self.outputProfiles = dict(zip(self.bamNames,
                                           ["f_profile.txt" for _ in self.bamNames]))


        # test parameters
        self.params = {
          "none": ['--use_secondary',
                   '--use_supplementary']
        }

    @classmethod
    def teardown_class(self):
        for name in self.bamNames:
            self.rmTestFile(name)

    @classmethod
    def rmTestFile(self, name):
        paths = [os.path.join(self.dataDir, self.outputProfiles[name])]
        print(paths);
        for path in paths:
            if os.path.exists(path):
                os.remove(path)
            else:
                sys.stderr.write("No file: %s\n" % path)


    def generate_profile(self, name, args):
        cmd = " ".join(["bamm profile -b",
                        self.bamFiles[name],
                        " ".join(args),
                        " > ",
                        os.path.join(self.dataDir, self.outputProfiles[name])])
        subprocess.check_call(cmd, shell=True)
        

    def assert_profile(self, out, bam):
        try:
            aln = pysam.AlignmentFile(bam, "rb")
        except:
            raise AssertionError('File of reads "%s" exists and is readable.' % expected)

        try:
            prf_out = open(out, "rb")
        except:
            raise AssertionError('Generated profile "%s" exists and is readable.' % out)

        count = 0;
        
        try:
            out_header = prf_out.next().split()
        except StopIteration:
            raise AssertionError('Generated profile "%s" is not an empty file.' % out) 

        assert_equals(out_header,
                      "line\tsupp\tsecondary\tmapQ\tmismatches\tmatches\tqLen\tpcId\tpcAln".split())

        while True:
            try:
                read = aln.next()
                count += 1
            except StopIteration:
                read = None

            try:
                out_vals = prf_out.next().split()
            except StopIteration:
                out_vals = None

            if read is None:
                assert_true(out_vals is None, 'Profile of "%s" contains expected number of reads.' %out)
                break
                
            assert_true(out_vals is not None, 'Profile of "%s" contains expected number of reads.' %out)

            assert_equals(len(out_vals), 9, 'Profile "%s" contains 9 fields' % out)
            assert_equals(out_vals[0], '%d'%count, '"Line" field of %d read from "%s" matches expected value.' % (count, out))
            assert_equals(out_vals[1], '%d'%read.is_supplementary, '"Supplementary" field of %d read of "%s" matches expected value.' % (count,out))
            assert_equals(out_vals[2], '%d'%read.is_secondary, '"Secondary" field of %d read of "%s" matches expected value.' % (count,out))
            assert_equals(out_vals[3], '%d'%read.mapping_quality, '"mapQ" field of %d read of "%s" matches expected value.' % (count,out))
            
            # mismatches (NM tag)
            cigar = read.get_cigar_stats()[0]
            mismatches = cigar[10]
            assert_equals(out_vals[4], '%d'%mismatches, '"mismatches" field of %d read of "%s" matches expected value.' % (count,out))
            
            # matches (CMATCH+CEQUAL+CDIFF)
            matches = cigar[0]+cigar[7]+cigar[8]
            assert_equals(out_vals[5], '%d'%matches, '"matches" field of %d read of "%s" matches expected value.' % (count,out))
            
            # qLen
            qLen = read.infer_query_length()
            assert_equals(out_vals[6], '%d'%qLen, '"qLen" field of %d read of "%s" matches expected value.' % (count,out))
            
            # pcId
            pcId = float(matches - mismatches) / matches;
            assert_almost_equals(float(out_vals[7]), pcId, '"pcId" field of %d read of "%s" matches expected value.' % (count,out))
            
            # pcAln
            pcAln = float(matches) / qLen;
            assert_almost_equals(float(out_vals[8]), pcAln, '"pcAln" field of %d read of "%s" matches expected value.' % (count,out))

    
            

    def testProfiler(self):
        for bamName in self.bamNames:
            for (testName, args) in self.params.iteritems():
                self.generate_profile(bamName, args)
                bam = self.bamFiles[bamName]
                out = os.path.join(self.dataDir, self.outputProfiles[bamName])
                self.assert_profile(out, bam)



def assert_almost_equals(a, b, msg=None):
    if msg==None:
        msg = "%s != %s" % (a, b)
    return assert_true(abs(a - b) < 0.001, msg)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
