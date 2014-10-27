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


import unittest

from bamm.cWrapper import *

from numpy import mean, median

class VerifyCoverageStats(unittest.TestCase):
    @classmethod
    def setup_class(self):
        """Setup class variables before any tests."""
        
        self.even_list = range(10)
        self.even_c_array = (c.c_uint32 * len(self.even_list))()
        self.even_c_array[:] = self.even_list
        
        self.odd_list = range(9)
        self.odd_c_array = (c.c_uint32 * len(self.odd_list))()
        self.odd_c_array[:] = self.odd_list
        
        self.empty_list = []
        self.empty_c_array = (c.c_uint32 * len(self.empty_list))()
        self.empty_c_array[:] = self.empty_list
        
        self.single_num_list = [10]
        self.single_num_c_array = (c.c_uint32 * len(self.single_num_list))()
        self.single_num_c_array[:] = self.single_num_list
        
        self.zero_list = [0]*10
        self.zero_c_array = (c.c_uint32 * len(self.zero_list))()
        self.zero_c_array[:] = self.zero_list
        
        BCT = BM_coverageType_C()
        BCT.type = CT.C_MEAN
        BCT.upperCut = 1
        BCT.lowerCut = 1
        
        self.pBCT = c.POINTER(BM_coverageType_C)
        self.pBCT = c.pointer(BCT)
        
        self.CW = CWrapper(unitTests=True)

        
    def testMedianEven(self):
        """Verify C implementation of median with even number of values."""

        rtn = self.CW._BM_median(self.even_c_array, len(self.even_list))
        self.assertEqual(rtn, median(self.even_list))
        
    def testMedianOdd(self):
        """Verify C implementation of median with odd number of values."""

        rtn = self.CW._BM_median(self.odd_c_array, len(self.odd_list))
        self.assertEqual(rtn, median(self.odd_list))
        
    def testMean(self):
        """Verify C implementation of mean."""

        rtn = self.CW._BM_mean(self.even_c_array, len(self.even_list))
        self.assertEqual(rtn, mean(self.even_list))
        
    def testMeanEmpty(self):
        """Verify C implementation of mean with empty array."""

        rtn = self.CW._BM_mean(self.empty_c_array, len(self.empty_list))
        self.assertEqual(rtn, c.c_float(float('nan')))
        
    def testMeanSingle(self):
        """Verify C implementation of mean when given array with a single element."""

        rtn = self.CW._BM_mean(self.single_num_c_array, len(self.single_num_list))  
        self.assertEqual(rtn, mean(self.single_num_list))
        
    def testMeanZero(self):
        """Verify C implementation of mean when given array of zeros."""

        rtn = self.CW._BM_mean(self.zero_c_array, len(self.zero_list))  
        self.assertEqual(rtn, mean(self.zero_list))
        
    def testCountCoverage(self):
        """Verify computation of count coverage."""

        coverage = self.CW._estimate_COUNT_Coverage(self.even_c_array, self.pBCT, len(self.even_list))
        self.assertEqual(coverage, sum(self.even_list))
        
        coverage = self.CW._estimate_COUNT_Coverage(self.odd_c_array, self.pBCT, len(self.odd_list))
        self.assertEqual(coverage, sum(self.odd_list))
        
    def testCountCoverageEmpty(self):
        """Verify computation of count coverage when given an empty array."""

        coverage = self.CW._estimate_COUNT_Coverage(self.empty_c_array, self.pBCT, len(self.empty_list))
        self.assertEqual(coverage, sum(self.empty_list))
        
    def testCountCoverageSingle(self):
        """Verify computation of count coverage when given an array with a single element."""

        coverage = self.CW._estimate_COUNT_Coverage(self.single_num_c_array, self.pBCT, len(self.single_num_list))
        self.assertEqual(coverage, sum(self.single_num_list))
        
    def testCountCoverageZero(self):
        """Verify computation of count coverage when given array of zeros."""

        coverage = self.CW._estimate_COUNT_Coverage(self.zero_c_array, self.pBCT, len(self.zero_list))
        self.assertEqual(coverage, sum(self.zero_list))
        
    def testCountMeanCoverage(self):
        """Verify computation of count mean coverage."""

        coverage = self.CW._estimate_C_MEAN_Coverage(self.even_c_array, self.pBCT, len(self.even_list))
        self.assertEqual(coverage, mean(self.even_list))
        
        coverage = self.CW._estimate_C_MEAN_Coverage(self.odd_c_array, self.pBCT, len(self.odd_list))
        self.assertEqual(coverage, mean(self.odd_list))
        
    def testCountMeanCoverageEmpty(self):
        """Verify computation of count coverage when given an empty array."""

        coverage = self.CW._estimate_C_MEAN_Coverage(self.empty_c_array, self.pBCT, len(self.empty_list))
        self.assertEqual(coverage, mean(self.empty_list))
        
    def testCountMeanCoverageSingle(self):
        """Verify computation of count coverage when given an array with a single element."""

        coverage = self.CW._estimate_C_MEAN_Coverage(self.single_num_c_array, self.pBCT, len(self.single_num_list))
        self.assertEqual(coverage, mean(self.single_num_list))
        
    def testCountMeanCoverageZero(self):
        """Verify computation of count coverage when given array of zeros."""

        coverage = self.CW._estimate_C_MEAN_Coverage(self.zero_c_array, self.pBCT, len(self.zero_list))
        self.assertEqual(coverage, mean(self.zero_list))

if __name__ == "__main__":
    unittest.main()
