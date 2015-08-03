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
import random
from math import isnan
from numpy import mean, median, std, array
from nose.tools import assert_equals, assert_true
import sys
from os.path import join
import ctypes as c

# local imports
from bamm.cWrapper import BM_coverageType_C, CT, CWrapper

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class TestCoverageStats:
    @classmethod
    def setup_class(self):
        """Setup class variables before any tests."""
        self.even_list = range(10)
        self.even_c_array = (c.c_uint32 * len(self.even_list))()
        self.even_c_array[:] = self.even_list

        self.odd_list = range(9)
        self.odd_c_array = (c.c_uint32 * len(self.odd_list))()
        self.odd_c_array[:] = self.odd_list

        self.unsorted_list = range(9)
        random.shuffle(self.unsorted_list)
        self.unsorted_c_array = (c.c_uint32 * len(self.unsorted_list))()
        self.unsorted_c_array[:] = self.unsorted_list

        self.empty_list = []
        self.empty_c_array = (c.c_uint32 * len(self.empty_list))()
        self.empty_c_array[:] = self.empty_list

        self.single_num_list = [10]
        self.single_num_c_array = (c.c_uint32 * len(self.single_num_list))()
        self.single_num_c_array[:] = self.single_num_list

        self.zero_list = [0]*10
        self.zero_c_array = (c.c_uint32 * len(self.zero_list))()
        self.zero_c_array[:] = self.zero_list

        self.BCT = BM_coverageType_C()
        self.BCT.type = CT.C_MEAN
        self.BCT.upperCut = 10.0
        self.BCT.lowerCut = 10.0

        self.pBCT = c.POINTER(BM_coverageType_C)
        self.pBCT = c.pointer(self.BCT)

        self.CW = CWrapper(UT=True)

    #########################################################################
    # Median
    #########################################################################
    def testMedianEven(self):
        """Verify C implementation of median with even number of values."""

        rtn = self.CW._BM_median(self.even_c_array, len(self.even_list))
        assert_equals(rtn, median(self.even_list))

    def testMedianOdd(self):
        """Verify C implementation of median with odd number of values."""

        rtn = self.CW._BM_median(self.odd_c_array, len(self.odd_list))
        assert_equals(rtn, median(self.odd_list))

    def testMedianEmpty(self):
        """Verify C implementation of median with empty array."""

        pass
        #rtn = self.CW._BM_median(self.empty_c_array, len(self.empty_list))
        #assert_true(isnan(rtn))

    def testMedianSingle(self):
        """Verify C implementation of median when given array with a single element."""

        rtn = self.CW._BM_median(self.single_num_c_array, len(self.single_num_list))
        assert_equals(rtn, median(self.single_num_list))

    def testMedianZero(self):
        """Verify C implementation of median when given array of zeros."""

        rtn = self.CW._BM_median(self.zero_c_array, len(self.zero_list))
        assert_equals(rtn, median(self.zero_list))

    #########################################################################
    # Mean
    #########################################################################
    def testMean(self):
        """Verify C implementation of mean."""

        rtn = self.CW._BM_mean(self.even_c_array, len(self.even_list))
        assert_equals(rtn, mean(self.even_list))

        rtn = self.CW._BM_mean(self.odd_c_array, len(self.odd_list))
        assert_equals(rtn, mean(self.odd_list))

    def testMeanEmpty(self):
        """Verify C implementation of mean with empty array."""

        rtn = self.CW._BM_mean(self.empty_c_array, len(self.empty_list))
        assert_true(isnan(rtn))

    def testMeanSingle(self):
        """Verify C implementation of mean when given array with a single element."""

        rtn = self.CW._BM_mean(self.single_num_c_array, len(self.single_num_list))
        assert_equals(rtn, mean(self.single_num_list))

    def testMeanZero(self):
        """Verify C implementation of mean when given array of zeros."""

        rtn = self.CW._BM_mean(self.zero_c_array, len(self.zero_list))
        assert_equals(rtn, mean(self.zero_list))

    #########################################################################
    # Count Coverage
    #########################################################################
    def testCountCoverage(self):
        """Verify computation of count coverage."""

        coverage = self.CW._estimate_COUNT_Coverage(self.even_c_array, self.pBCT, len(self.even_list))
        assert_equals(coverage, sum(self.even_list))

        coverage = self.CW._estimate_COUNT_Coverage(self.odd_c_array, self.pBCT, len(self.odd_list))
        assert_equals(coverage, sum(self.odd_list))

    def testCountCoverageEmpty(self):
        """Verify computation of count coverage when given an empty array."""

        coverage = self.CW._estimate_COUNT_Coverage(self.empty_c_array, self.pBCT, len(self.empty_list))
        assert_equals(coverage, 0)

    def testCountCoverageSingle(self):
        """Verify computation of count coverage when given an array with a single element."""

        coverage = self.CW._estimate_COUNT_Coverage(self.single_num_c_array, self.pBCT, len(self.single_num_list))
        assert_equals(coverage, sum(self.single_num_list))

    def testCountCoverageZero(self):
        """Verify computation of count coverage when given array of zeros."""

        coverage = self.CW._estimate_COUNT_Coverage(self.zero_c_array, self.pBCT, len(self.zero_list))
        assert_equals(coverage, sum(self.zero_list))

    #########################################################################
    # Count Mean Coverage
    #########################################################################
    def testCountMeanCoverage(self):
        """Verify computation of count mean coverage."""

        coverage = self.CW._estimate_C_MEAN_Coverage(self.even_c_array, self.pBCT, len(self.even_list))
        assert_equals(coverage, mean(self.even_list))

        coverage = self.CW._estimate_C_MEAN_Coverage(self.odd_c_array, self.pBCT, len(self.odd_list))
        assert_equals(coverage, mean(self.odd_list))

    def testCountMeanCoverageEmpty(self):
        """Verify computation of count mean coverage when given an empty array."""

        coverage = self.CW._estimate_C_MEAN_Coverage(self.empty_c_array, self.pBCT, len(self.empty_list))
        assert_true(isnan(coverage))

    def testCountMeanCoverageSingle(self):
        """Verify computation of count mean coverage when given an array with a single element."""

        coverage = self.CW._estimate_C_MEAN_Coverage(self.single_num_c_array, self.pBCT, len(self.single_num_list))
        assert_equals(coverage, mean(self.single_num_list))

    def testCountMeanCoverageZero(self):
        """Verify computation of count mean coverage when given array of zeros."""

        coverage = self.CW._estimate_C_MEAN_Coverage(self.zero_c_array, self.pBCT, len(self.zero_list))
        assert_equals(coverage, mean(self.zero_list))

    #########################################################################
    # Pileup Mean Coverage
    #########################################################################
    def testPileupMeanCoverage(self):
        """Verify computation of pileup mean coverage."""

        coverage = self.CW._estimate_P_MEAN_Coverage(self.even_c_array, self.pBCT, len(self.even_list))
        assert_equals(coverage, mean(self.even_list))

        coverage = self.CW._estimate_P_MEAN_Coverage(self.odd_c_array, self.pBCT, len(self.odd_list))
        assert_equals(coverage, mean(self.odd_list))

    def testPileupMeanCoverageEmpty(self):
        """Verify computation of pileup mean coverage when given an empty array."""

        coverage = self.CW._estimate_P_MEAN_Coverage(self.empty_c_array, self.pBCT, len(self.empty_list))
        assert_true(isnan(coverage))

    def testPileupMeanCoverageSingle(self):
        """Verify computation of pileup mean coverage when given an array with a single element."""

        coverage = self.CW._estimate_P_MEAN_Coverage(self.single_num_c_array, self.pBCT, len(self.single_num_list))
        assert_equals(coverage, mean(self.single_num_list))

    def testPileupMeanCoverageZero(self):
        """Verify computation of pileup mean coverage when given array of zeros."""

        coverage = self.CW._estimate_P_MEAN_Coverage(self.zero_c_array, self.pBCT, len(self.zero_list))
        assert_equals(coverage, mean(self.zero_list))

    #########################################################################
    # Pileup Median Coverage
    #########################################################################
    def testPileupMedianCoverage(self):
        """Verify computation of pileup median coverage."""

        coverage = self.CW._estimate_P_MEDIAN_Coverage(self.even_c_array, self.pBCT, len(self.even_list))
        assert_equals(coverage, median(self.even_list))

        coverage = self.CW._estimate_P_MEDIAN_Coverage(self.odd_c_array, self.pBCT, len(self.odd_list))
        assert_equals(coverage, median(self.odd_list))

    def testPileupMedianCoverageEmpty(self):
        """Verify computation of pileup median coverage when given an empty array."""

        coverage = self.CW._estimate_P_MEDIAN_Coverage(self.empty_c_array, self.pBCT, len(self.empty_list))
        assert_true(isnan(coverage))

    def testPileupMedianCoverageSingle(self):
        """Verify computation of pileup median coverage when given an array with a single element."""

        coverage = self.CW._estimate_P_MEDIAN_Coverage(self.single_num_c_array, self.pBCT, len(self.single_num_list))
        assert_equals(coverage, median(self.single_num_list))

    def testPileupMedianCoverageZero(self):
        """Verify computation of pileup median coverage when given array of zeros."""

        coverage = self.CW._estimate_P_MEDIAN_Coverage(self.zero_c_array, self.pBCT, len(self.zero_list))
        assert_equals(coverage, median(self.zero_list))

    #########################################################################
    # Pileup Outlier Coverage
    #########################################################################
    def testPileupMeanOutlierCoverage(self):
        """Verify computation of pileup outlier coverage."""

        BCT = BM_coverageType_C()
        BCT.type = CT.P_MEAN_OUTLIER
        BCT.upperCut = 1.0
        BCT.lowerCut = 1.0

        pBCT = c.POINTER(BM_coverageType_C)
        pBCT = c.pointer(BCT)

        coverage = self.CW._estimate_P_MEAN_OUTLIER_Coverage(self.even_c_array, pBCT, len(self.even_list))
        lower_bound = mean(self.even_list) - BCT.lowerCut * std(self.even_list)
        if lower_bound < 0:
            lower_bound = 0
        upper_bound = mean(self.even_list) + BCT.upperCut * std(self.even_list)
        valid_indices = [i for i, x in enumerate(self.even_list) if (x >= lower_bound and x <= upper_bound)]
        np_even_list = array(self.even_list)
        assert_equals(coverage, mean(np_even_list[valid_indices]))

        coverage = self.CW._estimate_P_MEAN_OUTLIER_Coverage(self.odd_c_array, pBCT, len(self.odd_list))
        lower_bound = mean(self.odd_list) - BCT.lowerCut * std(self.odd_list)
        if lower_bound < 0:
            lower_bound = 0
        upper_bound = mean(self.odd_list) + BCT.upperCut * std(self.odd_list)
        valid_indices = [i for i, x in enumerate(self.odd_list) if (x >= lower_bound and x <= upper_bound)]
        np_odd_list = array(self.odd_list)
        assert_equals(coverage, mean(np_odd_list[valid_indices]))

    def testPileupMeanOutlierCoverageEmpty(self):
        """Verify computation of pileup outlier coverage when given an empty array."""

        coverage = self.CW._estimate_P_MEAN_OUTLIER_Coverage(self.empty_c_array, self.pBCT, len(self.empty_list))
        assert_true(isnan(coverage))

    def testPileupMeanOutlierCoverageSingle(self):
        """Verify computation of pileup outlier coverage when given an array with a single element."""

        BCT = BM_coverageType_C()
        BCT.type = CT.P_MEAN_OUTLIER
        BCT.upperCut = 1.0
        BCT.lowerCut = 1.0

        pBCT = c.POINTER(BM_coverageType_C)
        pBCT = c.pointer(BCT)

        coverage = self.CW._estimate_P_MEAN_OUTLIER_Coverage(self.single_num_c_array, self.pBCT, len(self.single_num_list))
        assert_equals(coverage, self.single_num_list[0])

    def testPileupMeanOutlierCoverageZero(self):
        """Verify computation of pileup outlier coverage when given array of zeros."""

        coverage = self.CW._estimate_P_MEAN_OUTLIER_Coverage(self.zero_c_array, self.pBCT, len(self.zero_list))
        assert_equals(coverage, 0)

    def testPileupMeanOutlierCoverageZeroStd(self):
        """Verify computation of pileup outlier coverage when given a std of zero."""

        BCT = BM_coverageType_C()
        BCT.type = CT.P_MEAN_OUTLIER
        BCT.upperCut = 0.0
        BCT.lowerCut = 0.0

        pBCT = c.POINTER(BM_coverageType_C)
        pBCT = c.pointer(BCT)

        coverage = self.CW._estimate_P_MEAN_OUTLIER_Coverage(self.even_c_array, pBCT, len(self.even_list))
        assert_equals(coverage, 0)

    def testPileupMeanOutlierCoverageUnequalStd(self):
        """Verify computation of pileup outlier coverage when given unequal standard deviations."""

        BCT = BM_coverageType_C()
        BCT.type = CT.P_MEAN_OUTLIER
        BCT.upperCut = 1.0
        BCT.lowerCut = 2.5

        pBCT = c.POINTER(BM_coverageType_C)
        pBCT = c.pointer(BCT)

        coverage = self.CW._estimate_P_MEAN_OUTLIER_Coverage(self.even_c_array, pBCT, len(self.even_list))
        lower_bound = mean(self.even_list) - BCT.lowerCut * std(self.even_list)
        if lower_bound < 0:
            lower_bound = 0
        upper_bound = mean(self.even_list) + BCT.upperCut * std(self.even_list)
        valid_indices = [i for i, x in enumerate(self.even_list) if (x >= lower_bound and x <= upper_bound)]
        np_even_list = array(self.even_list)
        assert_equals(coverage, mean(np_even_list[valid_indices]))

    #########################################################################
    # Pileup Mean Trimmed Coverage
    #########################################################################
    def testPileupMeanTrimmedCoverage(self):
        """Verify computation of pileup mean trimmed coverage."""

        coverage = self.CW._estimate_P_MEAN_TRIMMED_Coverage(self.even_c_array, self.pBCT, len(self.even_list))
        trim_lower = int(self.BCT.lowerCut/100.0 * len(self.even_list)) + 1
        trim_upper = len(self.even_list) - int(self.BCT.upperCut/100.0 * len(self.even_list)) - 1
        assert_equals(coverage, mean(self.even_list[trim_lower:trim_upper]))

        coverage = self.CW._estimate_P_MEAN_TRIMMED_Coverage(self.odd_c_array, self.pBCT, len(self.odd_list))
        trim_lower = int(self.BCT.lowerCut/100.0 * len(self.odd_list)) + 1
        trim_upper = len(self.odd_list) - int(self.BCT.upperCut/100.0 * len(self.odd_list)) - 1
        assert_equals(coverage, mean(self.odd_list[trim_lower:trim_upper]))

    def testPileupMeanTrimmedCoverageUnsorted(self):
        """Verify computation of pileup mean trimmed coverage on an unsorted array."""

        coverage = self.CW._estimate_P_MEAN_TRIMMED_Coverage(self.unsorted_c_array, self.pBCT, len(self.unsorted_list))
        sorted_list = sorted(self.unsorted_list)
        trim_lower = int(self.BCT.lowerCut/100.0 * len(sorted_list)) + 1
        trim_upper = len(sorted_list) - int(self.BCT.upperCut/100.0 * len(sorted_list)) - 1
        assert_equals(coverage, mean(sorted_list[trim_lower:trim_upper]))

    def testPileupMeanTrimmedCoverageEmpty(self):
        """Verify computation of pileup mean trimmed coverage when given an empty array."""

        coverage = self.CW._estimate_P_MEAN_TRIMMED_Coverage(self.empty_c_array, self.pBCT, len(self.empty_list))
        assert_true(isnan(coverage))

    def testPileupMeanTrimmedCoverageSingle(self):
        """Verify computation of pileup mean trimmed coverage when given an array with a single element."""

        coverage = self.CW._estimate_P_MEAN_TRIMMED_Coverage(self.single_num_c_array, self.pBCT, len(self.single_num_list))
        assert_equals(coverage, self.single_num_list[0])

    def testPileupMeanTrimmedCoverageZero(self):
        """Verify computation of pileup mean trimmed coverage when given array of zeros."""

        coverage = self.CW._estimate_P_MEAN_TRIMMED_Coverage(self.zero_c_array, self.pBCT, len(self.zero_list))
        assert_equals(coverage, 0)

    def testPileupMeanTrimmedCoverageInvalidTrim(self):
        """Verify computation of pileup mean trimmed coverage when given an invalid trimming range."""

        BCT = BM_coverageType_C()
        BCT.type = CT.P_MEAN_TRIMMED
        BCT.upperCut = 60.0
        BCT.lowerCut = 60.0

        pBCT = c.POINTER(BM_coverageType_C)
        pBCT = c.pointer(BCT)

        coverage = self.CW._estimate_P_MEAN_TRIMMED_Coverage(self.even_c_array, pBCT, len(self.even_list))
        assert_true(isnan(coverage))

    def testPileupMeanTrimmedCoverageUnequalTrim(self):
        """Verify computation of pileup mean trimmed coverage when given an unequal trim percentages."""

        BCT = BM_coverageType_C()
        BCT.type = CT.P_MEAN_TRIMMED
        BCT.upperCut = 10.0
        BCT.lowerCut = 20.0

        pBCT = c.POINTER(BM_coverageType_C)
        pBCT = c.pointer(BCT)

        coverage = self.CW._estimate_P_MEAN_TRIMMED_Coverage(self.even_c_array, pBCT, len(self.even_list))
        trim_lower = int(BCT.lowerCut/100.0 * len(self.even_list)) + 1
        trim_upper = len(self.even_list) - int(BCT.upperCut/100.0 * len(self.even_list)) - 1
        assert_equals(coverage, mean(self.even_list[trim_lower:trim_upper]))

    #########################################################################
    # Pileup variance mode
    #########################################################################
    def testVarianceCoverage(self):
        """Verify computation of pileup variance coverage."""

        observed = self.CW._estimate_P_VARIANCE_Coverage(self.even_c_array, self.pBCT, len(self.even_list))
        expected = pow(std(self.even_list), 2)
        assert_equals(expected, observed)
 
    def testVarianceSingle(self):
        """Verify C implementation of variance when given array with a single element."""
 
        observed = self.CW._estimate_P_VARIANCE_Coverage(self.single_num_c_array, self.pBCT, len(self.single_num_list))
        expected = pow(std(self.single_num_list), 2)
        assert_equals(expected, observed)
 
    def testVarianceZero(self):
        """Verify C implementation of variance when given array of zeros."""
 
        observed = self.CW._estimate_P_VARIANCE_Coverage(self.zero_c_array, self.pBCT, len(self.zero_list))
        assert_equals(0.0, observed)
###############################################################################
###############################################################################
###############################################################################
###############################################################################
