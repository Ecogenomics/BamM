//#############################################################################
//
//   coverageEstimators.h
//
//   Methods for estimating coverage values from pileups
//
//   Copyright (C) Michael Imelfort
//
//   This library is free software; you can redistribute it and/or
//   modify it under the terms of the GNU Lesser General Public
//   License as published by the Free Software Foundation; either
//   version 3.0 of the License, or (at your option) any later version.
//
//   This library is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//   Lesser General Public License for more details.
//
//   You should have received a copy of the GNU Lesser General Public
//   License along with this library.
//
//#############################################################################

#ifndef BM_COV_ESTIMATOR
  #define BM_COV_ESTIMATOR

// system includes
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// there are several different definitions of coverage types
typedef enum {CT_NONE,           // do not calculate coverage
              CT_COUNT,          // read counts, unaffected by contig length
              CT_C_MEAN,         // read counts, divided by contig length
              CT_P_MEAN,         // mean pileup depth
              CT_P_MEDIAN,       // median pileup depth
              CT_P_MEAN_TRIMMED, // pileup mean trancated based on upper lower %
              CT_P_MEAN_OUTLIER, // pileup mean trancated based on distributions
              } CT;

/*! @typedef
 * @abstract Structure for storing information about a coverage type
 *
 * @field type                Identifies the method used to calculate coverage
 * @field range               Upper and lower bounds for coverage calculation
                              (if needed)
 */
 typedef struct BM_coverageType {
    CT type;
    float range;
 } BM_coverageType;

/*!
 * @abstract Estimate the coverage along a single contig given pileup
 *           information and a coverage type to calculate
 *
 * @param  coverageValues   array of floats to set (len == numBams)
 * @param  pileupValues     matrix of int pileup depths (numBams x contigLength)
 * @param  BM_coverageType  BM_coverageType struct with initialised values
 * @param  contigLength     Length of the contig being assesed
 * @param  numBams          the number of bams responsible for the pileup
 * @return void
 *
 * @discussion This function is essentially a wrapper that calls specific
 * coverage calculators based on the type of coverage supplied. It updates the
 * values in coverageValues before exiting.
 */
void estimateCoverages(float * coverageValues,
                       uint32_t ** pileupValues,
                       BM_coverageType * covType,
                       uint32_t contigLength,
                       uint32_t numBams
                       );


/*!
 * @abstract Estimate the raw read count coverage along a single contig
 *
 * @param  coverageValues   array of floats to set (len == numBams)
 * @param  readstarts       matrix of int read starts (numBams x contigLength)
 * @param  contigLength     Length of the contig being assesed
 * @param  numBams          the number of bams responsible for the pileup
 * @return void
 *
 * @discussion This function updates the values in coverageValues before exiting
*/
void estimate_RC_Coverage(float * coverageValues,
                          uint32_t ** readStarts,
                          uint32_t contigLength,
                          uint32_t numBams
                          );

/*!
 * @abstract Estimate the mean pileup coverage along a single contig given
 *           pileup information
 *
 * @param  coverageValues   array of floats to set (len == numBams)
 * @param  pileupValues     matrix of int pileup depths (numBams x contigLength)
 * @param  contigLength     Length of the contig being assesed
 * @param  numBams          the number of bams responsible for the pileup
 * @return void
 *
 * @discussion This function updates the values in coverageValues before exiting
*/
void estimate_PM_Coverage(float * coverageValues,
                          uint32_t ** pileupValues,
                          uint32_t contigLength,
                          uint32_t numBams
                          );


/*!
 * @abstract Estimate tuncated mean (outlier) pileup coverage along a single
 *           contig given pileup information and a symetric stdev cutoff
 *
 * @param  coverageValues   array of floats to set (len == numBams)
 * @param  pileupValues     matrix of int pileup depths (numBams x contigLength)
 * @param  stdevs           number of standard deviations to determine limits
 * @param  contigLength     Length of the contig being assesed
 * @param  numBams          the number of bams responsible for the pileup
 * @return void
 *
 * @discussion This function updates the values in coverageValues before exiting
 */
void estimate_PMO_Coverage(float * coverageValues,
                           uint32_t ** pileupValues,
                           float stdevs,
                           uint32_t contigLength,
                           uint32_t numBams
                           );

        /***********************
        ***   MATHY EXTRAS   ***
        ***********************/
/*!
 * @abstract Calculate the mean of an array values
 *
 * @param  values       array of (integer) values
 * @param  size         size of values array
 * @return float        the caluclated mean
 */
float BM_mean(uint32_t * values, uint32_t size);

/*!
 * @abstract Calculate the standard deviations of an array values
 *
 * @param  values       array of (integer) values
 * @param  size         size of values array
 * @param  m            mean of the values array
 * @return float        the caluclated standard deviation
 */
float BM_stdDev(uint32_t * values, uint32_t size, float m);

/*!
 * @abstract Calculate the standard deviations of an array values
 *
 * @param  values       array of (integer) values
 * @param  size         size of values array
 * @return float        the caluclated standard deviation
 *
 * @discussion Everything is 3 stdevs from the mean right?
*/
float BM_fakeStdDev(uint32_t * values, uint32_t size);


#ifdef __cplusplus
}
#endif

#endif // BM_COV_ESTIMATOR

