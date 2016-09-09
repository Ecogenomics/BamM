//#############################################################################
//
//   coverageEstimators.c
//
//   Methods for estimating coverage values from pileups
//
//   Copyright (C) Michael Imelfort, Ben Woodcroft
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

// system includes
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <libgen.h>

// local includes
#include "coverageEstimators.h"
#include "stats.h"

// need this huh?
#ifndef NAN
    #define NAN (0.0/0.0)
#endif

//****************** partial quicksort includes and defines *********/
#include "pqsort.h"
static inline int Uint32_tCmp(const uint32_t a, const uint32_t b) {
        return (a - b);
}
define_pqsort(my_uint32_t, uint32_t, Uint32_tCmp);
//*******************************************************************/

typedef float (*CovFunc)(uint32_t * data,
                         BM_coverageType * covType,
                         uint32_t contigLength);

void estimateCoverages(float * coverageValues,
                       uint32_t ** data,
                       BM_coverageType * covType,
                       uint32_t contigLength,
                       uint32_t numBams
) {
    CovFunc cf = 0;
    switch(covType->type) {
        case CT_COUNT:
            cf = estimate_COUNT_Coverage;
            break;
        case CT_C_MEAN:
            cf = estimate_C_MEAN_Coverage;
            break;
        case CT_P_MEAN:
            cf = estimate_P_MEAN_Coverage;
            break;
        case CT_P_MEDIAN:
            cf = estimate_P_MEDIAN_Coverage;
            break;
        case CT_P_MEAN_TRIMMED:
            cf = estimate_P_MEAN_TRIMMED_Coverage;
            break;
        case CT_P_MEAN_OUTLIER:
            cf = estimate_P_MEAN_OUTLIER_Coverage;
            break;
        case CT_P_VARIANCE:
            cf = estimate_P_VARIANCE_Coverage;
            break;
        case CT_MAPPED_COUNT:
            cf = estimate_MAPPED_COUNT_Coverage;
            break;
        case CT_MAPPED_MEAN:
            cf = estimate_MAPPED_MEAN_Coverage;
            break;
        case CT_MAPPED_MEAN_TRIMMED:
            cf = estimate_MAPPED_MEAN_TRIMMED_Coverage;
            break;
        case CT_NONE:
        default:
            return;
    }
    int b = 0;
    for(; b < numBams; ++b) {
        coverageValues[b] = (*cf)(data[b], covType, contigLength);
    }
}

//------------------------------------------------------------------------------
//
float estimate_COUNT_Coverage(uint32_t * readStarts,
                              BM_coverageType * covType,
                              uint32_t contigLength
) {
    int pos = 0;
    uint32_t rc_sum = 0;
    for(pos = 0; pos < contigLength; ++pos) {
        rc_sum += readStarts[pos];
    }
    return (float)(rc_sum);
}

//------------------------------------------------------------------------------
//
float estimate_C_MEAN_Coverage(uint32_t * readStarts,
                               BM_coverageType * covType,
                               uint32_t contigLength
) {
    return BM_mean(readStarts, contigLength);
}

//------------------------------------------------------------------------------
//
float estimate_P_MEAN_Coverage(uint32_t * pileupValues,
                               BM_coverageType * covType,
                               uint32_t contigLength
) {
    return BM_mean(pileupValues, contigLength);
}

//------------------------------------------------------------------------------
//
float estimate_P_MEDIAN_Coverage(uint32_t * pileupValues,
                                 BM_coverageType * covType,
                                 uint32_t contigLength
) {
    return BM_median(pileupValues, contigLength);
}

//------------------------------------------------------------------------------
//
float estimate_P_MEAN_TRIMMED_Coverage(uint32_t * pileupValues,
                                       BM_coverageType * covType,
                                       uint32_t contigLength
) {
    // Convert the top and bottom percentages to an absolute number
    // of pileups to remove and keep, respectively.
    uint32_t numToRemoveOffBottom = (uint32_t) (covType->lowerCut/100*contigLength);
    uint32_t numToRemoveOffTop = (uint32_t) (covType->upperCut/100*contigLength);

    uint32_t divisor = contigLength-numToRemoveOffBottom-numToRemoveOffTop;

    int pos = 0;
    uint32_t plp_sum;

    // to avoid dividing by zero or a negative number
    if (numToRemoveOffBottom+numToRemoveOffTop >= contigLength) {
        return NAN;
    } else {
        // usual sensible length contig
        // pqsort the coverage array around the upper and lower limits
        my_uint32_t_pqsort(pileupValues,
                           contigLength,
                           numToRemoveOffBottom,
                           contigLength-numToRemoveOffTop);

        // calculate the mean of those values between the limits
        plp_sum = 0;
        for(pos = numToRemoveOffBottom; \
            pos < contigLength-numToRemoveOffTop; \
            ++pos) {
            plp_sum += pileupValues[pos];
        }
        return (float)(plp_sum) / divisor;
    }
}

//------------------------------------------------------------------------------
//
float estimate_P_MEAN_OUTLIER_Coverage(uint32_t * pileupValues,
                                       BM_coverageType * covType,
                                       uint32_t contigLength
) {
    if (contigLength == 0)
        return NAN;

    int pos = 0;
    uint32_t plp_sum = 0;
    uint32_t drops = 0;
    // set the cut off at a stdev either side of the mean
    float m = BM_mean(pileupValues, contigLength);
    float std = BM_stdDev(pileupValues,
                          contigLength,
                          m);
    float lower_cut = ((m-(covType->lowerCut*std)) < 0) ? \
                        0 : (m-(covType->lowerCut*std));
    float upper_cut = m+(covType->upperCut*std);

    for(pos = 0; pos < contigLength; ++pos) {
        if((pileupValues[pos] <= upper_cut) &&
            (pileupValues[pos] >= lower_cut)) {
            // OK
            plp_sum += pileupValues[pos];
        } else {
            // DROP
            ++drops;
        }
    }

    float divisor = (float)(contigLength - drops);
    if (divisor == 0) {
        return 0.0;
    }
    else {
        return (float)(plp_sum) / divisor;
    }
}

//------------------------------------------------------------------------------
//
float estimate_P_VARIANCE_Coverage(uint32_t * pileupValues,
                                   BM_coverageType * covType,
                                   uint32_t contigLength
                                  ) {
    return BM_variance(pileupValues, contigLength, -1);
}

//------------------------------------------------------------------------------
//
float estimate_MAPPED_COUNT_Coverage(uint32_t * pileupValues,
                                     BM_coverageType * covType,
                                     uint32_t contigLength
                                    ) {
    int pos = 0;
    uint32_t rc_sum = 0;
    for(pos = 0; pos < contigLength; ++pos) {
        rc_sum += pileupValues[pos] > 0;
    }
    return (float)(rc_sum);
}

//------------------------------------------------------------------------------
//
float estimate_MAPPED_MEAN_Coverage(uint32_t * pileupValues,
                                    BM_coverageType * covType,
                                    uint32_t contigLength
                                   ) {
    return BM_nzmean(pileupValues, contigLength);
}

//------------------------------------------------------------------------------
//
float estimate_MAPPED_MEAN_TRIMMED_Coverage(uint32_t * pileupValues,
                                            BM_coverageType * covType,
                                            uint32_t contigLength
                                           ) {
    // Convert the top and bottom percentages to an absolute number
    // of pileups to remove and keep, respectively.
    uint32_t numToRemoveOffBottom = (uint32_t) (covType->lowerCut/100*contigLength);
    uint32_t numToRemoveOffTop = (uint32_t) (covType->upperCut/100*contigLength);

    uint32_t divisor = contigLength-numToRemoveOffBottom-numToRemoveOffTop;

    int pos = 0;
    uint32_t plp_sum;

    // to avoid dividing by zero or a negative number
    if (numToRemoveOffBottom+numToRemoveOffTop >= contigLength) {
        return NAN;
    } else {
        // usual sensible length contig
        // pqsort the coverage array around the upper and lower limits
        my_uint32_t_pqsort(pileupValues,
                           contigLength,
                           numToRemoveOffBottom,
                           contigLength-numToRemoveOffTop);

        // calculate the mean of those values between the limits
        plp_sum = 0;
        for(pos = numToRemoveOffBottom; \
            pos < contigLength-numToRemoveOffTop; \
            ++pos) {
            plp_sum += pileupValues[pos] > 0;
        }
        return (float)(plp_sum) / divisor;
    }
}
