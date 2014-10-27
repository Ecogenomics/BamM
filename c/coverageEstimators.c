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
        //usual sensible length contig
        // pqsort the coverage array around the upper and lower limits
        my_uint32_t_pqsort(pileupValues,
                           contigLength,
                           0,
                           numToRemoveOffBottom); //bottom
        my_uint32_t_pqsort(pileupValues,
                           contigLength,
                           contigLength-numToRemoveOffTop,
                           numToRemoveOffTop); //top

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

uint32_t BM_median(uint32_t * values,
                   uint32_t size)
{
    uint32_t low, high ;
    uint32_t median;
    uint32_t middle, ll, hh;
    low = 0 ; high = size-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return values[median] ;
        if (high == low + 1) { /* Two elements only */
            if (values[low] > values[high])
                ELEM_SWAP(values[low], values[high]) ;
            return values[median] ;
        }
        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if (values[middle] > values[high])
            ELEM_SWAP(values[middle], values[high]) ;
        if (values[low] > values[high])
            ELEM_SWAP(values[low], values[high]) ;
        if (values[middle] > values[low])
            ELEM_SWAP(values[middle], values[low]) ;
        /* Swap low item (now in position middle) into position (low+1) */
        ELEM_SWAP(values[middle], values[low+1]) ;
        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        for (;;) {
            do ll++; while (values[low] > values[ll]) ;
            do hh--; while (values[hh] > values[low]) ;
            if (hh < ll)
                break;
            ELEM_SWAP(values[ll], values[hh]) ;
        }
        /* Swap middle item (in position low) back into correct position */
        ELEM_SWAP(values[low], values[hh]) ;
        /* Re-set active partition */
        if (hh <= median)
            low = ll;
        if (hh >= median)
            high = hh - 1;
    }
    return values[median] ;
}

float BM_mean(uint32_t * values, uint32_t size) {
    uint32_t sum = 0;
    int i = 0;
    for( i = 0; i < size; ++i) {
        sum += *(values + i);
    }

    return (float)sum/(float)size;
}

float BM_stdDev(uint32_t * values, uint32_t size, float m) {
    float sum = 0;
    int i = 0;
    if(m == -1)
        m = BM_mean(values, size);
    for(i = 0; i < size; ++i) {
        sum += pow((float)*(values + i) - m, 2);
    }
    return sqrt(sum/(float)size);
}

float BM_fakeStdDev(uint32_t * values, uint32_t size) {
    // everything is 3 stdevs from the mean right?
    uint32_t max = 0, min = 1<<30;
    int i = 0;
    for( i = 0; i < size; ++i) {
        if (*(values + i) > max)
            max = *(values + i);
        else if (*(values + i) < min)
            min = *(values + i);
    }
    return (float)(max-min)/6;
}
