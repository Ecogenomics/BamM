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

void estimateCoverages(float * coverageValues,
                       uint32_t ** pileupValues,
                       BM_coverageType * covType,
                       uint32_t contigLength,
                       uint32_t numBams
) {
    switch(covType->type) {
        case CT_COUNT:
            // pilupValues are interpreted as readStarts
            estimate_COUNT_Coverage(coverageValues,
                                    pileupValues,
                                    contigLength,
                                    numBams);
            break;
        case CT_C_MEAN:
            // pilupValues are interpreted as readStarts
            estimate_C_MEAN_Coverage(coverageValues,
                                     pileupValues,
                                     contigLength,
                                     numBams);
            break;
        case CT_P_MEAN:
            estimate_P_MEAN_Coverage(coverageValues,
                                     pileupValues,
                                     contigLength,
                                     numBams);
            break;
        case CT_P_MEDIAN:
            estimate_P_MEDIAN_Coverage(coverageValues,
                                       pileupValues,
                                       contigLength,
                                       numBams);
            break;
        case CT_P_MEAN_TRIMMED:
            estimate_P_MEAN_TRIMMED_Coverage(coverageValues,
                                             pileupValues,
                                             covType->range,
                                             contigLength,
                                             numBams);
            break;
        case CT_P_MEAN_OUTLIER:
            estimate_P_MEAN_OUTLIER_Coverage(coverageValues,
                                             pileupValues,
                                             covType->range,
                                             contigLength,
                                             numBams);
            break;
        case CT_NONE:
        default:
            break;
    }
}

//------------------------------------------------------------------------------
//
void estimate_COUNT_Coverage(float * coverageValues,
                             uint32_t ** readStarts,
                             uint32_t contigLength,
                             uint32_t numBams
) {
    int pos = 0, b = 0;
    for(; b < numBams; ++b) {
        uint32_t rc_sum = 0;
        for(pos = 0; pos < contigLength; ++pos) {
            rc_sum += readStarts[b][pos];
        }
        coverageValues[b] = (float)(rc_sum);
    }
}

//------------------------------------------------------------------------------
//
void estimate_C_MEAN_Coverage(float * coverageValues,
                              uint32_t ** readStarts,
                              uint32_t contigLength,
                              uint32_t numBams
) {
    int b = 0;
    for(; b < numBams; ++b) {
        coverageValues[b] = BM_mean(readStarts[b], contigLength);
    }
}

//------------------------------------------------------------------------------
//
void estimate_P_MEAN_Coverage(float * coverageValues,
                              uint32_t ** pileupValues,
                              uint32_t contigLength,
                              uint32_t numBams
) {
    int b = 0;
    for(; b < numBams; ++b) {
        coverageValues[b] = BM_mean(pileupValues[b], contigLength);
    }
}

//------------------------------------------------------------------------------
//
void estimate_P_MEDIAN_Coverage(float * coverageValues,
                                uint32_t ** pileupValues,
                                uint32_t contigLength,
                                uint32_t numBams
) {
    int b = 0;
    for(; b < numBams; ++b) {
        coverageValues[b] = BM_median(pileupValues[b], contigLength);
    }
}

//------------------------------------------------------------------------------
//
void estimate_P_MEAN_TRIMMED_Coverage(float * coverageValues,
                                      uint32_t ** pileupValues,
                                      float percent,
                                      uint32_t contigLength,
                                      uint32_t numBams
) {}

//------------------------------------------------------------------------------
//
void estimate_P_MEAN_OUTLIER_Coverage(float * coverageValues,
                                      uint32_t ** pileupValues,
                                      float stdevs,
                                      uint32_t contigLength,
                                      uint32_t numBams
) {
    int pos = 0, b = 0;
    for(; b < numBams; ++b) {
        uint32_t plp_sum = 0;
        uint32_t drops = 0;
        // set the cut off at a stdev either side of the mean
        float m = BM_mean(pileupValues[b], contigLength);
        float std = BM_stdDev(pileupValues[b],
                              contigLength,
                              m);
        float lower_cut = ((m-(stdevs*std)) < 0) ? 0 : (m-(stdevs*std));
        float upper_cut = m+(stdevs*std);
        for(pos = 0; pos < contigLength; ++pos) {
            if((pileupValues[b][pos] <= upper_cut) &&
                (pileupValues[b][pos] >= lower_cut)) {
                // OK
                plp_sum += pileupValues[b][pos];
            } else {
                // DROP
                ++drops;
            }
        }
        coverageValues[b] = (float)(plp_sum) / (float)(contigLength - drops);
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
