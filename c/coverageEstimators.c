//#############################################################################
//
//   coverageEstimators.c
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
        case CT_NONE:
            break;
        case CT_COUNT:
            // pilupValues are interpreted as readStarts
            estimate_RC_Coverage(coverageValues,
                                 pileupValues,
                                 contigLength,
                                 numBams);
            break;
        case CT_C_MEAN:
            break;
        case CT_P_MEAN:
            estimate_PM_Coverage(coverageValues,
                                 pileupValues,
                                 contigLength,
                                 numBams);
            break;
        case CT_P_MEDIAN:
            break;
        case CT_P_MEAN_TRIMMED:
            break;
        case CT_P_MEAN_OUTLIER:
            estimate_PMO_Coverage(coverageValues,
                                  pileupValues,
                                  covType->range,
                                  contigLength,
                                  numBams);
            break;
    }
}


void estimate_RC_Coverage(float * coverageValues,
                          uint32_t ** readStarts,
                          uint32_t contigLength,
                          uint32_t numBams
) {}


void estimate_PM_Coverage(float * coverageValues,
                          uint32_t ** pileupValues,
                          uint32_t contigLength,
                          uint32_t numBams
) {
    int pos = 0, b = 0;
    for(; b < numBams; ++b) {
        uint32_t plp_sum = 0;
        for(pos = 0; pos < contigLength; ++pos) {
            plp_sum += pileupValues[b][pos];
        }
        coverageValues[b] = (float)(plp_sum) / (float)(contigLength);
    }
}

void estimate_PMO_Coverage(float * coverageValues,
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
