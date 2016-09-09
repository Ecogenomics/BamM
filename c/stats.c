//#############################################################################
//
//   stats.c
//
//   Curse you C!!!!!
//
//   Copyright (C) Michael Imelfort
//
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//#############################################################################

// system includes
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

// local includes
#include "stats.h"
#include "pqsort.h"

float BM_mean(uint32_t * values, uint32_t size)
{
    uint32_t sum = 0;
    int i = 0;
    for( i = 0; i < size; ++i) {
        sum += *(values + i);
    }
    return (float)sum/(float)size;
}

float BM_nzmean(uint32_t * values, uint32_t size)
{
    uint32_t sum = 0;
    int i = 0;
    for( i = 0; i < size; ++i) {
        sum += *(values + i) > 0;
    }
    return (float)sum/(float)size;
}

float BM_stdDev(uint32_t * values, uint32_t size, float m)
{
    return sqrt(BM_variance(values, size, m));
}

float BM_variance(uint32_t * values, uint32_t size, float m)
{
    float sum = 0;
    int i = 0;
    if(m == -1)
        m = BM_mean(values, size);
    for(i = 0; i < size; ++i) {
        sum += pow((float)*(values + i) - m, 2);
    }
    return sum/(float)size;
}

float BM_fakeStdDev(uint32_t * values, uint32_t size)
{
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

int cmpfunc_uint32(const void * a, const void * b)
{
    const uint32_t *x = a, *y = b;
    if(*x > *y)
        return 1;
    else
        return(*x < *y) ? -1: 0;
}

float BM_median(uint32_t * values,
                uint32_t size)
{
    if (size == 0)
        return NAN;

    qsort(values, size, sizeof(uint32_t), cmpfunc_uint32);

    if (size % 2) {
        // number is odd so return middle number
        return values[size/2];
    }

    // number is even so return mean of two middle numbers
    return (float)(values[size/2 - 1] + values[size/2]) / 2.0f;
}

/*
#define ELEM_SWAP(a,b) { register int t=(a);(a)=(b);(b)=t; }

uint32_t BM_median(uint32_t * values,
                   uint32_t size)
{
    uint32_t low, high ;
    uint32_t median;
    uint32_t middle, ll, hh;
    low = 0 ; high = size-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) // One element only
            return values[median] ;
        if (high == low + 1) { // Two elements only
            if (values[low] > values[high])
                ELEM_SWAP(values[low], values[high]) ;
            return values[median] ;
        }
        // Find median of low, middle and high items; swap into position low
        middle = (low + high) / 2;
        if (values[middle] > values[high])
            ELEM_SWAP(values[middle], values[high]) ;
        if (values[low] > values[high])
            ELEM_SWAP(values[low], values[high]) ;
        if (values[middle] > values[low])
            ELEM_SWAP(values[middle], values[low]) ;
        // Swap low item (now in position middle) into position (low+1)
        ELEM_SWAP(values[middle], values[low+1]) ;
        // Nibble from each end towards middle, swapping items when stuck
        ll = low + 1;
        hh = high;
        for (;;) {
            do ll++; while (values[low] > values[ll]) ;
            do hh--; while (values[hh] > values[low]) ;
            if (hh < ll)
                break;
            ELEM_SWAP(values[ll], values[hh]) ;
        }
        // Swap middle item (in position low) back into correct position
        ELEM_SWAP(values[low], values[hh]) ;
        // Re-set active partition
        if (hh <= median)
            low = ll;
        if (hh >= median)
            high = hh - 1;
    }
    return values[median] ;
}
*/
