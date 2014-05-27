//#############################################################################
//
//   stats.h
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

#include <stdio.h> 
#include <stdint.h>
#include <math.h> 

/*find and print the average*/ 
float BMM_mean(uint32_t * values, uint32_t size) 
{ 
    uint32_t sum = 0; 
    int i = 0; 
    for( i = 0; i < size; ++i) {
        sum += *(values + i); 
    }
    return (float)sum/(float)size; 
} 

float BMM_stdDev(uint32_t * values, uint32_t size, float m) 
{ 
    float sum = 0; 
    int i = 0; 
    if(m == -1)
        m = BMM_mean(values, size); 
    for(i = 0; i < size; ++i) {
        sum += pow((float)*(values + i) - m, 2); 
    }
    return sqrt(sum/(float)size); 
} 

float BMM_fakeStdDev(uint32_t * values, uint32_t size) 
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
