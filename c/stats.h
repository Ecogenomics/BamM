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

#ifndef BM_STATS_H
#define BM_STATS_H

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * @abstract Calculate the mean of an array values
 *
 * @param  values       array of (integer) values
 * @param  size         size of values array
 * @return float        the caluclated mean
 */
float BM_mean(uint32_t * values, uint32_t size);

/*!
 * @abstract Calculate the mean number of non-zero array values
 *
 * @param  values       array of (integer) values
 * @param  size         size of values array
 * @return float        the caluclated mean
 */
float BM_nzmean(uint32_t * values, uint32_t size);

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

int cmpfunc_uint32(const void * a, const void * b);

/*!
 * @abstract Calculate the median of an array values
 *
 * @param  values       array of (integer) values
 * @param  size         size of values array
 * @return float        the caluclated median
 *
 * MODIFIES THE ORDER OF VALUES IN values ARRAY!
 */
float BM_median(uint32_t * values, uint32_t size);

/*!
 * @abstract Calculate the variance of an array values
 *
 * @param values array of (integer) values
 * @param size size of values array
 * @param m mean of the values array
 * @return float the calculated variance
 */
float BM_variance(uint32_t * values, uint32_t size, float m);

#ifdef __cplusplus
}
#endif

#endif // BM_STATS_H
