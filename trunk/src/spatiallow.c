/*
 * IPL (Image Processing Library)
 *
 * Copyright (C) 2011 Andrea Gagliardi La Gala.
 * Email: andrea.lagala@gmail.com
 * Web: http://code.google.com/p/ipl/
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * IPL depends on and extends:
 *   - Leptonica library by Dan S. Bloomberg (http://www.leptonica.com/)
 *
 * and (optionally) integrates:
 *   - FFTW library by Matteo Frigo and Steven G. Johnson (http://www.fftw.org/)
 *   - fuzzylite C++ library by Juan Rada-Vilela (http://code.google.com/p/fuzzy-lite/)
 *
 * to achieve fast and efficient image filtering and manipulation within
 * the frequency domain (using the Discrete Fourier Transform and its inverse)
 * and by using Fuzzy Sets and Logic.
 *
 */

/*
 *  spatiallow.c
 *
 *   Low-level implementation routines for filters
 *       l_float32     calculateLocalMeanLow()
 *       l_float32     calculateLocalVarianceLow()
 */

#include "ipl.h"
#include <math.h>

/*-----------------------------------------------------------------------*
 *                      Low-level implementation                         *
 *-----------------------------------------------------------------------*/
/*!
 *  calculateLocalMeanLow()
 *
 *  Calculate the mean pixel intensity within a region.
 */
l_float32
calculateLocalMeanLow(l_uint32  *data,
					  l_int32    wpl,
					  l_int32    width,
					  l_int32    height,
					  l_int32    i,
					  l_int32    j)
{
	l_uint32 *line;
	l_int32   x, y;
	l_float32 mean = 0;
	
	PROCNAME("calculateLocalMeanLow");
	
	
	for (y = 0; y < height; y++) {
		line = data + ((i + y) * wpl);
		for (x = 0; x < width; x++) {
			mean += GET_DATA_BYTE(line, j + x);
		}
	}
	
	return (mean /= (width * height));
}

/*!
 *  calculateLocalVarianceLow()
 *
 *  Calculate the variance of pixel intensity within a region,
 *  which is given by the sum of the squares of the differences
 *  in the intensities between each pixel in the region and the
 *  regional mean value.
 */
l_float32
calculateLocalVarianceLow(l_uint32  *data,
						  l_int32    wpl,
						  l_int32    width,
						  l_int32    height,
						  l_int32    i,
						  l_int32    j,
						  l_float32  mean)
{
	l_uint32 *line;
	l_int32   val;
	l_int32   x, y;
	l_float32 var = 0;
	
	PROCNAME("calculateLocalVarianceLow");
	
	for (y = 0; y < height; y++) {
		line = data + ((i + y) * wpl);
		for (x = 0; x < width; x++) {
			val = GET_DATA_BYTE(line, j + x);
			var += pow(val - mean, 2);
		}
	}
	
	return (var /= (width * height));
}