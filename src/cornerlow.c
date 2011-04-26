/*
 * IPL (Image Processing Library)
 *
 * Copyright (C) 2011 Andrea Gagliardi La Gala.
 * Email: andrea@ameliemedia.com
 * Web: http://www.ameliemedia.com/
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
 *  cornerlow.c
 *
 *   Low-level implementation routines for corner points detection
 *       l_int32     multiplyLow()
 *       l_int32     calculateVarianceLow()
 *       l_float32   findLocalMaximaLow()
 *       l_int16     isLocalMaximaLow()
 */

#include "ipl.h"

/*-----------------------------------------------------------------------*
 *                      Low-level implementation                         *
 *-----------------------------------------------------------------------*/
/*!
 *  multiplyLow()
 *
 *  Multiply two pix blocks of equal dimension @wsize,
 *  where (i, j) is the origin of the first block and (k, l) is
 *  the origin of the second one.
 */
l_int32
multiplyLow(l_uint32 *data,
			l_int32   wpl,
			l_int32   d,
			l_int32   wsize,
			l_int32   i,
			l_int32   j,
			l_int32   k,
			l_int32   l)
{
	l_uint32 *line1, *line2;
	l_int32   val1, val2;
	l_int32   x, y, sum = 0;
	
	PROCNAME("multiplyLow");
	
	if (d == 1) {
		for (y = 0; y < wsize; y++) {
			
			line1 = data + ((i + y) * wpl);
			line2 = data + ((k + y) * wpl);
			
			for (x = 0; x < wsize; x++) {
				val1 = GET_DATA_BIT(line1, j + x);
				val2 = GET_DATA_BIT(line2, l + x);
				if (val1 != 0 && val2 != 0)
					sum += val1 * val2;
			}
		}
	}
	else if (d == 8) {
		for (y = 0; y < wsize; y++) {
			
			line1 = data + ((i + y) * wpl);
			line2 = data + ((k + y) * wpl);
			
			for (x = 0; x < wsize; x++) {
				val1 = GET_DATA_BYTE(line1, j + x);
				val2 = GET_DATA_BYTE(line2, l + x);
				if (val1 != 0 && val2 != 0)
					sum += val1 * val2;
			}
		}
	}
	
	return sum;
}

/*!
 *  calculateVarianceLow()
 *
 *  Calculate the variance of pixel intensity, which is given
 *  by the sum of the squares of the differences in the intensities
 *  between the neighborhood of a pixel under consideration and its
 *  adjacent regions.
 */
l_int32
calculateVarianceLow(l_uint32 *data,
					 l_int32   wpl,
					 l_int32   d,
					 l_int32   wsize,
					 l_int32   i,
					 l_int32   j,
					 l_int32   k,
					 l_int32   l)
{
	l_uint32 *srcline, *dstline;
	l_int32   srcval, dstval;
	l_int32   x, y, var = 0;
	
	PROCNAME("calculateVarianceLow");

	if (d == 1) {
		for (y = 0; y < wsize; y++) {
		
			srcline = data + ((i + y) * wpl);
			dstline = data + ((k + y) * wpl);
		
			for (x = 0; x < wsize; x++) {
				srcval = GET_DATA_BIT(srcline, j + x);
				dstval = GET_DATA_BIT(dstline, l + x);
				var += (srcval - dstval) * (srcval - dstval);
			}
		}
	}
	else if (d == 8) {
		for (y = 0; y < wsize; y++) {
			
			srcline = data + ((i + y) * wpl);
			dstline = data + ((k + y) * wpl);
			
			for (x = 0; x < wsize; x++) {
				srcval = GET_DATA_BYTE(srcline, j + x);
				dstval = GET_DATA_BYTE(dstline, l + x);
				var += (srcval - dstval) * (srcval - dstval);
			}
		}
	}
	
	return var;
}

/*!
 *  findLocalMaximaLow()
 *
 *  Return the local maxima value within the window with origin
 *  in (i, j) and of size @wsize.
 */
l_float32
findLocalMaximaLow(l_float32 *data,
				   l_int32    wpl,
				   l_int32    wsize,
				   l_int32    i,
				   l_int32    j)
{
	l_float32  *line;
	l_float32  val, maxval;
	l_int32    x, y;
	
	PROCNAME("findLocalMaximaLow");
	
	maxval = 0.;
	for (y = 0; y < wsize; y++) {
		line = data + ((i + y) * wpl);
		for (x = 0; x < wsize; x++) {
			val = line[j + x];
			if (val > maxval)
				maxval = val;
		}
	}
	
	return maxval;
}

/*!
 *  isLocalMaximaLow()
 *
 *  Return 1 if the val is a local maxima within the window with origin
 *  in (i, j) and of size @wsize.
 */
l_int16
isLocalMaximaLow(l_float32  val,
				 l_float32 *data,
				 l_int32    wpl,
				 l_int32    wsize,
				 l_int32    i,
				 l_int32    j)
{
	l_float32  *line;
	l_float32  currval;
	l_int32    x, y;
	
	PROCNAME("isLocalMaximaLow");
	
	for (y = 0; y < wsize; y++) {
		line = data + ((i + y) * wpl);
		for (x = 0; x < wsize; x++) {
			currval = line[j + x];
			if (val < currval)
				return 0;
		}
	}
	
	return 1;
}