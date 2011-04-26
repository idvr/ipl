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
 *  motionlow.c
 *
 *   Low-level implementation routines for motion analysis
 *       void     absThreholdWithDifferenceLow()
 */

#include "ipl.h"


/*-----------------------------------------------------------------------*
 *                      Low-level implementation                         *
 *-----------------------------------------------------------------------*/					 
/*!
 *  absThreholdWithDifferenceLow()
 *
 *  Finds the absolute value of the difference of each pixel,
 *  and threshold based on the specified value.
 *  The results are written into datad.  
 */
void
absThreholdWithDifferenceLow(l_uint32  *datad,
							 l_int32    w,
							 l_int32    h,
							 l_int32    wpld,
							 l_uint32  *datas1,
							 l_uint32  *datas2,
							 l_int32    wpls,
							 l_int32	thresh,
							 l_uint32  *datag,
							 l_int32    wplg)
{
	l_int32    i, j, val1, val2, diff;
	l_uint32  *lines1, *lines2, *lined, *lineg;
	
    PROCNAME("absThreholdWithDifferenceLow");
	
	for (i = 0; i < h; i++) {
		lines1 = datas1 + i * wpls;
		lines2 = datas2 + i * wpls;
		lined = datad + i * wpld;
		if (datag != NULL) lineg = datag + i * wplg;
		for (j = 0; j < w; j++) {
			val1 = GET_DATA_BYTE(lines1, j);
			val2 = GET_DATA_BYTE(lines2, j);
			diff = L_ABS(val1 - val2);
			if (diff > thresh)
				SET_DATA_BIT(lined, j);
			
			if (datag != NULL)
				SET_DATA_BYTE(lineg, j, diff);
		}
	}
	
    return;
}

