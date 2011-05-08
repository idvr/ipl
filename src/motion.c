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
 *  motion.c
 *
 *   Motion analysis
 *       PIX      *pixThresholdWithAbsDifference()
 */

#include "ipl.h"

/*!
 *  pixThresholdWithAbsDifference()
 *
 *      Input:  pixs1, pixs2  (both 8 bpp gray)
 *				threshold value
 *				pixgray <optional>
 *      Return: pixd (1 bpp), or null on error
 *
 *  Notes:
 *      (1) The depth of pixs1 and pixs2 must be equal.
 *      (2) Clips computation to the min size, aligning the UL corners
 *      (3) Computes the absolute value of the difference between
 *          each component value.
 *		(4) If the absolute value is less than the threshold value,
 *          the dest will be 0; otherwise, it will be 1
 */
PIX *
pixThresholdWithAbsDifference(PIX  *pixs1,
							  PIX  *pixs2,
							  l_int32 thresh,
							  PIX **pixgray)
{
	l_int32    w, h, w2, h2, d, wpls, wpld, wplg;
	l_uint32  *datas1, *datas2, *datad, *datag;
	PIX       *pixd;
	
    PROCNAME("pixThresholdWithAbsDifference");
	
    if (!pixs1)
        return (PIX *)ERROR_PTR("pixs1 not defined", procName, NULL);
    if (!pixs2)
        return (PIX *)ERROR_PTR("pixs2 not defined", procName, NULL);
    d = pixGetDepth(pixs1);
    if (d != pixGetDepth(pixs2))
        return (PIX *)ERROR_PTR("src1 and src2 depths unequal", procName, NULL);
    if (d != 8)
        return (PIX *)ERROR_PTR("depths not in {8}", procName, NULL);
	
    pixGetDimensions(pixs1, &w, &h, NULL);
    pixGetDimensions(pixs2, &w2, &h2, NULL);
    w = L_MIN(w, w2);
    h = L_MIN(h, h2);
    if ((pixd = pixCreate(w, h, 1)) == NULL)
        return (PIX *)ERROR_PTR("pixd not made", procName, NULL);
	
	if (pixgray != NULL) {
		if ((*pixgray = pixCreate(w, h, 8)) == NULL)
			return (PIX *)ERROR_PTR("pixgray not made", procName, NULL);
		datag = pixGetData(*pixgray);
		wplg = pixGetWpl(*pixgray);
	}
	else {
		datag = NULL;
		wplg = 0;
	}

    datas1 = pixGetData(pixs1);
    datas2 = pixGetData(pixs2);
    datad = pixGetData(pixd);
    wpls = pixGetWpl(pixs1);
    wpld = pixGetWpl(pixd);
	
    absThreholdWithDifferenceLow(datad, w, h, wpld, datas1, datas2, wpls, thresh, datag, wplg);
	
    return pixd;
}
