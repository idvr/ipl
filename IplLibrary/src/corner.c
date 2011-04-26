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
 *  corner.c
 *
 *   Corner points detection
 *       PTA      *pixCornerPointsMoravec()
 *       PTA      *fpixLocalMaxima()
 */

#include "ipl.h"

/*!
 *  pixCornerPointsMoravec()
 *
 *      Input:  (1, 8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              whsize (window half-width for measuring variance and local maxima)
 *              thresh (threshold value to determine the "cornerness" of a point)
 *              factor (factor used for subsampling of image)
 *
 *      Output: pfpixv (the measured variance per pixel)
 *
 *      Return: pta, or null on error
 *
 *  Notes:
 *      (1) Moravec algorithm calculates the minimum of the
 *			intensity variance across the 4 main directions of
 *			the neighborhood of square size (2 * whsize + 1) for
 *          each pixel in pix.
 *			If the variance value is above a specified threshold
 *          and it is a local maxima, then the pixel is considered
 *          a corner point.
 *		(2) The variance of intensity is given from the sum of
 *			the squares of the differences in the intensities
 *			between the neighborood and its adjacent regions.
 *      (3) The full width and height of the neighborhood is
 *          (2 * whsize + 1).
 *      (4) Require that w and h >= (2 * whsize + 1), where (w, h)
 *          are the dimensions of pixs.
 *      (5) If colormapped, remove to grayscale.
 */
PTA *
pixCornerPointsMoravec(PIX      *pixs,
					   l_int32   whsize,
					   l_int32   thresh,
					   l_int32   factor,
					   FPIX    **pfpixv)
{
	l_int32    i, j, k;
	l_int32    wsize, var, minvar;
	l_int32    w, h, d, wpl, wplv;
	l_uint32  *data, *line;
	l_float32 *datav, *linev;
	PIX       *pixt;
	FPIX      *fpixv;
	PTA       *locations, *corners;
	
	l_int16   shifts[] = {
			 1,  0,  1,  1,  0,  1, -1,  1,	/* ( 1,  0), ( 1,  1), ( 0,  1), (−1,  1) */
			-1,  0, -1, -1,  0, -1,  1, -1	/* (-1,  0), (-1, -1), ( 0, -1), ( 1, -1) */
		};
	
	PROCNAME("pixCornerPointsMoravec");
	
    if (!pixs)
        return (PTA *)ERROR_PTR("pixs not defined", procName, NULL);
	if (factor < 1)
        return (PTA *)ERROR_PTR("factor < 1", procName, NULL);
    
	/* Process the input image */
	if (pixGetColormap(pixs))
        pixt = pixRemoveColormap(pixs, REMOVE_CMAP_TO_GRAYSCALE);
    else if (pixGetDepth(pixs) == 32)
        pixt = pixConvertRGBToLuminance(pixs);
    else
        pixt = pixClone(pixs);
    
	/* Create a FPIX structure where to store
	   the computed variance values */
	pixGetDimensions(pixt, &w, &h, &d);
	wsize =  2 * whsize + 1;
    if (w < wsize || h < wsize)
        return (PTA *)ERROR_PTR("kernel too large", procName, NULL);
	if ((fpixv = fpixCreate(w, h)) == NULL)
        return (PTA *)ERROR_PTR("fpixv not made", procName, NULL);
	
	wplv = fpixGetWpl(fpixv);
	datav = fpixGetData(fpixv);
	wpl = pixGetWpl(pixt);
    data = pixGetData(pixt);
	
	/* Iterate through each pixel in the image */
	locations = ptaCreate((w * h) >> 8);	/* window area / 256 */
	for (i = whsize + 1; i < h - whsize - 1; i += factor) {
		linev = datav + i * wplv;
		for (j = whsize + 1; j < w - whsize - 1; j += factor) {
			
			/* Calculate the variance for each shift in the
			   4 main directions and select the minimum one */
			minvar = +1000000000;
			for (k = 0; k < 16; k += 2) {
				var = calculateVarianceLow(
						data, wpl, d, wsize,
						i - whsize, j - whsize,
						i - whsize + shifts[k + 1], j - whsize + shifts[k]
					);
				
				if (var < minvar)
					minvar = var;
			}
			
			/* Classify the pixel as a potential corner point if the minimum
			   of the variance values is above the specified threshold */
			linev[j] = minvar > thresh ? minvar : 0;
			if (minvar > thresh)
				ptaAddPt(locations, j, i);
		}
	}
	
	/* Find local maxima within the image
	   (sorted by decreasing variance) */
	corners = fpixLocalMaxima(fpixv, whsize, factor, locations);
	
	if (pfpixv)
		*pfpixv = fpixv;
	else
		fpixDestroy(&fpixv);
	
	/* Release the allocated resources */
	ptaDestroy(&locations);
	pixDestroy(&pixt);
	
	return corners;
}

/*!
 *  fpixLocalMaxima()
 *
 *      Input:  fpix
 *              whsize (window half-width for finding local maxima)
 *              factor (factor used for subsampling of image)
 *
 *      Return: pta, or null on error
 *
 *  Notes:
 *      (1) The full width and height of the neighborhood is
 *          (2 * whsize + 1).
 *      (2) Require that w and h >= (2 * whsize + 1), where (w, h)
 *          are the dimensions of fpix.
 */
PTA *
fpixLocalMaxima(FPIX    *fpix,
				l_int32  whsize,
				l_int32  factor,
				PTA     *locations)
{
	l_int32    i, j;
	l_int32    w, h, wsize, wpl;
	l_float32 *data, *line;
	l_float32  val;
	PTA       *points;
	
	PROCNAME("fpixLocalMaxima");
	
	if (!fpix)
        return (PTA *)ERROR_PTR("fpix not defined", procName, NULL);
	if (factor < 1)
        return (PTA *)ERROR_PTR("factor < 1", procName, NULL);
    
	wsize =  2 * whsize + 1;
	fpixGetDimensions(fpix, &w, &h);
    if (w < wsize || h < wsize) {
        return (PTA *)ERROR_PTR("kernel too large", procName, NULL);
    }
	
	points = ptaCreate(locations != NULL ? ptaGetCount(locations) : (w * h) >> 8);
	
	fpixGetDimensions(fpix, &w, &h);
	wpl = fpixGetWpl(fpix);
	data = fpixGetData(fpix);
	
	if (!locations) {
		/* Iterate through each pixel in the FPIX image */
		for (i = whsize + 1; i < h - whsize - 1; i += factor) {
			line = data + i * wpl;
			for (j = whsize + 1; j < w - whsize - 1; j += factor) {
				if ((val = line[j]) != 0) {
					/* If the pixel is a local maxima, then add it to
					   the list of returned points */
					if (isLocalMaximaLow(val, data, wpl, wsize, i - whsize, j - whsize) == 1)
						ptaAddPt(points, j, i);
				}
			}
		}
	}
	else {
		/* Iterate through each location in the FPIX image */
		l_int32 x, y;
		l_int32 count = ptaGetCount(locations);
		for (i = 0; i < count; i++) {
			ptaGetIPt(locations, i, &x, &y);
			val = *(data + y * w + x);
			/* If the pixel is a local maxima, then add it to
			   the list of returned points */
			if (isLocalMaximaLow(val, data, wpl, wsize, y - whsize, x - whsize) == 1)
				ptaAddPt(points, x, y);
		}
	}
	
	return points;
}