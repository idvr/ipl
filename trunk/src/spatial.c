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
 *  spatial.c
 *
 *   Statistical measures
 *       l_int32   pixGlobalStats()
 *
 *   Filters
 *       PIX      *pixAdaptiveMeanFilter()
 */

#include "ipl.h"

#include <math.h>

/*-----------------------------------------------------------------------*
 *                        Statistical measures                           *
 *-----------------------------------------------------------------------*/
/*!
 *  pixGlobalStats()
 *
 *      Input:  pixs   (8 bpp grayscale)
 *              &mean  (<optional return> pixs mean)
 *              &var   (<optional return> pixs variance)
 *              &std   (<optional return> pixs standard deviation)
 *      Return: 0 if OK; 1 on error
 */
l_int32
pixGlobalStats(PIX       *pixs,
			   l_float32 *mean,
			   l_float32 *var,
			   l_float32 *std)
{
	l_int32    w, h, d, i, j;
	l_int32    wpl;
	l_uint32  *data, *line;
	l_float32  m, v;
	
	PROCNAME("pixGlobalStats");
	
	if (!mean && !var && !std)
        return ERROR_INT("nothing to do", procName, 1);
	if (!pixs)
        return ERROR_INT("pixs not defined", procName, 1);
    pixGetDimensions(pixs, &w, &h, &d);
    if (d != 8)
        return ERROR_INT("pixs not 8 bpp", procName, 1);
	
	wpl = pixGetWpl(pixs);
	data = pixGetData(pixs);
	
	/* Calculate the global mean */
	m = 0.;
	for (i = 0; i < h; i++) {
		line = data + i * wpl;
		for (j = 0; j < w; j++)
			m += GET_DATA_BYTE(line, j);
	}
	m /= (w * h);
	
	/* Calculate the global variance */
	v = 0.;
	for (i = 0; i < h; i++) {
		line = data + i * wpl;
		for (j = 0; j < w; j++)
			v += pow((GET_DATA_BYTE(line, j) - m), 2);
	}
	v /= (w * h);
	
	if (mean)
		*mean = m;
	if (var)
		*var = v;
	if (std)
		*std = sqrt(v);
	
	return 0;
}

/*-----------------------------------------------------------------------*
 *                                Filters                                *
 *-----------------------------------------------------------------------*/
/*!
 *  pixAdaptiveMeanFilter()
 *
 *      Input:  pixs   (8 bpp grayscale)
 *              wc, hc (half width/height of convolution kernel)
 *              varn   (value of overall noise variance)
 *      Return: pixd (8 bpp, filtered image)
 *
 *  Notes:
 *      (1) The filter reduces gaussian noise, achieving results similar
 *          to the arithmetic and geometric mean filters but avoiding the
 *          considerable image blurring effect introduced by those filters.
 *      (2) The filter can be expressed mathematically by:
 *            f'(x, y) = g(x, y) - varN / varL * [ g(x, y) - meanL ]
 *          where:
 *            -- g(x, y) is the pixel at the center of local region S of
 *               width (2 * wc + 1) and height (2 * wh + 1)
 *            -- varN and varL are the overall noise variance (given in input)
 *               and local variance of S, respectively
 *            -- meanL is the local mean of S
 *      (3) Typically @varn is estimated by studying the PDFs produced by
 *          the camera or equipment sensors.
 */
PIX *
pixAdaptiveMeanFilter(PIX       *pixs,
					  l_int32    wc,
					  l_int32    hc,
					  l_float32  varn)
{
	l_int32    i, j, w, h, d, wplt, wpld, wincr, hincr;
	l_uint32   val;
	l_uint32  *datat, *datad, *linet, *lined;
	l_float32  norm, meanl, varl, ratio;
	PIX       *pixt, *pixd;
	
    PROCNAME("pixAdaptiveMeanFilter");
    
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
    pixGetDimensions(pixs, &w, &h, &d);
    if (d != 8)
        return (PIX *)ERROR_PTR("pixs not 8 bpp", procName, NULL);
    if (wc < 1 || hc < 1)
        return (PIX *)ERROR_PTR("wc and hc not >= 1", procName, NULL);
	
	/* Add wc to each side, and hc to top and bottom of the image,
	 * mirroring for accuracy and to avoid special-casing the boundary. */
    if ((pixt = pixAddMirroredBorder(pixs, wc, wc, hc, hc)) == NULL)
        return (PIX *)ERROR_PTR("pixt not made", procName, NULL);
	
	/* Place the filter center at (0, 0).  This is just a
	 * convenient location, because it allows us to perform
	 * the filtering over x:(0 ... w - 1) and y:(0 ... h - 1). */
	pixd = pixCreateTemplate(pixs);
	wplt = pixGetWpl(pixt);
    wpld = pixGetWpl(pixd);
	datat = pixGetData(pixt);
    datad = pixGetData(pixd);
	
    wincr = 2 * wc + 1;
    hincr = 2 * hc + 1;
	norm = 1.0 / (wincr * hincr);
	for (i = 0; i < h; i++) {
        linet = datat + (i + hc) * wplt;
		lined = datad + i * wpld;
        for (j = 0; j < w; j++) {
            /* Calculate mean intensity value */
			meanl = calculateLocalMeanLow(datat, wplt, wincr, hincr, i, j);
			/* Calculate local variance */
			varl = calculateLocalVarianceLow(datat, wplt, wincr, hincr, i, j, meanl);
			/* Account for special case in which varN is more than varL */
			ratio = (varn > varl) ? 1 : varn / varl;
			val = GET_DATA_BYTE(linet, j + wc);
			SET_DATA_BYTE(lined, j, (l_uint8) (val - ratio * (val - meanl)));
        } 
    }
	
	pixDestroy(&pixt);	
    return pixd;
}