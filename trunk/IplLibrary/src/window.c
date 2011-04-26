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
 *  window.c
 *
 *   Global window functions
 *       PIX      *pixWindowBartlettGlobal()
 *       PIX      *pixWindowHannGlobal()
 *       PIX      *pixWindowTukeyGlobal()
 *       PIX      *pixWindowBlackmanGlobal()
 *
 *   Local window functions
 *       PIX      *pixWindowBartlettLocal()
 *       PIX      *pixWindowHannLocal()
 *       PIX      *pixWindowTukeyLocal()
 *       PIX      *pixWindowBlackmanLocal()
 *
 *   Pix padding
 *       PIX      *pixAddPadding()
 *       PIX      *pixRemovePadding()
 *
 *  Window functions are used to gradually reduce the intensity values
 *  to zero around a center pixel.
 *
 *  One possible application of these functions is to smooth edges at
 *  the border of an image before transforming the image with the
 *  Discrete Fourier Transform.
 *  The DFT assumes that an image is periodic (that is, it repeats
 *  itself infinite times along the four main directions), therefore
 *  the left and right sides and top and bottom should coincide.
 *  This is seldom the case, and discontinuity in the periodicity of
 *  the image distorts the Fourier spectrum.
 *  The proper way to work around this issue is to pad the image by
 *  appending zeroes on the right and bottom sides. However the
 *  drawback is the doubling of the memory requirements.
 *  Window functions, even if less accurate than padding, can be used
 *  instead to minimize the discontinuities at the edges of the image.
 *
 *  Global window functions act on the entire image, having their
 *  center in the middle of the image.
 *  Local window functions acts only on a region of the image, whose
 *  center and radius are specified in the input parameters.
 *  
 */

#include "ipl.h"
#include <math.h>

/*------------------------------------------------------------------*
 *                      Global window functions                     *
 *------------------------------------------------------------------*/
/*!
 *  pixWindowBartlettGlobal()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              typeflag (0 to use Bartlett coefficients; or
 *                        L_USE_TRIANG to use triangular window)
 *
 *      Return: pixd (windowed, 8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) Together with the rectangle window function, the Bartlett
 *          window is the simplest and most computationally efficient
 *          function.
 *      (2) The Bartlett window is very similar to a triangular window.
 *          The different is that it always hs zero values at the
 *          endpoints, while the triangular window is non-zero at those
 *          points.
 */
PIX *
pixWindowBartlettGlobal(PIX      *pixs,
						l_int32   typeflag)
{
	l_int32    w, h, d;
	
	PROCNAME("pixWindowBartlettGlobal");
	
	pixGetDimensions(pixs, &w, &h, &d);	
	return pixWindowBartlettLocal(pixs,
								  typeflag,
								  (l_int32) w / 2,
								  (l_int32) h / 2,
								  (l_int32) L_MIN(w / 2, h / 2));
}

/*!
 *  pixWindowHannGlobal()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              typeflag (0 to use Hann coefficients; or
 *                        L_USE_HAMMING to use Hamming ones)
 *
 *      Return: pixd (windowed, 8 or 32 bpp), or null on error
 */
PIX *
pixWindowHannGlobal(PIX      *pixs,
					l_int32   typeflag)
{
	l_int32    w, h, d;

	PROCNAME("pixWindowHannGlobal");
	
	pixGetDimensions(pixs, &w, &h, &d);
	return pixWindowHannLocal(pixs,
							  typeflag,
							  (l_int32) w / 2,
							  (l_int32) h / 2,
							  (l_int32) L_MIN(w / 2, h / 2));
}

/*!
 *  pixWindowTukeyGlobal()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              alpha (fraction of the radius within which the
 *                     intensity values in pix are not changed)
 *
 *      Return: pixd (windowed, 8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) @alpha values can range between 0 and 1 (typically = 0.5).
 *          When @alpha is 0 the function becomes a rectangle window.
 *          When @alpha is 1 the function becomes a Hann window.
 *      (2) The Tukey window function is also known as the tapered
 *          cosine window.
 */
PIX *
pixWindowTukeyGlobal(PIX       *pixs,
					 l_float32  alpha)
{
	l_int32    w, h, d;
	
	PROCNAME("pixWindowTukeyGlobal");
	
	pixGetDimensions(pixs, &w, &h, &d);	
	return pixWindowTukeyLocal(pixs,
							   alpha,
							   (l_int32) w / 2,
							   (l_int32) h / 2,
							   (l_int32) L_MIN(w / 2, h / 2));
}

/*!
 *  pixWindowBlackmanGlobal()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              alpha (term of the Blackman window. By convention
 *                     it is set to 0.16)
 *
 *      Return: pixd (windowed, 8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) @alpha value is typically set to 0.16.
 */
PIX *
pixWindowBlackmanGlobal(PIX       *pixs,
						l_float32  alpha)
{
	l_int32    w, h, d;
	
	PROCNAME("pixWindowBlackmanGlobal");
	
	pixGetDimensions(pixs, &w, &h, &d);	
	return pixWindowBlackmanLocal(pixs,
								  alpha,
								  (l_int32) w / 2,
								  (l_int32) h / 2,
								  (l_int32) L_MIN(w / 2, h / 2));
}

/*------------------------------------------------------------------*
 *                      Local window functions                      *
 *------------------------------------------------------------------*/
/*!
 *  pixWindowBartlettLocal()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              typeflag (0 to use Bartlett coefficients; or
 *                        L_USE_TRIANG to use triangular window)
 *              x, y (the center of the circle of @r radius)
 *              r (the radius within which to apply the windowing function,
 *                 and outside which all pixels are set to 0)
 *
 *      Return: pixd (windowed, 8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) Together with the rectangle window function, the Bartlett
 *          window is the simplest and most computationally efficient
 *          function.
 *      (2) The Bartlett window is very similar to a triangular window.
 *          The different is that it always hs zero values at the
 *          endpoints, while the triangular window is non-zero at those
 *          points.
 */
PIX *
pixWindowBartlettLocal(PIX      *pixs,
					   l_int32   typeflag,
					   l_int32   x,
					   l_int32   y,
					   l_int32   r)
{
	l_int32    w, h, d;
	l_int32    i, j, xstart, ystart, xend, yend, wpl, val;
	l_uint32   uval;
	l_uint32  *datas, *datad, *lines, *lined;
	l_float32  n, dist, win;
	PIX       *pixt, *pixd;
	
	PROCNAME("pixWindowBartlettLocal");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	if (typeflag != 0 && typeflag != L_USE_TRIANG)
        return (PIX *)ERROR_PTR("invalid typeflag", procName, NULL);
	if (r < 0)
		return (PIX *)ERROR_PTR("invalid r", procName, NULL);
	
	/* Remove colormap if necessary */
	pixGetDimensions(pixs, &w, &h, &d);
    if ((d == 2 || d == 4 || d == 8) && pixGetColormap(pixs)) {
        L_WARNING("pix has colormap; removing", procName);
        pixt = pixRemoveColormap(pixs, REMOVE_CMAP_BASED_ON_SRC);
        d = pixGetDepth(pixt);
    }
    else
        pixt = pixClone(pixs);
	
    if (d != 8 && d != 32) {
        pixDestroy(&pixt);
        return (PIX *)ERROR_PTR("depth not 8 or 32 bpp", procName, NULL);
    }
	
	if ((pixd = pixCreateTemplate(pixt)) == NULL) {
		pixDestroy(&pixt);
		return (PIX *)ERROR_PTR("pixd not made", procName, NULL);
	}
		
	datas = pixGetData(pixt);
	datad = pixGetData(pixd);
	wpl = pixGetWpl(pixd);
	
	ystart = L_MAX(0, y - r);
	xstart = L_MAX(0, x - r);
	yend = L_MIN(h, y + r);
	xend = L_MIN(w, x + r);
	n = typeflag == L_USE_TRIANG ? r + 0.5 : r;
	for (i = ystart; i < yend; i++) {
		lines = datas + i * wpl;
		lined = datad + i * wpl;
		for (j = xstart; j < xend; j++) {
			dist = sqrt(pow(j - x, 2) + pow(i - y, 2));
			win = dist <= r ? win = 1 - (dist / n) : 0;
			if (d == 8) {
                val = GET_DATA_BYTE(lines, j);
                SET_DATA_BYTE(lined, j, val * win);
			}
			else { /* d == 32 */
                uval = GET_DATA_FOUR_BYTES(lines, j);
				uval =
					(l_int32) (((uval >> L_RED_SHIFT) & 0xff) * win) << L_RED_SHIFT |
					(l_int32) (((uval >> L_GREEN_SHIFT) & 0xff) * win) << L_GREEN_SHIFT |
					(l_int32) (((uval >> L_BLUE_SHIFT) & 0xff) * win) << L_BLUE_SHIFT;
				SET_DATA_FOUR_BYTES(lined, j, uval);
            }
		}
	}
	
	pixDestroy(&pixt);
	return pixd;
}

/*!
 *  pixWindowHannLocal()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              typeflag (0 to use Hann coefficients; or
 *                        L_USE_HAMMING to use Hamming ones)
 *              x, y (the center of the circle of @r radius)
 *              r (the radius within which to apply the windowing function,
 *                 and outside which all pixels are set to 0)
 *
 *      Return: pixd (windowed, 8 or 32 bpp), or null on error
 */
PIX *
pixWindowHannLocal(PIX      *pixs,
				   l_int32   typeflag,
				   l_int32   x,
				   l_int32   y,
				   l_int32   r)
{
	l_int32    w, h, d;
	l_int32    i, j, xstart, ystart, xend, yend, wpl, val;
	l_uint32   uval;
	l_uint32  *datas, *datad, *lines, *lined;
	l_float32  a, b, dist, win;
	PIX       *pixt, *pixd;
	
	PROCNAME("pixWindowHannLocal");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	if (typeflag != 0 && typeflag != L_USE_HAMMING)
        return (PIX *)ERROR_PTR("invalid typeflag", procName, NULL);
	if (r < 0)
		return (PIX *)ERROR_PTR("invalid r", procName, NULL);
	
	/* Remove colormap if necessary */
	pixGetDimensions(pixs, &w, &h, &d);
    if ((d == 2 || d == 4 || d == 8) && pixGetColormap(pixs)) {
        L_WARNING("pix has colormap; removing", procName);
        pixt = pixRemoveColormap(pixs, REMOVE_CMAP_BASED_ON_SRC);
        d = pixGetDepth(pixt);
    }
    else
        pixt = pixClone(pixs);
	
    if (d != 8 && d != 32) {
        pixDestroy(&pixt);
        return (PIX *)ERROR_PTR("depth not 8 or 32 bpp", procName, NULL);
    }
	
	if ((pixd = pixCreateTemplate(pixt)) == NULL) {
		pixDestroy(&pixt);
		return (PIX *)ERROR_PTR("pixd not made", procName, NULL);
	}
	
	datas = pixGetData(pixt);
	datad = pixGetData(pixd);
	wpl = pixGetWpl(pixd);
	
	/* Setup the coefficients to Hann or Hamming depending
	   on the specified @typeflag */
	a = typeflag == L_USE_HAMMING ? 0.54 : 0.5;
	b = 1 - a;
	
	ystart = L_MAX(0, y - r);
	xstart = L_MAX(0, x - r);
	yend = L_MIN(h, y + r);
	xend = L_MIN(w, x + r);
	for (i = ystart; i < yend; i++) {
		lines = datas + i * wpl;
		lined = datad + i * wpl;
		for (j = xstart; j < xend; j++) {
			dist = sqrt(pow(j - x, 2) + pow(i - y, 2));
			win = (dist <= r) ? win = a - b * cos(M_PI * (1 - (dist / r))) : 0;
			if (d == 8) {
                val = GET_DATA_BYTE(lines, j);
                SET_DATA_BYTE(lined, j, val * win);
			}
			else { /* d == 32 */
                uval = GET_DATA_FOUR_BYTES(lines, j);
				uval =
				(l_int32) (((uval >> L_RED_SHIFT) & 0xff) * win) << L_RED_SHIFT |
				(l_int32) (((uval >> L_GREEN_SHIFT) & 0xff) * win) << L_GREEN_SHIFT |
				(l_int32) (((uval >> L_BLUE_SHIFT) & 0xff) * win) << L_BLUE_SHIFT;
				SET_DATA_FOUR_BYTES(lined, j, uval);
            }
		}
	}
	
	pixDestroy(&pixt);
	return pixd;
}

/*!
 *  pixWindowTukeyLocal()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              alpha (fraction of the @r radius within which the
 *                     intensity values in pix are not changed)
 *              x, y (the center of the circle of @r radius)
 *              r (the radius within which to apply the windowing function,
 *                 and outside which all pixels are set to 0)
 *
 *      Return: pixd (windowed, 8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) @alpha values can range between 0 and 1 (typically = 0.5).
 *          When @alpha is 0 the function becomes a rectangle window.
 *          When @alpha is 1 the function becomes a Hann window.
 *      (2) The Tukey window function is also known as the tapered
 *          cosine window.
 */
PIX *
pixWindowTukeyLocal(PIX       *pixs,
					l_float32  alpha,
					l_int32    x,
					l_int32    y,
					l_int32    r)
{
	l_int32    w, h, d;
	l_int32    i, j, xstart, ystart, xend, yend, wpl, val;
	l_uint32   uval;
	l_uint32  *datas, *datad, *lines, *lined;
	l_float32  fract, dist, win;
	PIX       *pixt, *pixd;
	
	PROCNAME("pixWindowTukeyLocal");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	if (alpha < 0.0 || alpha > 1.0)
        return (PIX *)ERROR_PTR("invalid alpha", procName, NULL);
	if (r < 0)
		return (PIX *)ERROR_PTR("invalid r", procName, NULL);
	
	/* Remove colormap if necessary */
	pixGetDimensions(pixs, &w, &h, &d);
    if ((d == 2 || d == 4 || d == 8) && pixGetColormap(pixs)) {
        L_WARNING("pix has colormap; removing", procName);
        pixt = pixRemoveColormap(pixs, REMOVE_CMAP_BASED_ON_SRC);
        d = pixGetDepth(pixt);
    }
    else
        pixt = pixClone(pixs);
	
    if (d != 8 && d != 32) {
        pixDestroy(&pixt);
        return (PIX *)ERROR_PTR("depth not 8 or 32 bpp", procName, NULL);
    }
	
	if ((pixd = pixCreateTemplate(pixt)) == NULL) {
		pixDestroy(&pixt);
		return (PIX *)ERROR_PTR("pixd not made", procName, NULL);
	}
	
	datas = pixGetData(pixt);
	datad = pixGetData(pixd);
	wpl = pixGetWpl(pixd);
	
	ystart = L_MAX(0, y - r);
	xstart = L_MAX(0, x - r);
	yend = L_MIN(h, y + r);
	xend = L_MIN(w, x + r);
	fract = r * (1 - alpha);
	for (i = ystart; i < yend; i++) {
		lines = datas + i * wpl;
		lined = datad + i * wpl;
		for (j = xstart; j < xend; j++) {
			dist = sqrt(pow(j - x, 2) + pow(i - y, 2));
			if (dist <= fract)
				win = 1;
			else
				win = (dist <= r) ? win = 0.5 - 0.5 * cos(M_PI * (1 - (dist - fract) / (r - fract))) : 0;
			if (d == 8) {
                val = GET_DATA_BYTE(lines, j);
                SET_DATA_BYTE(lined, j, val * win);
			}
			else { /* d == 32 */
                uval = GET_DATA_FOUR_BYTES(lines, j);
				uval =
				(l_int32) (((uval >> L_RED_SHIFT) & 0xff) * win) << L_RED_SHIFT |
				(l_int32) (((uval >> L_GREEN_SHIFT) & 0xff) * win) << L_GREEN_SHIFT |
				(l_int32) (((uval >> L_BLUE_SHIFT) & 0xff) * win) << L_BLUE_SHIFT;
				SET_DATA_FOUR_BYTES(lined, j, uval);
            }
		}
	}
	
	pixDestroy(&pixt);
	return pixd;
}

/*!
 *  pixWindowBlackmanLocal()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              alpha (term of the Blackman window. By convention
 *                     it is set to 0.16)
 *              x, y (the center of the circle of @r radius)
 *              r (the radius within which to apply the windowing function,
 *                 and outside which all pixels are set to 0)
 *
 *      Return: pixd (windowed, 8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) @alpha value is typically set to 0.16.
 */
PIX *
pixWindowBlackmanLocal(PIX       *pixs,
					   l_float32  alpha,
					   l_int32    x,
					   l_int32    y,
					   l_int32    r)
{
	l_int32    w, h, d;
	l_int32    i, j, xstart, ystart, xend, yend, wpl, val;
	l_uint32   uval;
	l_uint32  *datas, *datad, *lines, *lined;
	l_float32  a, b, dist, win;
	PIX       *pixt, *pixd;
	
	PROCNAME("pixWindowBlackmanLocal");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	if (r < 0)
		return (PIX *)ERROR_PTR("invalid r", procName, NULL);
	
	/* Remove colormap if necessary */
	pixGetDimensions(pixs, &w, &h, &d);
    if ((d == 2 || d == 4 || d == 8) && pixGetColormap(pixs)) {
        L_WARNING("pix has colormap; removing", procName);
        pixt = pixRemoveColormap(pixs, REMOVE_CMAP_BASED_ON_SRC);
        d = pixGetDepth(pixt);
    }
    else
        pixt = pixClone(pixs);
	
    if (d != 8 && d != 32) {
        pixDestroy(&pixt);
        return (PIX *)ERROR_PTR("depth not 8 or 32 bpp", procName, NULL);
    }
	
	if ((pixd = pixCreateTemplate(pixt)) == NULL) {
		pixDestroy(&pixt);
		return (PIX *)ERROR_PTR("pixd not made", procName, NULL);
	}
	
	datas = pixGetData(pixt);
	datad = pixGetData(pixd);
	wpl = pixGetWpl(pixd);
	
	/* Setup the coefficients to Hann or Hamming depending
	 on the specified @typeflag */
	a = (1 - alpha) / 2;
	b = alpha / 2;
	
	ystart = L_MAX(0, y - r);
	xstart = L_MAX(0, x - r);
	yend = L_MIN(h, y + r);
	xend = L_MIN(w, x + r);
	for (i = ystart; i < yend; i++) {
		lines = datas + i * wpl;
		lined = datad + i * wpl;
		for (j = xstart; j < xend; j++) {
			win = 0;
			dist = sqrt(pow(j - x, 2) + pow(i - y, 2));
			if (dist <= r)
				win = a - 0.5 * cos(M_PI * (1 - (dist / r))) + b * cos(2 * M_PI * (1 - (dist / r)));
			if (d == 8) {
                val = GET_DATA_BYTE(lines, j);
                SET_DATA_BYTE(lined, j, val * win);
			}
			else { /* d == 32 */
                uval = GET_DATA_FOUR_BYTES(lines, j);
				uval =
				(l_int32) (((uval >> L_RED_SHIFT) & 0xff) * win) << L_RED_SHIFT |
				(l_int32) (((uval >> L_GREEN_SHIFT) & 0xff) * win) << L_GREEN_SHIFT |
				(l_int32) (((uval >> L_BLUE_SHIFT) & 0xff) * win) << L_BLUE_SHIFT;
				SET_DATA_FOUR_BYTES(lined, j, uval);
            }
		}
	}
	
	pixDestroy(&pixt);
	return pixd;
}

/*------------------------------------------------------------------*
 *                             Pix padding                          *
 *------------------------------------------------------------------*/
/*!
 *  pixAddPadding()
 *
 *      Input:  pix (all depths; colormap ok)
 *      Return: pixd (with the added padding), or null on error
 *
 *  Notes:
 *      (1) Pad pix by adding 0 valued pixels at the right and bottom.
 *          This is helpful to construct a periodic image for obtaining
 *          a proper Discrete Fourier Transform.
 *          The padded pix is twice in width and height than the
 *          original source image.
 */
PIX *
pixAddPadding(PIX       *pixs)
{
	l_int32    w, h, d;
	PIX       *pixd;
	
	PROCNAME("pixAddPadding");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
	pixGetDimensions(pixs, &w, &h, &d);
	pixd = pixAddBorderGeneral(pixs, 0, w, 0, h, 0x00);
	
	return pixd;
}

/*!
 *  pixRemovePadding()
 *
 *      Input:  pix (all depths; colormap ok)
 *      Return: pixd (with padding removed), or null on error
 *
 *  Notes:
 *      (1) Remove padding from pix by clipping pixels at the right
 *          and bottom.
 *          The output pix is half in width and height than the
 *          original source image.
 */
PIX *
pixRemovePadding(PIX       *pixs)
{
	l_int32    w, h, d;
	PIX       *pixd;
	
	PROCNAME("pixRemovePadding");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
	pixGetDimensions(pixs, &w, &h, &d);
	pixd = pixRemoveBorderGeneral(pixs, 0, w / 2, 0, h / 2);
	
	return pixd;
}