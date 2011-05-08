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
 *  noise.c
 *
 *   Additive noise
 *       PIX          *pixAddNoiseUniform()
 *       PIX          *pixAddNoiseGaussian()
 *       PIX          *pixAddNoiseRayleigh()
 *       PIX          *pixAddNoiseErlang()
 *       PIX          *pixAddNoiseExponential()
 *       PIX          *pixAddNoiseImpulseBipolar()
 *       PIX          *pixAddNoiseImpulseUnipolar()
 *       PIX          *pixAddNoiseImpulse()
 */

#include "ipl.h"
#include <math.h>

/*--------------------------------------------------------------------*
 *                           Additive noise                           *
 *--------------------------------------------------------------------*/
/*!
 *  pixAddNoiseUniform()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              min, max (lower and upper bounds of noise variation)
 *      Return: pixd (8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) Return pix with additive uniform noise.
 *      (2) If colormapped, remove to grayscale.
 */
PIX *
pixAddNoiseUniform(PIX       *pixs,
				   l_float32  min,
				   l_float32  max)
{
	l_float64   r;
	l_int32     w, h, d, i, j, val, wpl;
	l_uint32    uval;
	l_uint32   *datas, *datad, *lines, *lined;
	PIX        *pixt, *pixd;
	
    PROCNAME("pixAddNoiseUniform");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
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
	for (i = 0; i < h; i++) {
		lines = datas + i * wpl;
		lined = datad + i * wpl;
		for (j = 0; j < w; j++) {
			r = randomUniformRange(min, max);
			if (d == 8) {
				val = GET_DATA_BYTE(lines, j);
				SET_DATA_BYTE(lined, j, L_MAX(0, L_MIN(val + r, 255)));
			}
			else { /* d == 32 */
				uval = GET_DATA_FOUR_BYTES(lines, j);
				uval =
				(l_int32) L_MAX(0, L_MIN((((uval >> L_RED_SHIFT) & 0xff) + r), 255)) << L_RED_SHIFT |
				(l_int32) L_MAX(0, L_MIN((((uval >> L_GREEN_SHIFT) & 0xff) + r), 255)) << L_GREEN_SHIFT |
				(l_int32) L_MAX(0, L_MIN((((uval >> L_BLUE_SHIFT) & 0xff) + r), 255)) << L_BLUE_SHIFT;
				SET_DATA_FOUR_BYTES(lined, j, uval);
			}
		}
	}
	
	pixDestroy(&pixt);
	return pixd;
}

/*!
 *  pixAddNoiseGaussian()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              mean (mean value of the gaussian PDF)
 *              std  (standard deviation valie of the PDF)
 *      Return: pixd (8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) Return pix with additive Gaussian noise.
 *      (2) If colormapped, remove to grayscale.
 */
PIX *
pixAddNoiseGaussian(PIX       *pixs,
					l_float32  mean,
					l_float32  std)
{
	l_float64   r;
	l_int32     w, h, d, i, j, val, wpl;
	l_uint32    uval;
	l_uint32   *datas, *datad, *lines, *lined;
	PIX        *pixt, *pixd;
	
    PROCNAME("pixAddNoiseGaussian");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
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
	for (i = 0; i < h; i++) {
		lines = datas + i * wpl;
		lined = datad + i * wpl;
		for (j = 0; j < w; j++) {
			r = randomGaussianThreadSafe(mean, std);
			if (d == 8) {
				val = GET_DATA_BYTE(lines, j);
				SET_DATA_BYTE(lined, j, L_MAX(0, L_MIN(val + r, 255)));
			}
			else { /* d == 32 */
				uval = GET_DATA_FOUR_BYTES(lines, j);
				uval =
				(l_int32) L_MAX(0, L_MIN((((uval >> L_RED_SHIFT) & 0xff) + r), 255)) << L_RED_SHIFT |
				(l_int32) L_MAX(0, L_MIN((((uval >> L_GREEN_SHIFT) & 0xff) + r), 255)) << L_GREEN_SHIFT |
				(l_int32) L_MAX(0, L_MIN((((uval >> L_BLUE_SHIFT) & 0xff) + r), 255)) << L_BLUE_SHIFT;
				SET_DATA_FOUR_BYTES(lined, j, uval);
			}
		}
	}
	
	pixDestroy(&pixt);
	return pixd;
}

/*!
 *  pixAddNoiseRayleigh()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              sigma
 *              var (variation of additive noise)
 *      Return: pixd (8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) Return pix with additive Rayleigh noise.
 *      (2) If colormapped, remove to grayscale.
 */
PIX *
pixAddNoiseRayleigh(PIX       *pixs,
					l_float32  sigma,
					l_int32    var)
{
	l_float64   r;
	l_int32     w, h, d, i, j, val, wpl;
	l_uint32    uval;
	l_uint32   *datas, *datad, *lines, *lined;
	PIX        *pixt, *pixd;
	
    PROCNAME("pixAddNoiseRayleigh");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
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
	for (i = 0; i < h; i++) {
		lines = datas + i * wpl;
		lined = datad + i * wpl;
		for (j = 0; j < w; j++) {
			r = randomRayleigh(sigma) * var;
			if (d == 8) {
				val = GET_DATA_BYTE(lines, j);
				SET_DATA_BYTE(lined, j, L_MAX(0, L_MIN(val + r, 255)));
			}
			else { /* d == 32 */
				uval = GET_DATA_FOUR_BYTES(lines, j);
				uval =
				(l_int32) L_MAX(0, L_MIN((((uval >> L_RED_SHIFT) & 0xff) + r), 255)) << L_RED_SHIFT |
				(l_int32) L_MAX(0, L_MIN((((uval >> L_GREEN_SHIFT) & 0xff) + r), 255)) << L_GREEN_SHIFT |
				(l_int32) L_MAX(0, L_MIN((((uval >> L_BLUE_SHIFT) & 0xff) + r), 255)) << L_BLUE_SHIFT;
				SET_DATA_FOUR_BYTES(lined, j, uval);
			}
		}
	}
	
	pixDestroy(&pixt);
	return pixd;
}

/*!
 *  pixAddNoiseErlang()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              k   (shape of the Erlang distribution)
 *              var (variation of additive noise)
 *      Return: pixd (8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) Return pix with additive Erlang noise.
 *      (2) If colormapped, remove to grayscale.
 */
PIX *
pixAddNoiseErlang(PIX       *pixs,
				  l_int32    k,
				  l_int32    var)
{
	l_float64   r;
	l_int32     w, h, d, i, j, val, wpl;
	l_uint32    uval;
	l_uint32   *datas, *datad, *lines, *lined;
	PIX        *pixt, *pixd;
	
    PROCNAME("pixAddNoiseErlang");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
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
	for (i = 0; i < h; i++) {
		lines = datas + i * wpl;
		lined = datad + i * wpl;
		for (j = 0; j < w; j++) {
			r = randomErlang(k) * var;
			if (d == 8) {
				val = GET_DATA_BYTE(lines, j);
				SET_DATA_BYTE(lined, j, L_MAX(0, L_MIN(val + r, 255)));
			}
			else { /* d == 32 */
				uval = GET_DATA_FOUR_BYTES(lines, j);
				uval =
				(l_int32) L_MAX(0, L_MIN((((uval >> L_RED_SHIFT) & 0xff) + r), 255)) << L_RED_SHIFT |
				(l_int32) L_MAX(0, L_MIN((((uval >> L_GREEN_SHIFT) & 0xff) + r), 255)) << L_GREEN_SHIFT |
				(l_int32) L_MAX(0, L_MIN((((uval >> L_BLUE_SHIFT) & 0xff) + r), 255)) << L_BLUE_SHIFT;
				SET_DATA_FOUR_BYTES(lined, j, uval);
			}
		}
	}
	
	pixDestroy(&pixt);
	return pixd;
}

/*!
 *  pixAddNoiseExponential()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              lambda
 *              var (variation of additive noise)
 *      Return: pixd (8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) Return pix with additive exponential noise.
 *      (2) If colormapped, remove to grayscale.
 */
PIX *
pixAddNoiseExponential(PIX       *pixs,
					   l_float32  lambda,
					   l_int32    var)
{
	l_float64   r;
	l_int32     w, h, d, i, j, val, wpl;
	l_uint32    uval;
	l_uint32   *datas, *datad, *lines, *lined;
	PIX        *pixt, *pixd;
	
    PROCNAME("pixAddNoiseExponential");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
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
	for (i = 0; i < h; i++) {
		lines = datas + i * wpl;
		lined = datad + i * wpl;
		for (j = 0; j < w; j++) {
			r = randomExponential(lambda) * var;
			if (d == 8) {
				val = GET_DATA_BYTE(lines, j);
				SET_DATA_BYTE(lined, j, L_MAX(0, L_MIN(val + r, 255)));
			}
			else { /* d == 32 */
				uval = GET_DATA_FOUR_BYTES(lines, j);
				uval =
				(l_int32) L_MAX(0, L_MIN((((uval >> L_RED_SHIFT) & 0xff) + r), 255)) << L_RED_SHIFT |
				(l_int32) L_MAX(0, L_MIN((((uval >> L_GREEN_SHIFT) & 0xff) + r), 255)) << L_GREEN_SHIFT |
				(l_int32) L_MAX(0, L_MIN((((uval >> L_BLUE_SHIFT) & 0xff) + r), 255)) << L_BLUE_SHIFT;
				SET_DATA_FOUR_BYTES(lined, j, uval);
			}
		}
	}
	
	pixDestroy(&pixt);
	return pixd;
}

/*!
 *  pixAddNoiseImpulseBipolar()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              d   (density of noise occurrence)
 *      Return: pixd (8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) Return pix with additive bipolar impulse noise.
 *      (2) If colormapped, remove to grayscale.
 */
PIX *
pixAddNoiseImpulseBipolar(PIX       *pixs,
						  l_float32  d)
{
    PROCNAME("pixAddNoiseImpulseBipolar");
	
	d /= 2;
	return pixAddNoiseImpulse(pixs, d, d);
}

/*!
 *  pixAddNoiseImpulseUnipolar()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              noisetype (L_NOISE_SALT or L_NOISE_PEPPER)
 *              d   (density of noise occurrence)
 *      Return: pixd (8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) Return pix with additive unipolar impulse noise.
 *      (2) If colormapped, remove to grayscale.
 */
PIX *
pixAddNoiseImpulseUnipolar(PIX       *pixs,
						   l_int32    noisetype,
						   l_float32  d)
{
    PROCNAME("pixAddNoiseImpulseUnipolar");
	
	if (noisetype != L_NOISE_SALT && noisetype != L_NOISE_PEPPER)
        return (PIX *)ERROR_PTR("invalid noisetype", procName, NULL);
	
	if (noisetype == L_NOISE_SALT)
		return pixAddNoiseImpulse(pixs, d, 0.);
	else
		return pixAddNoiseImpulse(pixs, 0., d);
}

/*!
 *  pixAddNoiseImpulse()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              salt, pepper (probability of salt and pepper occurrence)
 *      Return: pixd (8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) Return pix with additive impulse noise.
 *      (2) If colormapped, remove to grayscale.
 */
PIX *
pixAddNoiseImpulse(PIX       *pixs,
				   l_float32  salt,
				   l_float32  pepper)
{
	l_float64   r;
	l_int32     w, h, d, i, j, val, wpl;
	l_uint32    uval;
	l_uint32   *datas, *datad, *lines, *lined;
	PIX        *pixt, *pixd;
	
    PROCNAME("pixAddNoiseImpulse");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
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
	for (i = 0; i < h; i++) {
		lines = datas + i * wpl;
		lined = datad + i * wpl;
		for (j = 0; j < w; j++) {
			r = randomImpulseBipolar(salt, pepper) * 255;
			if (d == 8) {
				val = GET_DATA_BYTE(lines, j);
				SET_DATA_BYTE(lined, j, L_MAX(0, L_MIN(val + r, 255)));
			}
			else { /* d == 32 */
				uval = GET_DATA_FOUR_BYTES(lines, j);
				uval =
				(l_int32) L_MAX(0, L_MIN((((uval >> L_RED_SHIFT) & 0xff) + r), 255)) << L_RED_SHIFT |
				(l_int32) L_MAX(0, L_MIN((((uval >> L_GREEN_SHIFT) & 0xff) + r), 255)) << L_GREEN_SHIFT |
				(l_int32) L_MAX(0, L_MIN((((uval >> L_BLUE_SHIFT) & 0xff) + r), 255)) << L_BLUE_SHIFT;
				SET_DATA_FOUR_BYTES(lined, j, uval);
			}
		}
	}
	
	pixDestroy(&pixt);
	return pixd;
}
