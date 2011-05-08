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
 *  frequency1.c
 *
 *   Low-level filter creation functions
 *       DPIX          *filterCreateIdeal()
 *       DPIX          *filterCreateButterworth()
 *       DPIX          *filterCreateGaussian()
 *       DPIX          *filterCreateIdealSelective()
 *       DPIX          *filterCreateButterworthSelective()
 *       DPIX          *filterCreateGaussianSelective()
 *
 *   Display filters
 *       PIX           *pixDisplayFilter()
 */

#include "ipl.h"
#include <math.h>

/*--------------------------------------------------------------------*
 *               Low-level filter creation functions                  *
 *--------------------------------------------------------------------*/
/*!
 *  filterCreateIdeal()
 *
 *      Input:  w, h (size of the image to be filtered)
 *              radius (circle radius)
 *              typeflag (L_LO_PASS or L_HI_PASS)
 *      Return: filter, or null on error
 *
 *  Notes:
 *      (1) Implements an ideal lowpass or hipass filter, depending
 *          on the value of @typeflag, where the lowpass version
 *          smoothes the image while the hipass one sharpens it.
 *          If @typeflag == L_LO_PASS the filter passes without
 *          attenuation all frequencies within a circle of radius
 *          @radius from the center of the image and cuts off all
 *          frequencies outside this circle.
 *          Otherwise, if @typeflag == L_HI_PASS the filter behaves
 *          as the opposite, setting to zero all frequencies inside
 *          the circle while passing, without attenuation, all
 *          frequencies outside the circle.
 *      (2) Usually, the filtered image is characterized by ringing
 *          and therefore the ideal filter is not very practical.
 *      (3) The returned filter is of size [w / 2 + 1][h].
 */
DPIX *
filterCreateIdeal(l_int32    w,
				  l_int32    h,
				  l_float32  radius,
				  l_int32    typeflag)
{
	l_int32	       hw, hh, i, j;
	l_float32      dist;
	l_float64     *data, *line;
	l_int32        wpl;
	DPIX          *dpix;
	
	PROCNAME("filterCreateIdeal");
	
	if (typeflag != L_LO_PASS && typeflag != L_HI_PASS)
        return (DPIX *)ERROR_PTR("invalid typeflag", procName, NULL);
	
	if ((dpix = dpixCreate(w / 2 + 1, h)) == NULL)
		return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	/* Create the filter */
	hw = w / 2, hh = h / 2;
	data = dpixGetData(dpix);
	wpl = dpixGetWpl(dpix);
	for (i = 0; i < h; i++) {
		line = data + i * wpl;
		for (j = 0; j < w / 2 + 1; j++) {
			dist = sqrt(pow(j - hw, 2) + pow(i - hh, 2));
			if ((typeflag == L_LO_PASS && dist > radius) ||
				(typeflag == L_HI_PASS && dist <= radius))
				line[j] = 0;
			else
				line[j] = 1;
		}
	}
	
	return dpix;
}

/*!
 *  filterCreateButterworth()
 *
 *      Input:  w, h (size of the image to be filtered)
 *              radius (circle radius)
 *              order (order of the filter function)
 *              typeflag (L_LO_PASS or L_HI_PASS)
 *      Return: filter, or null on error
 *
 *  Notes:
 *      (1) Implements a Butterworth lowpass or hipass filter, depending
 *          on the value of @typeflag, where the lowpass version
 *          smoothes the image while the hipass one sharpens it.
 *          If @typeflag == L_LO_PASS the filter passes all frequencies
 *          within a circle of radius @radius from the center of the
 *          image and cuts off all frequencies outside this circle.
 *          Otherwise, if @typeflag == L_HI_PASS the filter behaves
 *          as the opposite, setting to zero all frequencies inside
 *          the circle while passing all frequencies outside the circle.
 *      (2) Unlike the ideal filter, the Butterworth transfer function
 *          does not have a sharp discontinuity that gives a clear
 *          cutoff between passed and filtered frequencies, which
 *          accounts for no strong ringing effects in the final image.
 *      (3) The higher the @order the more the filter exhibits
 *          characteristics similar to those of the ideal filter.
 *          Usually, Butterworth filters of @order 2 or 4 are a good
 *          compromise between effective lowpass or hipass filtering
 *          and ringing.
 *      (4) The returned filter is of size [w / 2 + 1][h].
 */
DPIX *
filterCreateButterworth(l_int32    w,
						l_int32    h,
						l_float32  radius,
						l_int32    order,
						l_int32    typeflag)
{
	l_int32	       hw, hh, i, j;
	l_float32      dist;
	l_float64     *data, *line;
	l_int32        wpl;
	DPIX          *dpix;
	
	PROCNAME("filterCreateButterworth");
	
	if (typeflag != L_LO_PASS && typeflag != L_HI_PASS)
        return (DPIX *)ERROR_PTR("invalid typeflag", procName, NULL);
	
	if ((dpix = dpixCreate(w / 2 + 1, h)) == NULL)
		return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	/* Create the filter */
	hw = w / 2, hh = h / 2;
	data = dpixGetData(dpix);
	wpl = dpixGetWpl(dpix);
	for (i = 0; i < h; i++) {
		line = data + i * wpl;
		for (j = 0; j < w / 2 + 1; j++) {
			dist = sqrt(pow(j - hw, 2) + pow(i - hh, 2));
			if (typeflag == L_LO_PASS)
				line[j] = 1 / (1 + pow(dist / radius, 2 * order));
			else /* typeflag == L_HI_PASS */
				line[j] = 1 / (1 + pow(radius / dist, 2 * order));
		}
	}
	
	return dpix;
}

/*!
 *  filterCreateGaussian()
 *
 *      Input:  w, h (size of the image to be filtered)
 *              radius (circle radius)
 *              typeflag (L_LO_PASS or L_HI_PASS)
 *      Return: filter, or null on error
 *
 *  Notes:
 *      (1) Implements a Gaussian lowpass or hipass filter, depending
 *          on the value of @typeflag, where the lowpass version
 *          smoothes the image while the hipass one sharpens it.
 *          If @typeflag == L_LO_PASS the filter passes all frequencies
 *          within a circle of radius @radius from the center of the
 *          image and cuts off all frequencies outside this circle.
 *          Otherwise, if @typeflag == L_HI_PASS the filter behaves
 *          as the opposite, setting to zero all frequencies inside
 *          the circle while passing all frequencies outside the circle.
 *      (2) Like the Butterworth transfer function, the Gaussian filter
 *          does not have a sharp discontinuity that gives a clear
 *          cutoff between passed and filtered frequencies.
 *      (3) The results between the Gaussian and Butterworth filters
 *          are quite comparable, but the Gaussian does not exhibit
 *          ringing. This is an important characteristic, especially
 *          in situations (e.g. medical imaging) in which any type of
 *          artifact is unacceptable. In cases where tight control of
 *          the transition between low and high frequencies about the
 *          cutoff frequency are needed, then the Butterworth filter
 *          presents a more suitable choice. The price of this additional
 *          control over the filter profile is the possibility of ringing.
 *      (4) The returned filter is of size [w / 2 + 1][h].
 */
DPIX *
filterCreateGaussian(l_int32    w,
					 l_int32    h,
					 l_float32  radius,
					 l_int32    typeflag)
{
	l_int32	       hw, hh, i, j;
	l_float32      dist;
	l_float64     *data, *line;
	l_int32        wpl;
	DPIX          *dpix;
	
	PROCNAME("filterCreateGaussian");
	
	if (typeflag != L_LO_PASS && typeflag != L_HI_PASS)
        return (DPIX *)ERROR_PTR("invalid typeflag", procName, NULL);
	
	if ((dpix = dpixCreate(w / 2 + 1, h)) == NULL)
		return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	/* Create the filter */
	hw = w / 2, hh = h / 2;
	data = dpixGetData(dpix);
	wpl = dpixGetWpl(dpix);
	for (i = 0; i < h; i++) {
		line = data + i * wpl;
		for (j = 0; j < w / 2 + 1; j++) {
			dist = sqrt(pow(j - hw, 2) + pow(i - hh, 2));
			if (typeflag == L_LO_PASS)
				line[j] = exp(-pow(dist, 2) / (2 * pow(radius, 2)));
			else /* typeflag == L_HI_PASS */
				line[j] = 1 - exp(-pow(dist, 2) / (2 * pow(radius, 2)));
		}
	}
	
	return dpix;
}

/*!
 *  filterCreateIdealSelective()
 *
 *      Input:  w, h (size of the image to be filtered)
 *              radius (circle radius)
 *              bandwidth (width of the band)
 *              typeflag (L_BAND_REJECT or L_BAND_PASS)
 *      Return: filter, or null on error
 *
 *  Notes:
 *      (1) Implements an ideal bandreject or bandpass filter,
 *          depending on the value of @typeflag.
 *      (2) The same concepts and arguments of filterCreateIdeal()
 *          apply, the difference being that the latter operates over the
 *          entire image frequency rectangle while the selective filter
 *          processes only specific bands of frequencies of radius @radius
 *          and width @bandwidth.
 *      (3) The returned filter is of size [w / 2 + 1][h].
 */
DPIX *
filterCreateIdealSelective(l_int32    w,
						   l_int32    h,
						   l_float32  radius,
						   l_float32  bandwidth,
						   l_int32    typeflag)
{
	l_int32	       hw, hh, i, j;
	l_float32      dist, prof;
	l_float64     *data, *line;
	l_int32        wpl;
	DPIX          *dpix;
	
	PROCNAME("filterCreateIdealSelective");
	
	if (typeflag != L_BAND_REJECT && typeflag != L_BAND_PASS)
        return (DPIX *)ERROR_PTR("invalid typeflag", procName, NULL);
	
	if ((dpix = dpixCreate(w / 2 + 1, h)) == NULL)
		return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	/* Create the filter */
	hw = w / 2, hh = h / 2;
	data = dpixGetData(dpix);
	wpl = dpixGetWpl(dpix);
	for (i = 0; i < h; i++) {
		line = data + i * wpl;
		for (j = 0; j < w / 2 + 1; j++) {
			dist = sqrt(pow(j - hw, 2) + pow(i - hh, 2));
			prof = 1;
			if ((dist >= radius - bandwidth / 2) && (dist <= radius + bandwidth / 2))
				prof = 0;
			if (typeflag == L_BAND_PASS)
				prof = 1 - prof;
			line[j] = prof;
		}
	}
	
	return dpix;
}

/*!
 *  filterCreateButterworthSelective()
 *
 *      Input:  w, h (size of the image to be filtered)
 *              radius (circle radius)
 *              bandwidth (width of the band)
 *              order (order of the filter function)
 *              typeflag (L_BAND_REJECT or L_BAND_PASS)
 *      Return: filter, or null on error
 *
 *  Notes:
 *      (1) Implements a Butterworth bandreject or bandpass filter,
 *          depending on the value of @typeflag.
 *      (2) The same concepts and arguments of filterCreateButterworth()
 *          apply, the difference being that the latter operates over the
 *          entire image frequency rectangle while the selective filter
 *          processes only specific bands of frequencies of radius @radius
 *          and width @bandwidth.
 *      (3) The returned filter is of size [w / 2 + 1][h].
 */
DPIX *
filterCreateButterworthSelective(l_int32    w,
								 l_int32    h,
								 l_float32  radius,
								 l_float32  bandwidth,
								 l_int32    order,
								 l_int32    typeflag)
{
	l_int32	       hw, hh, i, j;
	l_float32      dist, prof;
	l_float64     *data, *line;
	l_int32        wpl;
	DPIX          *dpix;
	
	PROCNAME("filterCreateButterworthSelective");
	
	if (typeflag != L_BAND_REJECT && typeflag != L_BAND_PASS)
        return (DPIX *)ERROR_PTR("invalid typeflag", procName, NULL);
	
	if ((dpix = dpixCreate(w / 2 + 1, h)) == NULL)
		return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	/* Create the filter */
	hw = w / 2, hh = h / 2;
	data = dpixGetData(dpix);
	wpl = dpixGetWpl(dpix);
	for (i = 0; i < h; i++) {
		line = data + i * wpl;
		for (j = 0; j < w / 2 + 1; j++) {
			dist = sqrt(pow(j - hw, 2) + pow(i - hh, 2));
			prof = 1 / (1 + pow(dist * bandwidth / (dist * dist - radius * radius), 2 * order));
			if (typeflag == L_BAND_PASS)
				prof = 1 - prof;
			line[j] = prof;
		}
	}
	
	return dpix;
}

/*!
 *  filterCreateGaussianSelective()
 *
 *      Input:  w, h (size of the image to be filtered)
 *              radius (circle radius)
 *              bandwidth (width of the band)
 *              typeflag (L_BAND_REJECT or L_BAND_PASS)
 *      Return: filter, or null on error
 *
 *  Notes:
 *      (1) Implements a Gaussian bandreject or bandpass filter,
 *          depending on the value of @typeflag.
 *      (2) The same concepts and arguments of filterCreateGaussian()
 *          apply, the difference being that the latter operates over the
 *          entire image frequency rectangle while the selective filter
 *          processes only specific bands of frequencies of radius @radius
 *          and width @bandwidth.
 *      (3) The returned filter is of size [w / 2 + 1][h].
 */
DPIX *
filterCreateGaussianSelective(l_int32    w,
							  l_int32    h,
							  l_float32  radius,
							  l_float32  bandwidth,
							  l_int32    typeflag)
{
	l_int32	       hw, hh, i, j;
	l_float32      dist, prof;
	l_float64     *data, *line;
	l_int32        wpl;
	DPIX          *dpix;
	
	PROCNAME("filterCreateGaussianSelective");
	
	if (typeflag != L_BAND_REJECT && typeflag != L_BAND_PASS)
        return (DPIX *)ERROR_PTR("invalid typeflag", procName, NULL);
	
	if ((dpix = dpixCreate(w / 2 + 1, h)) == NULL)
		return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	/* Create the filter */
	hw = w / 2, hh = h / 2;
	data = dpixGetData(dpix);
	wpl = dpixGetWpl(dpix);
	for (i = 0; i < h; i++) {
		line = data + i * wpl;
		for (j = 0; j < w / 2 + 1; j++) {
			dist = sqrt(pow(j - hw, 2) + pow(i - hh, 2));
			prof = 1 - exp(-pow((dist * dist - radius * radius) / (dist * bandwidth), 2));
			if (typeflag == L_BAND_PASS)
				prof = 1 - prof;
			line[j] = prof;
		}
	}
	
	return dpix;
}

/*--------------------------------------------------------------------*
 *                          Display filters                           *
 *--------------------------------------------------------------------*/
/*!
 *  pixDisplayFilter()
 *
 *      Input:  dpix (filter)
 *      Return: pixd (8 bpp), or null on error
 *
 *  Notes:
 *      (1) Assumes that the filter is half size and symmetric.
 */
PIX *
pixDisplayFilter(DPIX    *dpix)
{
	l_int32	 w, h;
	PIX     *pixf, *pixd;
	DPIX    *dpixf;
	
	PROCNAME("pixDisplayFilter");
	
	if (!dpix)
        return (PIX *)ERROR_PTR("dpix not defined", procName, NULL);
	dpixGetDimensions(dpix, &w, &h);
	if ((pixd = pixCreate((w - 1) * 2, h, 8)) == NULL)
		return (PIX *)ERROR_PTR("pixd not made", procName, NULL);
	
	dpixf = dpixMaxDynamicRange(NULL, dpix, 255, L_HI_TO_WHITE);
	pixf = dpixConvertToPix(dpixf, 8, L_CLIP_TO_ZERO, 0, 0);
	
	/* Copy the left half of filter to pixd */
	pixRasterop(pixd, 0, 0, w, h, PIX_SRC, pixf, 0, 0);
	
	/* Copy the right half of filter to pixd */
	pixFlipLR(pixf, pixf);
	pixRasterop(pixd, w, 0, w, h, PIX_SRC, pixf, 0, 0);
	
	dpixDestroy(&dpixf);
	pixDestroy(&pixf);
	
	return pixd;
}