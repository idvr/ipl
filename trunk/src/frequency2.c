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
 *  frequency2.c
 *
 *   Top-level lowpass and hipass standard filters
 *       PIX           *pixIdealFilter()
 *       PIX           *pixButterworthFilter()
 *       PIX           *pixGaussianFilter()
 *
 *   Top-level bandreject and bandpass selective filters
 *       PIX           *pixIdealSelectiveFilter()
 *       PIX           *pixButterworthSelectiveFilter()
 *       PIX           *pixGaussianSelectiveFilter()
 *
 *   Filter application
 *       PIX           *pixApplyFilter()
 *       PIX           *pixApplyFilterGray()
 *
 */

#include "ipl.h"

#ifdef HAVE_CONFIG_H
#include "config_auto.h"
#endif  /* HAVE_CONFIG_H */

/* --------------------------------------------*/
#if  HAVE_LIBFFTW3
/* --------------------------------------------*/

#include <math.h>
#include <complex.h>
#include <fftw3.h>

fftw_complex *
pixDFT(PIX     *pixs,
	   l_int32  shiftflag);

PIX *
pixInverseDFT(fftw_complex *dft,
			  l_int32       w,
			  l_int32       h,
			  l_int32       shiftflag,
			  l_int32       outflag);

/*--------------------------------------------------------------------*
 *          Top-level lowpass and hipass standard filters             *
 *--------------------------------------------------------------------*/
/*!
 *  pixIdealFilter()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              radius (circle radius)
 *              typeflag (L_LO_PASS or L_HI_PASS)
 *              outflag (L_CLIP_TO_ZERO, L_TAKE_ABSVAL,
 *                       L_THRESH_NEG_TO_BLACK or L_THRESH_NEG_TO_WHITE)
 *              dpixf (<optional return> filter used)
 *      Return: filtered pix (8 or 32 bpp), or null on error
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
 */
PIX *
pixIdealFilter(PIX        *pixs,
			   l_float32   radius,
			   l_int32     typeflag,
			   l_int32     outflag,
			   DPIX      **dpixf)
{
	DPIX    *dpix;
	PIX     *pixd;
	
	PROCNAME("pixIdealFilter");
	
	if (!pixs)
		return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
	dpix = filterCreateIdeal(pixGetWidth(pixs), pixGetHeight(pixs), radius, typeflag);
	if (dpix == NULL)
		return (PIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	pixd = pixApplyFilter(pixs, dpix, outflag);
	
	if (dpixf)
		*dpixf = dpix;
	else
		dpixDestroy(&dpix);
	
	return pixd;
}

/*!
 *  pixButterworthFilter()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              radius (circle radius)
 *              order (order of the filter function)
 *              typeflag (L_LO_PASS or L_HI_PASS)
 *              outflag (L_CLIP_TO_ZERO, L_TAKE_ABSVAL,
 *                       L_THRESH_NEG_TO_BLACK or L_THRESH_NEG_TO_WHITE)
 *              dpixf (<optional return> filter used)
 *      Return: filtered pix (8 or 32 bpp), or null on error
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
 */
PIX *
pixButterworthFilter(PIX        *pixs,
					 l_float32   radius,
					 l_int32     order,
					 l_int32     typeflag,
					 l_int32     outflag,
					 DPIX      **dpixf)
{
	DPIX    *dpix;
	PIX     *pixd;
	
	PROCNAME("pixButterworthFilter");
		
	if (!pixs)
		return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
	dpix = filterCreateButterworth(pixGetWidth(pixs), pixGetHeight(pixs), radius, order, typeflag);
	if (dpix == NULL)
		return (PIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	pixd = pixApplyFilter(pixs, dpix, outflag);
	
	if (dpixf)
		*dpixf = dpix;
	else
		dpixDestroy(&dpix);

	return pixd;
}

/*!
 *  pixGaussianFilter()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              radius (circle radius)
 *              typeflag (L_LO_PASS or L_HI_PASS)
 *              outflag (L_CLIP_TO_ZERO, L_TAKE_ABSVAL,
 *                       L_THRESH_NEG_TO_BLACK or L_THRESH_NEG_TO_WHITE)
 *              dpixf (<optional return> filter used)
 *      Return: filtered pix, or null on error
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
 */
PIX *
pixGaussianFilter(PIX        *pixs,
				  l_float32   radius,
				  l_int32     typeflag,
				  l_int32     outflag,
				  DPIX      **dpixf)
{
	DPIX    *dpix;
	PIX     *pixd;
	
	PROCNAME("pixGaussianFilter");
	
	if (!pixs)
		return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
	dpix = filterCreateGaussian(pixGetWidth(pixs), pixGetHeight(pixs), radius, typeflag);
	if (dpix == NULL)
		return (PIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	pixd = pixApplyFilter(pixs, dpix, outflag);
	
	if (dpixf)
		*dpixf = dpix;
	else
		dpixDestroy(&dpix);
	
	return pixd;
}

/*--------------------------------------------------------------------*
 *        Top-level bandreject and bandpass selective filters         *
 *--------------------------------------------------------------------*/
/*!
 *  pixIdealSelectiveFilter()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              radius (circle radius)
 *              bandwidth (width of the band)
 *              typeflag (L_BAND_REJECT or L_BAND_PASS)
 *              outflag (L_CLIP_TO_ZERO, L_TAKE_ABSVAL,
 *                       L_THRESH_NEG_TO_BLACK or L_THRESH_NEG_TO_WHITE)
 *              dpixf (<optional return> filter used)
 *      Return: filtered pix, or null on error
 *
 *  Notes:
 *      (1) Implements an ideal bandreject or bandpass filter, depending
 *          on the value of @typeflag.
 *      (2) The same concepts and arguments of pixIdealFilter() apply,
 *          the difference being that the latter operates over the entire
 *          frequency rectangle while the selective filter processes only
 *          specific bands of frequencies of radius @radius and width
 *          @bandwidth.
 */
PIX *
pixIdealSelectiveFilter(PIX        *pixs,
						l_float32   radius,
						l_float32   bandwidth,
						l_int32     typeflag,
						l_int32     outflag,
						DPIX      **dpixf)
{
	DPIX    *dpix;
	PIX     *pixd;
	
	PROCNAME("pixIdealSelectiveFilter");
	
	if (!pixs)
		return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
	dpix = filterCreateIdealSelective(pixGetWidth(pixs), pixGetHeight(pixs), radius, bandwidth, typeflag);
	if (dpix == NULL)
		return (PIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	pixd = pixApplyFilter(pixs, dpix, outflag);
	
	if (dpixf)
		*dpixf = dpix;
	else
		dpixDestroy(&dpix);
	
	return pixd;
}

/*!
 *  pixButterworthSelectiveFilter()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              radius (circle radius)
 *              bandwidth (width of the band)
 *              order (order of the filter function)
 *              typeflag (L_BAND_REJECT or L_BAND_PASS)
 *              outflag (L_CLIP_TO_ZERO, L_TAKE_ABSVAL,
 *                       L_THRESH_NEG_TO_BLACK or L_THRESH_NEG_TO_WHITE)
 *              dpixf (<optional return> filter used)
 *      Return: filtered pix, or null on error
 *
 *  Notes:
 *      (1) Implements a Butterworth bandreject or bandpass filter,
 *          depending on the value of @typeflag.
 *      (2) The same concepts and arguments of pixButterworthFilter() apply,
 *          the difference being that the latter operates over the entire
 *          frequency rectangle while the selective filter processes only
 *          specific bands of frequencies of radius @radius and width
 *          @bandwidth.
 */
PIX *
pixButterworthSelectiveFilter(PIX        *pixs,
							  l_float32   radius,
							  l_float32   bandwidth,
							  l_int32     order,
							  l_int32     typeflag,
							  l_int32     outflag,
							  DPIX      **dpixf)
{
	DPIX    *dpix;
	PIX     *pixd;
	
	PROCNAME("pixButterworthSelectiveFilter");
	
	if (!pixs)
		return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
	dpix = filterCreateButterworthSelective(pixGetWidth(pixs), pixGetHeight(pixs), radius, bandwidth, order, typeflag);
	if (dpix == NULL)
		return (PIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	pixd = pixApplyFilter(pixs, dpix, outflag);
	
	if (dpixf)
		*dpixf = dpix;
	else
		dpixDestroy(&dpix);
	
	return pixd;
}

/*!
 *  pixGaussianSelectiveFilter()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              radius (circle radius)
 *              bandwidth (width of the band)
 *              typeflag (L_BAND_REJECT or L_BAND_PASS)
 *              outflag (L_CLIP_TO_ZERO, L_TAKE_ABSVAL,
 *                       L_THRESH_NEG_TO_BLACK or L_THRESH_NEG_TO_WHITE)
 *              dpixf (<optional return> filter used)
 *      Return: filtered pix, or null on error
 *
 *  Notes:
 *      (1) Implements a Gaussian bandreject or bandpass filter,
 *          depending on the value of @typeflag.
 *      (2) The same concepts and arguments of pixGaussianFilter() apply,
 *          the difference being that the latter operates over the entire
 *          frequency rectangle while the selective filter processes only
 *          specific bands of frequencies of radius @radius and width
 *          @bandwidth.
 */
PIX *
pixGaussianSelectiveFilter(PIX        *pixs,
						   l_float32   radius,
						   l_float32   bandwidth,
						   l_int32     typeflag,
						   l_int32     outflag,
						   DPIX      **dpixf)
{
	DPIX    *dpix;
	PIX     *pixd;
	
	PROCNAME("pixGaussianSelectiveFilter");
	
	if (!pixs)
		return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	
	dpix = filterCreateGaussianSelective(pixGetWidth(pixs), pixGetHeight(pixs), radius, bandwidth, typeflag);
	if (dpix == NULL)
		return (PIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	pixd = pixApplyFilter(pixs, dpix, outflag);
	
	if (dpixf)
		*dpixf = dpix;
	else
		dpixDestroy(&dpix);
	
	return pixd;
}

/*--------------------------------------------------------------------*
 *                         Filter application                         *
 *--------------------------------------------------------------------*/
/*!
 *  pixApplyFilter()
 *
 *      Input:  pix (8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              dpix (filter)
 *              outflag (L_CLIP_TO_ZERO, L_TAKE_ABSVAL,
 *                       L_THRESH_NEG_TO_BLACK or L_THRESH_NEG_TO_WHITE)
 *      Return: pixd, or null on error
 *
 *  Notes:
 *      (1) Apply the input filter to pix by multiplying it for the
 *          Fourier transform of pix. The inverse Fourier transform
 *          is then used on the result to return the filtered image.
 *      (2) If pix is 32 bpp RGB, the filter is applied to each color
 *          channel separately.
 *      (3) If colormapped, remove to grayscale.
 */
PIX *
pixApplyFilter(PIX     *pixs,
			   DPIX    *dpix,
			   l_int32  outflag)
{
	l_int32	 w, h, d;
	PIX     *pixt, *pixd, *pixr, *pixrc, *pixg, *pixgc, *pixb, *pixbc;
	
	PROCNAME("pixApplyFilter");
	
    if (!pixs && !dpix)
        return (PIX *)ERROR_PTR("pixs or dpix not defined", procName, NULL);
	
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
	
	if (d == 8) {
		pixd = pixApplyFilterGray(pixt, dpix, outflag);
	}
	else { /* d == 32 */
		pixr = pixGetRGBComponent(pixt, COLOR_RED);
        pixrc = pixApplyFilterGray(pixr, dpix, outflag);
        pixDestroy(&pixr);
        pixg = pixGetRGBComponent(pixt, COLOR_GREEN);
        pixgc = pixApplyFilterGray(pixg, dpix, outflag);
        pixDestroy(&pixg);
        pixb = pixGetRGBComponent(pixt, COLOR_BLUE);
        pixbc = pixApplyFilterGray(pixb, dpix, outflag);
        pixDestroy(&pixb);
        pixd = pixCreateRGBImage(pixrc, pixgc, pixbc);
        pixDestroy(&pixrc);
        pixDestroy(&pixgc);
        pixDestroy(&pixbc);
	}
	
	pixDestroy(&pixt);
	return pixd;
}

/*!
 *  pixApplyFilterGray()
 *
 *      Input:  pix  (8 bpp)
 *              dpix (filter)
 *              outflag (L_CLIP_TO_ZERO, L_TAKE_ABSVAL,
 *                       L_THRESH_NEG_TO_BLACK or L_THRESH_NEG_TO_WHITE)
 *      Return: pixd (8 bpp), or null on error
 *
 *  Notes:
 *      (1) Apply the input filter to pix by multiplying it for the
 *          Fourier transform of pix. The inverse Fourier transform
 *          is then used on the result to return the filtered image.
 */
PIX *
pixApplyFilterGray(PIX     *pixs,
				   DPIX    *dpix,
				   l_int32  outflag)
{
	l_int32	       w, h, d;
	l_int32        i, j, k;
	fftw_complex  *output;
	l_float64     *data, *line;
	l_int32        wpl;
	PIX           *pixd;
	
	PROCNAME("pixApplyFilterGray");
	
	if (!pixs && !dpix)
        return (PIX *)ERROR_PTR("pixs or dpix not defined", procName, NULL);
	pixGetDimensions(pixs, &w, &h, &d);
	if (dpixGetWidth(dpix) > w / 2 + 1 || dpixGetHeight(dpix) > h)
		return (PIX *)ERROR_PTR("dpix is smaller than pix", procName, NULL);
	if (d != 8)
        return (PIX *)ERROR_PTR("pixs not 8 bpp", procName, NULL);
	if (outflag != L_CLIP_TO_ZERO && outflag != L_TAKE_ABSVAL &&
		outflag != L_THRESH_NEG_TO_BLACK && outflag != L_THRESH_NEG_TO_WHITE)
        return (PIX *)ERROR_PTR("invalid outflag", procName, NULL);
	
	/* Calculate the DFT of pixs */
	if ((output = pixDFT(pixs, L_WITH_SHIFTING)) == NULL)
        return (PIX *)ERROR_PTR("pixs DFT not computed", procName, NULL);
	
	/* Filter the DFT results */
	data = dpixGetData(dpix);
	wpl = dpixGetWpl(dpix);
	for (i = 0, k = 0; i < h; i++) {
		line = data + i * wpl;
		for (j = 0; j < w / 2 + 1; j++, k++) {
			output[k] *= line[j];
		}
	}
	
	/* Compute the inverse of the DFT */
	pixd = pixInverseDFT(output, w, h, L_WITH_SHIFTING, outflag);
	
	/* Release the allocated resources */
	fftw_free(output);
	
	return pixd;
}

/* --------------------------------------------*/
#endif  /* HAVE_LIBFFTW3 */
/* --------------------------------------------*/