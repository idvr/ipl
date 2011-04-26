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
 *  spectra.c
 *
 *   Top-level functions
 *       DPIX          *pixGetFourierSpectrumGray()
 *       DPIX          *pixGetPowerSpectrumGray()
 *       DPIX          *pixGetPhaseAngleGray()
 *       DPIX          *pixGetFourierSpectrumRGB()
 *       DPIX          *pixGetPowerSpectrumRGB()
 *       DPIX          *pixGetPhaseAngleRGB()
 *
 *   Spectra computation
 *       l_int32        pixGetDFTSpectraGray()
 *       l_int32        pixGetDFTSpectraRGB()
 *
 *   Phase correlation
 *       l_int32        pixPhaseCorrelation()
 *
 *   Display utilities
 *       PIX           *pixDisplaySpectrum()
 *       PIX           *pixConvertSpectrumToPix()
 *
 *   Pix reconstruction from spectra
 *       PIX           *pixRestoreFromSpectraGray()
 *       PIX           *pixRestoreFromSpectraRGB() 
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

DPIX *
dpixInverseDFT(fftw_complex *dft,
			   l_int32       w,
			   l_int32       h);

/*--------------------------------------------------------------------*
 *                       Top-level functions                          *
 *--------------------------------------------------------------------*/
/*!
 *  pixGetFourierSpectrumGray()
 *
 *      Input:  pix (1 or 8 bpp; or 2, 4 or 8 bpp with colormap)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: dpix, or null on error
 *
 *  Notes:
 *      (1) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the spectra to the center of pix.
 */
DPIX *
pixGetFourierSpectrumGray(PIX      *pixs,
						  l_int32   shiftflag)
{
	DPIX *dpix;
	
	PROCNAME("pixGetFourierSpectrumGray");
	
	if (pixGetDFTSpectraGray(pixs, &dpix, NULL, NULL, shiftflag) != 0)
		return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	return dpix;
}

/*!
 *  pixGetPowerSpectrumGray()
 *
 *      Input:  pix (1 or 8 bpp; or 2, 4 or 8 bpp with colormap)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: dpix, or null on error
 *
 *  Notes:
 *      (1) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the spectra to the center of pix.
 */
DPIX *
pixGetPowerSpectrumGray(PIX      *pixs,
						l_int32   shiftflag)
{
	DPIX *dpix;
	
	PROCNAME("pixGetPowerSpectrumGray");
	
	if (pixGetDFTSpectraGray(pixs, NULL, NULL, &dpix, shiftflag) != 0)
		return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	return dpix;
}

/*!
 *  pixGetPhaseAngleGray()
 *
 *      Input:  pix (1 or 8 bpp; or 2, 4 or 8 bpp with colormap)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: dpix, or null on error
 *
 *  Notes:
 *      (1) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the spectra to the center of pix.
 */
DPIX *
pixGetPhaseAngleGray(PIX      *pixs,
					 l_int32   shiftflag)
{
	DPIX *dpix;
	
	PROCNAME("pixGetPhaseAngleGray");
	
	if (pixGetDFTSpectraGray(pixs, NULL, &dpix, NULL, shiftflag) != 0)
		return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	return dpix;
}

/*!
 *  pixGetFourierSpectrumRGB()
 *
 *      Input:  pix (32 bpp)
 *              color  (one of {COLOR_RED, COLOR_GREEN, COLOR_BLUE})
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: dpix, or null on error
 *
 *  Notes:
 *      (1) Calculate the Fourier spectrum of the RGB component of pix
 *          specified by @color.
 *      (2) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the spectra to the center of pix.
 */
DPIX *
pixGetFourierSpectrumRGB(PIX      *pixs,
						 l_int32   color,
						 l_int32   shiftflag)
{
	DPIX *dpix;
	
	PROCNAME("pixGetFourierSpectrumRGB");
	
	if (pixGetDFTSpectraRGB(pixs, color, &dpix, NULL, NULL, shiftflag) != 0)
		return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	return dpix;
}

/*!
 *  pixGetPowerSpectrumRGB()
 *
 *      Input:  pix (32 bpp)
 *              color  (one of {COLOR_RED, COLOR_GREEN, COLOR_BLUE})
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: dpix, or null on error
 *
 *  Notes:
 *      (1) Calculate the power spectrum of the RGB component of pix
 *          specified by @color.
 *      (2) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the spectra to the center of pix.
 */
DPIX *
pixGetPowerSpectrumRGB(PIX      *pixs,
					   l_int32   color,
					   l_int32   shiftflag)
{
	DPIX *dpix;
	
	PROCNAME("pixGetPowerSpectrumRGB");
	
	if (pixGetDFTSpectraRGB(pixs, color, NULL, NULL, &dpix, shiftflag) != 0)
		return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	return dpix;
}

/*!
 *  pixGetPhaseAngleRGB()
 *
 *      Input:  pix (32 bpp)
 *              color  (one of {COLOR_RED, COLOR_GREEN, COLOR_BLUE})
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: dpix, or null on error
 *
 *  Notes:
 *      (1) Calculate the phase angle of the RGB component of pix
 *          specified by @color.
 *      (2) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the spectra to the center of pix.
 */
DPIX *
pixGetPhaseAngleRGB(PIX      *pixs,
					l_int32   color,
					l_int32   shiftflag)
{
	DPIX *dpix;
	
	PROCNAME("pixGetPhaseAngleRGB");
	
	if (pixGetDFTSpectraRGB(pixs, color, NULL, &dpix, NULL, shiftflag) != 0)
		return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	return dpix;
}

/*--------------------------------------------------------------------*
 *                       Spectra computation                          *
 *--------------------------------------------------------------------*/
/*!
 *  pixGetDFTSpectraGray()
 *
 *      Input:  pix (1 or 8 bpp; or 2, 4 or 8 bpp with colormap)
 *              &spectrum (<optional return> Fourier spectrum)
 *              &phaseAngle (<optional return> phase angle)
 *              &powerSpectrum (<optional return> power spectrum)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: 0 if OK; 1 on error
 *
 *  Notes:
 *      (1) Calculate the Fourier spectrum, phase angle and power
 *          spectrum of pix.
 *      (2) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the spectra to the center of pix.
 */
l_int32
pixGetDFTSpectraGray(PIX      *pixs,
					 DPIX    **spectrum,
					 DPIX    **phaseAngle,
					 DPIX    **powerSpectrum,
					 l_int32   shiftflag)
{
	l_int32	       w, h, d, wpl;
	l_int32        i, j, k;
	l_float64     *datas, *dataa, *datap;
	l_float64     *lines, *linea, *linep, *linehs, *lineha, *linehp;
	l_float64      power;
	fftw_complex  *output;
	PIX           *pixt;
	
	PROCNAME("pixGetDFTSpectraGray");
	
    if (!pixs)
        return ERROR_INT("pixs not defined", procName, 1);
	if (!spectrum && !phaseAngle && !powerSpectrum)
        return ERROR_INT("nothing to do", procName, 1);
	if (shiftflag != L_NO_SHIFTING && shiftflag != L_WITH_SHIFTING)
        return ERROR_INT("invalid shiftflag", procName, 1);
	
	/* Remove colormap if necessary */
	pixGetDimensions(pixs, &w, &h, &d);
    if ((d == 2 || d == 4 || d == 8) && pixGetColormap(pixs)) {
        L_WARNING("pix has colormap; removing", procName);
        pixt = pixRemoveColormap(pixs, REMOVE_CMAP_BASED_ON_SRC);
        d = pixGetDepth(pixt);
    }
    else
        pixt = pixClone(pixs);
	
    if (d != 1 && d != 8) {
        pixDestroy(&pixt);
        return ERROR_INT("depth not 1 or 8 bpp", procName, 1);
    }
	
	/* Calculate the DFT of pixs */
	if ((output = pixDFT(pixt, shiftflag)) == NULL)
        return ERROR_INT("pixt DFT not computed", procName, 1);
	
	if (spectrum) {
		*spectrum = dpixCreate(w, h);
		datas = dpixGetData(*spectrum);
		wpl = dpixGetWpl(*spectrum);
	}
	if (phaseAngle) {
		*phaseAngle = dpixCreate(w, h);
		dataa = dpixGetData(*phaseAngle);
		wpl = dpixGetWpl(*phaseAngle);
	}
	if (powerSpectrum) {
		*powerSpectrum = dpixCreate(w, h);
		datap = dpixGetData(*powerSpectrum);
		wpl = dpixGetWpl(*powerSpectrum);
	}
	
	/* Calculate the frequency spectrum, phase angle and power spectrum */
	for (i = 0, k = 0; i < h; i++) {
		if (spectrum) lines = datas + i * wpl;
		if (phaseAngle) linea = dataa + i * wpl;
		if (powerSpectrum) linep = datap + i * wpl;
		for (j = 0; j < w / 2 + 1; j++, k++) {
			if (spectrum || powerSpectrum)
				power = pow(creal(output[k]), 2) + pow(cimag(output[k]), 2);
			if (spectrum)
				lines[j] = sqrt(power);
			if (phaseAngle)
				linea[j] = atan2(cimag(output[k]), creal(output[k]));
			if (powerSpectrum)
				linep[j] = power;
		}
	}
	/* Compute the second half of the spectra */
	for (i = 0, k = 0; i < h; i++) {
		if (spectrum) lines = datas + i * wpl, linehs = datas + ((h - i) * wpl);
		if (phaseAngle) linea = dataa + i * wpl, lineha = dataa + ((h - i) * wpl);
		if (powerSpectrum) linep = datap + i * wpl, linehp = datap + ((h - i) * wpl);
		for (j = 1; j < w / 2 + (w % 2); j++, k++) {
			if (spectrum)
				lines[w - j] = i > 0 ? linehs[j] : lines[j];
			if (phaseAngle)
				linea[w - j] = i > 0 ? -lineha[j] : -linea[j];
			if (powerSpectrum)
				linep[w - j] = i > 0 ? linehp[j] : linep[j];
		}
	}
	
	/* Release the allocated resources */
	fftw_free(output);
	pixDestroy(&pixt);
	
	return 0;
}

/*!
 *  pixGetDFTSpectraRGB()
 *
 *      Input:  pix (32 bpp)
 *              color  (one of {COLOR_RED, COLOR_GREEN, COLOR_BLUE})
 *              &spectrum (<optional return> Fourier spectrum)
 *              &phaseAngle (<optional return> phase angle)
 *              &powerSpectrum (<optional return> power spectrum)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: 0 if OK; 1 on error
 *
 *  Notes:
 *      (1) Calculate the Fourier spectrum, phase angle and power
 *          spectrum of the RGB component of pix specified by @color.
 *      (2) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the spectra to the center of pix.
 */
l_int32
pixGetDFTSpectraRGB(PIX      *pixs,
					l_int32   color,
					DPIX    **spectrum,
					DPIX    **phaseAngle,
					DPIX    **powerSpectrum,
					l_int32   shiftflag)
{
	PIX     *pixc;
	
	PROCNAME("pixGetDFTSpectraRGB");
	
    if (pixGetDepth(pixs) != 32)
        return ERROR_INT("pixs not 32 bpp", procName, 1);
    if (color != COLOR_RED && color != COLOR_GREEN && color != COLOR_BLUE)
        return ERROR_INT("invalid color", procName, 1);
	
	/* Extract the selected color channel and get its spectra */
	pixc = pixGetRGBComponent(pixs, color);
	if (pixGetDFTSpectraGray(pixc, spectrum, phaseAngle, powerSpectrum, shiftflag) != 0) {
		pixDestroy(&pixc);
		return ERROR_INT("pixc not made", procName, 1);
	}
	
	pixDestroy(&pixc);
	return 0;
}

/*--------------------------------------------------------------------*
 *                         Phase correlation                          *
 *--------------------------------------------------------------------*/
/*!
 *  pixPhaseCorrelation()
 *
 *      Input:  pixr, pixs (1, 8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              &peak (<optional return> phase correlation peak)
 *              &xloc (<optional return> x location of the peak)
 *              &yloc (<optional return> y location of the peak)
 *      Return: 0 if OK; 1 on error
 *
 *  Notes:
 *      (1) Phase correlation is a method of image registration,
 *          and uses a fast frequency-domain approach to estimate the
 *          relative translative offset between two similar images.
 *      (2) The reference and input images must have same width
 *          and height.
 *      (3) Determine the location of the highest value (peak) in @r,
 *          defined as the phase cross-correlated 2D array:
 *          r = inverse Fourier (F * G / |F * G|)
 *          where F is the Fourier transform of the reference image and
 *          G is the complex conjugate of the transform of the input.
 *          The peak in @r corresponds to the object translation movement
 *          from the reference to the input images. A value of 1.0  at
 *          location (0, 0) means that the images are identical.
 *      (4) If colormapped, remove to grayscale.
 */
l_int32
pixPhaseCorrelation(PIX       *pixr,
					PIX       *pixs,
					l_float64 *ppeak,
					l_int32   *pxloc,
					l_int32   *pyloc)
{
	l_int32	       w, h, d;
	l_int32        i, j, k;
	l_float64      cr, ci, r;
	fftw_complex  *outputr, *outputs, *outputd;
	PIX           *pixt1, *pixt2;
	DPIX          *dpix;
	
	PROCNAME("pixPhaseCorrelation");
	
    if (!pixr && !pixs)
        return ERROR_INT("pixr or pixs not defined", procName, 1);
	if (!pixSizesEqual(pixr, pixs))
        return ERROR_INT("pixr and pixs unequal size", procName, 1);
	if (!ppeak && !pxloc && !pyloc)
        return ERROR_INT("nothing to do", procName, 1);
	
	/* Process the reference image */
	if (pixGetColormap(pixr))
        pixt1 = pixRemoveColormap(pixr, REMOVE_CMAP_TO_GRAYSCALE);
    else if (pixGetDepth(pixr) == 32)
        pixt1 = pixConvertRGBToLuminance(pixr);
    else
        pixt1 = pixClone(pixr);
	
	/* Process the input image */
	if (pixGetColormap(pixs))
        pixt2 = pixRemoveColormap(pixs, REMOVE_CMAP_TO_GRAYSCALE);
    else if (pixGetDepth(pixs) == 32)
        pixt2 = pixConvertRGBToLuminance(pixs);
    else
        pixt2 = pixClone(pixs);
    
	/* Calculate the DFT of pixr and pixs */
	if ((outputr = pixDFT(pixt1, L_NO_SHIFTING)) == NULL)
		return ERROR_INT("outputr not made", procName, 1);
	if ((outputs = pixDFT(pixt2, L_NO_SHIFTING)) == NULL) {
		fftw_free(outputr);
		return ERROR_INT("outputs not made", procName, 1);
	}
	pixGetDimensions(pixr, &w, &h, &d);
	outputd = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * h * (w / 2 + 1));
	if (outputd == NULL) {
		fftw_free(outputr);
		fftw_free(outputs);
		return ERROR_INT("outputd not made", procName, 1);
	}
	
	/* Calculate the cross-power spectrum */
	for (i = 0, k = 0; i < h; i++) {
		for (j = 0; j < w / 2 + 1; j++, k++) {
			cr = creal(outputs[k]) * creal(outputr[k]) - cimag(outputs[k]) * (-cimag(outputr[k]));
			ci = creal(outputs[k]) * (-cimag(outputr[k])) + cimag(outputs[k]) * creal(outputr[k]);
			r = sqrt(pow(cr, 2.) + pow(ci, 2.));
			outputd[k] = (cr / r) + I * (ci / r);
		}
	}
	
	/* Compute the inverse DFT of the cross-power spectrum
	    and find its peak */
	dpix = dpixInverseDFT(outputd, w, h);
	dpixGetMax(dpix, ppeak, pxloc, pyloc);
	
	if (*pxloc >= w / 2)
		*pxloc -= w;
	if (*pyloc >= h / 2)
		*pyloc -= h;
		
	/* Release the allocated resources */
	fftw_free(outputr);
	fftw_free(outputs);
	fftw_free(outputd);
	dpixDestroy(&dpix);
	pixDestroy(&pixt1);
	pixDestroy(&pixt2);
}

/*--------------------------------------------------------------------*
 *                         Display utilities                          *
 *--------------------------------------------------------------------*/
/*!
 *  pixDisplaySpectrum()
 *
 *      Input:  pix (1, 8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              scale (to expand or compress the intensity of the spectrum)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *              hivals (L_HI_TO_WHITE or L_HI_TO_BLACK)
 *      Return: pixd (8 or 32 bpp), or null on error
 *
 *  Notes:
 *      (1) Typically @scale is set to 1.0.
 *      (2) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the spectra to the center of pix.
 *      (3) Set @hivals to L_HI_TO_WHITE to convert higher intensities to
 *          white, or L_HI_TO_BLACK to convert them to black instead.
 */
PIX *
pixDisplaySpectrum(PIX       *pixs,
				   l_float32  scale,
				   l_int32    shiftflag,
				   l_int32    hivals)
{
	l_int32  w, h, d;
	DPIX    *dpix;
	PIX     *pixt, *pixd, *pixr, *pixg, *pixb;
	
	PROCNAME("pixDisplaySpectrum");
	
	if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	if (shiftflag != L_NO_SHIFTING && shiftflag != L_WITH_SHIFTING)
        return (PIX *)ERROR_PTR("invalid shiftflag", procName, NULL);
	if (hivals != L_HI_TO_WHITE && hivals != L_HI_TO_BLACK)
        return (PIX *)ERROR_PTR("invalid hivals", procName, NULL);
	
	/* Remove colormap if necessary */
	pixGetDimensions(pixs, &w, &h, &d);
    if ((d == 2 || d == 4 || d == 8) && pixGetColormap(pixs)) {
        L_WARNING("pix has colormap; removing", procName);
        pixt = pixRemoveColormap(pixs, REMOVE_CMAP_BASED_ON_SRC);
        d = pixGetDepth(pixt);
    }
    else
        pixt = pixClone(pixs);
	
    if (d != 1 && d != 8 && d != 32) {
        pixDestroy(&pixt);
        return (PIX *)ERROR_PTR("depth not 1, 8 or 32 bpp", procName, NULL);
    }
	
	if (d == 1 || d == 8) {
		dpix = pixGetFourierSpectrumGray(pixt, shiftflag);
		pixd = pixConvertSpectrumToPix(dpix, scale, shiftflag, hivals);
		dpixDestroy(&dpix);
	}
	else { /* d == 32 */
		dpix = pixGetFourierSpectrumRGB(pixt, COLOR_RED, shiftflag);
		pixr = pixConvertSpectrumToPix(dpix, scale, shiftflag, hivals);
		dpixDestroy(&dpix);
		dpix = pixGetFourierSpectrumRGB(pixt, COLOR_GREEN, shiftflag);
		pixg = pixConvertSpectrumToPix(dpix, scale, shiftflag, hivals);
		dpixDestroy(&dpix);
		dpix = pixGetFourierSpectrumRGB(pixt, COLOR_BLUE, shiftflag);
		pixb = pixConvertSpectrumToPix(dpix, scale, shiftflag, hivals);
		dpixDestroy(&dpix);
		pixd = pixCreateRGBImage(pixr, pixg, pixb);
        pixDestroy(&pixr);
        pixDestroy(&pixg);
        pixDestroy(&pixb);
	}
	
	pixDestroy(&pixt);
	return pixd;
}

/*!
 *  pixConvertSpectrumToPix()
 *
 *      Input:  dpix
 *              scale (to expand or compress the intensity of the spectrum)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *              hivals (L_HI_TO_WHITE or L_HI_TO_BLACK)
 *      Return: pixd (8 bpp), or null on error
 *
 *  Notes:
 *      (1) Typically @scale is set to 1.0.
 *      (2) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the spectra to the center of pix.
 *      (3) Set @hivals to L_HI_TO_WHITE to convert higher intensities to
 *          white, or L_HI_TO_BLACK to convert them to black instead.
 */
PIX *
pixConvertSpectrumToPix(DPIX      *dpix,
						l_float32  scale,
						l_int32    shiftflag,
						l_int32    hivals)
{
	l_int32    w, h;
	l_int32    i, j, wpl;
	l_float64  max, c;
	l_float64 *data, *line;
	DPIX      *dpixc;
	PIX       *pixd;

	PROCNAME("pixConvertSpectrumToPix");
	
	if (!dpix)
        return (PIX *)ERROR_PTR("dpix not defined", procName, NULL);
	if (shiftflag != L_NO_SHIFTING && shiftflag != L_WITH_SHIFTING)
        return (PIX *)ERROR_PTR("invalid shiftflag", procName, NULL);
	if (hivals != L_HI_TO_WHITE && hivals != L_HI_TO_BLACK)
        return (PIX *)ERROR_PTR("invalid hivals", procName, NULL);
	
	/* Make a working copy of the DPix */
	if ((dpixc = dpixCopy(NULL, dpix)) == NULL)
        return (PIX *)ERROR_PTR("dpixc not made", procName, NULL);
	
	wpl = dpixGetWpl(dpixc);
	data = dpixGetData(dpixc);
	dpixGetDimensions(dpixc, &w, &h);
	for (i = 0; i < h; i++) {
		line = data + i * wpl;
		for (j = 0; j < w; j++) {
			line[j] = scale * log(1 + line[j]);
		}
	}
	
	/* Rescale intensity values and convert to Pix */
	dpixMaxDynamicRange(dpixc, dpixc, 255, hivals);
	pixd = dpixConvertToPix(dpixc, 8, L_CLIP_TO_ZERO, 0, shiftflag);
	
	dpixDestroy(&dpixc);
	
	return pixd;
}

/*--------------------------------------------------------------------*
 *                Pix reconstruction from spectra                     *
 *--------------------------------------------------------------------*/
/*!
 *  pixRestoreFromSpectraGray()
 *
 *      Input:  spectrum (the Fourier spectrum of the original image)
 *              phaseAngle (the phase angle of the original image)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: pixd (8 bpp), or null on error
 *
 *  Notes:
 *      (1) Restore an image starting from the Fourier spectrum and
 *          the phase angle.
 *          The real part of the complex numbers is computed by
 *          spectrum[k] * cos(phaseAngle[k]).
 *          The imaginary part of the complex numbers is computed by
 *          spectrum[k] * sin(phaseAngle[k]).
 *          Once the complex numbers are computed, the inverse DFT
 *          is applied to generate pixd.
 *      (2) The spectrum and phaseAngle DPIX must have same width
 *          and height.
 *      (3) Set @shiftflag to L_WITH_SHIFTING if the spectrum and phase angle
 *          where shifted to move the DC to the center of pix.
 *      (4) Set @padflag to L_WITH_PADDING if the original image was padded
 *          before computing its spectrum and phase angle.
 */
PIX *
pixRestoreFromSpectraGray(DPIX    *spectrum,
						  DPIX    *phaseAngle,
						  l_int32  shiftflag)
{
	l_int32       w, h;
	l_int32       i, j, k, wpl;
	l_float64    *datas, *datap, *lines, *linep;
	PIX          *pixd;
	fftw_complex *output;
	
	PROCNAME("pixRestoreFromSpectraGray");
	
	if (shiftflag != L_NO_SHIFTING && shiftflag != L_WITH_SHIFTING)
        return (PIX *)ERROR_PTR("invalid shiftflag", procName, NULL);
	if (!spectrum || !phaseAngle)
        return (PIX *)ERROR_PTR("spectrum or phaseAngle not defined", procName, NULL);
	if (!dpixSizesEqual(spectrum, phaseAngle))
        return (PIX *)ERROR_PTR("spectrum and phaseAngle unequal size", procName, NULL);
	dpixGetDimensions(spectrum, &w, &h);
	if ((output = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * h * (w / 2 + 1))) == NULL)
		return (PIX *)ERROR_PTR("output not made", procName, NULL);
	
	/* Compute the complex array values */
	wpl = dpixGetWpl(spectrum);
	datas = dpixGetData(spectrum);
	datap = dpixGetData(phaseAngle);
	for (i = 0, k = 0; i < h; i++) {
		lines = datas + i * wpl;
		linep = datap + i * wpl;
		for (j = 0; j < w / 2 + 1; j++, k++) {
			output[k] = (lines[j] * cos(linep[j])) + I * (lines[j] * sin(linep[j]));
		}
	}
	
	/* Compute the inverse of the DFT */
	pixd = pixInverseDFT(output, w, h, shiftflag, L_CLIP_TO_ZERO);
	
	fftw_free(output);
	return pixd;
}

/* --------------------------------------------*/
#endif  /* HAVE_LIBFFTW3 */
/* --------------------------------------------*/