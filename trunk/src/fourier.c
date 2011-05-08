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
 *  fourier.c
 *
 *   DFT testing
 *       PIX           *pixTestDFT()
 *       PIX           *pixTestDFTGray()
 *
 *   Low-level implementation of Discrete Fourier Transform 
 *       fftw_complex  *pixDFT()
 *       PIX           *pixInverseDFT()
 *       DPIX          *dpixInverseDFT()
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
 *                            DFT testing                             *
 *--------------------------------------------------------------------*/
/*!
 *  pixTestDFT()
 *
 *      Input:  pix (1, 8 or 32 bpp; or 2, 4 or 8 bpp with colormap)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: pixd, or null on error
 *
 *  Notes:
 *      (1) This is only a testing utility routine that performs
 *          DFT and inverse DFT in sequence on pixs. If everything
 *          works correctly, then pixd yields an image equal to pixs.
 *      (2) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the DFT to the center of pix.
 */
PIX *
pixTestDFT(PIX     *pixs,
		   l_int32  shiftflag)
{
	l_int32	 w, h, d;
	PIX     *pixt, *pixd, *pixr, *pixrc, *pixg, *pixgc, *pixb, *pixbc;
	
	PROCNAME("pixTestDFT");
	
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
	
    if (d != 1 && d != 8 && d != 32) {
        pixDestroy(&pixt);
        return (PIX *)ERROR_PTR("depth not 1, 8 or 32 bpp", procName, NULL);
    }
	
	if (d == 1 || d == 8) {
		pixd = pixTestDFTGray(pixt, shiftflag);
	}
	else { /* d == 32 */
		pixr = pixGetRGBComponent(pixt, COLOR_RED);
        pixrc = pixTestDFTGray(pixr, shiftflag);
        pixDestroy(&pixr);
        pixg = pixGetRGBComponent(pixt, COLOR_GREEN);
        pixgc = pixTestDFTGray(pixg, shiftflag);
        pixDestroy(&pixg);
        pixb = pixGetRGBComponent(pixt, COLOR_BLUE);
        pixbc = pixTestDFTGray(pixb, shiftflag);
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
 *  pixTestDFTGray()
 *
 *      Input:  pix (1 or 8 bpp)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: pixd, or null on error
 *
 *  Notes:
 *      (1) This is only a testing utility routine that performs
 *          DFT and inverse DFT in sequence on pixs. If everything
 *          works correctly, then pixd yields an image equal to pixs.
 *      (2) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the DFT to the center of pix.
 */
PIX *
pixTestDFTGray(PIX     *pixs,
			   l_int32  shiftflag)
{
	l_int32	       w, h, d;
	PIX           *pixd;
	fftw_complex  *output;
	
	PROCNAME("pixTestDFTGray");
	
    if (!pixs)
        return (PIX *)ERROR_PTR("pixs not defined", procName, NULL);
	if (shiftflag != L_NO_SHIFTING && shiftflag != L_WITH_SHIFTING)
        return (PIX *)ERROR_PTR("invalid shiftflag", procName, NULL);
	
	pixGetDimensions(pixs, &w, &h, &d);
	if (d != 1 && d != 8)
        return (PIX *)ERROR_PTR("pixs not 1 bpp or 8 bpp", procName, NULL);
	
	/* Calculate the DFT of pixs */
	if ((output = pixDFT(pixs, shiftflag)) == NULL)
        return (PIX *)ERROR_PTR("pixs DFT not computed", procName, NULL);
		
	/* Compute the inverse of the DFT */
	pixd = pixInverseDFT(output, w, h, shiftflag, L_CLIP_TO_ZERO);
	
	/* Release the allocated resources */
	fftw_free(output);
	
	return pixd;
}

/*--------------------------------------------------------------------*
 *      Low-level implementation of Discrete Fourier Transform        *
 *--------------------------------------------------------------------*/
/*!
 *  pixDFT()
 *
 *      Input:  pix (1 bpp or 8 bpp)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: complex array, or null on error
 *
 *  Notes:
 *      (1) The complex array returned has size (pixs->h) * (pixs->w / 2 + 1).
 *          This is to save space, given the fact the other half of the
 *          transform can be calculated by the complex conjugate.
 *      (2) By default, the DC of the DFT is in the top left corner (0, 0).
 *          Set @shiftflag to L_WITH_SHIFTING to move the DC to the center.
 *      (3) It is the responsibility of the caller to release the allocated
 *          complex array by invoking fftw_free().
 */
fftw_complex *
pixDFT(PIX      *pixs,
	   l_int32   shiftflag)
{
	l_int32	       w, h, d;
	l_int32        i, j, k;
	DPIX          *dpix;
	fftw_complex  *output;
	fftw_plan      plan;
	
	PROCNAME("pixDFT");
	
    if (!pixs)
        return (fftw_complex *)ERROR_PTR("pixs not defined", procName, NULL);
	if (shiftflag != L_NO_SHIFTING && shiftflag != L_WITH_SHIFTING)
        return (fftw_complex *)ERROR_PTR("invalid shiftflag", procName, NULL);
	
	pixGetDimensions(pixs, &w, &h, &d);
	if (d != 1 && d != 8)
        return (fftw_complex *)ERROR_PTR("pixs not 1 bpp or 8 bpp", procName, NULL);
	
	/* Convert Pix to a DPix that can be fed to the FFTW library */
	if ((dpix = pixConvertToDPix(pixs, 1, shiftflag)) == NULL)
        return (fftw_complex *)ERROR_PTR("dpix not made", procName, NULL);
	
	/* Compute the DFT of the DPix */
	output = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * h * (w / 2 + 1));
	plan = fftw_plan_dft_r2c_2d(h, w, (double *) dpixGetData(dpix), output, FFTW_ESTIMATE);
	fftw_execute(plan);
	
	dpixDestroy(&dpix);
	fftw_destroy_plan(plan);
	
	return output;
}

/*!
 *  pixInverseDFT()
 *
 *      Input:  dft
 *              w, h (image size)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *              outflag (L_CLIP_TO_ZERO, L_TAKE_ABSVAL,
 *                       L_THRESH_NEG_TO_BLACK or L_THRESH_NEG_TO_WHITE)
 *      Return: pixd (8bpp), or null on error
 *
 *  Notes:
 *      (1) Set @shiftflag to L_WITH_SHIFTING if the DC was moved to the
 *          center of the image during the DFT computation. Reshifting it
 *          will move the DC back to the top left corner (0, 0).
 *      (2) It is the responsibility of the caller to release the allocated
 *          complex array by invoking fftw_free().
 */
PIX *
pixInverseDFT(fftw_complex *dft,
			  l_int32       w,
			  l_int32       h,
			  l_int32       shiftflag,
			  l_int32       outflag)
{
	PIX          *pixd;
	DPIX         *dpix;
	
	PROCNAME("pixInverseDFT");
	
	dpix = dpixInverseDFT(dft, w, h);
	
	/* Convert DPix to a Pix */
	pixd = dpixConvertToPix(dpix, 8, outflag, 0, shiftflag);
	dpixDestroy(&dpix);
	
	return pixd;
}

/*!
 *  dpixInverseDFT()
 *
 *      Input:  dft
 *              w, h (image size)
 *      Return: dpix (unnormalized), or null on error
 */
DPIX *
dpixInverseDFT(fftw_complex *dft,
			   l_int32       w,
			   l_int32       h)
{
	DPIX         *dpix;
	fftw_plan    plan;
	
	PROCNAME("dpixInverseDFT");
	
    if (!dft)
        return (DPIX *)ERROR_PTR("dft not defined", procName, NULL);
	if ((dpix = dpixCreate(w, h)) == NULL)
        return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
	/* Compute the inverse DFT, storing the results into DPix */
	plan = fftw_plan_dft_c2r_2d(h, w, dft, (double *) dpixGetData(dpix), FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	
	dpixNormalize(dpix, dpix);
	
	return dpix;
}

/* --------------------------------------------*/
#endif  /* HAVE_LIBFFTW3 */
/* --------------------------------------------*/