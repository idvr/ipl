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
 *  dpix.c
 *
 *   DPix general-purpose utilities
 *       l_int32        dpixSizesEqual()
 *       l_int32        dpixGetWidth()
 *       l_int32        dpixGetHeight()
 *
 *   DPix  <-->  Pix conversions
 *       DPIX          *pixConvertToDPix()
 *       PIX           *dpixConvertToPix()
 *
 *   DPix normalization
 *       DPIX          *dpixNormalize()
 *
 *   DPix scaling for maximum dynamic range
 *       DPIX          *dpixMaxDynamicRange()
 *
 *   DPix min/max value
 *       l_int32        dpixGetMin()
 *       l_int32        dpixGetMax()
 *       l_int32        dpixGetMinMax()
 *
 *   DPix convolution for mean square (with 1 bpp handling)
 *       DPIX          *dpixMeanSquareAccum()
 */

#include "ipl.h"
#include <math.h>

/*-------------------------------------------------------------------------*
 *                    DPix general-purpose utilities                       *
 *-------------------------------------------------------------------------*/
/*!
 *  dpixSizesEqual()
 *
 *      Input:  two dpix
 *      Return: 1 if the two dpix have same {h, w}; 0 otherwise.
 */
l_int32
dpixSizesEqual(DPIX  *dpix1,
			   DPIX  *dpix2)
{
    PROCNAME("dpixSizesEqual");
	
    if (!dpix1 || !dpix2)
        return ERROR_INT("dpix1 and dpix2 not both defined", procName, 0);
	
    if (dpix1 == dpix2)
        return 1;
	
    if ((dpixGetWidth(dpix1) != dpixGetWidth(dpix2)) ||
        (dpixGetHeight(dpix1) != dpixGetHeight(dpix2)))
        return 0;
    else
        return 1;
}

l_int32
dpixGetWidth(DPIX  *dpix)
{
    PROCNAME("dpixGetWidth");
	
    if (!dpix)
        return ERROR_INT("dpix not defined", procName, UNDEF);
	
    return dpix->w;
}

l_int32
dpixGetHeight(DPIX  *dpix)
{
    PROCNAME("dpixGetHeight");
	
    if (!dpix)
        return ERROR_INT("dpix not defined", procName, UNDEF);
	
    return dpix->h;
}

/*--------------------------------------------------------------------*
 *                     DPix  <-->  Pix conversions                    *
 *--------------------------------------------------------------------*/
/*!
 *  pixConvertToDPix()
 *
 *      Input:  pix (1, 2, 4, 8, 16 or 32 bpp)
 *              ncomps (number of components: 3 for RGB, 1 otherwise)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: dpix, or null on error
 *
 *  Notes:
 *      (1) If colormapped, remove to grayscale.
 *      (2) If 32 bpp and @ncomps == 3, this is RGB; convert to luminance.
 *          In all other cases the src image is treated as having a single
 *          component of pixel values.
 *      (3) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the DFT to the center of pix.
 */
DPIX *
pixConvertToDPix(PIX     *pixs,
                 l_int32  ncomps,
				 l_int32  shiftflag)
{
	l_int32     w, h, d;
	l_int32     i, j, val, wplt, wpld;
	l_uint32    uval;
	l_uint32   *datat, *linet;
	l_float64  *datad, *lined;
	PIX        *pixt;
	DPIX       *dpixd;
	
    PROCNAME("pixConvertToDPix");
	
    if (!pixs)
        return (DPIX *)ERROR_PTR("pixs not defined", procName, NULL);
	if (shiftflag != L_NO_SHIFTING && shiftflag != L_WITH_SHIFTING)
        return (DPIX *)ERROR_PTR("invalid shiftflag", procName, NULL);
	
    if (pixGetColormap(pixs))
        pixt = pixRemoveColormap(pixs, REMOVE_CMAP_TO_GRAYSCALE);
    else if (pixGetDepth(pixs) == 32 && ncomps == 3)
        pixt = pixConvertRGBToLuminance(pixs);
    else
        pixt = pixClone(pixs);
	
    pixGetDimensions(pixt, &w, &h, &d);
    if ((dpixd = dpixCreate(w, h)) == NULL)
        return (DPIX *)ERROR_PTR("dpixd not made", procName, NULL);
    datat = pixGetData(pixt);
    wplt = pixGetWpl(pixt);
    datad = dpixGetData(dpixd);
    wpld = dpixGetWpl(dpixd);
    for (i = 0; i < h; i++) {
        linet = datat + i * wplt;
        lined = datad + i * wpld;
        if (d == 1) {
            for (j = 0; j < w; j++) {
                val = GET_DATA_BIT(linet, j);
                lined[j] = shiftflag ? (l_float64)(val * pow(-1, i + j)) : (l_float64)val;
            }
        }
        else if (d == 2) {
            for (j = 0; j < w; j++) {
                val = GET_DATA_DIBIT(linet, j);
                lined[j] = shiftflag ? (l_float64)(val * pow(-1, i + j)) : (l_float64)val;
            }
        }
        else if (d == 4) {
            for (j = 0; j < w; j++) {
                val = GET_DATA_QBIT(linet, j);
                lined[j] = shiftflag ? (l_float64)(val * pow(-1, i + j)) : (l_float64)val;
            }
        }
        else if (d == 8) {
            for (j = 0; j < w; j++) {
                val = GET_DATA_BYTE(linet, j);
                lined[j] = shiftflag ? (l_float64)(val * pow(-1, i + j)) : (l_float64)val;
            }
        }
        else if (d == 16) {
            for (j = 0; j < w; j++) {
                val = GET_DATA_TWO_BYTES(linet, j);
                lined[j] = shiftflag ? (l_float64)(val * pow(-1, i + j)) : (l_float64)val;
            }
        }
        else if (d == 32) {
            for (j = 0; j < w; j++) {
                uval = GET_DATA_FOUR_BYTES(linet, j);
                lined[j] = shiftflag ? (l_float64)(uval * pow(-1, i + j)) : (l_float64)uval;
            }
        }
    }
	
    pixDestroy(&pixt);
    return dpixd;
}

/*!
 *  dpixConvertToPix()
 *
 *      Input:  dpixs
 *              outdepth (0, 1, 8, 16 or 32 bpp)
 *              negvals (L_CLIP_TO_ZERO, L_TAKE_ABSVAL,
 *                       L_THRESH_NEG_TO_BLACK or L_THRESH_NEG_TO_WHITE)
 *              errorflag (1 to output error stats; 0 otherwise)
 *              shiftflag (L_NO_SHIFTING or L_WITH_SHIFTING)
 *      Return: pixd, or null on error
 *
 *  Notes:
 *      (1) Use @outdepth = 0 to programmatically determine the
 *          output depth.  If no values are greater than 255,
 *          it will set outdepth = 8; otherwise to 16 or 32.
 *      (2) Because we are converting a float to an unsigned int
 *          with a specified dynamic range (8, 16 or 32 bits), errors
 *          can occur.  If errorflag == TRUE, output the number
 *          of values out of range, both negative and positive.
 *      (3) If a pixel value is positive and out of range, clip to
 *          the maximum value represented at the outdepth of 8, 16
 *          or 32 bits.
 *      (4) Set @shiftflag to L_WITH_SHIFTING to move the DC of the
 *          the DFT to the center of pix.
 *          When @shiftflag == L_WITH_SHIFTING, pixd is multiplied
 *          by pow(-1, x + y).
 */
PIX *
dpixConvertToPix(DPIX    *dpixs,
                 l_int32  outdepth,
                 l_int32  negvals,
                 l_int32  errorflag,
				 l_int32  shiftflag)
{
	l_int32     w, h, wpls, wpld, maxval;
	l_int32     i, j;
	l_uint32    vald;
	l_float64   val;
	l_float64  *datas, *lines;
	l_uint32   *datad, *lined;
	PIX        *pixd;
	
    PROCNAME("dpixConvertToPix");
	
    if (!dpixs)
        return (PIX *)ERROR_PTR("dpixs not defined", procName, NULL);
	if (negvals != L_CLIP_TO_ZERO && negvals != L_TAKE_ABSVAL &&
		negvals != L_THRESH_NEG_TO_BLACK && negvals != L_THRESH_NEG_TO_WHITE)
        return (PIX *)ERROR_PTR("invalid negvals", procName, NULL);
    if (outdepth != 0 && outdepth != 8 && outdepth != 16 && outdepth != 32)
        return (PIX *)ERROR_PTR("outdepth not in {0,8,16,32}", procName, NULL);
	if (shiftflag != L_NO_SHIFTING && shiftflag != L_WITH_SHIFTING)
        return (PIX *)ERROR_PTR("invalid shiftflag", procName, NULL);
	
    dpixGetDimensions(dpixs, &w, &h);
    datas = dpixGetData(dpixs);
    wpls = dpixGetWpl(dpixs);
	
	/* Adaptive determination of output depth */
    if (outdepth == 0) {
		outdepth = 8;
		for (i = 0; i < h; i++) {
			lines = datas + i * wpls;
			for (j = 0; j < w; j++) {
				val = lines[j];
				if (val > 65535.5) {
					outdepth = 32;
					break;
				}
				if (val > 255.5)
					outdepth = 16;
			}
			if (outdepth == 32) break;
		}
    }
    maxval = (1 << outdepth) - 1;
	
	/* Gather statistics if @errorflag = TRUE */
    if (errorflag) {
        l_int32  negs = 0;
        l_int32  overvals = 0;
        for (i = 0; i < h; i++) {
            lines = datas + i * wpls;
            for (j = 0; j < w; j++) {
                val = lines[j];
                if (val < 0.0)
                    negs++;
                else if (val > maxval)
                    overvals++;
            }
        }
        if (negs > 0)
            L_ERROR_INT("Number of negative values: %d", procName, negs);
        if (overvals > 0)
            L_ERROR_INT("Number of too-large values: %d", procName, overvals);
    }
	
	/* Make the pix and convert the data */
    if ((pixd = pixCreate(w, h, outdepth)) == NULL)
        return (PIX *)ERROR_PTR("pixd not made", procName, NULL);
    datad = pixGetData(pixd);
    wpld = pixGetWpl(pixd);
    for (i = 0; i < h; i++) {
        lines = datas + i * wpls;
        lined = datad + i * wpld;
        for (j = 0; j < w; j++) {
            val = lines[j];
			if (shiftflag)
				val *= pow(-1, i + j);
            if (val > 0.0) {
				if (negvals == L_THRESH_NEG_TO_BLACK || negvals == L_THRESH_NEG_TO_WHITE)
					vald = negvals == L_THRESH_NEG_TO_BLACK ? maxval : 0;
				else
					vald = (l_uint32)(val + 0.5);
			}
            else { /* val <= 0.0 */
				if (negvals == L_THRESH_NEG_TO_BLACK || negvals == L_THRESH_NEG_TO_WHITE)
					vald = negvals == L_THRESH_NEG_TO_BLACK ? 0 : maxval;
                else if (negvals == L_CLIP_TO_ZERO)
                    vald = 0;
                else
                    vald = (l_uint32)(-val + 0.5);
            }
            if (vald > maxval)
                vald = maxval;
            if (outdepth == 8)
                SET_DATA_BYTE(lined, j, vald);
            else if (outdepth == 16)
                SET_DATA_TWO_BYTES(lined, j, vald);
            else  /* outdepth == 32 */
                SET_DATA_FOUR_BYTES(lined, j, vald);
        }
    }
	
    return pixd;
}

/*--------------------------------------------------------------------*
 *                 DPix normalization and shifting                    *
 *--------------------------------------------------------------------*/
/*!
 *  dpixNormalize()
 *
 *      Input:  dpixd  (<optional>; this can be null, equal to dpixs,
 *                     or different from dpixs)
 *              dpixs
 *      Return: dpixd, or null on error
 *
 *  Notes:
 *      (1) Normalize dpixs by dividing each value by the size of the
 *          image. This is to compensate for the FFTW library returning
 *          unnormalized DFT values.
 *      (2) There are 3 cases:
 *           (a) dpixd == null,   ~src --> new dpixd
 *           (b) dpixd == dpixs,  ~src --> src  (in-place)
 *           (c) dpixd != dpixs,  ~src --> input dpixd
 *      (3) For clarity, if the case is known, use these patterns:
 *           (a) dpixd = dpixNormalize(NULL, dpixs);
 *           (b) dpixNormalize(dpixs, dpixs);
 *           (c) dpixNormalize(dpixd, dpixs);
 */
DPIX *
dpixNormalize(DPIX *dpixd,
			  DPIX *dpixs)
{
	l_int32     w, h, n, i, j, wpl;
	l_float64  *data, *line;
	
    PROCNAME("dpixNormalize");
	
    if (!dpixs)
        return (DPIX *)ERROR_PTR("dpixs not defined", procName, NULL);
	
	/* Prepare dpixd for in-place operation */
    if ((dpixd = dpixCopy(dpixd, dpixs)) == NULL)
		return (DPIX *)ERROR_PTR("dpixd not made", procName, NULL);
	
    dpixGetDimensions(dpixd, &w, &h);
    data = dpixGetData(dpixd);
    wpl = dpixGetWpl(dpixd);
	n = w * h;
	
	for (i = 0; i < h; i++) {
		line = data + i * wpl;
		for (j = 0; j < w; j++) {
			line[j] /= n;
		}
	}
	
	return dpixd;
}

/*--------------------------------------------------------------------*
 *             DPix scaling for maximum dynamic range                 *
 *--------------------------------------------------------------------*/
/*!
 *  dpixMaxDynamicRange()
 *
 *      Input:  dpixd  (<optional>; this can be null, equal to dpixs,
 *                     or different from dpixs)
 *              dpixs
 *              k      (maximum intensity value)
 *              hivals (L_HI_TO_WHITE or L_HI_TO_BLACK)
 *      Return: dpixd, or null on error
 *
 *  Notes:
 *      (1) Linearly rescale dpix values to maximize the dynamic range
 *          within range [0, k]
 *      (2) There are 3 cases:
 *           (a) dpixd == null,   ~src --> new dpixd
 *           (b) dpixd == dpixs,  ~src --> src  (in-place)
 *           (c) dpixd != dpixs,  ~src --> input dpixd
 *      (3) For clarity, if the case is known, use these patterns:
 *           (a) dpixd = dpixNormalize(NULL, dpixs);
 *           (b) dpixNormalize(dpixs, dpixs);
 *           (c) dpixNormalize(dpixd, dpixs);
 *      (4) Set @hivals to L_HI_TO_WHITE to convert higher intensities to
 *          white, or L_HI_TO_BLACK to convert them to black instead.
 */
DPIX *
dpixMaxDynamicRange(DPIX    *dpixd,
					DPIX    *dpixs,
					l_int32  k,
					l_int32  hivals)
{
	l_int32     w, h, n, i, j, wpl;
	l_float64  *data, *line;
	l_float64   minval, maxval;
	
    PROCNAME("dpixMaxDynamicRange");
	
    if (!dpixs)
        return (DPIX *)ERROR_PTR("dpixs not defined", procName, NULL);
	
	if (hivals != L_HI_TO_WHITE && hivals != L_HI_TO_BLACK)
        return (DPIX *)ERROR_PTR("invalid hivals", procName, NULL);
	
	/* Prepare dpixd for in-place operation */
    if ((dpixd = dpixCopy(dpixd, dpixs)) == NULL)
		return (DPIX *)ERROR_PTR("dpixd not made", procName, NULL);
	
    dpixGetDimensions(dpixd, &w, &h);
    data = dpixGetData(dpixd);
    wpl = dpixGetWpl(dpixd);

	dpixGetMinMax(dpixd, &minval, &maxval);
	maxval -= minval;
	
	/* Prevent division by 0 */
	if (!maxval)
		maxval += 1.0;
	
	for (i = 0; i < h; i++) {
		line = data + i * wpl;
		for (j = 0; j < w; j++) {
			if (hivals == L_HI_TO_WHITE)
				line[j] = k * ((line[j] - minval) / maxval);
			else
				line[j] = k - (k * ((line[j] - minval) / maxval));
		}
	}
	
	return dpixd;
}

/*--------------------------------------------------------------------*
 *                           Min/max value                            *
 *--------------------------------------------------------------------*/
/*!
 *  dpixGetMin()
 *
 *      Input:  dpix
 *              &minval (<optional return> min value)
 *              &xminloc (<optional return> x location of min)
 *              &yminloc (<optional return> y location of min)
 *      Return: 0 if OK; 1 on error
 */
l_int32
dpixGetMin(DPIX       *dpix,
           l_float64  *pminval,
           l_int32    *pxminloc,
           l_int32    *pyminloc)
{
	l_int32     i, j, w, h, wpl, xminloc, yminloc;
	l_float64  *data, *line;
	l_float64   minval;
	
    PROCNAME("dpixGetMin");
	
    if (!pminval && !pxminloc && !pyminloc)
        return ERROR_INT("nothing to do", procName, 1);
    if (pminval) *pminval = 0.0;
    if (pxminloc) *pxminloc = 0;
    if (pyminloc) *pyminloc = 0;
    if (!dpix)
        return ERROR_INT("dpix not defined", procName, 1);
	
    minval = +1.0e40;
    xminloc = 0;
    yminloc = 0;
    dpixGetDimensions(dpix, &w, &h);
    data = dpixGetData(dpix);
    wpl = dpixGetWpl(dpix);
    for (i = 0; i < h; i++) {
        line = data + i * wpl;
        for (j = 0; j < w; j++) {
            if (line[j] < minval) {
                minval = line[j];
                xminloc = j;
                yminloc = i;
            }
        }
    }
	
    if (pminval) *pminval = minval;
    if (pxminloc) *pxminloc = xminloc;
    if (pyminloc) *pyminloc = yminloc;
    return 0;
}

/*!
 *  dpixGetMax()
 *
 *      Input:  dpix
 *              &maxval (<optional return> max value)
 *              &xmaxloc (<optional return> x location of max)
 *              &ymaxloc (<optional return> y location of max)
 *      Return: 0 if OK; 1 on error
 */
l_int32
dpixGetMax(DPIX       *dpix,
           l_float64  *pmaxval,
           l_int32    *pxmaxloc,
           l_int32    *pymaxloc)
{
	l_int32     i, j, w, h, wpl, xmaxloc, ymaxloc;
	l_float64  *data, *line;
	l_float64   maxval;
	
    PROCNAME("dpixGetMax");
	
    if (!pmaxval && !pxmaxloc && !pymaxloc)
        return ERROR_INT("nothing to do", procName, 1);
    if (pmaxval) *pmaxval = 0.0;
    if (pxmaxloc) *pxmaxloc = 0;
    if (pymaxloc) *pymaxloc = 0;
    if (!dpix)
        return ERROR_INT("dpix not defined", procName, 1);
	
    maxval = -1.0e40;
    xmaxloc = 0;
    ymaxloc = 0;
    dpixGetDimensions(dpix, &w, &h);
    data = dpixGetData(dpix);
    wpl = dpixGetWpl(dpix);
    for (i = 0; i < h; i++) {
        line = data + i * wpl;
        for (j = 0; j < w; j++) {
            if (line[j] > maxval) {
                maxval = line[j];
                xmaxloc = j;
                ymaxloc = i;
            }
        }
    }
	
    if (pmaxval) *pmaxval = maxval;
    if (pxmaxloc) *pxmaxloc = xmaxloc;
    if (pymaxloc) *pymaxloc = ymaxloc;
    return 0;
}

/*!
 *  dpixGetMinMax()
 *
 *      Input:  dpix
 *              &minval (<optional return> min value)
 *              &maxval (<optional return> max value)
 *      Return: 0 if OK; 1 on error
 */
l_int32
dpixGetMinMax(DPIX       *dpix,
			  l_float64  *pminval,
			  l_float64  *pmaxval)
{
	l_int32     i, j, w, h, wpl;
	l_float64  *data, *line;
	l_float64   minval, maxval;
	
    PROCNAME("dpixGetMinMax");
	
    if (!pminval && !pmaxval)
        return ERROR_INT("nothing to do", procName, 1);
    if (pminval) *pminval = 0.0;
	if (pmaxval) *pmaxval = 0.0;
    if (!dpix)
        return ERROR_INT("dpix not defined", procName, 1);
	
    minval = +1.0e40;
    maxval = -1.0e40;
    dpixGetDimensions(dpix, &w, &h);
    data = dpixGetData(dpix);
    wpl = dpixGetWpl(dpix);
    for (i = 0; i < h; i++) {
        line = data + i * wpl;
        for (j = 0; j < w; j++) {
            if (line[j] < minval) minval = line[j];
			if (line[j] > maxval) maxval = line[j];
        }
    }
	
    if (pminval) *pminval = minval;
    if (pmaxval) *pmaxval = maxval;
    return 0;
}

/*--------------------------------------------------------------------*
 *       DPix convolution for mean square (with 1 bpp handling)       *
 *--------------------------------------------------------------------*/
/*!
 *  dpixMeanSquareAccum()
 *
 *      Input:  pixs (1 bpp or 8 bpp grayscale)
 *      Return: dpix (64 bit array), or null on error
 *
 *  Notes:
 *      (1) This is an extension to the standard pixMeanSquareAccum()
 *          implementation provided by Leptonica, to handle 1bpp binary pix
 *          transparently.
 *      (1) Similar to pixBlockconvAccum(), this computes the
 *          sum of the squares of the pixel values in such a way
 *          that the value at (i,j) is the sum of all squares in
 *          the rectangle from the origin to (i,j).
 *      (2) The general recursion relation (v are squared pixel values) is
 *            a(i,j) = v(i,j) + a(i-1, j) + a(i, j-1) - a(i-1, j-1)
 *          For the first line, this reduces to the special case
 *            a(i,j) = v(i,j) + a(i, j-1)
 *          For the first column, the special case is
 *            a(i,j) = v(i,j) + a(i-1, j)
 */
DPIX *
dpixMeanSquareAccum(PIX  *pixs)
{
	l_int32     i, j, w, h, d, wpl, wpls, val;
	l_uint32   *datas, *lines;
	l_float64  *data, *line, *linep;
	DPIX       *dpix;
	
    PROCNAME("dpixMeanSquareAccum");
	
    if (!pixs)
        return (DPIX *)ERROR_PTR("pixs not defined", procName, NULL);
    pixGetDimensions(pixs, &w, &h, &d);
    if (d != 1 && d != 8)
        return (DPIX *)ERROR_PTR("pixs not 1 bpp or 8 bpp", procName, NULL);
    if ((dpix = dpixCreate(w, h)) ==  NULL)
        return (DPIX *)ERROR_PTR("dpix not made", procName, NULL);
	
    datas = pixGetData(pixs);
    wpls = pixGetWpl(pixs);
    data = dpixGetData(dpix);
    wpl = dpixGetWpl(dpix);
	
    lines = datas;
    line = data;
    for (j = 0; j < w; j++) {   /* first line */
        val = d == 1 ? GET_DATA_BIT(lines, j) : GET_DATA_BYTE(lines, j);
        if (j == 0)
            line[0] = val * val;
        else
            line[j] = line[j - 1] + val * val;
    }
	
	/* Do the other lines */
    for (i = 1; i < h; i++) {
        lines = datas + i * wpls;
        line = data + i * wpl;  /* current dest line */
        linep = line - wpl;;  /* prev dest line */
        for (j = 0; j < w; j++) {
            val = d == 1 ? GET_DATA_BIT(lines, j) : GET_DATA_BYTE(lines, j);
            if (j == 0)
                line[0] = linep[0] + val * val;
            else
                line[j] = line[j - 1] + linep[j] - linep[j - 1] + val * val;
        }
    }
	
    return dpix;
}
