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
 *  ipl.h
 */

#ifndef __IPL_H__
#define __IPL_H__

#include <leptonica/allheaders.h>

/*------------------------------------------------------------------------*
 *  Defines and includes differ for Unix and Windows.  Also for Windows,  *
 *  differentiate between conditionals based on platform and compiler.    *
 *      For platforms:                                                    *
 *          _WIN32       =>     Windows, 32- or 64-bit                    *
 *          _WIN64       =>     Windows, 64-bit only                      *
 *          __CYGWIN__   =>     Cygwin                                    *
 *      For compilers:                                                    *
 *          __GNUC__     =>     gcc                                       *
 *          _MSC_VER     =>     msvc                                      *
 *------------------------------------------------------------------------*/

/* Windows specifics */

#ifdef _WIN32

/* DLL EXPORT/IMPORT */
#ifdef IPLLIB_EXPORTS
#define IPL_DLL __declspec(dllexport)
#elif defined(IPLLIB_IMPORTS)
#define IPL_DLL __declspec(dllimport)
#else
#define IPL_DLL
#endif

#else  /* non-WINDOWS-SPECIFICS */
#include <stdint.h>
#define IPL_DLL
#endif  /* _WIN32 */

/*-------------------------------------------------------------------------*
 *             Non-standard math definitiions and declarations             *
 *-------------------------------------------------------------------------*/
#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* pi */
#endif

/*-------------------------------------------------------------------------*
 *          Handling negative values in conversion to unsigned int         *
 *-------------------------------------------------------------------------*/
enum {
	/* Inherited from Leptonica:                                           */
    /* L_CLIP_TO_ZERO = 1,        Clip negative values to 0                */
    /* L_TAKE_ABSVAL = 2          Convert to positive using L_ABS()        */
	
	/* Extended by IPL:                                                    */
	L_THRESH_NEG_TO_BLACK = 3, /* Set to black all negative values and to
								  white all positive values.
								  The output is a binary image.            */
	L_THRESH_NEG_TO_WHITE = 4  /* Set to black all positive values and to
								  white all negative values.
								  The output is a binary image.            */
};

/*-------------------------------------------------------------------------*
 *                Dynamic range intensity values in DPix                   *
 *-------------------------------------------------------------------------*/
enum {
	L_HI_TO_WHITE = 1,         /* Set higher intensities to white          */
	L_HI_TO_BLACK = 2          /* Set higher intensities to black          */
};

/*-------------------------------------------------------------------------*
 *           Coefficient values to use in windowing functions              *
 *-------------------------------------------------------------------------*/
enum {
	L_USE_TRIANG = 1,          /* Used to differentiate between Bartlett
								  and the triangular window functions.     */
	L_USE_HAMMING = 2          /* Used to differentiate between Hann
								  and the Hamming window functions.        */
};

/*-------------------------------------------------------------------------*
 *                Types of filters in the frequency domain                 *
 *-------------------------------------------------------------------------*/
enum {
	L_NO_SHIFTING = 0,
	L_WITH_SHIFTING = 1        /* Shift the DC of the DFT to the center */
};

enum {
    L_LO_PASS = 1,
    L_HI_PASS = 2
};

enum {
	L_BAND_REJECT = 1,
    L_BAND_PASS = 2
};

/*-------------------------------------------------------------------------*
 *                          Types of impulse noise                         *
 *-------------------------------------------------------------------------*/
enum {
	L_NOISE_SALT = 1,
	L_NOISE_PEPPER = 2
};

/*------------------------------------------------------------------------*
 *  Exported routines
 *------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif  /* __cplusplus */
	
	/* corner.c */
	IPL_DLL extern PTA * pixCornerPointsMoravec(PIX *pixs, l_int32 whsize, l_int32 thresh, l_int32 factor, FPIX **pfpixv);
	IPL_DLL extern PTA * fpixLocalMaxima(FPIX *fpix, l_int32 whsize, l_int32 factor, PTA *locations);
	
	/* cornerlow.c */
	IPL_DLL extern l_int32 multiplyLow(l_uint32 *data, l_int32 wpl, l_int32 d, l_int32 wsize, l_int32 i, l_int32 j, l_int32 k, l_int32 l);
	IPL_DLL extern l_int32 calculateVarianceLow(l_uint32 *data, l_int32 wpl, l_int32 d, l_int32 wsize, l_int32 i, l_int32 j, l_int32 k, l_int32 l);
	IPL_DLL extern l_float32 findLocalMaximaLow(l_float32 *data, l_int32 wpl, l_int32 wsize, l_int32 i, l_int32 j);
	IPL_DLL extern l_int16 isLocalMaximaLow(l_float32 val, l_float32 *data, l_int32 wpl, l_int32 wsize, l_int32 i, l_int32 j);

	/* dpix.c */
	IPL_DLL extern l_int32 dpixSizesEqual(DPIX *dpix1, DPIX *dpix2);
	IPL_DLL extern l_int32 dpixGetWidth(DPIX *dpix);
	IPL_DLL extern l_int32 dpixGetHeight(DPIX *dpix);
	IPL_DLL extern DPIX * pixConvertToDPix(PIX *pixs, l_int32 ncomps, l_int32 shiftflag);
	IPL_DLL extern PIX * dpixConvertToPix(DPIX *dpixs, l_int32 outdepth, l_int32 negvals, l_int32 errorflag, l_int32 shiftflag);
	IPL_DLL extern DPIX * dpixNormalize(DPIX *dpixd, DPIX *dpixs);
	IPL_DLL extern DPIX * dpixMaxDynamicRange(DPIX *dpixd, DPIX *dpixs, l_int32 k, l_int32 hivals);
	IPL_DLL extern l_int32 dpixGetMin(DPIX *dpix, l_float64 *pminval, l_int32 *pxminloc, l_int32 *pyminloc);
	IPL_DLL extern l_int32 dpixGetMax(DPIX *dpix, l_float64 *pmaxval, l_int32 *pxmaxloc, l_int32 *pymaxloc);
	IPL_DLL extern l_int32 dpixGetMinMax(DPIX *dpix, l_float64 *pminval, l_float64 *pmaxval);
	IPL_DLL extern DPIX * dpixMeanSquareAccum(PIX *pixs);
	
	/* fourier.c - fourierstub.c: require FFTW3 library */
	IPL_DLL extern PIX * pixTestDFT(PIX *pixs, l_int32 shiftflag);
	IPL_DLL extern PIX * pixTestDFTGray(PIX *pixs, l_int32 shiftflag);
	
	/* frequency1.c */
	IPL_DLL extern DPIX * filterCreateIdeal(l_int32 w, l_int32 h, l_float32 radius, l_int32 typeflag);
	IPL_DLL extern DPIX * filterCreateButterworth(l_int32 w, l_int32 h, l_float32 radius, l_int32 order, l_int32 typeflag);
	IPL_DLL extern DPIX * filterCreateGaussian(l_int32 w, l_int32 h, l_float32 radius, l_int32 typeflag);
	IPL_DLL extern DPIX * filterCreateIdealSelective(l_int32 w, l_int32 h, l_float32 radius, l_float32 bandwidth, l_int32 typeflag);
	IPL_DLL extern DPIX * filterCreateButterworthSelective(l_int32 w, l_int32 h, l_float32 radius, l_float32 bandwidth, l_int32 order, l_int32 typeflag);
	IPL_DLL extern DPIX * filterCreateGaussianSelective(l_int32 w, l_int32 h, l_float32 radius, l_float32 bandwidth, l_int32 typeflag);
	IPL_DLL extern PIX * pixDisplayFilter(DPIX *dpix);
	
	/* frequency2.c - frequency2stub.c: require FFTW3 library */
	IPL_DLL extern PIX * pixIdealFilter(PIX *pixs, l_float32 radius, l_int32 typeflag, l_int32 outflag, DPIX **dpixf);
	IPL_DLL extern PIX * pixButterworthFilter(PIX *pixs, l_float32 radius, l_int32 order, l_int32 typeflag, l_int32 outflag, DPIX **dpixf);
	IPL_DLL extern PIX * pixGaussianFilter(PIX *pixs, l_float32 radius, l_int32 typeflag, l_int32 outflag, DPIX **dpixf);
	IPL_DLL extern PIX * pixIdealSelectiveFilter(PIX *pixs, l_float32 radius, l_float32 bandwidth, l_int32 typeflag, l_int32 outflag, DPIX **dpixf);
	IPL_DLL extern PIX * pixButterworthSelectiveFilter(PIX *pixs, l_float32 radius, l_float32 bandwidth, l_int32 order, l_int32 typeflag, l_int32 outflag, DPIX **dpixf);
	IPL_DLL extern PIX * pixGaussianSelectiveFilter(PIX *pixs, l_float32 radius, l_float32 bandwidth, l_int32 typeflag, l_int32 outflag, DPIX **dpixf);
	IPL_DLL extern PIX * pixApplyFilter(PIX *pixs, DPIX *dpix, l_int32 outflag);
	IPL_DLL extern PIX * pixApplyFilterGray(PIX *pixs, DPIX *dpix, l_int32 outflag);
	
	/* motion.c */
	IPL_DLL extern PIX * pixThresholdWithAbsDifference(PIX *pixs1, PIX *pixs2, l_int32 thresh, PIX **pixgray);
	
	/* motionlow.c */
	IPL_DLL extern void absThreholdWithDifferenceLow(l_uint32 *datad, l_int32 w, l_int32 h, l_int32 wpld, l_uint32 *datas1, l_uint32 *datas2, l_int32 wpls, l_int32 thresh, l_uint32 *datag, l_int32 wplg);
	
	/* noise.c */
	IPL_DLL extern PIX * pixAddNoiseUniform(PIX *pixs, l_float32 min, l_float32 max);
	IPL_DLL extern PIX * pixAddNoiseGaussian(PIX *pixs, l_float32 mean, l_float32 std);
	IPL_DLL extern PIX * pixAddNoiseRayleigh(PIX *pixs, l_float32 sigma, l_int32 var);
	IPL_DLL extern PIX * pixAddNoiseErlang(PIX *pixs, l_int32 k, l_int32 var);
	IPL_DLL extern PIX * pixAddNoiseExponential(PIX *pixs, l_float32 lambda, l_int32 var);
	IPL_DLL extern PIX * pixAddNoiseImpulseBipolar(PIX *pixs, l_float32 d);
	IPL_DLL extern PIX * pixAddNoiseImpulseUnipolar(PIX *pixs, l_int32 noisetype, l_float32 d);
	IPL_DLL extern PIX * pixAddNoiseImpulse(PIX *pixs, l_float32 salt, l_float32 pepper);
	
	/* pixafunc.c */
	IPL_DLL extern PIX * pixInnerSelectBySize(PIX *pixs, l_int32 width, l_int32 height, l_int32 connectivity, l_int32 type, l_int32 relation, l_int32 *pchanged);
	
	/* random.c */
	IPL_DLL extern l_float64 randomUniform();
	IPL_DLL extern l_float64 randomUniformRange(l_float32 min, l_float32 max);
	IPL_DLL extern l_float64 randomGaussian(l_float32 mean, l_float32 std);
	IPL_DLL extern l_float64 randomGaussianThreadSafe(l_float32 mean, l_float32 std);
	IPL_DLL extern l_float64 randomRayleigh(l_float32 sigma);
	IPL_DLL extern l_float64 randomErlang(l_int32 k);
	IPL_DLL extern l_float64 randomExponential(l_float32 lambda);
	IPL_DLL extern l_int8 randomImpulseBipolar(l_float32 white, l_float32 black);
	IPL_DLL extern l_int8 randomImpulseUnipolarWhite(l_float32 d);
	IPL_DLL extern l_int8 randomImpulseUnipolarBlack(l_float32 d);
	
	/* spatial.c */
	IPL_DLL extern l_int32 pixGlobalStats(PIX *pixs, l_float32 *mean, l_float32 *var, l_float32 *std);
	IPL_DLL extern PIX * pixAdaptiveMeanFilter(PIX *pixs, l_int32 wc, l_int32 hc, l_float32 varn);
	
	/* spatiallow.c */
	IPL_DLL extern l_float32 calculateLocalMeanLow(l_uint32 *data, l_int32 wpl, l_int32 width, l_int32 height, l_int32 i, l_int32 j);
	IPL_DLL extern l_float32 calculateLocalVarianceLow(l_uint32 *data, l_int32 wpl, l_int32 width, l_int32 height, l_int32 i, l_int32 j, l_float32  mean);
	
	/* spectra.c - spectrastub.c: require FFTW3 library */
	IPL_DLL extern DPIX * pixGetFourierSpectrumGray(PIX *pixs, l_int32 shiftflag);
	IPL_DLL extern DPIX * pixGetPowerSpectrumGray(PIX *pixs, l_int32 shiftflag);
	IPL_DLL extern DPIX * pixGetPhaseAngleGray(PIX *pixs, l_int32 shiftflag);
	IPL_DLL extern DPIX * pixGetFourierSpectrumRGB(PIX *pixs, l_int32 color, l_int32 shiftflag);
	IPL_DLL extern DPIX * pixGetPowerSpectrumRGB(PIX *pixs, l_int32 color, l_int32 shiftflag);
	IPL_DLL extern DPIX * pixGetPhaseAngleRGB(PIX *pixs, l_int32 color, l_int32 shiftflag);
	IPL_DLL extern l_int32 pixGetDFTSpectraGray(PIX *pixs, DPIX **spectrum, DPIX **phaseAngle, DPIX **powerSpectrum, l_int32 shiftflag);
	IPL_DLL extern l_int32 pixGetDFTSpectraRGB(PIX *pixs, l_int32 color, DPIX **spectrum, DPIX **phaseAngle, DPIX **powerSpectrum, l_int32 shiftflag);
	IPL_DLL extern l_int32 pixPhaseCorrelation(PIX *pixr, PIX *pixs, l_float64 *ppeak, l_int32 *pxloc, l_int32 *pyloc);
	IPL_DLL extern PIX * pixDisplaySpectrum(PIX *pixs, l_float32 scale, l_int32 shiftflag, l_int32 hivals);
	IPL_DLL extern PIX * pixConvertSpectrumToPix(DPIX *dpix, l_float32 scale, l_int32 shiftflag, l_int32 hivals);
	IPL_DLL extern PIX * pixDisplayDFTWithColor(PIX *pixs, l_float32 scale, l_int32 shiftflag);
	IPL_DLL extern PIX * pixRestoreFromSpectraGray(DPIX *spectrum, DPIX *phaseAngle, l_int32 shiftflag);
	
	/* window.c */
	IPL_DLL extern PIX * pixWindowBartlettGlobal(PIX *pixs, l_int32 typeflag);
	IPL_DLL extern PIX * pixWindowHannGlobal(PIX *pixs, l_int32 typeflag);
	IPL_DLL extern PIX * pixWindowTukeyGlobal(PIX *pixs, l_float32 alpha);
	IPL_DLL extern PIX * pixWindowBlackmanGlobal(PIX *pixs, l_float32 alpha);
	IPL_DLL extern PIX * pixWindowBartlettLocal(PIX *pixs, l_int32 typeflag, l_int32 x, l_int32 y, l_int32 r);
	IPL_DLL extern PIX * pixWindowHannLocal(PIX *pixs, l_int32 typeflag, l_int32 x, l_int32 y, l_int32 r);
	IPL_DLL extern PIX * pixWindowTukeyLocal(PIX *pixs, l_float32 alpha, l_int32 x, l_int32 y, l_int32 r);
	IPL_DLL extern PIX * pixWindowBlackmanLocal(PIX *pixs, l_float32 alpha, l_int32 x, l_int32 y, l_int32 r);
	IPL_DLL extern PIX * pixAddPadding(PIX *pixs);
	IPL_DLL extern PIX * pixRemovePadding(PIX *pixs);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* __IPL_H__ */