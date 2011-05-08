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
 *  spectrastub.c
 *
 *   Stubs for spectra.c functions
 */

#include "ipl.h"

#ifdef HAVE_CONFIG_H
#include "config_auto.h"
#endif  /* HAVE_CONFIG_H */

/* --------------------------------------------*/
#if  !HAVE_LIBFFTW3
/* --------------------------------------------*/

DPIX *
pixGetFourierSpectrumGray(PIX     *pixs,
						  l_int32  shiftflag)
{
	return (DPIX *)ERROR_PTR("function not present", "pixGetFourierSpectrumGray", NULL);
}

DPIX *
pixGetPowerSpectrumGray(PIX     *pixs,
						l_int32  shiftflag)
{
	return (DPIX *)ERROR_PTR("function not present", "pixGetPowerSpectrumGray", NULL);
}

DPIX *
pixGetPhaseAngleGray(PIX     *pixs,
					 l_int32  shiftflag)
{
	return (DPIX *)ERROR_PTR("function not present", "pixGetPhaseAngleGray", NULL);
}

DPIX *
pixGetFourierSpectrumRGB(PIX     *pixs,
						 l_int32  color,
						 l_int32  shiftflag)
{
	return (DPIX *)ERROR_PTR("function not present", "pixGetFourierSpectrumRGB", NULL);
}

DPIX *
pixGetPowerSpectrumRGB(PIX     *pixs,
					   l_int32  color,
					   l_int32  shiftflag)
{
	return (DPIX *)ERROR_PTR("function not present", "pixGetPowerSpectrumRGB", NULL);
}

DPIX *
pixGetPhaseAngleRGB(PIX     *pixs,
					l_int32  color,
					l_int32  shiftflag)
{
	return (DPIX *)ERROR_PTR("function not present", "pixGetPhaseAngleRGB", NULL);
}

l_int32
pixGetDFTSpectraGray(PIX      *pixs,
					 DPIX    **spectrum,
					 DPIX    **phaseAngle,
					 DPIX    **powerSpectrum,
					 l_int32   shiftflag)
{
	return ERROR_INT("function not present", "pixGetDFTSpectraGray", 1);
}

l_int32
pixGetDFTSpectraRGB(PIX      *pixs,
					l_int32   color,
					DPIX    **spectrum,
					DPIX    **phaseAngle,
					DPIX    **powerSpectrum,
					l_int32   shiftflag)
{
	return ERROR_INT("function not present", "pixGetDFTSpectraRGB", 1);
}

l_int32
pixPhaseCorrelation(PIX       *pixr,
					PIX       *pixs,
					l_float64 *ppeak,
					l_int32   *pxloc,
					l_int32   *pyloc)
{
	return ERROR_INT("function not present", "pixPhaseCorrelation", 1);
}

PIX *
pixDisplaySpectrum(PIX       *pixs,
				   l_float32  scale,
				   l_int32    shiftflag,
				   l_int32    hivals)
{
	return (PIX *)ERROR_PTR("function not present", "pixGetFourierSpectrumGray", NULL);
}

PIX *
pixConvertSpectrumToPix(DPIX      *dpix,
						l_float32  scale,
						l_int32    shiftflag,
						l_int32    hivals)
{
	return (PIX *)ERROR_PTR("function not present", "pixConvertSpectrumToPix", NULL);
}

PIX *
pixDisplayDFTWithColor(PIX       *pixs,
					   l_float32  scale,
					   l_int32    shiftflag)
{
	return (PIX *)ERROR_PTR("function not present", "pixDisplayDFTWithColor", NULL);
}

PIX *
pixRestoreFromSpectraGray(DPIX    *spectrum,
						  DPIX    *phaseAngle,
						  l_int32  shiftflag)
{
	return (PIX *)ERROR_PTR("function not present", "pixRestoreFromSpectraGray", NULL);
}

/* --------------------------------------------*/
#endif  /* !HAVE_LIBFFTW3 */
/* --------------------------------------------*/