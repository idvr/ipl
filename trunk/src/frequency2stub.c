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
 *  frequency2stub.c
 *
 *   Stubs for frequency2.c functions
 */

#include "ipl.h"

#ifdef HAVE_CONFIG_H
#include "config_auto.h"
#endif  /* HAVE_CONFIG_H */

/* --------------------------------------------*/
#if  !HAVE_LIBFFTW3
/* --------------------------------------------*/

PIX *
pixIdealFilter(PIX        *pixs,
			   l_float32   radius,
			   l_int32     typeflag,
			   l_int32     outflag,
			   DPIX      **dpixf)
{
	return (PIX *)ERROR_PTR("function not present", "pixIdealFilter", NULL);
}

PIX *
pixButterworthFilter(PIX        *pixs,
					 l_float32   radius,
					 l_int32     order,
					 l_int32     typeflag,
					 l_int32     outflag,
					 DPIX      **dpixf)
{
	return (PIX *)ERROR_PTR("function not present", "pixButterworthFilter", NULL);
}

PIX *
pixGaussianFilter(PIX        *pixs,
				  l_float32   radius,
				  l_int32     typeflag,
				  l_int32     outflag,
				  DPIX       **dpixf)
{
	return (PIX *)ERROR_PTR("function not present", "pixGaussianFilter", NULL);
}

PIX *
pixIdealSelectiveFilter(PIX        *pixs,
						l_float32   radius,
						l_float32   bandwidth,
						l_int32     typeflag,
						l_int32     outflag,
						DPIX      **dpixf)
{
	return (PIX *)ERROR_PTR("function not present", "pixIdealSelectiveFilter", NULL);
}

PIX *
pixButterworthSelectiveFilter(PIX        *pixs,
							  l_float32   radius,
							  l_float32   bandwidth,
							  l_int32     order,
							  l_int32     typeflag,
							  l_int32     outflag,
							  DPIX      **dpixf)
{
	return (PIX *)ERROR_PTR("function not present", "pixButterworthSelectiveFilter", NULL);
}

PIX *
pixGaussianSelectiveFilter(PIX        *pixs,
						   l_float32   radius,
						   l_float32   bandwidth,
						   l_int32     typeflag,
						   l_int32     outflag,
						   DPIX      **dpixf)
{
	return (PIX *)ERROR_PTR("function not present", "pixGaussianSelectiveFilter", NULL);
}

PIX *
pixApplyFilter(PIX     *pixs,
			   DPIX    *dpix,
			   l_int32  outflag)
{
	return (PIX *)ERROR_PTR("function not present", "pixApplyFilter", NULL);
}

PIX *
pixApplyFilterGray(PIX     *pixs,
				   DPIX    *dpix,
				   l_int32  outflag)
{
	return (PIX *)ERROR_PTR("function not present", "pixApplyFilterGray", NULL);
}

/* --------------------------------------------*/
#endif  /* !HAVE_LIBFFTW3 */
/* --------------------------------------------*/