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
 *  random.c
 *
 *   Random number generators
 *      l_float64      randomUniform()
 *      l_float64      randomUniformRange()
 *      l_float64      randomGaussian()
 *      l_float64      randomGaussianThreadSafe()
 *      l_float64      randomRayleigh()
 *      l_float64      randomErlang()
 *      l_float64      randomExponential()
 *      l_int8         randomImpulseBipolar()
 *      l_int8         randomImpulseUnipolarWhite()
 *      l_int8         randomImpulseUnipolarBlack()
 */

#include "ipl.h"
#include <math.h>

/*--------------------------------------------------------------------*
 *                      Random number generators                      *
 *--------------------------------------------------------------------*/
/*!
 *  randomUniform()
 *
 *      Return: random number, in the range [0..1)
 *
 *  Notes:
 *      (1) Generate a random number with uniform distribution.
 */
l_float64
randomUniform()
{
	l_float64 val;
	
	PROCNAME("randomUniform");
	
	val = rand() / (RAND_MAX + 1.);

	return val;
}

/*!
 *  randomUniformRange()
 *
 *      Input:  min, max (lower and upper bounds)
 *      Return: random number, in the range [min..max)
 *
 *  Notes:
 *      (1) Generate a random number with uniform distribution.
 */
l_float64
randomUniformRange(l_float32 min, l_float32 max)
{
	l_float64 val;
	
	PROCNAME("randomUniformRange");
	
	val = min + randomUniform() * (max - min);
	
	return val;
}

/*!
 *  randomGaussian()
 *
 *      Input:  mean (mean value of the gaussian PDF)
 *              std  (standard deviation valie of the PDF)
 *      Return: random number, in the range [0..1)
 *
 *  Notes:
 *      (1) Generate a random number with Gaussian distribution.
 *      (2) The Box–Muller transform is applied to generate pairs
 *          of independent standard normally distributed random
 *          numbers, given a source of uniformly distributed random
 *          numbers.
 *      (3) This method is not thread-safe since its implementation
 *          in C depends on static variables. If you need to ensure
 *          thread-safety, use randomGaussianThreadSafe() instead.
 */
l_float64
randomGaussian(l_float32 mean,
			   l_float32 std)
{
	static l_int8 cached = 0;
	static l_float64 v, t;
	l_float64 u, val, r;
	
	PROCNAME("randomGaussian");
	
	if (cached)
		val = v * t;
	else {
		do {
			u = 2. * randomUniform() - 1.;
			v = 2. * randomUniform() - 1.;
			r = u * u + v * v;
		} while (r >= 1.);
		
		t = sqrt(-2. * log(r) / r);
		val = u * t;
	}
	
	cached = 1 - cached;
	
	return mean + val * std;
}

/*!
 *  randomGaussianThreadSafe()
 *
 *      Input:  mean (mean value of the gaussian PDF)
 *              std  (standard deviation valie of the PDF)
 *      Return: random number, in the range [0..1)
 *
 *  Notes:
 *      (1) Generate a random number with Gaussian distribution.
 *      (2) This method is thread-safe unlike the Box–Muller
 *          transform used by randomGaussian().
 *      (3) Three random numbers, with values in range [-1..1)
 *          and uniform distribution, are added together giving
 *          a normal distribution with mean = 0 and standard
 *          deviation = 1.
 */
l_float64
randomGaussianThreadSafe(l_float32 mean,
						 l_float32 std)
{
	l_float64 a, b, c, val;
	
	PROCNAME("randomGaussianThreadSafe");
	
	a = (randomUniform() * 2 - 1);
	b = (randomUniform() * 2 - 1);
	c = (randomUniform() * 2 - 1);
	val = mean + ((a + b + c) * std);
	
	return val;
}

/*!
 *  randomRayleigh()
 *
 *      Input:  sigma
 *      Return: random number, in the range [0..1)
 *
 *  Notes:
 *      (1) Generate a random number with Rayleigh distribution.
 */
l_float64
randomRayleigh(l_float32 sigma)
{
	PROCNAME("randomRayleigh");
	
	return sigma * sqrt(-1. * log(randomUniform()));
}

/*!
 *  randomErlang()
 *
 *      Input:  k (shape of the Erlang distribution)
 *
 *      Return: random number, in the range [0..1)
 *
 *  Notes:
 *      (1) Generate a random number with Erlang distribution,
 *          shape = @k and scale = 1.
 */
l_float64
randomErlang(l_int32 k)
{
	l_int32   i;
	l_float64 val;
	
	PROCNAME("randomErlang");
	
	val = 0;
	for (i = 0; i < k; i++)
		val += -log(randomUniform());
	
	return val;
}

/*!
 *  randomExponential()
 *
 *      Input:  lambda
 *
 *      Return: random number, in the range [0..1)
 *
 *  Notes:
 *      (1) Generate a random number with exponential distribution.
 */
l_float64
randomExponential(l_float32 lambda)
{
	PROCNAME("randomExponential");
	
	return -lambda * log(randomUniform());
}

/*!
 *  randomImpulseBipolar()
 *
 *      Input:  white (probability of white impulse)
 *              black (probability of black impulse)
 *
 *      Return: random number, with value -1, 0 or 1.
 *
 *  Notes:
 *      (1) Generate a bipolar impulse random number.
 */
l_int8
randomImpulseBipolar(l_float32 white,
					 l_float32 black)
{
	l_float64 val;
	
	PROCNAME("randomImpulseBipolar");
	
	val = randomUniform();
	if (val > 1 - white)
		return 1;
	else if (val < black)
		return -1;
	
	return 0;
}

/*!
 *  randomImpulseUnipolarWhite()
 *
 *      Input:  d (density or probability value)
 *      Return: random number, with value 0 or 1.
 *
 *  Notes:
 *      (1) Generate a white unipolar impulse random number.
 */
l_int8
randomImpulseUnipolarWhite(l_float32 d)
{
	PROCNAME("randomImpulseUnipolarWhite");
	
	return randomImpulseBipolar(d, 0.);
}

/*!
 *  randomImpulseUnipolarBlack()
 *
 *      Input:  d (density or probability value)
 *      Return: random number, with value 0 or -1.
 *
 *  Notes:
 *      (1) Generate a black unipolar impulse random number.
 */
l_int8
randomImpulseUnipolarBlack(l_float32 d)
{
	PROCNAME("randomImpulseUnipolarBlack");
	
	return randomImpulseBipolar(0., d);
}