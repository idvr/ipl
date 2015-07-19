# Current features and expected roadmap of IPL #

Unless otherwise stated, all the following functions and utilities work on both grayscale and color images.

### Rel. 0.1 - Apr 26, 2011 ###

Image processing in the frequency domain:
  * 2D Discrete Fourier Transform and its inverse
  * Ideal, Butterworth and Gaussian lowpass and hipass filters
  * Ideal, Butterworth and Gaussian bandpass and bandreject filters
  * Spectrum, phase angle and power spectrum calculation and display
  * Phase correlation of two images

Spatial filters:
  * Adaptive noise reduction filter

Windowing functions:
  * Zero padding
  * Bartlett window
  * Hann and Hamming windows
  * Tukey window
  * Blackman window

Random distributions and noise generation:
  * Uniform noise
  * Gaussian noise
  * Rayleigh noise
  * Erlang noise
  * Exponential noise
  * Impulse (salt-and-pepper) noise

Corner detectors:
  * Moravec interest operator

Motion analysis:
  * Absolute difference of two images


## Roadmap ##

### Rel. 0.2 - Jun 1, 2011 (tentative) ###

Image processing in the frequency domain:
  * Polar and log-polar images
  * Rotation and scale invariant phase correlation

Spatial filters:
  * Geometric mean filter
  * Harmonic mean filter
  * Contra-harmonic mean filter
  * Midpoint filter
  * Alpha-trimmed mean filter
  * Adaptive median filter
  * Olympic filter

Image sharpening:
  * Fast laplacian
  * Whole body bone scan enhancement

Fuzzy logic transformations:
  * Fuzzy contrast enhancement

### Rel. 0.3 - Jul 1, 2011 (tentative) ###

Image processing in the frequency domain:
  * Homomorphic filtering
  * Notch filters
  * Inverse and Wiener filters
  * Cross-correlation and auto-correlation
  * Laplacian in the frequency domain

Fuzzy logic transformations:
  * Fuzzy image segmentation