# IPL #
IPL (Image Processing Library) is a collection of C routines for digital image processing, providing algorithms for image filtering, enhancement and segmentation, motion analysis and features detection.

IPL is developed and maintained by Andrea Gagliardi La Gala.


## Applicability ##
IPL's objective is to support research in a variety of domains and industries, including:
  * Medical Imaging (enhancement of Gamma-Ray, X-Ray, Magnetic Resonances and ultrasound images, such as bone and CT scans)
  * Astronomical Imaging
  * Forensic Imaging (enhancement of fingerprints and their ridges)
  * Visual inspection of manufactured goods
  * Surveillance systems

_Disclaimer: IPL is not approved or endorsed by any Organizations operating in the industries mentioned above, not even for clinical or diagnostic usage, and should not be relied upon for any diagnostic decisions. It is provided 'as is' and any risks associated with its misuse or lack of quality lie completely with the user. The only purpose of IPL is to support and enhance research in the enlisted areas._


## Licensing ##
IPL is released under the terms of the GNU General Public License version 3, which lets you freely redistribute it and/or modify it, as long as you preserve the credit information and license your derived work under the identical terms.

Contact us at andrea _dot_ lagala _at_ gmail _dot_ com if you need to explore commercial licensing options.


## External Libraries ##
IPL depends on and extends:
  * _Leptonica_ library by Dan S. Bloomberg (http://www.leptonica.com/)

and (optionally) integrates:
  * _FFTW_ library by Matteo Frigo and Steven G. Johnson (http://www.fftw.org/)
  * _fuzzylite_ C++ library by Juan Rada-Vilela (http://code.google.com/p/fuzzy-lite/)

to achieve fast and efficient image filtering and manipulation within the frequency domain (using the Discrete Fourier Transform and its inverse) and by using Fuzzy Logic (based on the theory of Fuzzy Sets first devised by Lotfi Zadeh).


## Build Instructions ##
IPL can be compiled by any modern C99 compiler and runs on all major UNIX-based platforms. To build it:
  * install _Leptonica 1.68_ (see instructions on [Leptonica web pages](http://www.leptonica.com/))
  * optionally install _FFTW 3.2.2_ (see instructions on  [FFTW web pages](http://www.fftw.org/))
  * optionally install _fuzzylite 1.03_ (see instructions on [fuzzylite web pages](http://code.google.com/p/fuzzy-lite/))
  * type in console:
```
./configure [--prefix install-dir]
make
make install
```


## Tutorials ##
Watch out for the IPL tutorials published on the [Wiki pages](http://code.google.com/p/ipl/w/list).