# Fingerprint enhancement and phase correlation #

This tutorial is intended to get you started with IPL and to quickly showcase the functionalities provided for image enhancement and phase correlation.

In image processing, [phase correlation](http://en.wikipedia.org/wiki/Phase_correlation) is a method of image registration, and uses a fast frequency-domain approach to estimate the relative translative offset between two similar images.

For example, we can make use of phase correlation to find out if two fingerprint images are actually the same. In fact, if we acquire a reference image of a person's fingerprint, we should be able to quickly discover the degree of similarity of subsequent scans of that same fingerprint, provided that the person is likely to touch the fingerprint sensor every time in a slightly different position.

## Highpass image filtering ##

For the sake of simplicity, let's assume in this demo that we have acquired the following reference fingeprint image:

![http://ipl.googlecode.com/svn/images/fingerprints/original.png](http://ipl.googlecode.com/svn/images/fingerprints/original.png)

_(Original image courtesy of the U.S. National Institute of Standards and Technology)_

As a preliminary step, we can enhance the image details by reducing the effect of smudges and highlight the fingerprint ridges. One way to achieve this is by operating in the frequency domain and apply a [Butterworth](http://en.wikipedia.org/wiki/Butterworth_filter) highpass filter to the image.

```
/* 1. Read, pad and display fingerprint pix */
pixt = pixRead("images/original.png");
pixs = pixAddPadding(pixt);
pixWrite("padded.png", pixs, IFF_PNG); 
pixDestroy(&pixt);

/* 2. Display spectrum of padded pix */
pixd = pixDisplaySpectrum(pixs, 1, L_WITH_SHIFTING, L_HI_TO_WHITE);
pixWrite("spectrum.png", pixs, IFF_PNG); 
pixDestroy(&pixd);

/* 3. Filter via hipass Butterworth filter, and display filter */
pixt = pixButterworthFilter(pixs, 50, 4, L_HI_PASS, L_CLIP_TO_ZERO, &dpix);
pixd = pixDisplayFilter(dpix);
pixWrite("filter.png", pixd, IFF_PNG); 
dpixDestroy(&dpix);
pixDestroy(&pixd);
pixDestroy(&pixs);

/* 4. Remove pix padding and display filtered pix */
pixs = pixRemovePadding(pixt);
pixDestroy(&pixt);
pixWrite("hipass.png", pixs, IFF_PNG); 
pixDestroy(&pixs);
```

The code above will generate and display the following images:

![http://ipl.googlecode.com/svn/images/fingerprints/padded.png](http://ipl.googlecode.com/svn/images/fingerprints/padded.png)
![http://ipl.googlecode.com/svn/images/fingerprints/spectrum.png](http://ipl.googlecode.com/svn/images/fingerprints/spectrum.png)
![http://ipl.googlecode.com/svn/images/fingerprints/filter.png](http://ipl.googlecode.com/svn/images/fingerprints/filter.png)
![http://ipl.googlecode.com/svn/images/fingerprints/hipass.png](http://ipl.googlecode.com/svn/images/fingerprints/hipass.png)

The first step is to pad the reference image with zeroes at its right and bottom, thus doubling its size. This is necessary if we are to avoid inaccuracies in the calculation of the Fourier spectrum, since the Fourier transform assumes that we deal with periodic images (we could avoid doubling the size of the image with padding by using one of the IPL windowing functions, which are the subject of a different tutorial).

Then we display the spectrum and filter the image by a Butterworth highpass filter of order 4 with a cutoff frequency of 50. As expected, the filter sharpens the image but it also darkens the gray tones since the DC term (shifted in the middle of the spectrum) is reduced to zero.

We also show the spatial representation of the Butterworth filter applied as well as the resulting filtered image of the fingerprint (after padding has been removed).

In this case, it is useful to enhance details of interest by thresholding the filtered image. If in our code above we replace the statement:

`pixt = pixButterworthFilter(pixs, 50, 4, L_HI_PASS, L_CLIP_TO_ZERO, &dpix);`

with:

`pixt = pixButterworthFilter(pixs, 50, 4, L_HI_PASS, L_THRESH_NEG_TO_BLACK, &dpix);`

we set to black all negative values and to white all positive values in the filtered image.

You can see the results in the image below:

![http://ipl.googlecode.com/svn/images/fingerprints/original.png](http://ipl.googlecode.com/svn/images/fingerprints/original.png)
![http://ipl.googlecode.com/svn/images/fingerprints/threshold.png](http://ipl.googlecode.com/svn/images/fingerprints/threshold.png)

Compare the thresholded image to the reference one and observe how ridges are clearer and the effect of smudges are reduced considerably.

## Image alignment via phase correlation ##

Now that we have a clearer picture of the reference fingerprint image, we can examine it to identify its _minutiae_ and other features, thus determining the region of interest (ROI) that uniquely distinguishes that fingerprint from others belonging to different persons.

For simplicity, we assume that the ROI has been determined by the square with the top-left corner at coordinates 338, 312 (on the _x_ and _y_ planes, respectively), of width 160 pixels.

We can display a red ROI box with the following code:

```
box = boxCreate(338, 312, 160, 160);
pixRenderBoxArb(pixs, box, 7, 0xFF, 0x00, 0x00);
pixWrite("hipass.png", pixs, IFF_PNG); 
boxDestroy(&box);
```

If the image is grayscale, we can convert it from 8-bit to 32-bit depth by a call to the `pixConvertTo32()` function.

![http://ipl.googlecode.com/svn/images/fingerprints/roi.png](http://ipl.googlecode.com/svn/images/fingerprints/roi.png)

We can use the phase correlation method of registration to compare and align our reference image to another image of the same fingerprint, characterized by a translative offset.

The images below show the reference fingerprint and newly acquired one:

![http://ipl.googlecode.com/svn/images/fingerprints/original.png](http://ipl.googlecode.com/svn/images/fingerprints/original.png)
![http://ipl.googlecode.com/svn/images/fingerprints/scan.png](http://ipl.googlecode.com/svn/images/fingerprints/scan.png)

We filter and threshold the new fingerprint image following the same process described above. Then we compare the two images by phase correlation using the following code:

```
pixPhaseCorrelation(pixfp1, pixfp2, &peak, &x, &y);
printf("phase correlation: peak = %f at %d, %d\n", peak, x, y);
```

where `pixfp1` is the reference image and `pixfp2` is the newly acquired one (both images have been filtered and thresholded).

IPL is going to output something like:

`phase correlation: peak = 0.172683 at -226, -162`

indicating that the ROI in the new image is shifted 226 pixels on the right and 162 pixels on the top, compared to the reference fingerprint. This offset can be used to align the two images.

In fact, we can visualize the results by displaying a box around the ROI as determined in the new image:

```
box = boxCreate(338 + x, 312 + y, 160, 160);
pixRenderBoxArb(pixfp2, box, 7, 0xFF, 0x00, 0x00);
pixWrite("correl.png", pixfp2, IFF_PNG); 
boxDestroy(&box);
```

thus obtaining:

![http://ipl.googlecode.com/svn/images/fingerprints/roi.png](http://ipl.googlecode.com/svn/images/fingerprints/roi.png)
![http://ipl.googlecode.com/svn/images/fingerprints/correl.png](http://ipl.googlecode.com/svn/images/fingerprints/correl.png)

## Conclusions ##

In this tutorial we have shown how to filter images in the frequency domain (using, in our particular example, a Butterworth highpass filter for sharpening and thresholding) and align two images using the phase correlation registration method.

Although the example was kept brief for reasons of simplicity, it should be sufficient to give you a flavor of the functionalities provided by IPL and to get you started with its usage.