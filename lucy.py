""" I want to do lucy-richardson in python!!!
"""

import numpy as np
import scipy
from scipy import stats
from scipy.ndimage.filters import convolve
from astropy.io import fits


def convolute(input,input2):
    fourier_space_convolution = (np.fft.rfft2(input,input.shape))*(np.fft.rfft2(input2,input2.shape))
    convert_realspace_convolution = np.fft.fftshift(np.fft.irfft2(fourier_space_convolution,fourier_space_convolution.shape))

    return convert_realspace_convolution


def deconvolve( raw_image, psf, niter=20):
    """ do it """

    # Initial estimate is arbitrary
    lucy_image = np.ones( raw_image.shape ) * .5
    psf /= psf.max()
    psf_inverse = psf[::-1]

    for i in range( niter ):
        tmp_image = raw_image / convolve(lucy_image, psf, mode='mirror')
        tmp_image[ np.isnan(tmp_image) ] = 0
        lucy_image *= convolve(tmp_image, psf_inverse, mode='mirror')
        tmp_image[ np.isnan(tmp_image) ] = 0
        chisq = scipy.nansum((lucy_image - raw_image)**2 / (lucy_image)) / (raw_image.size-1)

    return lucy_image
