""" I want to do lucy-richardson in python!!!
"""

import numpy as np
import scipy
from scipy import stats
from scipy.ndimage.filters import convolve
from astropy.io import fits

import matplotlib.pyplot as plt

import pdb

def pad_psf(psf, spectrum):
    out_psf = np.zeros( spectrum.shape )
    start = len(spectrum)/2 - len(psf)/2
    end = start + len(psf)
    out_psf[start:end] = psf

    return out_psf

def rl_standard(raw_image, psf, niter):
    """ Standerd lucy-richardson convolution

    arXiv 2002 Lauer

    """

    psf_inverse = psf[::-1]
    lucy = np.ones( raw_image.shape )

    for i in xrange( niter ):
        estimate = convolve(lucy, psf, mode='mirror')
        estimate[ np.isnan(estimate) ] = 0

        correction = convolve(raw_image/estimate, psf_inverse, mode='mirror')
        correction[ np.isnan(correction) ] = 0
        print 'Correction:',correction.mean()
        
        lucy *= correction
        print 'Means:', raw_image.mean(), lucy.mean()
        chisq = scipy.nansum((lucy - raw_image)**2 / (lucy)) / (raw_image.size-1)
        print chisq

    return lucy

def sample_noise( spectrum ):
    samples = [spectrum[start:start+300].std() for start in range(0, len(spectrum), 300) ]

    return np.mean( samples )

def rl_damped(raw, psf, niter=2, damped=True):
    """ working on it"""

    psf /= psf.sum()
    lucy = np.ones(raw.shape) * raw.mean()

    #plt.ion()
    #plt.figure()
    #plt.plot(raw)
    #plt.axhline(y=0, lw=2, color='black')
    for i in xrange(niter):
        if damped:
            print "dampening"
            lucy_temp = convolve( lucy, psf, mode='mirror')
            ratio = dampen(lucy_temp, raw)
        else:
            ratio = raw / convolve(lucy, psf, mode='mirror')

        ratio[ np.isnan(ratio) ] = 0

        top = convolve( ratio, psf, mode='mirror')
        top[ np.isnan(top) ] = 0

        lucy = lucy * (top / psf.sum())
        #plt.plot( lucy )
        print 'iteration', i, lucy.mean()
        print

    #raw_input('Done')
    return lucy

def u_factor(lucy, raw_image, T=None):
    """ Equation 7 
    http://spider.ipac.caltech.edu/staff/fmasci/home/astro_refs/DampledLR94.pdf 

    """
    assert np.all( lucy > 0 ), 'Negative values'

    T = T or 3 * sample_noise(raw_image)
    print 'Using {} for T'.format(T)

    first = (-2.0 / T**2)
    ratio = lucy / raw_image
    ratio[ np.isnan(ratio) ] = 0

    logarithm = np.log( ratio )
    logarithm[ np.isnan(logarithm) ] = 0

    second = (raw_image * logarithm - lucy + raw_image)
              
    factor = first * second

    #factor = (-2.0 / T**2) * (raw_image * np.log( lucy / raw_image) - lucy + raw_image)
    factor[ np.isnan(factor) ] = 0
    print 'Factor_pre',np.median(factor[factor>0]), factor[factor>0].min()
    factor = np.where(factor > 1, 1, factor)
    print 'Factor=',np.median([factor>0]), factor[factor>0].min()
    return factor

def dampen(lucy, raw, N=3):

    first = u_factor(lucy, raw)**(N-1)
    first[ np.isnan(first) ] = 0
    print first.mean()

    second = (N - (N-1) * u_factor(lucy, raw))
    second[ np.isnan(second) ] = 0
    print second.mean()

    third = (raw - lucy)/lucy
    third[ np.isnan(third) ] = 0
    print third.mean()

    return 1 + first * second * third
