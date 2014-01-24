""" I want to do lucy-richardson in python!!!
"""

import numpy as np
import scipy
from scipy import stats
from scipy.ndimage.filters import convolve
try:
    from astropy.io import fits
except ImportError:
    import pyfits as fits

import matplotlib.pyplot as plt

import pdb

def pad_psf(psf, spectrum):
    out_psf = np.zeros( spectrum.shape )
    start = len(spectrum)/2 - len(psf)/2
    end = start + len(psf)
    out_psf[start:end] = psf

    return out_psf

def rl_fft(raw_image, psf, niter, k=1, con_var=None):
    """ Implementing the one i got from Jerry

    """

    calc_chisq = lambda a, b, c, d: np.sum((a - b)**2 / (a + c)**2 / (d-1))
    
    conversion =  raw_image.mean() / 10
    raw_image /= conversion
    
    lucy = np.ones(raw_image.shape)
    ratio = k * np.ones(raw_image.shape)
    fft_psf = np.fft.fft(psf)
    
    con_var = sample_noise(raw_image)
    print "using: ", con_var

    norm = np.fft.ifft(np.fft.fft(ratio) * np.conj(fft_psf))
    #index = np.where(norm <= 1E-3 * norm.max())
    #norm[index] = 1
    #raw_image[index] = 0

    fft_conv = fft_psf * np.fft.fft(lucy)
    lucy_conv = np.fft.ifft(fft_conv)

    chisq = calc_chisq(lucy_conv, raw_image, con_var, raw_image.size)
    print "initial Chisq: {}".format(chisq)

    #plt.figure()
    #plt.plot(raw_image)

    for iteration in range(niter):
        ratio = k * (raw_image + con_var) / (lucy_conv + con_var)
        fft_srat = np.fft.fft(ratio) * np.conj(fft_psf)

        lucy *= np.fft.ifft(fft_srat) / norm
        print lucy.max(), lucy.mean(), lucy.min()
        fft_conv = fft_psf * np.fft.fft(lucy)
        lucy_conv = np.fft.ifft(fft_conv)
        size = lucy.size
        #plt.plot(lucy[range(size/2,size)+range(0,size/2)])
        chisq = calc_chisq(lucy_conv, raw_image, con_var, raw_image.size)
        print "Iteration {} Chisq: {}".format(iteration, chisq)

    #pdb.set_trace()
    #raw_input('done')

    #--Why?!
    lucy = lucy[range(size/2,size)+range(0,size/2)]
    return lucy * conversion

def rl_standard(raw_image, psf, niter):
    """ Standerd lucy-richardson convolution

    arXiv 2002 Lauer

    """
   
    psf /= psf.sum()
    psf_inverse = psf[::-1]
    lucy = np.ones( raw_image.shape ) * raw_image.mean()

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

    return np.median( samples )

def rl_damped(raw, psf, niter=2, damped=True, N=3, T=None, multiplier=1):
    """ working on it"""

    #psf /= psf.sum()

    conversion = raw.mean() / 10
    raw /= conversion
    lucy = np.ones(raw.shape) * raw.mean()

    #plt.ion()
    #plt.figure()
    #plt.plot(raw)
    #plt.axhline(y=0, lw=2, color='black')
    for i in xrange(niter):
        if damped:
            print "dampening"
            lucy_temp = convolve( lucy, psf, mode='mirror')
            ratio = dampen(lucy_temp, raw, N, T, multiplier)
        else:
            ratio = raw / convolve(lucy, psf, mode='mirror')

        ratio[ np.isnan(ratio) ] = 0

        top = convolve( ratio, psf, mode='mirror')
        top[ np.isnan(top) ] = 0

        lucy = lucy * (top / psf.sum())
        #plt.plot( lucy )
        print 'iteration', i, lucy.mean(), raw.mean()
        print

    #raw_input('Done')
    return lucy * conversion

def u_factor(lucy, raw_image, T=None, multiplier=1):
    """ Equation 7 
    http://spider.ipac.caltech.edu/staff/fmasci/home/astro_refs/DampledLR94.pdf 

    """
    assert np.all( lucy > 0 ), 'Negative values'

    T = T or multiplier * sample_noise(raw_image)
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

def dampen(lucy, raw, N=3, T=None, multiplier=1):

    first = u_factor(lucy, raw, T, multiplier)**(N-1)
    first[ np.isnan(first) ] = 0
    print first.mean()

    second = (N - (N-1) * u_factor(lucy, raw, T, multiplier))
    second[ np.isnan(second) ] = 0
    print second.mean()

    third = (raw - lucy)/lucy
    third[ np.isnan(third) ] = 0
    print third.mean()

    return 1 + first * second * third
