#===============================================================================
# Computes the Norm, Match and Overlaps.
# SM 9/ 2016
# TODO: Add the possibility of computing match between time domian and
#       fourier domain waveforms.
#===============================================================================

import numpy as np
from numpy import sqrt, sin, cos, pi, exp, arctan2, arccos, array
from lal import MTSUN_SI

import pycbc
from pycbc import waveform as WF
from pycbc.psd import from_txt
from pycbc.filter import sigma, overlap, match
from pycbc.types import TimeSeries, FrequencySeries, zeros
from pycbc.waveform import get_fd_waveform, get_td_waveform

import STF2_psd_cache as psd_cache
import STF2_waveform as wave

def norm(H, psd, f_low, f_cut):
    return pycbc.filter.sigma(H, psd=psd, low_frequency_cutoff=f_low, \
    high_frequency_cutoff=f_cut)

def overlap(H1, H2, psd, f_low, f_cut, normalized=True):
    return pycbc.filter.overlap(H1, H2, psd=psd, low_frequency_cutoff=f_low, \
    high_frequency_cutoff=f_cut, normalized=normalized)

def match(H1, H2, psd, f_low, f_cut):
    return pycbc.filter.match(H1, H2, psd=psd, low_frequency_cutoff=f_low, \
    high_frequency_cutoff=f_cut)


def compute_overlap(**options):
    
    psd_choice = 'HPZD' # TODO: Include this in global dictionary?
    psd = psd_cache.load_psd(psd_choice,options['F_MAX'], options['DEL_F'])

    nsamples = int(options['F_MAX']/options['DEL_F']) + 1
    
    SIDEBAND = [None, 2, 1, 0, -1, -2]
    H = []
    
    for i, band in enumerate(SIDEBAND):
        options['BAND'] = band
        H.append(wave.generate_template(**options))
    
    options['APPROX'] = 'SpinTaylorT2'
    T2   = wave.generate_template(**options)
    
    #TODO: Change these to add any overlaps of your choice.    
    OLVP = np.zeros(9)
    
    #SNR
    OLVP[0]  = norm(T2, psd, options['F_MIN'], options['F_MAX'])  
    OLVP[1]  = norm(H[0], psd, options['F_MIN'], options['F_MAX'])      
    
    #Overlaps
    OLVP[2]  = overlap(T2, H[0], psd, options['F_MIN'], options['F_MAX'])
    OLVP[3]  = overlap(T2, H[1], psd, options['F_MIN'], options['F_MAX'])
    OLVP[4]  = overlap(T2, H[2], psd, options['F_MIN'], options['F_MAX'])
    OLVP[5]  = overlap(T2, H[3], psd, options['F_MIN'], options['F_MAX'])
    OLVP[6]  = overlap(T2, H[4], psd, options['F_MIN'], options['F_MAX'])
    OLVP[7]  = overlap(T2, H[5], psd, options['F_MIN'], options['F_MAX'])
    OLVP[7]  = overlap(T2, H[1] + H[3], psd, options['F_MIN'], options['F_MAX'])

    return OLVP