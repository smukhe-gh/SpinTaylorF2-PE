#=======================================================================
# Computes the OLVPs over the parameter space.
# Returns: Overlaps over thetaJ, Chi, eta and Kappa
# SM 17/ 16
#=======================================================================

import numpy as np
from STF2_overlaps import compute_overlap
from time import localtime, strftime
from joblib import Parallel, delayed
import os
import argparse

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

#=======================================================================
# BEGIN CONTROL PANEL
#=======================================================================

options = {
    'APPROX' : 'SpinTaylorF2',
    'DEL_F'  : 1./256.,
    'F_MIN'  : 20.,
    'F_INJ'  : 20.,
    'F_MAX'  : 2000.,
    'BAND'   : None,
    }

#=======================================================================
# END CONTROL PANEL
#=======================================================================


parser = argparse.ArgumentParser(description='Code to check SNR distribution')
parser.add_argument('-tag','--tag', help='FileTag',required=True)
parser.add_argument('-N','--N', help='N',type=int,required=True)
parser.add_argument('-procs','--procs', help='N',type=int,required=True)

args  = parser.parse_args()
tag   = args.tag
N     = args.N
procs = args.procs

# Choose to compute at random points
print "==> Choosing random points for computation"
m1_vec     = np.random.uniform(low=12.0, high=100.,  size=N)
m2_vec     = np.random.uniform(low=1.2,  high=2.0,   size=N)
chi1_vec   = np.random.uniform(low=0.,   high=1.,    size=N)
kappa_vec  = np.random.uniform(low=-0.5, high=1.,    size=N)
alpha0_vec = np.random.uniform(low=0.,   high=2*pi,  size=N)
thetaJ_vec = np.random.uniform(low=0.,   high=pi,    size=N)
psiJ_vec   = np.random.uniform(low=0.,   high=2*pi,  size=N)
phi0_vec   = np.random.uniform(low=0.,   high=2*pi,  size=N)

def norm(H, psd, f_low, f_cut):
    return pycbc.filter.sigma(H, psd=psd, low_frequency_cutoff=f_low, \
    high_frequency_cutoff=f_cut)

def compute_SNR(**options):
    psd_choice = 'HPZD'
    psd = psd_cache.load_psd(psd_choice,options['F_MAX'], options['DEL_F'])
    nsamples = int(options['F_MAX']/options['DEL_F']) + 1
    SNR  = norm(wave.generate_template(**options), psd, \
        options['F_MIN'], options['F_MAX'])
    return SNR

def generate_SNR(mm):
    options['M1']       = m1_vec[mm]
    options['M2']       = 1.4    #m2_vec[mm]
    options['CHI1']     = 0.001  #chi1_vec[mm]
    options['KAPPA']    = 1.000  #kappa_vec[mm]
    options['ALPHA0']   = 0.001  #alpha0_vec[mm]
    options['THETAJ']   = 1.570  #thetaJ_vec[mm]
    options['PSIJ']     = 0.001  #psiJ_vec[mm]
    options['PHI0']     = 0.001  #phi0_vec[mm]
    options['SIDEBAND'] = None

    SNR = compute_SNR(**options)

    return SNR

INDEX      = np.arange(0, len(kappa_vec), 1)
RESULTS    = Parallel(n_jobs=procs, verbose=5)(
             map(delayed(generate_SNR), INDEX))

RESULTS    = np.array(RESULTS)

FILENAME = ('./immediate/SNR_distribution_N%r_%s'%(N, str(tag)))

np.savez(FILENAME,
         SNR        = RESULTS,
         M1_VEC     = m1_vec,
         M2_VEC     = m2_vec,
         CHI1_VEC   = chi1_vec,
         KAPPA_VEC  = kappa_vec,
         ALPHA0_VEC = alpha0_vec,
         THETAJ_VEC = thetaJ_vec,
         PSIJ_VEC   = psiJ_vec)
