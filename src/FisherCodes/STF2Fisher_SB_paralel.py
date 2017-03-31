#======================================================================
# Soham M 3/2017
# Code to compute Fisher determinant for 
#   > randomly sampled points in 7D parameter space 
#     (maximized over tC and phi)
#   > For different sidebands (None, 2, and 0)
#   > Can use more than one processes
#======================================================================

import argparse
__author__ = 'soham.m'

import numpy as np
from numpy import pi
import STF2FisherMatrix as STF2_FM
from joblib import Parallel, delayed
import sys
import os

parser = argparse.ArgumentParser(description='Code to compute Fisher Matrix at Random points in Parallel')
parser.add_argument('-tag','--tag', help='FileTag',required=True)
parser.add_argument('-N','--N', help='N',type=int,required=True)

args = parser.parse_args()
tag =  args.tag
N = args.N

kappa_vec  = np.random.uniform(low=-1., high=1.,   size=N)
alpha0_vec = np.random.uniform(low=0.,  high=2*pi, size=N)
thetaJ_vec = np.random.uniform(low=0.,  high=pi,   size=N)
psiJ_vec   = np.random.uniform(low=0.,  high=2*pi, size=N)

Fisher_data = np.zeros([N,7,7])
Fisher_det  = np.zeros(N)
Err_flag    = np.zeros(N)
wf_params   = {'phi0': 0.001, 'tC': 1. }

def COMPUTE_FISHER(mm):
  wf_params['m1']     = 15.4
  wf_params['m2']     = 1.4
  wf_params['chi1']   = 0.9
  wf_params['kappa']  = kappa_vec[mm]
  wf_params['alpha0'] = alpha0_vec[mm]
  wf_params['thetaJ'] = thetaJ_vec[mm]
  wf_params['psiJ']   = psiJ_vec[mm]

  # Compute for different sidebands here. 
  # FIXME: This is a very crude hacky way to do this. Fix this.

  wf_params['sideband'] = None
  print 80*"="
  print "Computing Fisher DET for sideband: %r \t Iteration: %r" %(wf_params['sideband'], mm)
  print 80*"="
  PROJ_FISHER_None, FISHER_DET_None, ERR_FLAG_NONE = STF2_FM.FisherMatrix(**wf_params)
  
  wf_params['sideband'] = 2
  print 80*"="
  print "Computing Fisher DET for sideband: %r \t Iteration: %r" %(wf_params['sideband'], mm)
  print 80*"="
  PROJ_FISHER_M2, FISHER_DET_M2, ERR_FLAG_M2       = STF2_FM.FisherMatrix(**wf_params)
  
  wf_params['sideband'] = 0
  print 80*"="
  print "Computing Fisher DET for sideband: %r \t Iteration: %r" %(wf_params['sideband'], mm)
  print 80*"="
  PROJ_FISHER_M0, FISHER_DET_M0, ERR_FLAG_M0       = STF2_FM.FisherMatrix( **wf_params)

  return np.array([mm, FISHER_DET_None, FISHER_DET_M2, FISHER_DET_M0])
  
INDEX      = np.arange(0, N, 1)
RESULTS    = Parallel(n_jobs=2, verbose=5)(
             map(delayed(COMPUTE_FISHER), INDEX))

RESULTS    = np.array(RESULTS)
FISHER_DET = RESULTS.reshape((N,4))

if not os.path.exists("./immediate"):
        os.makedirs("./immediate")

FILENAME = ('./immediate/FISHER_SB_parallel_N%r_%s'%(N, str(tag)))

np.savez(FILENAME,
         FISHER_DET = FISHER_DET,
         KAPPA_VEC  = kappa_vec,
         THETAJ_VEC = thetaJ_vec)

