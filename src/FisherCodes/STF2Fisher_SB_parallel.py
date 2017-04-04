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
parser.add_argument('-procs','--procs', help='N',type=int, required=True)

args = parser.parse_args()
tag =  args.tag
N = args.N
procs =  args.procs

kappa_vec  = np.random.uniform(low=-1., high=1.,   size=N)
alpha0_vec = np.random.uniform(low=0.,  high=2*pi, size=N)
thetaJ_vec = np.random.uniform(low=0.,  high=pi,   size=N)
psiJ_vec   = np.random.uniform(low=0.,  high=2*pi, size=N)

Fisher_data = np.zeros([N,7,7])
Fisher_det  = np.zeros(N)
Err_flag    = np.zeros(N)
wf_params   = {'phi0': 0.001, 'tC': 1. }

def CHECK_BOUNDARY(m1, m2, chi1, thetaJ, kappa):

  eta = m1*m2/(m1+m2)/(m1+m2)
  K_AC = (145.83*chi1 - 155.92)*eta**2.0 - (1.10*chi1 + 0.16)*eta + 0.08*chi1 + 0.50
  K_CB = (49.98*chi1 - 57.13)*eta**2.0 + (1.96*chi1 - 2.39)*eta - 0.09*chi1 + 0.82

  C_AC = (-6.23*chi1 + 10.47)*eta + 0.01*chi1 + 0.72
  C_CB = (-8.37*chi1 + 10.62)*eta + 0.10*chi1 + 0.30

  RESA = K_AC - C_AC*(thetaJ - np.pi/2)**2.0 - kappa
  RESB = K_CB - C_CB*(thetaJ - np.pi/2)**2.0 - kappa

  if RESA > 0:  # REGION A
      return True
  else:
    return False

def COMPUTE_FISHER(mm):
  wf_params['m1']     = 15.4
  wf_params['m2']     = 1.4
  wf_params['chi1']   = 0.9
  wf_params['kappa']  = kappa_vec[mm]
  wf_params['alpha0'] = alpha0_vec[mm]
  wf_params['thetaJ'] = thetaJ_vec[mm]


  # Check if the point is inside the boundary. Switch to control how points are chosen
  if(1):
    FLAG = CHECK_BOUNDARY(wf_params['m1'], wf_params['m2'], wf_params['chi1'], wf_params['thetaJ'], wf_params['kappa'])
    wf_params['psiJ']   = 0.001
  else:
    FLAG = True
    wf_params['psiJ'] = psiJ_vec[mm]

  # Compute for different sidebands here.
  if FLAG==True:
    wf_params['sideband'] = None
    PROJ_FISHER_None, FISHER_DET_None, ERR_FLAG_NONE = STF2_FM.FisherMatrix(**wf_params)

    wf_params['sideband'] = 2
    PROJ_FISHER_M2, FISHER_DET_M2, ERR_FLAG_M2       = STF2_FM.FisherMatrix(**wf_params)

    wf_params['sideband'] = 0
    PROJ_FISHER_M0, FISHER_DET_M0, ERR_FLAG_M0       = STF2_FM.FisherMatrix( **wf_params)

    return np.array([mm, kappa_vec[mm], thetaJ_vec[mm], FISHER_DET_None, FISHER_DET_M2, FISHER_DET_M0])

  else:
    print "Point outside our region of interest. Setting VALS to zero."
    return np.array([mm, kappa_vec[mm], thetaJ_vec[mm], np.nan, np.nan, np.nan])

INDEX      = np.arange(0, N, 1)
RESULTS    = Parallel(n_jobs=procs, verbose=5)(
             map(delayed(COMPUTE_FISHER), INDEX))

RESULTS    = np.array(RESULTS)
FISHER_DET = RESULTS.reshape((N, 6))

if not os.path.exists("./immediate"):
        os.makedirs("./immediate")

FILENAME = ('./immediate/FISHER_SB_parallel_N%r_%s'%(N, str(tag)))

np.savez(FILENAME,
         FISHER_DET = FISHER_DET,
         KAPPA_VEC  = kappa_vec,
         ALPHA0_VEC = alpha0_vec,
         THETAJ_VEC = thetaJ_vec,
         PSIJ_VEC   = psiJ_vec,
         WF_PARAMS  = wf_params)

