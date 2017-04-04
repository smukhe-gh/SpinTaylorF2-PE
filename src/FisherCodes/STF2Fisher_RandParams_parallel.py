#!/usr/bin/python
import argparse
__author__ = 'soham.m'

import numpy as np
from numpy import pi
import STF2FisherMatrix as STF2_FM
from joblib import Parallel, delayed
import sys

parser = argparse.ArgumentParser(description='Code randomize the parameters \
								and returns a npz file of fisher matrices')
parser.add_argument('-tag','--tag', help='FileTag',required=True)
parser.add_argument('-N','--N', help='N',type=int,required=True)
parser.add_argument('-procs','--procs', help='N',type=int,required=True)

args  = parser.parse_args()
tag   =  args.tag
N     = args.N
procs = args.procs

if(1):
	# Choose to compute at random points
	print "==> Choosing random points for computation"
	chi1_vec   = np.random.uniform(low=0.,   high=1.,   size=N)
	kappa_vec  = np.random.uniform(low=-0.5, high=1.,   size=N)
	alpha0_vec = np.random.uniform(low=0.,   high=2*pi, size=N)
	thetaJ_vec = np.random.uniform(low=0.,   high=pi,   size=N)
	psiJ_vec   = np.random.uniform(low=0.,   high=2*pi, size=N)
else: 
	# Choose points from Haris's computaion.
	print "==> Loading points from Haris's data."
	DATA = np.load("./immediate/fisher_data_4.npz")
	chi1_vec   = DATA["param1_vec"]
	kappa_vec  = DATA["param2_vec"]
	alpha0_vec = DATA["param3_vec"]
	thetaJ_vec = DATA["param4_vec"]
	psiJ_vec   = DATA["param5_vec"]

wf_params = {'phi0' : 0.001,
             'tC'   : 1. }

def COMPUTE_FISHER(mm):

	wf_params['m1']     = 10.0
	wf_params['m2']     = 1.4

	wf_params['chi1']   = chi1_vec[mm]
	wf_params['kappa']  = kappa_vec[mm]
	wf_params['alpha0'] = alpha0_vec[mm]
	wf_params['thetaJ'] = thetaJ_vec[mm]
	wf_params['psiJ']   = psiJ_vec[mm]

	proj_fisher, DET, err_flag = STF2_FM.FisherMatrix(**wf_params)
	return DET

INDEX      = np.arange(0, len(kappa_vec), 1)
RESULTS    = Parallel(n_jobs=procs, verbose=5)(
             map(delayed(COMPUTE_FISHER), INDEX))

RESULTS    = np.array(RESULTS)
FISHER_DET = RESULTS

FILENAME = ('./immediate/FISHER_SB_None_parallel_N%r_%s'%(N, str(tag)))

np.savez(filename,
         FISHER_DET = FISHER_DET,
		 CHI1_VEC 	= chi1_vec,
		 KAPPA_VEC	= kappa_vec,
		 ALPHA0_VEC = alpha0_vec,
		 THETAJ_VEC = thetaJ_vec,
		 PSIJ_VEC	= psiJ_vec,
		 M1         = wf_params['m1'],
		 M2         = wf_params['m2'])

