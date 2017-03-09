#!/usr/bin/python
import argparse
__author__ = 'soham.m'

import numpy as np
from numpy import pi
import STF2FisherMatrix as STF2_FM
from joblib import Parallel, delayed
import sys

parser = argparse.ArgumentParser(description='Code randomize the parameters \
		(m1,m2,chi1,kappa,alpha0,thetaJ,psiJ) and returns a npz file of fisher matrices')
parser.add_argument('-tag','--tag', help='FileTag',required=True)
parser.add_argument('-N','--N', help='N',type=int,required=True)

args = parser.parse_args()
tag =  args.tag
N = args.N

epsilon = 0.001

KAPPA  = np.linspace(-0.5 + epsilon, 1. - epsilon, N)
THETAJ = np.linspace(epsilon, np.pi - + epsilon,   N)

INDEX = []
ZEROS = np.zeros((N, N))

for index, value in np.ndenumerate(ZEROS):
	INDEX.append("%r.%r"%(index[0], index[1]))

wf_params   = {'phi0': 0.001,
              'tC': 1. }

def COMPUTE_FISHER(string):
	IND = string.split(".")
	M, N   = int(IND[0]), int(IND[1])

	wf_params['m1']     = 20
	wf_params['m2']     = 1.4
	wf_params['chi1']   = 0.8
	wf_params['alpha0'] = 0.001
	wf_params['psiJ']   = 0.001

	wf_params['kappa']  = KAPPA[M]
	wf_params['thetaJ'] = THETAJ[N]

	proj_fisher, DET, err_flag = STF2_FM.FisherMatrix(**wf_params)
	return DET

RESULTS = Parallel(n_jobs=-1, verbose=5, backend="threading")(
             map(delayed(COMPUTE_FISHER), INDEX))

RESULTS    = np.array(RESULTS)
FISHER_DET = RESULTS.reshape((N,N))

filename = ('./immediate/FISHER_parallel' + str(tag))

np.savez(filename,
         FISHER_DET = FISHER_DET,
		 OPTIONS    = wf_params)

