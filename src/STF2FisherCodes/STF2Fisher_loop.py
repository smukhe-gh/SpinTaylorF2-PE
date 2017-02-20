#!/usr/bin/python
import argparse
__author__ = 'soham.m'

import numpy as np
from numpy import pi
import STF2FisherMatrix as STF2_FM
import sys

parser = argparse.ArgumentParser(description='Code randomize the parameters \
		(m1,m2,chi1,kappa,alpha0,thetaJ,psiJ) and returns a npz file of fisher matrices')
parser.add_argument('-tag','--tag', help='FileTag',required=True)
parser.add_argument('-N','--N', help='N',type=int,required=True)

args = parser.parse_args()
tag =  args.tag
N = args.N

epsilon = 0.001

kappa_vec  = np.linspace(-0.5 + epsilon, 1. - epsilon, N)
thetaJ_vec = np.linspace(epsilon, np.pi - + epsilon,   N)

Fisher_det  = np.zeros((N, N))
wf_params   = {'phi0': 0.001, 
              'tC': 1. }
			 
for i, _kappa in enumerate(kappa_vec):
	for j, _thetaJ in enumerate(kappa_vec):
	
	  wf_params['m1']     = 20
	  wf_params['m2']     = 1.4
	  wf_params['chi1']   = 0.8
	  wf_params['kappa']  = _kappa
	  wf_params['alpha0'] = 0.001
	  wf_params['thetaJ'] = _thetaJ
	  wf_params['psiJ']   = 0.001
	
	  print('M1 = %1.2f \t M2 = %1.2f \t CHI1 = %1.2f \t kappa = %1.2f \t alpha0=%1.2f \t thetaJ=%1.2f \t psiJ=%1.2f' \
	  				%(wf_params['m1'], wf_params['m2'], wf_params['chi1'],  wf_params['kappa'], wf_params['alpha0'], wf_params['thetaJ'], wf_params['psiJ']))
  
	  proj_fisher, det, err_flag = STF2_FM.FisherMatrix(**wf_params)
	  Fisher_det[i, j] = det #rows -> kappa, columns -> thetaJ

  
filename = ('./immediate/FISHER_' + str(tag))
np.savez(filename,
         FISHER_DET = Fisher_det,
		 OPTIONS    = wf_params) 

