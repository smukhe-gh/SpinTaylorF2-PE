#=======================================================================
# Computes the Fisher Matrix over the parameter space.
# Returns: Fisher Matrix, Det, Err flag over thetaJ, Chi, eta and Kappa
__author__ = 'haris.k'
#=======================================================================

import argparse


import numpy as np
from numpy import pi
import STF2Fisher_compute_fishermatrix as STF2_FM
import sys

parser = argparse.ArgumentParser(description='Code randomize the  parameters (m1,m2,chi1,kappa,alpha0,thetaJ,psiJ) and returns a npz file of fisher matrices')
parser.add_argument('-tag','--tag', help='FileTag',required=True)
parser.add_argument('-N','--N', help='N',type=int,required=True)

args = parser.parse_args()
tag =  args.tag
N = args.N

m1_vec     = np.random.uniform(low=2., high=16., size=N)
m2_vec     = np.random.uniform(low=1., high=3., size=N)
chi1_vec   = np.random.uniform(low=0., high=1., size=N)
kappa_vec  = np.random.uniform(low=-1., high=1., size=N)
alpha0_vec = np.random.uniform(low=0., high=2*pi, size=N)
thetaJ_vec = np.random.uniform(low=0., high=pi, size=N)
psiJ_vec   = np.random.uniform(low=0., high=2*pi, size=N)


Fisher_data = np.zeros([N,7,7])
Fisher_det  = np.zeros(N)
Err_flag    = np.zeros(N)

wf_params   = {'phi0': 0.001, 'tC': 1., 'sideband'; None }

for mm in xrange(N):
  wf_params['m1']     = m1_vec[mm]
  wf_params['m2']     = m2_vec[mm]
  wf_params['chi1']   = chi1_vec[mm]
  wf_params['kappa']  = kappa_vec[mm]
  wf_params['alpha0'] = alpha0_vec[mm]
  wf_params['thetaJ'] = thetaJ_vec[mm]
  wf_params['psiJ']   = psiJ_vec[mm]
  print('m1=%f,m2=%f,chi1=%f,kappa=%f,alpha0=%f,thetaJ=%f,psiJ=%f'%(m1_vec[mm],m2_vec[mm],chi1_vec[mm],kappa_vec[mm],alpha0_vec[mm],thetaJ_vec[mm],psiJ_vec[mm]))
  proj_fisher,det,err_flag = STF2_FM.FisherMatrix(**wf_params)
  Fisher_data[mm,:,:] = proj_fisher
  Fisher_det[mm]= det
  Err_flag[mm] = err_flag
  
filename = ('../output/datasets/fisher_data_' + str(tag))
np.savez(filename,
         Fisher_data= Fisher_data,
         Fisher_det = Fisher_det,
         Err_flag = Err_flag,
         m1_vec = m1_vec,
         m2_vec = m2_vec,
         chi1_vec = chi1_vec,
         kappa_vec = kappa_vec,
         alpha0_vec = alpha0_vec,
         thetaJ_vec = thetaJ_vec,
         psiJ_vec = psiJ_vec) 

