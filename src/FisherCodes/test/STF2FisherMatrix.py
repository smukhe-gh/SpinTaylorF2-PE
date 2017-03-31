#!/usr/bin/python
import argparse
__author__ = 'haris.k'

from numpy import load, sqrt, exp, pi, cos, sin, array, arctan2, arccos
from numpy import linalg
import numpy as np
import pycbc
from pycbc.psd import from_txt
import psd_cache
import precessing_wf
import sys

# load psd and scale by dynamic range factor
f_max = 2098.
delta_f = 1/256.
f_inj = 20.
f_low = 30.

wf_options = {'approximant': 'SpinTaylorF2',
              'phase_order': 7,
              'spin_order': 6,
              'amplitude_order': 0 }

psd_choice = 'HPZD'
psd = psd_cache.load_psd(psd_choice, f_max, delta_f)

# Allocate the waveform generator and fisher generator
wf = precessing_wf.waveform(f_inj, f_max, delta_f, **wf_options)
fisher_gen = precessing_wf.fisher(psd, wf, f_low, f_high = None) 

#Parameter List
deriv_lst = ['chi1',
             'kappa',
             'alpha0',
             'thetaJ',
             'psiJ',
             'm1',
             'm2',
             'phi0',
             'tC']

wf_derivs = {'thetaJ': 2.e-5,
             'psiJ': 2.e-5,
             'kappa': 1.e-6,
             'alpha0': 2.e-5,
             'chi1': 1.e-6,
             'm1': 2.e-7,
             'm2': 2.e-7,
             'phi0': 2.e-5,
             'tC': 1.e-6 }

def FisherMatrix(**wf_params):
  err_flag = 0
  for key, val in fisher_gen.test_derivs(wf_params, wf_derivs, deriv_lst).items():
     if val < 0.9999:
       print "Warning, derivative", key, "has accuracy only", val
       err_flag = 1

  fisher_mat = fisher_gen.calc_matrix(wf_params, wf_derivs, deriv_lst)
  Inv_fisher = linalg.inv(fisher_mat)
  proj_fisher= linalg.inv(Inv_fisher[:7,:7]) #Marginalized fisher over tC and phi0
  fisher_det = linalg.det(proj_fisher)
  
  return proj_fisher,fisher_det,err_flag


