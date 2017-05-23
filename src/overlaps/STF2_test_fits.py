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

#=======================================================================
# BEGIN CONTROL PANEL
#=======================================================================

options = {

    'M1'     : 14.0,
    'M2'     : 1.4,    
    'CHI1'   : 0.8,

    'KAPPA'  : 0.5,
    'THETAJ' : 1.5708,

    'ALPHA0' : 0.001,
    'PSIJ'   : 0.001,
    'PHI0'   : 0.001,

    'DEL_F'  : 1./256.,
    'F_MIN'  : 20.,
    'F_INJ'  : 20.,
    'F_MAX'  : 2000.,
    'BAND'   : None,

    'APPROX' : 'SpinTaylorF2',
    'OUTPUT_DIR'     : "output-%s" %strftime("%Y_%m_%d_%H_%M_%S", localtime())
    }

#=======================================================================
# END CONTROL PANEL
#=======================================================================

parser = argparse.ArgumentParser(description='Code to check OVLP variation')
parser.add_argument('-tag','--tag', help='FileTag',required=True)
parser.add_argument('-N','--N', help='N',type=int,required=True)
parser.add_argument('-procs','--procs', help='N',type=int,required=True)

args  = parser.parse_args()
tag   = args.tag
N     = args.N
procs = args.procs

# Choose grid size
M1_VEC     = np.linspace(12.0, 100.,  N)
M2_VEC     = np.linspace(1.2,  2.0,   N)
CHI1_VEC   = np.linspace(0.,   1.,    N)

KAPPA_VEC  = np.linspace(-0.5, 1.-0.001,    N)
THETAJ_VEC = np.linspace(0.001,   np.pi-0.001,    N)

ALPHA0_VEC = np.linspace(0.,   2*np.pi,  N)
PSIJ_VEC   = np.linspace(0.,   2*np.pi,  N)
PHI0_VEC   = np.linspace(0.,   2*np.pi,  N)

def generate_OVLP(string):
    # pull out indices to compute at.
    ix, iy = int(string.split(".")[0]), int(string.split(".")[1])

    # compute over spin-precession parameter space.
    options['KAPPA']  = KAPPA_VEC[ix]
    options['THETAJ'] = THETAJ_VEC[iy]

    # call function
    OVLP = compute_overlap(**options)
    return np.array([OVLP[3], OVLP[5]])

INDEX = []
for _i in range(N):
    for _k in range(N):
        INDEX.append("%r.%r"%(_i, _k))

RESULTS = Parallel(n_jobs=procs, verbose=5)(
             map(delayed(generate_OVLP), INDEX))

OVLP20 = np.array(RESULTS)
OVLP_2 = np.reshape(OVLP20[:,0], (N,N))
OVLP_0 = np.reshape(OVLP20[:,1], (N,N))

FILENAME = ('./immediate/OVLP_N%r_%s'%(N, str(tag)))

np.savez(FILENAME,
         OVLP_2     = OVLP_2,
         OVLP_0     = OVLP_0,
         OPTIONS    = options,
         KAPPA_VEC  = KAPPA_VEC,
         THETAJ_VEC = THETAJ_VEC)

