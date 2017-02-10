#=======================================================================
# Computes the Fisher Matrix over the parameter space.
# Returns: Fisher Matrix, Det, Err flag over thetaJ, Chi, eta and Kappa
# SM 17/ 16
#=======================================================================

import numpy as np
import os
import STF2Fisher_compute_matrix as CM
from time import localtime, strftime
from joblib import Parallel, delayed

#=======================================================================
# BEGIN CONTROL PANEL
#=======================================================================

options = {
    'sideband': None,
    'thetaJ'  : 0.01,
    'psiJ'    : 0.01,
    'kappa'   : 0.5,
    'alpha0'  : 0.001,
    'chi1'    : 0.8,
    'm1'      : 12.,
    'm2'      : 1.4,
    'phi0'    : 0.01,
    'tC': 1.,
}


N = 10  # number of points in thetaJ kappa space
M = 2  # number of points in eta, chi1 space 

#=======================================================================
# END CONTROL PANEL
#=======================================================================

def M1_M2_to_MCHRIP_ETA(mass1, mass2):
    eta = mass1*mass2/(mass1+mass2)/(mass1+mass2)
    mc = (mass1+mass2) * pow(eta, 3./5.)
    return (mc, eta)

def generate_GRID(**options):

    DET = np.zeros((N, N))
    MCHIRP, ETA = M1_M2_to_MCHRIP_ETA(options['m1'], options['m2'])

    V_THETAJ = np.linspace(0.2, np.pi - 0.2, N)
    V_KAPPA  = np.linspace(-0.5, 0.8, N)

    for _thetaJ in xrange(N):
        for _kappa in xrange(N):

            options['thetaJ'] = V_THETAJ[_thetaJ]
            options['kappa']  = V_KAPPA[_kappa]
            
            print  "thetaJ: %1.2f \t kappa: %1.2f \t M1: %1.2f \t chi1: %1.2f" %(options['thetaJ'], options['kappa'], options['m1'], options['chi1'])
            print 60*"-"

            fisher_matrix,fisher_det, Err_flag = CM.FisherMatrix(**options)
            DET[_thetaJ][_kappa] = fisher_det


    save_params = {
    'OUTPUT_DIR': "output-%s" %strftime("%Y_%m_%d_%H_%M_%S", localtime()),
    'FILE_NAME' : "fisherdet_eta_%s_chi1_%s_N_%r.npz" %('{:.2f}'.format(ETA),\
    '{:.2f}'.format(options['chi1']), N)}

    if not os.path.exists("../output/datasets/%s" %save_params['OUTPUT_DIR']):
        os.makedirs("../output/datasets/%s" %save_params['OUTPUT_DIR'])

    np.savez("../output/datasets/%s/%s" %(save_params["OUTPUT_DIR"], save_params["FILE_NAME"]),
        DATE         = strftime("%Y-%m-%d %H:%M:%S", localtime()),
        FDET         = DET,
        THETAJ       = V_THETAJ,
        KAPPA        = V_KAPPA,
        CHI1         = options['chi1'],
        ETA          = ETA,
        OPTIONS      = options)

    return None

def parallel_GRID(_MASS, _CHI1, **options):

    options['m1']         = _MASS
    options['chi1']       = _CHI1
    generate_GRID(**options)

V_MASS1   = np.linspace(12, 20, M)
V_CHI1    = np.linspace(0.5, 0.9, M)

Parallel(n_jobs=-2, verbose=5)(delayed(parallel_GRID)(_MASS = V_MASS1[m], \
 _CHI1 = V_CHI1[c], **options) for m in xrange(M) for c in xrange(M))

