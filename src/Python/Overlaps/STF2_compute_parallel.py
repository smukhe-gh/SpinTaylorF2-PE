#=======================================================================
# Computes the OLVPs over the parameter space.
# Returns: Overlaps over thetaJ, Chi, eta and Kappa
# SM 17/ 16
# TODO: Add support for SpinTaylorF2
#=======================================================================

import numpy as np
from STF2_overlaps import compute_overlap
from numpy import sqrt
from time import localtime, strftime
from joblib import Parallel, delayed
import STF2_vis_overlaps as vs
import STF2_vis_grid as vsg
import os

#=======================================================================
# BEGIN CONTROL PANEL
#=======================================================================

options = {

    'ALPHA0' : 0.001,
    'KAPPA'  : 0.5,
    'CHI1'   : 0.5,
    'PSIJ'   : 0.001,
    'M1'     : 10.0,
    'PHI0'   : 0.001,
    'M2'     : 1.4,
    'THETAJ' : 0.7853981633974483,
    'APPROX' : 'SpinTaylorF2',
    'DEL_F'  : 1./256.,
    'F_MIN'  : 20.,
    'F_INJ'  : 20.,
    'F_MAX'  : 2000.,
    'BAND'   : None,

    'N'      : 2,
    'M'      : 2,

    'V_MASS1_RANGE' : [2.4, 50.0],
    'V_CHI1_RANGE'  : [0.5, 1.00],

    'V_KAPPA_RANGE'  : [-0.500, 0.999],
    'V_THETAJ_RANGE' : [0.001, 3.14],

    'OUTPUT_DIR'     : "output-%s" %strftime("%Y_%m_%d_%H_%M_%S", localtime()),
    'GENERATE_PLOTS' : 1}

#=======================================================================
# END CONTROL PANEL
#=======================================================================
def M1_M2_to_MCHRIP_ETA(mass1, mass2):
    eta = mass1*mass2/(mass1+mass2)/(mass1+mass2)
    mc = (mass1+mass2) * pow(eta, 3./5.)
    return (mc, eta)

def generate_GRID(**options):

    OVLP        = np.zeros((options['N'], options['N'], 9))
    MCHIRP, ETA = M1_M2_to_MCHRIP_ETA(options['M1'], options['M2'])

    V_THETAJ = np.linspace(options['V_THETAJ_RANGE'][0], \
    options['V_THETAJ_RANGE'][1], options['N'])
    V_KAPPA  = np.linspace(options['V_KAPPA_RANGE'][0],  \
    options['V_KAPPA_RANGE'][1],  options['N'])

    for _thetaJ in xrange(options['N']):
        for _kappa in xrange(options['N']):

            options['THETAJ'] = V_THETAJ[_thetaJ]
            options['KAPPA']  = V_KAPPA[_kappa]

            OVLP[_thetaJ][_kappa][:] = compute_overlap(**options)

    filename = "overlaps_eta_%s_chi1_%s_N_%r.npz" %('{:.2f}'.format(ETA),\
    '{:.2f}'.format(options['CHI1']), options['N'])

    np.savez("../../../output/datasets/%s/%s" %(options['OUTPUT_DIR'], filename),
        DATE         = strftime("%Y-%m-%d %H:%M:%S", localtime()),
        SNR_0F       = OVLP[:, :, 0],
        SNR_02       = OVLP[:, :, 1],
        SNR_00       = OVLP[:, :, 2],
        OLVP_0F_P2   = OVLP[:, :, 3],
        OLVP_0F_P1   = OVLP[:, :, 4],
        OLVP_0F_P0   = OVLP[:, :, 5],
        OLVP_0F_M1   = OVLP[:, :, 6],
        OLVP_0F_M2   = OVLP[:, :, 7],
        OLVP_0F_P2P0 = OVLP[:, :, 8],
        THETAJ       = V_THETAJ,
        KAPPA        = V_KAPPA,
        CHI1         = options['CHI1'],
        ETA          = ETA,
        OPTIONS      = options)

    return None

def parallel_GRID(_MASS, _CHI1, **options):

    if not os.path.exists("../../../output/datasets/%s" %options['OUTPUT_DIR']):
        os.makedirs("../../../output/datasets/%s" %options['OUTPUT_DIR'])

    options['M1']         = _MASS
    options['CHI1']       = _CHI1

    generate_GRID(**options)

V_MASS1   = np.linspace(options['V_MASS1_RANGE'][0], \
options['V_MASS1_RANGE'][1], options['M'])
V_CHI1    = np.linspace(options['V_CHI1_RANGE'][0], \
options['V_CHI1_RANGE'][1], options['M'])

Parallel(n_jobs=-1, verbose=5)(delayed(parallel_GRID)(_MASS = V_MASS1[m], \
 _CHI1 = V_CHI1[c], **options) for m in xrange(options['M']) for \
 c in xrange(options['M']))

if options['GENERATE_PLOTS'] == 1:
    print "\n[Generating plots]"
    vs.visualize_OVLP(options['OUTPUT_DIR'])
    vsg.visualize_OLVP_grid(options['OUTPUT_DIR'])

