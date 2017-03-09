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

    'M'      : 200,

    'V_MASS1_RANGE' : [12.4, 100.0],
    'V_CHI1_RANGE'  : [0.20, 0.93],

    'V_KAPPA_RANGE'  : [-0.500, 0.999],
    'V_THETAJ_RANGE' : [0.001, 3.14],

    'OUTPUT_DIR'     : "output-%s" %strftime("%Y_%m_%d_%H_%M_%S", localtime())
    }

#=======================================================================
# END CONTROL PANEL
#=======================================================================

def generate_SNR(**options):

	V_MASS1   = np.linspace(options['V_MASS1_RANGE'][0], \
	options['V_MASS1_RANGE'][1], options['M'])
	SNR       = np.zeros((options['M'], 9))

 	for _M1 in xrange(options['M']):
 	    options['M1'] = V_MASS1[_M1]
            options['CHI1'] = 0.01
            print "M1 : \t", V_MASS1[_M1]
            SNR[_M1][:] = compute_overlap(**options)

        filename = "SNR_M1_%s_chi1_%s_N_%r.npz" %('{:.2f}'.format(options["M1"]),\
                '{:.2f}'.format(options['CHI1']), options['M'])

        np.savez("../../output/datasets/%s/%s" %(options['OUTPUT_DIR'], filename),
                DATE         = strftime("%Y-%m-%d %H:%M:%S", localtime()),
                SNR_0F       = SNR[:, 0],
                SNR_02       = SNR[:, 1],
                SNR_00       = SNR[:, 2],
                M1_VEC       = V_MASS1,
                OPTIONS      = options)

        return None


if not os.path.exists("../../output/datasets/%s" %options['OUTPUT_DIR']):
        os.makedirs("../../output/datasets/%s" %options['OUTPUT_DIR'])


generate_SNR(**options)
