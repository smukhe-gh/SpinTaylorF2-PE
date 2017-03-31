#==============================================================================
# Code to test the fisher matrix computation.
# MARS
#==============================================================================
import STF2FisherMatrix
import numpy as np

wf_params = {
        'sideband': 0,
	'thetaJ' : 2.96015263036,
	'psiJ'   : 3.86315672024,
	'kappa'  : 0.0860697805227,
	'alpha0' : 4.82620100795,
	'chi1'   : 0.882077210724,
	'm1'     : 10.,
	'm2'     : 1.4,
	'phi0'   : 0.001,	# Check with Haris
	'tC'     : 1. 		# Check with Haris
}

fisher_matrix, fisher_det, Err_flag = STF2FisherMatrix.FisherMatrix(**wf_params)

print "Log10[Sqrt[Abs[DET]]] : ", np.log10(np.sqrt(np.abs(fisher_det)))


