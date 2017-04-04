#==============================================================================
# Code to test the fisher matrix computation.
# MARS
#==============================================================================
import STF2FisherMatrix
import numpy as np
import os

wf_params = {
        # 'sideband' : None, #[You need to be in STF2Fisher_parallel branch to use this option]
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

print "LAL source"
os.system("which lalapps_version")
print 60*"-"
#print "Sideband = ", wf_params["sideband"]
print "Log10[Sqrt[DET]] : %r DET: %r" %(np.log10(np.sqrt(fisher_det)), fisher_det)
print "Expected output"
print 60*"-"
if(1):
    print "Log10[Sqrt[DET]] : %r DET: %r" %(0.91234778852914733, 66.787559955346467)
    print "Residual: ", np.log10(np.sqrt(fisher_det)) - 0.91234778852914733
