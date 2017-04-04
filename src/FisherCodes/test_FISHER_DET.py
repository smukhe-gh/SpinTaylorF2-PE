#==============================================================================
# Code to test the fisher matrix computation.
# MARS
#==============================================================================
import STF2FisherMatrix
import numpy as np
import os

wf_params = {
        'sideband' : None,
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

sidebands = [None, 2, 0]


print "\n==> LAL source"
os.system("which lalapps_version")
print 60*"-"

for _sideband in sidebands:
	wf_params["sideband"] = _sideband
	print "==> Sideband = ", wf_params["sideband"]
	fisher_matrix, fisher_det, Err_flag = STF2FisherMatrix.FisherMatrix(**wf_params)
	print "DET: %r" %(fisher_det)
	if wf_params["sideband"] == None:
	    print "==> Expected output"
	    print "DET: %r" %(66.787559955346467)
	    print "Residual: ", fisher_det - 66.787559955346467

	print 60*"-"
