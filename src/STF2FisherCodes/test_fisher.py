#==============================================================================
# Code to test the fisher matrix computation. 
# Expected value for fisher_det : 2.33456435772
#==============================================================================
import STF2FisherMatrix
import numpy as np

wf_params = {
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
 
print "Sqrt[DET] : ", np.sqrt(fisher_det)
print "Expected  : ", 2.33456435772

np.savez("./fisher_test_wf_params_P256.npz",
	FISHER_MATRIX = fisher_matrix,
	FISHER_DET    = fisher_det,
	ERR_FLAG      = Err_flag,
	OPTIONS       = wf_params)
