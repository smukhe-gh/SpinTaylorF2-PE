#==============================================================================
# Code to test the coordinate transformations
# MARS
#==============================================================================

import STF2_to_lal_coords as tlc
import precessing_wf as pw

wf_params = {
	'thetaJ' : 2.96015263036,  
	'psiJ'   : 3.86315672024,  
	'kappa'  : 0.0860697805227,  
	'alpha0' : 4.82620100795,  
	'chi1'   : 0.882077210724,    
	'm1'     : 10.,  
	'm2'     : 1.4,  
	'phi0'   : 0.001,	# Check with Haris
	'tC'     : 1., 		# Check with Haris
	'f_inj'  : 20.		# Check with Haris
}

 
tlc.to_lal_coords(wf_params["m1"], wf_params["m2"], wf_params["chi1"], wf_params["kappa"], \
	wf_params["thetaJ"], wf_params["psiJ"], wf_params["alpha0"], wf_params["f_inj"])

pw.to_lal_coords(wf_params["m1"], wf_params["m2"], wf_params["chi1"], wf_params["kappa"], \
	wf_params["thetaJ"], wf_params["psiJ"], wf_params["alpha0"], wf_params["f_inj"])