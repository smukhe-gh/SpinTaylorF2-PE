Code to calculate the Fisher Matrix in 9D parameter space 
including phi (constant phase) and tC (constant time). 
While calculating the determinant one projects out 
the last two dimensions (phi and tC), 
which is equivalent to marginalizing Fisher over time and phase. 

The files include include

A function _FisherMatrix_ (defined in in file __STF2FisherMatrix.py__), 
which takes waveform parameters (m1,m2, chi1,kappa, alpha0,thetaJ,psiJ) 
and retun 7x7 Fisher Matrix (projected out phi and tC) and its determinant as output.

__usage__

	import STF2FisherMatrix
	wf_params = {'thetaJ': 0.01, 'psiJ': 0.01, 'kappa': 0.5, 'alpha0': 0.001, 'chi1': 0.2, 'm1': 5., 'm2': 1.4, 'phi0': 0.01, 'tC': 1.}
	fisher_matrix,fisher_det, Err_flag = STF2FisherMatrix.FisherMatrix(**wf_params)

A script __STF2Fisher_RandParams.py__, which randomize the 7 waveform parameters 
and produce fisher matrices, determinants and the parameters in an npz file.  
THe n^th element of arrays in npz file returns the Fisher matrix and its 
determinant corresponding to n^th point in  parameter space.
