#======================================================================
# Code to explore the bounds on H0/H2 using psiJ--thetaJ--kappa values.
# SM 1/2017
#======================================================================

"""
Here's the data we have: SNR_2/SNR_0 [chi1][thetaJ][kappa][psiJ][eta]
In this array all the values that are finite, belong to the SNR range 0.75 to 1.25
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import glob

def HRATIO_PSIJ(THETAJ, PSIJ, _thetaJ, _psiJ, A0 = 1, F = 1):
	thetaJ = THETAJ[_thetaJ]
	psiJ   = PSIJ[_psiJ]

	H0 = (3.0/4.0)*A0*F*np.cos(2.0*psiJ)*np.sin(thetaJ)**2.0
	H2 = A0*F*np.power(np.power((1.0 + np.cos(thetaJ)**2.0)/2.0, 2.0)*np.cos(2*psiJ)**2.0 + (np.cos(thetaJ)*np.sin(2*psiJ))**2.0, 0.5)
	return H0/H2

MATRIX = np.load("./temporary_datasets/SNR_RATIO_GENERIC_PSIJ.npz")

SNR_RATIO  = MATRIX["SNR_RATIO_GENERIC_PSIJ_GRID"]
H0H2_RATIO = np.zeros(np.shape(SNR_RATIO))


PSIJ   = np.linspace(0, 2*np.pi, 50) 
THETAJ = np.linspace(0, np.pi, 50)
GRIDH0H2 = np.zeros((50, 50))

for index_thetaJ, _thetaJ in enumerate(THETAJ):
	for index_psiJ, _psiJ in enumerate(PSIJ):
		GRIDH0H2[index_thetaJ][index_psiJ] = HRATIO_PSIJ(THETAJ, PSIJ, index_thetaJ, index_psiJ)
print "Finished computing H0/H2."

for index, value in np.ndenumerate(SNR_RATIO):
	if np.abs(value - 1.0) > 0.25:
		_thetaJ = index[1]
		_psiJ   = index[3]
		H0H2_RATIO[index] = GRIDH0H2[_thetaJ][_psiJ]
	else:
		H0H2_RATIO[index] = np.nan

print "Finished constructing H0H2 bounds matrix"
