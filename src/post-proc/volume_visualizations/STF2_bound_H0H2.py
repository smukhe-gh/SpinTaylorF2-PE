#======================================================================
# Code to explore the bounds on H0/H2 using psiJ--thetaJ--kappa values.
# SM 1/2017
#======================================================================

"""
What we are trying to do:

H0/H2 varies from -1.5 to 1.5 analytically. But
the actual range of H0/H2 might be smaller. 

Therefore, we are using DATA for the SNR ratio 
SNR_2/SNR_0 that has been computed over various 
values of : psiJ, thetaJ, kappa
keeping   : chi1, eta as *constant*.

Using this SNR ratio, we first--for each kappa and chi--find 
the points in the psiJ--thetaJ parameter space where this ratio 
is close to 1. This gives us a set of psiJ--thetaJ values. 

We compute H0/H2 at these values and see what are the global 
(i.e. over kappa and eta) maximum and mininum values.

These currently turn out of be: 
Max H0/H2  0.382914276499
Min H0/H2  -0.0125547027957

However, this is not a complete study. We also need the variation in chi1
to comment on the entire range. Currently chi1 = 0.5
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import interpolate
import glob

def HRATIO(thetaJ, psiJ, A0 = 1, F = 1):
	#TODO: Check this again!
	H0 = (3.0/4.0)*A0*F*np.cos(2.0*psiJ)*np.sin(thetaJ)**2.0
	H2 = A0*F*np.power(np.power((1.0 + np.cos(thetaJ)**2.0)/2.0, 2.0)*np.cos(2*psiJ)**2.0 + (np.cos(thetaJ)*np.sin(2*psiJ))**2.0, 0.5)
	return H0/H2

LIST = []
TJ = []
PJ = []
FILES = glob.glob("/Users/apple/Documents/Mars/output/datasets/output-2017_01_10_16_34_51/overlaps*")

for _file in FILES:
	DATA  = np.load(_file)
	LIST.append([0,0])
	
	print "ETA", DATA["ETA"]
	print "------------------------------"
	
	SNR_2 = DATA["SNR_02"]
	SNR_0 = DATA["SNR_00"]

	SRATIO = SNR_2/SNR_0

	"""
	Here we are only choosing the 
	points that are between 0.75 and 1.25
	"""

	for index, value in np.ndenumerate(SRATIO):
		if np.abs(value - 1.0) > 0.25:
			SRATIO[index] = np.nan
		else:
			SRATIO[index] = value

	for _kappa in range(30):
		SLICE = SRATIO[:, :, _kappa]
		
		INDEX = np.argwhere(np.isfinite(SLICE))
		
		PSIJ    = np.linspace(0.2, 0.8, 30)
		THETAJ  = np.linspace(0.0, 3.14, 30)

		R_PSIJ   = PSIJ[INDEX[:, 0]] 
		R_THETAJ = THETAJ[INDEX[:, 1]] 



		if len(R_THETAJ) > 2 and len(R_PSIJ) > 2: 

			TJ.append(np.amax(R_THETAJ))
			TJ.append(np.amin(R_THETAJ))

			PJ.append(np.amax(R_PSIJ))
			PJ.append(np.amin(R_PSIJ))

			MATRIX = np.zeros((len(R_PSIJ), len(R_THETAJ)))
			for i, _psiJ in enumerate(R_PSIJ):
				for j, _thetaJ in enumerate(R_THETAJ):
					MATRIX[i,j] = HRATIO(_thetaJ, _psiJ) 

			LIST.append([np.amax(MATRIX), np.amin(MATRIX)])

LIST = np.array(LIST)
TJ = np.array(TJ)
PJ = np.array(PJ)

print 'Max H0/H2 ', np.amax(LIST)
print 'Min H0/H2 ', np.amin(LIST)

print 'Max thetaJ ', np.amax(TJ)
print 'Min thetaJ ',   np.amin(TJ)

print 'Max psiJ ', np.amax(PJ)
print 'Min psiJ ',   np.amin(PJ)























