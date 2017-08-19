#======================================================================
# Code to plot the bounds on chi and kappa.
# SM 12/2016
#======================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import interpolate
import glob

goldenratio = 2. / (1 + 5**.5)
matplotlib.rcParams.update({
        "font.size": 18.0,
        "axes.titlesize": 18.0,
        "axes.labelsize": 18.0,
        "xtick.labelsize": 18.0,
        "ytick.labelsize": 18.0,
        "legend.fontsize": 18.0,
        "figure.figsize": (10.3, 10.3*goldenratio),
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "text.usetex": True
})


def HRATIO(CUBE, _thetaJ, psiJ = 0.001, A0 = 1, F = 1):
	thetaJ = CUBE["THETAJ"][_thetaJ]

	H0 = (3.0/4.0)*A0*F*np.cos(2.0*psiJ)*np.sin(thetaJ)**2.0
	H2 = A0*F*np.power(np.power((1.0 + np.cos(thetaJ)**2.0)/2.0, 2.0)*np.cos(2*psiJ)**2.0 + (np.cos(thetaJ)*np.sin(2*psiJ))**2.0, 0.5)
	# return H0/H2
	"""
	Changing things here to
	multply the SNR ratio with a factor
	between -1.5 to 1.5.
	"""
	H2H0 = np.linspace(-1.5, 1.5, 50)
	return H2H0[_thetaJ]

#======================================================================
# Control Panel
#======================================================================

BOUNDS = np.zeros((2, 50))
SNR    = [0.8, 1.2]
COLORS = ['magenta', 'red', 'green', 'blue']
files = glob.glob("/Users/apple/Documents/Mars/output/datasets/output-2017_01_10_16_34_51/overlaps*")

select_plot = "KAPPA"
select_plot = "CHI"

#======================================================================


DATA = np.load(files[2])

SNR_2 = DATA["SNR_02"]
SNR_0 = DATA["SNR_00"]
SRATIO = SNR_2/SNR_0

for index, value in np.ndenumerate(SRATIO):
	if np.abs(value - 1.0) > 0.25:
		SRATIO[index] = np.nan
		# SRATIO[index] = value*(HRATIO(DATA, index[1]))
	else:
		SRATIO[index] = value*(HRATIO(DATA, index[1]))

print "Finished modidying the array."

CHI_BOUNDS   = np.zeros(50)
KAPPA_BOUNDS = np.zeros(50)


for _chi in range(50):

	"""
	Explicitly mention what you're choosing here.
	If SRATIO[i, :, :] then kappa (x), thetaJ (y)
	If SRATIO[:, i, :] then kappa (x), chi1 (y)

	Accordingly you'll have to change some things
	here.
	"""

	SLICE = SRATIO[:, _chi, :]
	CHI1  = np.linspace(0.2, 0.8, 50)
	INDEX = np.argwhere(np.isfinite(SLICE))



	if INDEX.size > 0:
		print _chi
		MAX = np.amax(SLICE[np.isfinite(SLICE)])
		MIN = np.amin(SLICE[np.isfinite(SLICE)])

		"""
		KAPPA_U = SLICE[:, 0]
		CHI_L   = SLICE[-1]

		CL = np.argwhere(np.isfinite(CHI_L))
		KU = np.argwhere(np.isfinite(KAPPA_U))

		CHI_BOUNDS[_chi]   = DATA["CHI1"][np.amin(INDEX[:, 0])]
		KAPPA_BOUNDS[_chi] = DATA["KAPPA"][np.amax(INDEX[:, 1])]
		"""

		plt.contourf(DATA["CHI"], DATA["KAPPA"], SLICE, vmin=MIN, vmax=MAX)

		"""
		if CL.size > 0:
			plt.axvspan(DATA["KAPPA"][CL[0]], DATA["KAPPA"][np.amax(INDEX[:, 1])], alpha =0.1, color='red')
		if KU.size > 0:
			plt.axhspan(DATA["CHI1"][KU[-1]],DATA["CHI1"][np.amin(INDEX[:, 0])], alpha =0.1, color='blue')
		"""

		plt.xlabel(r"$\kappa$")
		plt.ylabel(r"$\theta_{J}$")
		plt.title(r"$\psi_{J}=%1.2f$"%CHI1[_chi])
		plt.colorbar()
		plt.tight_layout()
		plt.savefig("./bounds/SNR2SNR0_thetaJ_kappa_slice_psiJ_%1.2f.pdf"%(CHI1))
		plt.close()

#=======================================================================

# for _thetaJ in range(50):

# 	"""
# 	Explicitly mention what you're choosing here.
# 	If SRATIO[i, :, :] then kappa (x), thetaJ (y)
# 	If SRATIO[:, i, :] then kappa (x), chi1 (y)

# 	Accordingly you'll have to change some things
# 	here.
# 	"""

# 	SLICE = SRATIO[:, _thetaJ, :]
# 	THETAJ  = DATA["THETAJ"][_thetaJ]
# 	INDEX = np.argwhere(np.isfinite(SLICE))

#         print DATA["THETAJ"][38]

# 	if INDEX.size > 0:

# 		print 'MAX = %r \t MIN = %r'%(np.amax(SLICE[np.isfinite(SLICE)]), np.amin(SLICE[np.isfinite(SLICE)]))

# 		"""
# 		KAPPA_U = SLICE[:, 0]
# 		CHI_L   = SLICE[-1]

# 		CL = np.argwhere(np.isfinite(CHI_L))
# 		KU = np.argwhere(np.isfinite(KAPPA_U))

# 		CHI_BOUNDS[_chi]   = DATA["CHI1"][np.amin(INDEX[:, 0])]
# 		KAPPA_BOUNDS[_chi] = DATA["KAPPA"][np.amax(INDEX[:, 1])]
# 		"""

# 		# plt.contourf(DATA["KAPPA"], DATA["CHI1"], SLICE)

# 		"""
# 		if CL.size > 0:
# 			plt.axvspan(DATA["KAPPA"][CL[0]], DATA["KAPPA"][np.amax(INDEX[:, 1])], alpha =0.1, color='red')
# 		if KU.size > 0:
# 			plt.axhspan(DATA["CHI1"][KU[-1]],DATA["CHI1"][np.amin(INDEX[:, 0])], alpha =0.1, color='blue')
# 		"""

# 		# plt.xlabel(r"$\kappa$")
# 		# plt.ylabel(r"$\chi_{1}$")
# 		# plt.title(r"$\theta_{J}=%1.2f$"%THETAJ)
# 		# plt.colorbar()
# 		# plt.tight_layout()
# 		# plt.savefig("./bounds/chi_kappa_slice_thetaJ_%1.2f.pdf"%(THETAJ))
# 		# plt.close()
