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
	# return H2/H0
	
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
files = glob.glob("../../output/datasets/output-2016_10_19_19_24_50/overlaps*")
select_plot = "KAPPA"
select_plot = "CHI"

#======================================================================

for k, file in enumerate(files):
	if not (k+1)%2 == 0:	
		for j, snr in enumerate(SNR):
			DATA = np.load(file)
			print DATA["ETA"]
			SNR_2 = DATA["SNR_02"]
			SNR_0 = DATA["SNR_00"]
			SRATIO = SNR_2/SNR_0

			for index, value in np.ndenumerate(SRATIO):
				if np.abs(value - snr) > 0.1:
					SRATIO[index] = np.nan
				else:
					SRATIO[index] = value*HRATIO(DATA, index[1])

			CHI_BOUNDS   = np.zeros(50)
			KAPPA_BOUNDS = np.zeros(50)

			for _thetaJ in range(50):
				SLICE = SRATIO[:, _thetaJ, :]
				thetaJ = DATA["THETAJ"][_thetaJ]

				INDEX = np.argwhere(np.isfinite(SLICE))
				if INDEX.size > 0:

					CHI_BOUNDS[_thetaJ]   = DATA["CHI1"][np.amin(INDEX[:, 0])]
					KAPPA_BOUNDS[_thetaJ] = DATA["KAPPA"][np.amax(INDEX[:, 1])]

					"""
					Uncomment here if you want to plot each slice.
					"""
					# plt.axhline(DATA["CHI1"][np.amin(INDEX[:, 0])], color='red')
					# plt.axvline(DATA["KAPPA"][np.amax(INDEX[:, 1])])
					# plt.contourf(DATA["KAPPA"], DATA["CHI1"], SLICE)

					# plt.xlabel(r"$\kappa$")
					# plt.ylabel(r"$\chi_{1}$")
					# plt.title(r"$\theta_{J}=%1.2f$"%np.linspace(-1.5, 1.5, 50)[_thetaJ])
					# plt.xlim(-0.5, 1)
					# plt.ylim(0.2, 0.8)
					# plt.colorbar()
					# plt.savefig("./bounds/thetaJ_%1.2f.pdf"%(thetaJ))
					# plt.close()

			H2H0_RANGE = np.linspace(-1.5, 1.5, 50)

			if select_plot == "KAPPA":
				X = H2H0_RANGE[np.nonzero(KAPPA_BOUNDS)]
				Y = KAPPA_BOUNDS[np.nonzero(KAPPA_BOUNDS)]
				
				DEGREE = 2
				C = np.polyfit(X, Y, DEGREE)

				XF = H2H0_RANGE
				YF = 0

				for DG in range(DEGREE+1):
					YF = YF+ C[DG]*XF**(DEGREE-DG)

				BOUNDS[j] = YF

				AB = np.nonzero(KAPPA_BOUNDS)[0][0]
				BC = np.nonzero(KAPPA_BOUNDS)[0][-1]

				XR = H2H0_RANGE[AB: BC]
				BR0 = BOUNDS[0][AB: BC]
				BR1 = BOUNDS[1][AB: BC]


			elif select_plot == "CHI":
				X = H2H0_RANGE[np.nonzero(CHI_BOUNDS)]
				Y = CHI_BOUNDS[np.nonzero(CHI_BOUNDS)]
				
				DEGREE = 2
				C = np.polyfit(X, Y, DEGREE)

				XF = H2H0_RANGE
				YF = 0

				for DG in range(DEGREE+1):
					YF = YF+ C[DG]*XF**(DEGREE-DG)

				BOUNDS[j] = YF

				AB = np.nonzero(KAPPA_BOUNDS)[0][0]
				BC = np.nonzero(KAPPA_BOUNDS)[0][-1]

				XR = H2H0_RANGE[AB: BC]
				BR0 = BOUNDS[0][AB: BC]
				BR1 = BOUNDS[1][AB: BC]

		if select_plot == "KAPPA":
			# plt.plot(XF, BOUNDS[0], alpha = 0.8, linewidth=0.8, color=COLORS[k])
			# plt.plot(XF, BOUNDS[1], alpha = 0.8, linewidth=0.8, color=COLORS[k], label=r'$\eta=%1.2f$'%DATA["ETA"])
			plt.plot(XR, BR0, alpha = 0.8, linewidth=1.2, color=COLORS[k])
			plt.plot(XR, BR1, alpha = 0.8, linewidth=1.2, color=COLORS[k], label=r'$\eta=%1.2f$'%DATA["ETA"])
			plt.fill_between(XF, BOUNDS[0], BOUNDS[1], facecolor=COLORS[k], alpha=0.1)
		elif select_plot == "CHI":
			# plt.plot(XF, BOUNDS[0], alpha = 0.8, linewidth=1.2, color=COLORS[k+1])
			# plt.plot(XF, BOUNDS[1], alpha = 0.8, linewidth=1.2, color=COLORS[k+1], label=r'$\eta=%1.2f$'%DATA["ETA"])
			plt.plot(XR, BR0, alpha = 0.8, linewidth=1.2, color=COLORS[k+1])
			plt.plot(XR, BR1, alpha = 0.8, linewidth=1.2, color=COLORS[k+1], label=r'$\eta=%1.2f$'%DATA["ETA"])
			plt.fill_between(XF, BOUNDS[0], BOUNDS[1], facecolor=COLORS[k+1], alpha=0.1)

plt.legend(frameon=False)

if select_plot == "KAPPA":
	plt.axvline([0], linewidth=4, color='white')
	plt.xlim(-1.5, 1.5)
	plt.ylim(-0.5, 1.0)
	plt.xlabel(r'$H0/H2$')
	plt.ylabel(r'$\kappa_{\rm UB}$')
	plt.savefig("./test_kappa.pdf")
elif select_plot == "CHI":
	plt.axvline([0], linewidth=4, color='white')
	plt.xlim(-1.5, 1.5)
	plt.ylim(0.2, 0.8)
	plt.xlabel(r'$H0/H2$')
	plt.ylabel(r'$\chi_{\rm LB}$')
	plt.savefig("./test_chi.pdf")



