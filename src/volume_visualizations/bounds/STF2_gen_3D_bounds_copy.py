#======================================================================
# Code to plot the bounds on chi and kappa. 
# SM 12/2016
#======================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import interpolate
import matplotlib.gridspec as gridspec
import glob

goldenratio = 2. / (1 + 5**.5)
matplotlib.rcParams.update({
        "font.size": 20.0,
        "axes.titlesize": 20.0,
        "axes.labelsize": 20.0,
        "xtick.labelsize": 20.0,
        "ytick.labelsize": 20.0,
        "legend.fontsize": 20.0,
        "figure.figsize": (20.3, 10.3*goldenratio),
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "text.usetex": True
})

fig = plt.figure()
GS1 = gridspec.GridSpec(1, 2)
ax1 = fig.add_subplot(GS1[0])
ax2 = fig.add_subplot(GS1[1])

def HRATIO(CUBE, _thetaJ, psiJ = 0.001, A0 = 1, F = 1):
	thetaJ = CUBE["THETAJ"][_thetaJ]

	H0 = (3.0/4.0)*A0*F*np.cos(2.0*psiJ)*np.sin(thetaJ)**2.0
	H2 = A0*F*np.power(np.power((1.0 + np.cos(thetaJ)**2.0)/2.0, 2.0)*np.cos(2*psiJ)**2.0 + (np.cos(thetaJ)*np.sin(2*psiJ))**2.0, 0.5)
	# return H2/H0

	# Changing things here to multply the SNR ratio with a factor
	# between -1.5 to 1.5.

	H2H0 = np.linspace(-1.5, 1.5, 50)
	return H2H0[_thetaJ]

def plot_slice(DATA, SLICE, _thetaJ):

	plt.axhline(DATA["CHI1"][np.amin(INDEX[:, 0])], color='red')
	plt.axvline(DATA["KAPPA"][np.amax(INDEX[:, 1])])
	plt.contourf(DATA["KAPPA"], DATA["CHI1"], SLICE)

	plt.xlabel(r"$\kappa$")
	plt.ylabel(r"$\chi_{1}$")
	plt.title(r"$\theta_{J}=%1.2f$"%np.linspace(-1.5, 1.5, 50)[_thetaJ])
	plt.set_xlim(-0.5, 1)
	plt.set_ylim(0.2, 0.8)
	plt.colorbar()
	plt.savefig("./bounds/thetaJ_%1.2f.pdf"%(thetaJ))
	plt.close()

	pass

#======================================================================
# Control Panel [1]
#======================================================================

BOUNDS = np.zeros((2, 50))
SNR    = [0.8, 1.2]
COLORS = ['green', 'red', 'blue', 'magenta']
files = glob.glob("../../output/datasets/output-2017_01_16_00_36_16/overlaps*")
select_plot = "CHI"
THRESHOLD = 7.0
#======================================================================

for k, file in enumerate(files):
	if not (k+1)%2 == 0:	
	# if not k > 1:
		for j, snr in enumerate(SNR):
			DATA = np.load(file)
			print "Working on CHI1 \t eta=%1.2f"%(DATA["ETA"])
			SNR_2 = DATA["SNR_02"]
			SNR_0 = DATA["SNR_00"]
			SRATIO = SNR_2/SNR_0

			for index, value in np.ndenumerate(SRATIO):
				if np.abs(value - snr) > 0.1 or SNR_0[index] < THRESHOLD or SNR_2[index] < THRESHOLD:
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


			H2H0_RANGE = np.linspace(-1.5, 1.5, 50)

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

		# ax1.plot(XF, BOUNDS[0], alpha = 0.8, linewidth=1.2, color=COLORS[k])
		# ax1.plot(XF, BOUNDS[1], alpha = 0.8, linewidth=1.2, color=COLORS[k], label=r'$\eta=%1.2f$'%DATA["ETA"])
		ax1.plot(XR, BR0, alpha = 0.8, linewidth=1.2, color=COLORS[k])
		ax1.plot(XR, BR1, alpha = 0.8, linewidth=1.2, color=COLORS[k], label=r'$\eta=%1.2f$'%DATA["ETA"])
		ax1.fill_between(XF, BOUNDS[0], BOUNDS[1], facecolor=COLORS[k], alpha=0.1)

ax1.legend(frameon=False, loc='upper left')
ax1.axvline([0], linewidth=4, color='white')
ax1.set_xlim(0, 1.5)
ax1.set_ylim(0.0, 1.0)
ax1.set_xlabel(r'$H_{0}/H_{2}$')
ax1.set_ylabel(r'$\chi_{\rm LB}$')
ax1.set_title(r'Lower bound on $\chi_{1}$')


#======================================================================
# Control Panel [2]
#======================================================================

BOUNDS = np.zeros((2, 50))
SNR    = [0.8, 1.2]
files = glob.glob("../../output/datasets/output-2017_01_16_00_36_16/overlaps*")
select_plot = "KAPPA"

#======================================================================

for k, file in enumerate(files):
	if not (k+1)%2 == 0:	
	# if not k > 1:
		for j, snr in enumerate(SNR):
			DATA = np.load(file)
			print "Working on KAPPA \t eta=%1.2f"%(DATA["ETA"])
			SNR_2 = DATA["SNR_02"]
			SNR_0 = DATA["SNR_00"]
			SRATIO = SNR_2/SNR_0

			for index, value in np.ndenumerate(SRATIO):
				if np.abs(value - snr) > 0.1 or SNR_0[index] < THRESHOLD or SNR_2[index] < THRESHOLD:
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


			H2H0_RANGE = np.linspace(-1.5, 1.5, 50)

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

		# ax2.plot(XF, BOUNDS[0], alpha = 0.8, linewidth=0.8, color=COLORS[k])
		# ax2.plot(XF, BOUNDS[1], alpha = 0.8, linewidth=0.8, color=COLORS[k], label=r'$\eta=%1.2f$'%DATA["ETA"])
		ax2.plot(XR, BR0, alpha = 0.8, linewidth=1.2, color=COLORS[k])
		ax2.plot(XR, BR1, alpha = 0.8, linewidth=1.2, color=COLORS[k], label=r'$\eta=%1.2f$'%DATA["ETA"])
		ax2.fill_between(XF, BOUNDS[0], BOUNDS[1], facecolor=COLORS[k], alpha=0.1)


ax2.legend(frameon=False)
ax2.axvline([0], linewidth=4, color='white')
ax2.set_xlim(0, 1.5)
ax2.set_ylim(-0.5, 1.0)
ax2.set_xlabel(r'$H_{0}/H_{2}$')
ax2.set_ylabel(r'$\kappa_{\rm UB}$')
ax2.set_title(r'Upper bound on $\kappa$')

GS1.tight_layout(fig)
plt.savefig('./immediate/chi_kappa_bounds_copy.pdf')
