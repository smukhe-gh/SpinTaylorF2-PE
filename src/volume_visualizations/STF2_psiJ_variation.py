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
        "font.size": 18.0,
        "axes.titlesize": 18.0,
        "axes.labelsize": 18.0,
        "xtick.labelsize": 18.0,
        "ytick.labelsize": 18.0,
        "legend.fontsize": 18.0,
        "figure.figsize": (15.3, 15.3*goldenratio),
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "text.usetex": True
})


#======================================================================
# Control Panel
#======================================================================

BOUNDS = np.zeros((2, 50))
SNR    = [0.8, 1.2]
COLORS = ['magenta', 'red', 'green', 'blue']
files = glob.glob("/Users/apple/Documents/Mars/output/datasets/output-2017_01_10_16_34_51/overlaps*")

#======================================================================

DATA = np.load(files[2])

SNR_2 = DATA["SNR_02"]
SNR_0 = DATA["SNR_00"]

fig = plt.figure()
GS1 = gridspec.GridSpec(2, 3)

PSIJ =  np.linspace(0.2, 0.8, 30)

#for SNR_2
ax1 = fig.add_subplot(GS1[0, 0])
ax2 = fig.add_subplot(GS1[0, 1])
ax3 = fig.add_subplot(GS1[0, 2])

#for SNR_0
ax4 = fig.add_subplot(GS1[1, 0])
ax5 = fig.add_subplot(GS1[1, 1])
ax6 = fig.add_subplot(GS1[1, 2])

MAX = np.amax(np.vstack((SNR_2[5, :, :], SNR_2[15, :, :], SNR_2[25, :, :])))
MIN = np.amin(np.vstack((SNR_2[5, :, :], SNR_2[15, :, :], SNR_2[25, :, :])))
print MAX, MIN

ax1.contourf(DATA["KAPPA"], DATA["THETAJ"], SNR_2[5, :, :], vmin=MIN, vmax=MAX)
ax1.set_xlabel(r"$\kappa$")
ax1.set_ylabel(r"$\theta_{J}$")
ax1.set_title(r"$\rm SNR_{2},~\psi_{J}=%1.2f$"%PSIJ[5])
# ax1.set_colorbar()

ax2.contourf(DATA["KAPPA"], DATA["THETAJ"], SNR_2[15, :, :], vmin=MIN, vmax=MAX)
ax2.set_xlabel(r"$\kappa$")
ax2.set_ylabel(r"$\theta_{J}$")
ax2.set_title(r"$\rm SNR_{2},~\psi_{J}=%1.2f$"%PSIJ[15])
# ax2.set_colorbar()

ax3.contourf(DATA["KAPPA"], DATA["THETAJ"], SNR_2[25, :, :], vmin=MIN, vmax=MAX)
ax3.set_xlabel(r"$\kappa$")
ax3.set_ylabel(r"$\theta_{J}$")
ax3.set_title( r"$\rm SNR_{2},~\psi_{J}=%1.2f$"%PSIJ[25])
# ax3.set_colorbar()

MAX = np.amax(np.vstack((SNR_0[5, :, :], SNR_0[15, :, :], SNR_0[25, :, :])))
MIN = np.amin(np.vstack((SNR_0[5, :, :], SNR_0[15, :, :], SNR_0[25, :, :])))
print MAX, MIN

ax4.contourf(DATA["KAPPA"], DATA["THETAJ"], SNR_0[5, :, :], vmin=MIN, vmax=MAX)
ax4.set_xlabel(r"$\kappa$")
ax4.set_ylabel(r"$\theta_{J}$")
ax4.set_title(r"$\rm SNR_{0},~\psi_{J}=%1.2f$"%PSIJ[5])
# ax4.set_colorbar()

ax5.contourf(DATA["KAPPA"], DATA["THETAJ"], SNR_0[15, :, :], vmin=MIN, vmax=MAX)
ax5.set_xlabel(r"$\kappa$")
ax5.set_ylabel(r"$\theta_{J}$")
ax5.set_title(r"$\rm SNR_{0},~\psi_{J}=%1.2f$"%PSIJ[15])
# ax5.set_colorbar()

ax6.contourf(DATA["KAPPA"], DATA["THETAJ"], SNR_0[25, :, :], vmin=MIN, vmax=MAX)
ax6.set_xlabel(r"$\kappa$")
ax6.set_ylabel(r"$\theta_{J}$")
ax6.set_title(r"$\rm SNR_{0},~\psi_{J}=%1.2f$"%PSIJ[25])
# ax6.set_colorbar()

GS1.tight_layout(fig)
plt.savefig('./test_PSIJ.pdf')

# for _chi in range(30):

# 	"""
# 	Explicitly mention what you're choosing here.
# 	If SRATIO[i, :, :] then kappa (x), thetaJ (y)
# 	If SRATIO[:, i, :] then kappa (x), chi1 (y)

# 	Accordingly you'll have to change some things
# 	here.
# 	"""

# 	SLICE = SRATIO[_chi, :, :]
# 	CHI1  = DATA["CHI1"][_chi]
# 	INDEX = np.argwhere(np.isfinite(SLICE))



# 	if INDEX.size > 0:
# 		print _chi
# 		MAX = np.amax(SLICE[np.isfinite(SLICE)])
# 		MIN = np.amin(SLICE[np.isfinite(SLICE)])

# 		"""
# 		KAPPA_U = SLICE[:, 0]
# 		CHI_L   = SLICE[-1]

# 		CL = np.argwhere(np.isfinite(CHI_L))
# 		KU = np.argwhere(np.isfinite(KAPPA_U))

# 		CHI_BOUNDS[_chi]   = DATA["CHI1"][np.amin(INDEX[:, 0])]
# 		KAPPA_BOUNDS[_chi] = DATA["KAPPA"][np.amax(INDEX[:, 1])]
# 		"""

# 		plt.contourf(DATA["KAPPA"], DATA["THETAJ"], SLICE, vmin=MIN, vmax=MAX)

# 		"""
# 		if CL.size > 0:
# 			plt.axvspan(DATA["KAPPA"][CL[0]], DATA["KAPPA"][np.amax(INDEX[:, 1])], alpha =0.1, color='red')
# 		if KU.size > 0:
# 			plt.axhspan(DATA["CHI1"][KU[-1]],DATA["CHI1"][np.amin(INDEX[:, 0])], alpha =0.1, color='blue')
# 		"""

# 		plt.set_xlabel(r"$\kappa$")
# 		plt.set_ylabel(r"$\theta_{J}$")
# 		plt.set_title(r"$\psi_{J}=%1.2f$"%CHI1)
# 		plt.colorbar()
# 		plt.tight_layout()
# 		plt.savefig("./bounds/SNR2_thetaJ_kappa_slice_psiJ_%1.2f.pdf"%(CHI1))
# 		plt.close()

# #=======================================================================

# # for _thetaJ in range(50):

# # 	"""
# # 	Explicitly mention what you're choosing here.
# # 	If SRATIO[i, :, :] then kappa (x), thetaJ (y)
# # 	If SRATIO[:, i, :] then kappa (x), chi1 (y)

# # 	Accordingly you'll have to change some things
# # 	here.
# # 	"""

# # 	SLICE = SRATIO[:, _thetaJ, :]
# # 	THETAJ  = DATA["THETAJ"][_thetaJ]
# # 	INDEX = np.argwhere(np.isfinite(SLICE))

# #         print DATA["THETAJ"][38]

# # 	if INDEX.size > 0:

# # 		print 'MAX = %r \t MIN = %r'%(np.amax(SLICE[np.isfinite(SLICE)]), np.amin(SLICE[np.isfinite(SLICE)]))

# # 		"""
# # 		KAPPA_U = SLICE[:, 0]
# # 		CHI_L   = SLICE[-1]

# # 		CL = np.argwhere(np.isfinite(CHI_L))
# # 		KU = np.argwhere(np.isfinite(KAPPA_U))

# # 		CHI_BOUNDS[_chi]   = DATA["CHI1"][np.amin(INDEX[:, 0])]
# # 		KAPPA_BOUNDS[_chi] = DATA["KAPPA"][np.amax(INDEX[:, 1])]
# # 		"""

# # 		# plt.contourf(DATA["KAPPA"], DATA["CHI1"], SLICE)

# # 		"""
# # 		if CL.size > 0:
# # 			plt.axvspan(DATA["KAPPA"][CL[0]], DATA["KAPPA"][np.amax(INDEX[:, 1])], alpha =0.1, color='red')
# # 		if KU.size > 0:
# # 			plt.axhspan(DATA["CHI1"][KU[-1]],DATA["CHI1"][np.amin(INDEX[:, 0])], alpha =0.1, color='blue')
# # 		"""

# # 		# plt.xlabel(r"$\kappa$")
# # 		# plt.set_ylabel(r"$\chi_{1}$")
# # 		# plt.title(r"$\theta_{J}=%1.2f$"%THETAJ)
# # 		# plt.colorbar()
# # 		# plt.tight_layout()
# # 		# plt.savefig("./bounds/chi_kappa_slice_thetaJ_%1.2f.pdf"%(THETAJ))
# # 		# plt.close()
