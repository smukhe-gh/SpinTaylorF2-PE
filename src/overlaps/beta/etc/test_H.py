import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
goldenratio = 2. / (1 + 5**.5)
# matplotlib.rcParams.update({
#         "font.size": 14.0,
#         "axes.titlesize": 14.0,
#         "axes.labelsize": 14.0,
#         "xtick.labelsize": 14.0,
#         "ytick.labelsize": 14.0,
#         "legend.fontsize": 14.0,
#         "figure.figsize": (12.3, 12.3*goldenratio),
#         "figure.dpi": 300,
#         "savefig.dpi": 300,
#         "text.usetex": True
# })

A0 = 1.0
F   = 1.0

def H_ratio(theta_J, psi_J):

	#Check with Mathematica
	H2 = A0*F*np.sqrt(np.power((1.0 + np.cos(theta_J)**2.0)/2.0, 2.0)*np.cos(2*psi_J)**2.0 + \
		np.power(np.cos(theta_J)*np.sin(2*psi_J),2))
	H0 = (3.0/4.0)*A0*F*np.cos(2.0*psi_J)*np.sin(theta_J)**2.0

	return H0/H2


THETA_J = np.linspace(0, np.pi, 100)
PSI_J   = np.linspace(0, 2*np.pi, 200)

RATIO = np.zeros((100, 200))

for i, _theta_J in enumerate(THETA_J):
	for j, _psi_J in enumerate(PSI_J):
		RATIO[i][j] = H_ratio(_theta_J, _psi_J)


# plt.figure()
# ax = plt.gca()
# im = ax.imshow(RATIO)

# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.05)

# ax.set_xlabel(r'$\psi_{J}$')
# ax.set_ylabel(r'$\theta_J$')
# plt.colorbar(im, cax=cax)
# plt.savefig("./H_ratio.pdf")

DATA = np.load("./output-2016_10_04_23_22_25/overlaps_eta_0.08_chi1_0.80_N_100.npz")
SNR_0F  = DATA["SNR_0F"]
SNR_00  = DATA["SNR_00"]
SNR_02  = DATA["SNR_02"]
THETA_J = DATA["THETAJ"]
PSI_J   = 0.001

RATIO = SNR_00/SNR_02

H_RATIO = np.zeros((100, 100))
for i, _theta_J in enumerate(THETA_J):
	H_RATIO[i] = (1.0/H_ratio(_theta_J, PSI_J))*RATIO[i]


plt.imshow(np.flipud(H_RATIO))
plt.colorbar()
plt.xlabel(r'$\kappa$')
plt.ylabel(r'$\theta_J$')
plt.show()

# plt.imshow(SNR_0F)
# plt.colorbar()
# plt.plot(SNR_0F[50])



