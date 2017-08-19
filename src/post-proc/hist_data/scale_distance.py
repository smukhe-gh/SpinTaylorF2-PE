#===============================================================
# Code to spread out the sources in uniform volume.
# Soham 08 2017.
#===============================================================

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata

goldenratio = 2. / (1 + 5**.5)
matplotlib.rcParams.update({
     "font.size": 14.0,
     "axes.titlesize": 14.0,
     "axes.labelsize": 14.0,
     "xtick.labelsize": 14.0,
     "ytick.labelsize": 14.0,
     "legend.fontsize": 14.0,
     "figure.figsize": (15, 10),
     "figure.dpi": 300,
     "savefig.dpi": 300,
     "text.usetex": True
})

def plot_SNR_histogram(SNR_matrix):
	fig, axes = plt.subplots(2, 3)	#create subplot grid of 5x3
	fig.tight_layout()
	fig.subplots_adjust(left=0.2, bottom=0.2, right=None, top=None, wspace=0.28, hspace=0.28)
	
	for loc in range(np.shape(SNR_matrix)[0]):
		im = axes.flat[loc].hist(SNR_matrix[loc], bins=100, color='r')
		axes.flat[loc].set_xlim(0, 30)
		axes.flat[loc].set_xlabel(r"$%s$"%SNR_labels[loc])

	plt.savefig("test_SNR.png", bbox_inches='tight')
	plt.close()

def plot_SNR_contours(SNR_matrix):
	fig, axes = plt.subplots(2, 3)	#create subplot grid of 5x3
	fig.tight_layout()
	fig.subplots_adjust(left=0.2, bottom=0.2, right=None, top=None, wspace=0.28, hspace=0.28)
	
	for loc in range(np.shape(SNR_matrix)[0]):
		im = axes.flat[loc].scatter(data["m1"], data["chi1"], SNR_matrix[loc])
		# axes.flat[loc].set_xlim(0, 30)
		# axes.flat[loc].set_xlabel(r"$%s$"%SNR_labels[loc])

	plt.savefig("test_SNR.png", bbox_inches='tight')
	plt.close()


filename ='../../../datasets/Hist_data_random_cosine_angles_m2_1.4_ThetaJKappaChiAlphaPsiJM1_N_50000.npz'
data     = np.load(filename)

# define parameter space keys, labels, and limits
params = ["kappa", "thetaJ", "psiJ", "chi1", "alpha0", "m1"]
labels = [r"$\kappa$", r"$\theta_J$", r"$\psi_J$", r"$\chi_1$", r"$\alpha_0$", r"$m_1$"]
lims_L = [-0.5, 0.0, 0., 0., 0., 0.0, 3.0]
lims_U = [1.0, np.pi, np.pi, 1.0, np.pi, 50.0]

# define data labels and their dictionary keys.
overlap = ["Ov2", "Ov1", "Ov0", "Ovm1", "Ovm2"] 
overlap_labels = ["O_2", "O_1", "O_0", "O_{-1}", "O_{-2}"] 
SNR = ["SNR", "SNR2", "SNR1", "SNR0", "SNRm1", "SNRm2"]
SNR_codes  = [42, 2, 1, 0, -1, -2] 
SNR_labels = ["\\rm{Total~SNR}", "\\rm{SNR}_2", "\\rm{SNR}_1", "\\rm{SNR}_0", "\\rm{SNR}_{-1}", "\\rm{SNR}_{-2}"] 

# construct the scaling factors. 
N = len(data["SNR"])
u = np.random.rand(N);
R_max = 1600
R_min = 100
V_min = np.power(R_min,3)/3
V_max = np.power(R_max,3)/3
c = V_max - V_min
R = np.power(3*(c*u + V_min), 1./3 )

SNR_matrix = np.zeros((len(data["SNR"]), len(SNR))).T

# rescaling all the SNRs and populating a SNR_matrix...
for _i, vec in enumerate(SNR):
	if(1):
		SNR_matrix[_i]   = data[vec]*(400.0/R)
	else:
		SNR_matrix[_i]   = data[vec]

np.savez("../../../datasets/Hist_data_random_cosine_angles_m2_1.4_ThetaJKappaChiAlphaPsiJM1_N_50000_scaled.npz",
	R = R,
	SNR   = SNR_matrix[0],
	SNR2  = SNR_matrix[1],
	SNR1  = SNR_matrix[2],
	SNR0  = SNR_matrix[3],
	SNRm1 = SNR_matrix[4],
	SNRm2 = SNR_matrix[5],
)

#---------------------------------------------------------------
# Compute the probabilities of the all the combinations you like.
#---------------------------------------------------------------

import itertools

how_many_sidebands = 3
SNR_threshold = 5.0
combinations = np.array(list(itertools.combinations(np.arange(1,6,1), how_many_sidebands)))

gl = []
for _r in combinations:
	ll = []
	for _c in _r:
		ll.append(np.where(SNR_matrix[_c] > SNR_threshold))
	res = reduce(np.intersect1d, ll)
	gl.append(res)

for _i, _r in enumerate(combinations):
	for _c in _r:
		print "%r \t"%SNR_codes[_c],
	print " :\t %r" %(1.0*len(gl[_i])/len(R))








