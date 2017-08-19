import numpy as np
import matplotlib
import matplotlib.pyplot as plt

goldenratio = 2. / (1 + 5**.5)
matplotlib.rcParams.update({
     "font.size": 10.0,
     "axes.titlesize": 10.0,
     "axes.labelsize": 10.0,
     "xtick.labelsize": 10.0,
     "ytick.labelsize": 10.0,
     "legend.fontsize": 10.0,
     "figure.figsize": (10, 15),
     "figure.dpi": 300,
     "savefig.dpi": 300,
     "text.usetex": True
})

filename='../../datasets/Hist_data_random_cosine_angles_m2_1.4_ThetaJKappaChiAlphaPsiJM1_N_50000.npz'
#filename='Hist_data_random_angles_m2_1.4_ThetaJKappaChiAlphaPsiM1_N_50000.npz'
data   = np.load(filename)

# define parameter space keys, labels, and limits
params = ["kappa", "thetaJ", "psiJ", "chi1", "alpha0", "m1"]
labels = [r"$\kappa$", r"$\theta_J$", r"$\psi_J$", r"$\chi_1$", r"$\alpha_0$", r"$m_1$"]
lims_L = [-0.5, 0.0, 0., 0., 0., 0.0, 3.0]
lims_U = [1.0, np.pi, np.pi, 1.0, np.pi, 50.0]

# define data labels and their dictionary keys.
overlap = ["Ov2", "Ov1", "Ov0", "Ovm1", "Ovm2"] 
overlap_labels = ["O_2", "O_1", "O_0", "O_{-1}", "O_{-2}"] 
SNR = ["SNR", "SNR2", "SNR1", "SNR0", "SNRm1", "SNRm2"] 
SNR_labels = ["Total SNR", "SNR_2", "SNR_1", "SNR_0", "SNR_{-1}", "SNR_{-2}"] 



# main function definition
def loop_2_modes():
	fig, axes = plt.subplots(5, 3)	#create subplot grid of 5x3
	fig.tight_layout()
	fig.subplots_adjust(left=0.2, bottom=0.2, right=None, top=None, wspace=0.28, hspace=0.28)
	
	# XXX: If you need to change what you're plotting, you need to change in several places here.
	SNR_range = [8.0, np.amax(data['SNR'])]
	
	ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNR1'] > 8.0), np.logical_and(data['SNR0'] < 8.0, np.logical_and(data['SNRm1'] < 8.0, data['SNRm2'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNR0'] > 8.0), np.logical_and(data['SNR1'] < 8.0, np.logical_and(data['SNRm1'] < 8.0, data['SNRm2'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNRm1'] > 8.0), np.logical_and(data['SNR1'] < 8.0, np.logical_and(data['SNR0'] < 8.0, data['SNRm2'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNRm2'] > 8.0), np.logical_and(data['SNR1'] < 8.0, np.logical_and(data['SNRm1'] < 8.0, data['SNR0'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR1'] > 8.0, data['SNR0'] > 8.0), np.logical_and(data['SNR2'] < 8.0, np.logical_and(data['SNRm1'] < 8.0, data['SNRm2'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR1'] > 8.0, data['SNRm1'] > 8.0), np.logical_and(data['SNR2'] < 8.0, np.logical_and(data['SNR0'] < 8.0, data['SNRm2'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR1'] > 8.0, data['SNRm2'] > 8.0), np.logical_and(data['SNR2'] < 8.0, np.logical_and(data['SNRm1'] < 8.0, data['SNR0'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR0'] > 8.0, data['SNRm1'] > 8.0), np.logical_and(data['SNR2'] < 8.0, np.logical_and(data['SNR1'] < 8.0, data['SNRm2'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR0'] > 8.0, data['SNRm2'] > 8.0), np.logical_and(data['SNR1'] < 8.0, np.logical_and(data['SNRm1'] < 8.0, data['SNR2'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNRm1'] > 8.0, data['SNRm2'] > 8.0), np.logical_and(data['SNR1'] < 8.0, np.logical_and(data['SNR2'] < 8.0, data['SNR0'] < 8.0)))
	x = np.linspace(0,1,100)
	loc  = 0 

	for i, par1 in enumerate(params):	#loop over first param
		for j, par2 in enumerate(params):	#loop over second param
			if (i < j):	#choose only upper diagonal terms
				im = axes.flat[loc].scatter(data[params[i]][ind2], data[params[j]][ind2], c=data['SNR'][ind2], 
					s=10, alpha=1.0, edgecolor='face', cmap='viridis',
					vmin = SNR_range[0], vmax = SNR_range[1])
				axes.flat[loc].set_xlim(lims_L[i], lims_U[i])
				axes.flat[loc].set_ylim(lims_L[j], lims_U[j])
				axes.flat[loc].set_xlabel(labels[i])	# note that I'm passing a list element as a whole.
				axes.flat[loc].set_ylabel(labels[j])	# otherwise it typesets incorrectly.
				loc += 1 # keep increasing the axes.flat index so that it plots sequentially

	fig.colorbar(im, ax=axes.ravel().tolist(), aspect=100)	#common colorbar
	plt.savefig("SNR_all_params_2_modes_only_2_1.pdf", bbox_inches='tight')
plt.close()

#	return None

# call function [pretty slow]
#loop_2_modes()

def loop_3_modes():
	fig, axes = plt.subplots(5, 3)	#create subplot grid of 5x3
	fig.tight_layout()
	fig.subplots_adjust(left=0.2, bottom=0.2, right=None, top=None, wspace=0.28, hspace=0.28)
	
	# XXX: If you need to change what you're plotting, you need to change in several places here.
	SNR_range = [8.0, np.amax(data["SNR"])]
	
	ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNR1'] > 8.0), np.logical_and(data['SNR0'] > 8.0, np.logical_and(data['SNRm1'] < 8.0, data['SNRm2'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNR1'] > 8.0), np.logical_and(data['SNR0'] < 8.0, np.logical_and(data['SNRm1'] > 8.0, data['SNRm2'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNR1'] > 8.0), np.logical_and(data['SNR0'] < 8.0, np.logical_and(data['SNRm1'] < 8.0, data['SNRm2'] > 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNR1'] < 8.0), np.logical_and(data['SNR0'] > 8.0, np.logical_and(data['SNRm1'] > 8.0, data['SNRm2'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNR1'] < 8.0), np.logical_and(data['SNR0'] > 8.0, np.logical_and(data['SNRm1'] < 8.0, data['SNRm2'] > 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNR1'] < 8.0), np.logical_and(data['SNR0'] < 8.0, np.logical_and(data['SNRm1'] > 8.0, data['SNRm2'] > 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] < 8.0, data['SNR1'] > 8.0), np.logical_and(data['SNR0'] > 8.0, np.logical_and(data['SNRm1'] > 8.0, data['SNRm2'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] < 8.0, data['SNR1'] > 8.0), np.logical_and(data['SNR0'] > 8.0, np.logical_and(data['SNRm1'] < 8.0, data['SNRm2'] > 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] < 8.0, data['SNR1'] > 8.0), np.logical_and(data['SNR0'] < 8.0, np.logical_and(data['SNRm1'] > 8.0, data['SNRm2'] > 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] < 8.0, data['SNR1'] < 8.0), np.logical_and(data['SNR0'] > 8.0, np.logical_and(data['SNRm1'] > 8.0, data['SNRm2'] > 8.0)))
	x = np.linspace(0,1,100)
	loc  = 0 

	for i, par1 in enumerate(params):	#loop over first param
		for j, par2 in enumerate(params):	#loop over second param
			if (i < j):	#choose only upper diagonal terms
				im = axes.flat[loc].scatter(data[params[i]][ind2], data[params[j]][ind2], c=data["SNR"][ind2], 
					s=10, alpha=1.0, edgecolor='face', cmap='viridis',
					vmin = SNR_range[0], vmax = SNR_range[1])
				axes.flat[loc].set_xlim(lims_L[i], lims_U[i])
				axes.flat[loc].set_ylim(lims_L[j], lims_U[j])
				axes.flat[loc].set_xlabel(labels[i])	# note that I'm passing a list element as a whole.
				axes.flat[loc].set_ylabel(labels[j])	# otherwise it typesets incorrectly.
				loc += 1 # keep increasing the axes.flat index so that it plots sequentially

	fig.colorbar(im, ax=axes.ravel().tolist(), aspect=100)	#common colorbar
	plt.savefig("SNR_all_params_3_modes_only_2_1_0.pdf", bbox_inches='tight')
plt.close()

def loop_4_modes():
	fig, axes = plt.subplots(5, 3)	#create subplot grid of 5x3
	fig.tight_layout()
	fig.subplots_adjust(left=0.2, bottom=0.2, right=None, top=None, wspace=0.28, hspace=0.28)
	
	# XXX: If you need to change what you're plotting, you need to change in several places here.
	SNR_range = [8.0, np.amax(data["SNR"])]
	ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNR1'] > 8.0), np.logical_and(data['SNR0'] > 8.0, np.logical_and(data['SNRm1'] > 8.0, data['SNRm2'] < 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNR1'] > 8.0), np.logical_and(data['SNR0'] > 8.0, np.logical_and(data['SNRm1'] < 8.0, data['SNRm2'] > 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNR1'] > 8.0), np.logical_and(data['SNR0'] < 8.0, np.logical_and(data['SNRm1'] > 8.0, data['SNRm2'] > 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] > 8.0, data['SNR1'] < 8.0), np.logical_and(data['SNR0'] > 8.0, np.logical_and(data['SNRm1'] > 8.0, data['SNRm2'] > 8.0)))
	#ind2 = np.logical_and(np.logical_and(data['SNR2'] < 8.0, data['SNR1'] > 8.0), np.logical_and(data['SNR0'] > 8.0, np.logical_and(data['SNRm1'] > 8.0, data['SNRm2'] > 8.0)))
	
	x = np.linspace(0,1,100)
	loc  = 0 

	for i, par1 in enumerate(params):	#loop over first param
		for j, par2 in enumerate(params):	#loop over second param
			if (i < j):	#choose only upper diagonal terms
				im = axes.flat[loc].scatter(data[params[i]][ind2], data[params[j]][ind2], c=data["SNR"][ind2], 
					s=10, alpha=1.0, edgecolor='face', cmap='viridis',
					vmin = SNR_range[0], vmax = SNR_range[1])
				axes.flat[loc].set_xlim(lims_L[i], lims_U[i])
				axes.flat[loc].set_ylim(lims_L[j], lims_U[j])
				axes.flat[loc].set_xlabel(labels[i])	# note that I'm passing a list element as a whole.
				axes.flat[loc].set_ylabel(labels[j])	# otherwise it typesets incorrectly.
				loc += 1 # keep increasing the axes.flat index so that it plots sequentially

	fig.colorbar(im, ax=axes.ravel().tolist(), aspect=100)	#common colorbar
	plt.savefig("SNR_all_params_4_modes_only_2_1_0_m1.pdf", bbox_inches='tight')
plt.close()
