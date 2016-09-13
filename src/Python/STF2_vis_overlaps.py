#=======================================================================
# Plots the OLVPs over the parameter space.
# SM 13/ 16
# TODO: Add support for SpinTaylorF2
#=======================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import glob

goldenratio = 2. / (1 + 5**.5)
matplotlib.rcParams.update({
        "font.size": 8.0,
        "axes.titlesize": 8.0,
        "axes.labelsize": 8.0,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "legend.fontsize": 8.0,
        "figure.figsize": (20.3, 20.3*goldenratio),
        "figure.dpi": 300,
      # "subplots.left": 0.2,
      # "subplots.right": 0.75,
      # "subplots.bottom": 0.15,
      # "subplots.top": 0.75,
        "savefig.dpi": 600,
        "text.usetex": True
})

files = glob.glob("../../output/overlaps_eta_*_chi1_*_N_*.npz")

if not os.path.exists("../../output"):
    os.makedirs("../../output")

if not os.path.exists("../../output/plots"):
    os.makedirs("../../output/plots")

for file in files:
    data = np.load(file)

    plt.cm = plt.get_cmap('viridis')
    fig = plt.figure()
    fig.suptitle('SpinTaylorF2: Overlaps, '+ r'$\chi_{1}=%1.2f$'%data['CHI1'] + r' $\eta=%1.2f$'%data['ETA'])
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=18)

    plt.subplot(2,4,1)
    plt.ylabel(r'$\kappa$')
    plt.xlabel(r'$\theta_J$')
    plt.contourf(data['THETAJ'],data['KAPPA'],data['OLVP_0F_P2'])
    plt.clim(-0.02,1.)
    plt.grid()
    plt.title('Overlap (m = 2)')
    plt.colorbar()

    plt.subplot(2,4,2)
    plt.ylabel(r'$\kappa$')
    plt.xlabel(r'$\theta_J$')
    plt.contourf(data['THETAJ'],data['KAPPA'],data['OLVP_0F_P1'])
    plt.clim(-0.02,1.)
    plt.grid()
    plt.title('Overlap (m = 1)')
    plt.colorbar()

    plt.subplot(2,4,3)
    plt.ylabel(r'$\kappa$')
    plt.xlabel(r'$\theta_J$')
    plt.contourf(data['THETAJ'],data['KAPPA'],data['OLVP_0F_P0'])
    plt.clim(-0.02,1.)
    plt.grid()
    plt.title('Overlap (m = 0)')
    plt.colorbar()

    plt.subplot(2,4,4)
    plt.ylabel(r'$\kappa$')
    plt.xlabel(r'$\theta_J$')
    plt.contourf(data['THETAJ'],data['KAPPA'],data['OLVP_0F_M1'])
    plt.clim(-0.02,1.)
    plt.grid()
    plt.title('Overlap (m = -1)')
    plt.colorbar()

    plt.subplot(2,4,5)
    plt.ylabel(r'$\kappa$')
    plt.xlabel(r'$\theta_J$')
    plt.contourf(data['THETAJ'],data['KAPPA'],data['OLVP_0F_M2'])
    plt.clim(-0.02,1.)
    plt.grid()
    plt.title('Overlap (m = -2)')
    plt.colorbar()

    plt.subplot(2,4,6)
    plt.ylabel(r'$\kappa$')
    plt.xlabel(r'$\theta_J$')
    plt.contourf(data['THETAJ'],data['KAPPA'],data['SNR_00'])
    plt.grid()
    plt.title('SNR ($m=0$)')
    plt.colorbar()

    plt.subplot(2,4,7)
    plt.ylabel(r'$\kappa$')
    plt.xlabel(r'$\theta_J$')
    plt.contourf(data['THETAJ'],data['KAPPA'],data['SNR_02'])
    plt.grid()
    plt.title('SNR ($m=2$)')
    plt.colorbar()

    plt.subplot(2,4,8)
    plt.ylabel(r'$\kappa$')
    plt.xlabel(r'$\theta_J$')
    plt.contourf(data['THETAJ'],data['KAPPA'],data['SNR_0F'])
    plt.grid()
    plt.title('SNR (Full waveform)')
    plt.colorbar()

    plt.savefig('../../output/plots/' + 'OVLPS_CHI1_%1.2f_'%data['CHI1'] + 'ETA_%1.2f'%data['ETA'] + '.pdf')


