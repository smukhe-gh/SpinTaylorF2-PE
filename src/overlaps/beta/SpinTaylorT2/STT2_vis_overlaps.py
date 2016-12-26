#=======================================================================
# Plots the OLVPs over the parameter space.
# SM 15/16
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
        "figure.figsize": (20.3, 25.3*goldenratio),
        "figure.dpi": 300,
      # "subplots.left": 0.2,
      # "subplots.right": 0.75,
      # "subplots.bottom": 0.15,
      # "subplots.top": 0.75,
        "savefig.dpi": 600,
        "text.usetex": True
})

def visualize_OVLP(output_dir):

    files = glob.glob("../../../output/datasets/%s/overlaps*" %output_dir)

    if not os.path.exists("../../../output/plots/%s"%output_dir):
        os.makedirs("../../../output/plots/%s"%output_dir)

    for file in files:

        data = np.load(file)
        
        cmax = np.amax(np.vstack((data['OLVP_T2_0F'], data['OLVP_T2_P2P0'], data['OLVP_T2_P2'], data['OLVP_T2_P1'], data['OLVP_T2_P0'], data['OLVP_T2_M1'], data['OLVP_T2_M2'])))
        cmin = np.amin(np.vstack((data['OLVP_T2_0F'], data['OLVP_T2_P2P0'], data['OLVP_T2_P2'], data['OLVP_T2_P1'], data['OLVP_T2_P0'], data['OLVP_T2_M1'], data['OLVP_T2_M2'])))        

        plt.cm = plt.get_cmap('viridis')
        fig = plt.figure()
        fig.suptitle('Overlaps + SNR (SpinTaylorT2), '+ r'$\chi_{1}=%1.2f$'%data['CHI1'] + r' $\eta=%1.2f$'%data['ETA'])
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif',size=18)

        plt.subplot(3,3,1)
        plt.xlabel(r'$\kappa$')
        plt.ylabel(r'$\theta_J$')
        plt.contourf(data['KAPPA'],data['THETAJ'],data['OLVP_T2_0F'])
        plt.grid()
        plt.title(r'$O$' + 'F2')
        plt.clim(cmin, cmax)
        plt.colorbar()


        plt.subplot(3,3,2)
        plt.xlabel(r'$\kappa$')
        plt.ylabel(r'$\theta_J$')
        plt.contourf(data['KAPPA'],data['THETAJ'],data['OLVP_T2_P2'])
        plt.grid()
        plt.title(r'$O$' + ' (m = 2)')
        plt.clim(cmin, cmax)
        plt.colorbar()


        plt.subplot(3,3,3)
        plt.ylabel(r'$\kappa$')
        plt.ylabel(r'$\theta_J$')
        plt.contourf(data['KAPPA'],data['THETAJ'],data['OLVP_T2_P1'])
        plt.grid()
        plt.title(r'$O$' + ' (m = 1)')
        plt.clim(cmin, cmax)
        plt.colorbar()


        plt.subplot(3,3,4)
        plt.xlabel(r'$\kappa$')
        plt.ylabel(r'$\theta_J$')
        plt.contourf(data['KAPPA'],data['THETAJ'],data['OLVP_T2_P0'])
        plt.grid()
        plt.title(r'$O$' + ' (m = 0)')
        plt.clim(cmin, cmax)
        plt.colorbar()


        plt.subplot(3,3,5)
        plt.xlabel(r'$\kappa$')
        plt.ylabel(r'$\theta_J$')
        plt.contourf(data['KAPPA'],data['THETAJ'],data['OLVP_T2_M1'])
        plt.grid()
        plt.title(r'$O$' + ' (m = -1)')
        plt.clim(cmin, cmax)
        plt.colorbar()


        plt.subplot(3,3,6)
        plt.xlabel(r'$\kappa$')
        plt.ylabel(r'$\theta_J$')
        plt.contourf(data['KAPPA'],data['THETAJ'],data['OLVP_T2_M2'])
        plt.grid()
        plt.title(r'$O$' + ' ($m=-2$)')
        plt.clim(cmin, cmax)
        plt.colorbar()


        plt.subplot(3,3,7)
        plt.xlabel(r'$\kappa$')
        plt.ylabel(r'$\theta_J$')
        plt.contourf(data['KAPPA'],data['THETAJ'],data['OLVP_T2_P2P0'])
        plt.grid()
        plt.title(r'$O$' + ' ($m=0 + m=2$)')
        plt.clim(cmin, cmax)
        plt.colorbar()

        smax = np.amax(np.vstack((data['SNR_0F'], data['SNR_T2'])))
        smin = np.amin(np.vstack((data['SNR_0F'], data['SNR_T2'])))

        plt.subplot(3,3,8)
        plt.xlabel(r'$\kappa$')
        plt.ylabel(r'$\theta_J$')
        plt.contourf(data['KAPPA'],data['THETAJ'],data['SNR_0F'])
        plt.grid()
        plt.title('SNR F2')
        plt.clim(smin, smax)
        plt.colorbar()

        plt.subplot(3,3,9)
        plt.xlabel(r'$\kappa$')
        plt.ylabel(r'$\theta_J$')
        plt.contourf(data['KAPPA'],data['THETAJ'],data['SNR_T2'])
        plt.grid()
        plt.title('SNR T2')
        plt.clim(smin, smax)
        plt.colorbar()

        plt.savefig('../../../output/plots/%s/'%(output_dir) + 'OVLP_CHI1_%1.2f_'%data['CHI1'] + 'ETA_%1.2f'%data['ETA'] + '.pdf')
        plt.close()

    return None
