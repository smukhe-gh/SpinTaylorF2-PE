#=======================================================================
# Plots the grid of OLVPs of m=2 and m=0 over the parameter space.
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

def visualize_OLVP_grid(output_dir):

    files = glob.glob("../../output/datasets/%s/overlaps*" %output_dir)
    
    if not os.path.exists("../../output/plots/%s"%output_dir):
        os.makedirs("../../output/plots/%s"%output_dir)
    
   
    fig = plt.figure(1)
    plt.cm = plt.get_cmap('viridis')
    fig.suptitle('Overlap (m=0)')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=18)

    n = int(np.sqrt(len(files)))
    
    for i, file in enumerate(files):

        data = np.load(file)

        plt.subplot(n, n, i+1)
        plt.xlabel(r'$\kappa$')
        plt.ylabel(r'$\theta_J$')
        plt.contourf(data['THETAJ'],data['KAPPA'],data['OLVP_0F_P0'])
        plt.colorbar()
        plt.grid()
        plt.title(r'$\chi_{1}=%1.2f$'%data['CHI1'] + r' $\eta=%1.2f$'%data['ETA'])


    plt.savefig('../../output/plots/%s/'%(output_dir) + 'OVLP_GRID_m_0' + '.pdf')
    plt.close()
    
    fig = plt.figure(2)
    plt.cm = plt.get_cmap('viridis')
    fig.suptitle('Overlap (m=2)')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=18)

    for i, file in enumerate(files):

        data = np.load(file)

        plt.subplot(n, n, i+1)
        plt.xlabel(r'$\kappa$')
        plt.ylabel(r'$\theta_J$')
        plt.contourf(data['THETAJ'],data['KAPPA'],data['OLVP_0F_P2'])
        plt.colorbar()
        plt.grid()
        plt.title(r'$\chi_{1}=%1.2f$'%data['CHI1'] + r' $\eta=%1.2f$'%data['ETA'])


    plt.savefig('../../output/plots/%s/'%(output_dir) + 'OVLP_GRID_m_2' + '.pdf')
    plt.close()
    
    return None