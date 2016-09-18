#=======================================================================
# Plots the grid of OLVPs of m=2 and m=0 over the parameter space.
# SM 15/16
#=======================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import glob
from joblib import Parallel, delayed

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
        "savefig.dpi": 600,
        "text.usetex": True
})

def generate_GRID(OVLP, output_dir, files):
    
    fig = plt.figure()
    plt.cm = plt.get_cmap('viridis')
    fig.suptitle('Overlap (m=%r)'%int(OVLP[-1]))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=18)

    gridsize = int(np.sqrt(len(files)))
    cmin = 0.5
    cmax = 0.5
    
    for i, file in enumerate(files):
        data = np.load(file)
        if np.amax(data[OVLP]) > cmax:
            cmax = np.amax(data[OVLP])
        if np.amin(data[OVLP]) < cmin:
            cmin = np.amax(data[OVLP])

    for index, file in enumerate(files):
        data = np.load(file)
        plt.subplot(gridsize, gridsize, index+1)
        plt.xlabel(r'$\kappa$')
        plt.ylabel(r'$\theta_J$')
        plt.contourf(data['KAPPA'],data['THETAJ'],data[OVLP])
        plt.clim(cmin, cmax)
        plt.colorbar()
        plt.grid()
        plt.title(r'$\chi_{1}=%1.2f$'%data['CHI1'] + r' $\eta=%1.2f$'%data['ETA'])
        
    plt.savefig('../../output/plots/%s/'%(output_dir) + 'OVLP_GRID_m_%r'%int(OVLP[-1]) + '.pdf')
    plt.close()

    return None

def visualize_OLVP_grid(output_dir):

    files = set(glob.glob("../../output/datasets/%s/overlaps*" %output_dir))
    files = sorted(files, key=lambda item: (int(item.partition(' ')[0])
                                   if item[0].isdigit() else float('inf'), item))
    if not os.path.exists("../../output/plots/%s"%output_dir):
        os.makedirs("../../output/plots/%s"%output_dir)

    PLOTS = ['OLVP_0F_P0', 'OLVP_0F_P2']
    
    for OLVP in PLOTS:
        generate_GRID(OLVP, output_dir, files)
        
    # Parallel(n_jobs=-2, verbose=5)(delayed(generate_GRID)(OVLP = OVLP, output_dir = output_dir, files = files) for OVLP in PLOTS)
    
    return None

visualize_OLVP_grid("output-2016_09_18_11_45_17")