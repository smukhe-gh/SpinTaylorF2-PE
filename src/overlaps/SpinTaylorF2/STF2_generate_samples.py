#=======================================================================
# Plots the region of overlap where the difference in m2 and m0 is less
# than 0.3
# SM 15/16
#=======================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import glob

goldenratio = 2. / (1 + 5**.5)
matplotlib.rcParams.update({
        "font.size": 10.0,
        "axes.titlesize": 10.0,
        "axes.labelsize": 10.0,
        "xtick.labelsize": 10.0,
        "ytick.labelsize": 10.0,
        "legend.fontsize": 10.0,
        "figure.figsize": (20.3, 25.3*goldenratio),
        "figure.dpi": 300,
        "savefig.dpi": 600,
        "text.usetex": True
})

def visualize_masked_OLVP_grid(output_dir):
    files = set(glob.glob("../../../output/datasets/%s/overlaps*" %output_dir))

    files = sorted(files, key=lambda item: (int(item.partition(' ')[0])
                                   if item[0].isdigit() else float('inf'), item))
    if not os.path.exists("../../../output/plots/%s"%output_dir):
        os.makedirs("../../../output/plots/%s"%output_dir)

    fig = plt.figure(1)
    plt.cm = plt.get_cmap('viridis')
    fig.suptitle('SpinTaylorF2: ' + r'$O_{2} - O_{0} < 0.3$')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=18)

    n = int(np.sqrt(len(files)))

    for i, file in enumerate(files):
        if i < 6:
            data = np.load(file)
            REGION = data['OLVP_MASK']
            plt.subplot(2, 3, i+1)
            plt.xlabel(r'$\kappa$')
            plt.ylabel(r'$\theta_J$')
            plt.contourf(data['KAPPA'],data['THETAJ'], REGION)
            plt.colorbar()
            plt.grid()
            plt.title(r'$\chi_{1}=%1.2f$'%data['CHI1'] + r' $\eta=%1.2f$'%data['ETA'])

    plt.savefig('../../../output/plots/%s/'%(output_dir) + 'OVLP_MASK_GRID' + '.pdf')
    plt.close()

    return None
