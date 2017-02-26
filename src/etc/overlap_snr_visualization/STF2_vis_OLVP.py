#=======================================================================
# Plots the grid of OLVPs of m=2 and m=0 over the parameter space.
# SM 15/16
#=======================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob

goldenratio = 2. / (1 + 5**.5)
mpl.rcParams.update({
        "font.size": 32.0,
        "axes.titlesize": 32.0,
        "axes.labelsize": 32.0,
        "xtick.labelsize": 32.0,
        "ytick.labelsize": 32.0,
        "legend.fontsize": 32.0,
        "figure.figsize": (10.3, 12.3),
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "text.usetex": True
})

def visualize_OLVP_grid(output_dir, STRING):

    files = glob.glob("../../../output/datasets/%s/overlaps*" %output_dir)
    
    fig, axes = plt.subplots(nrows=3, ncols=3)
    fig.suptitle(r'$\rm %s$'%STRING)

    DAT = []
    for i, file in enumerate(files):

        data = np.load(file)
        DAT.append(data[STRING])

    cmin = np.amin(np.array(DAT))
    cmax = np.amax(np.array(DAT))

    files_ordered = files[6:9] + files[3:6] + files[0:3]
    FITS = []
    
    for i, ax in enumerate(axes.flat):
        data = np.load(files_ordered[i])

    #     x = np.array([0, 10, 20, 30, 40, 49])
    #     ax.set_xticks(x)
    #     ax.set_xticklabels([r"$%.2f$"%data['KAPPA'][i] for i in x])

    #     x = np.array([49, 40, 30, 20, 10, 0])
    #     ax.set_yticks(x)
    #     ax.set_yticklabels([r"$%.2f$"%data['THETAJ'][i] for i in x])
    #     ax.set_xlabel(r"$\kappa$")
    #     ax.set_ylabel(r"$\theta_{J}$")
    #     ax.grid()
    #     ax.set_title(r'$\chi_{1}=%1.2f$'%data['CHI1'] + r' $\eta=%1.2f$'%data['ETA'])
    #     im = ax.imshow(np.flipud(data[STRING]), cmap='gray_r', vmin=cmin, vmax=cmax)
        
    # fig.colorbar(im, ax=axes.ravel().tolist())
    # plt.savefig('./immediate/test.pdf')
    # plt.close()

    return None

visualize_OLVP_grid('output-2016_10_19_19_24_50', "SNR_0F")
