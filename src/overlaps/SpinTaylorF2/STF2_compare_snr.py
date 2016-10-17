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

def visualize_OLVP_grid(output_dir):

    files = set(glob.glob("../../../output/datasets/%s/overlaps*" %output_dir))

    files = sorted(files, key=lambda item: (int(item.partition(' ')[0])
                                   if item[0].isdigit() else float('inf'), item))
    if not os.path.exists("../../../output/plots/%s"%output_dir):
        os.makedirs("../../../output/plots/%s"%output_dir)

    #================================================================#
    #SNR 0F

    fig, axes = plt.subplots(nrows=3, ncols=3)
    plt.cm = plt.get_cmap('gray_r')
    fig.suptitle('SpinTaylorF2 SNR/SNR(max)')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=18)

    cmin = 0.5
    cmax = 0.5

    for i, file in enumerate(files):

        data = np.load(file)

        if np.amax(data['SNR_0F']) > cmax:
            cmax = np.amax(data['SNR_0F'])
        if np.amin(data['SNR_0F']) < cmin:
            cmin = np.amin(data['SNR_0F'])

    files_ordered = files[6:9] + files[3:6] + files[0:3]

    for i, ax in enumerate(axes.flat):
        data = np.load(files_ordered[i])

        x = np.array([0, 19, 39, 59, 79, 99])
        ax.set_xticks(x)
        ax.set_xticklabels([r"$%.2f$"%data['KAPPA'][i] for i in x])

        x = np.array([1, 20, 40, 60, 80, 100])
        ax.set_yticks(x)
        ax.set_yticklabels([r"$%.2f$"%data['THETAJ'][100-i] for i in x])

        ax.set_xlabel(r"$\kappa$")
        ax.set_ylabel(r"$\theta_{J}$")
        ax.grid()

        ax.set_title(r'$\chi_{1}=%1.2f$'%data['CHI1'] + r' $\eta=%1.2f$'%data['ETA'])
        im = ax.imshow(np.abs(np.flipud(data['SNR_0F']/cmax)),cmap='gray_r', vmin=0, vmax=1)

    fig.colorbar(im, ax=axes.ravel().tolist())
    plt.savefig('../../../output/plots/%s/'%(output_dir) + 'CM_OVLP_GRID_0F_SNR' + '.pdf')
    plt.close()

    return None

visualize_OLVP_grid('output-2016_10_04_23_22_25')
