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

from scipy import ndimage as ndi
from skimage import feature
from scipy.signal import argrelextrema
from scipy.optimize import curve_fit

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

def clean_inneredge(IS):
    """
    Fucking complicated. FIX IT!
    """
    IS_X = IS[:,0]
    L = int(len(IS_X)/2)
    A = argrelextrema(IS_X[:L], np.less)[0][-1]
    B = argrelextrema(IS_X[L:], np.less)[0][0]
    C = argrelextrema(IS_X[L:], np.greater)[0][0]
    D = argrelextrema(IS_X[:L], np.greater)[0][0]

    if B < C:
        IS = IS[A:L+B]
    elif B > C:
        IS = IS[A:L+C-1]
    elif A > D:
        IS = IS[A:L+B]
    elif A < D:
        IS = IS[D+1:L+B]

    return IS

def extract_edges(REGION):
    #TODO: Play with these values. Ask Archana about this.
    REGION[np.isnan(REGION)] = -1
    EDGES = np.flipud(feature.canny(REGION, sigma=5))
    X, Y = np.shape(EDGES)
    IS = []
    OS = []

    for j in range(Y):
        PROFILE = np.array(np.where(EDGES[j] == True)[0])
        if np.shape(PROFILE)[0] >= 2:
            IS.append([PROFILE[0],j])
            OS.append([PROFILE[-1],j])

    OS = np.array(OS)
    IS = np.array(IS)
    IS = clean_inneredge(IS)

    return IS, OS, EDGES

def quadratic(x, a, b):
    return - a*(x**2) + b

def visualize_masked_OLVP_grid(output_dir):
    files = set(glob.glob("../../../output/datasets/%s/overlaps*" %output_dir))

    files = sorted(files, key=lambda item: (int(item.partition(' ')[0])
                                   if item[0].isdigit() else float('inf'), item))
    if not os.path.exists("../../../output/plots/%s"%output_dir):
        os.makedirs("../../../output/plots/%s"%output_dir)

    fig = plt.figure(1)
    fig.suptitle('SpinTaylorF2: ' + r'$O_{2} - O_{0} < 0.3$')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=18)

    n = int(np.sqrt(len(files)))

    for i, file in enumerate(files):
        if i < 2*n:
            data = np.load(file)
            IS, OS, EDGES = extract_edges(data['OLVP_MASK'])
            plt.subplot(2, n, i+1)
            KAPPA  = data['KAPPA']
            THETAJ = data['THETAJ']

            plt.plot(KAPPA[IS[:,0]], THETAJ[IS[:,1]],'b-')
            plt.plot(KAPPA[OS[:,0]], THETAJ[OS[:,1]],'r-')
            plt.xlim(-0.500, 0.999)
            plt.ylim(0.001, 3.14)

            plt.title(r'$\chi_{1}=%1.2f$'%data['CHI1'] + r' $\eta=%1.2f$'%data['ETA'])

    plt.savefig('../../../output/plots/%s/'%(output_dir) + 'OVLP_GRID_EDGES' + '.pdf')
    plt.close()


    return None

visualize_masked_OLVP_grid('output-2016_09_30_20_10_03')
