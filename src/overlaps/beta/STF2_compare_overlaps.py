#=======================================================================
# Plots the grid of OLVPs of m=2 and m=0 over the parameter space.
# SM 15/16
#=======================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob
from scipy.signal import argrelextrema

import pandas
# For statistics. Requires statsmodels 5.0 or more
from statsmodels.formula.api import ols
# Analysis of Variance (ANOVA) on linear models
from statsmodels.stats.anova import anova_lm

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
    
    # def get_points(data):

    #     slices = np.arange(0,100,1)
    #     points = []        
    #     for slice in slices:
    #         somedata = data[STRING][:,slice]
    #         index = argrelextrema(somedata, np.less)[0]
    #         if len(index) == 1 and index > 10 and index < 90:
    #             points.append(np.array([slice, index]))

        
    #     points = np.array(points)
    #     if len(points) < 10:
    #         return None

    #     _kappa = data['KAPPA']
    #     _thetaJ  = data['THETAJ']

    #     img = np.polyfit(points[:,0], points[:,1], 1)
    #     # x = points[:,0]
    #     x = np.arange(0,100,1)
    #     y = points[:,1]
    #     z = img[0]*x + img[1]

    #     coeffs = np.polyfit(_kappa[points[:,0]], _thetaJ[points[:,1]], 1)

    #     return x, y, z, coeffs

    fig, axes = plt.subplots(nrows=3, ncols=3)
    fig.suptitle('SpinTaylorF2 Overlap (m=2)')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=18)

    n = int(np.sqrt(len(files)))

    cmin = 0.5
    cmax = 0.5
    
    #set what you want to plot here. 
    # STRING = 'SNR_00'

    DAT = []
    for i, file in enumerate(files):
        data = np.load(file)
        if np.amax(data['SNR_02']) > cmax:
            cmax = np.amax(data['SNR_02'])
        if np.amin(data['SNR_02']) < cmin:
            cmin = np.amin(data['SNR_02'])


    files_ordered = files[6:9] + files[3:6] + files[0:3]
    FITS = []
    
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
     
        # if get_points(data) != None:
        #     x, y, z, coeffs = get_points(data)            
        #     ax.plot(x, 99 - z, 'm-', linewidth=0.8)
        #     FITS.append([data['CHI1'], data['ETA'], coeffs[0], coeffs[1]])

        im = ax.imshow(np.flipud(data['SNR_02']), cmap='gray_r', vmin=cmin, vmax=cmax)
        
    fig.colorbar(im, ax=axes.ravel().tolist())
    plt.savefig('../../../output/plots/%s/'%(output_dir) + 'CM_OVLP_GRID_RATIO' + '.pdf')
    plt.close()

    # FITS = np.array(FITS)
    # _chi = FITS[:,0]
    # _eta = FITS[:,1]
    # _a   = FITS[:,2]
    # _b   = FITS[:,3]

    # print 60*"-"
    # print "A fit"
    # print 60*"-"

    # data = pandas.DataFrame({'x': _chi, 'y': _eta, 'z': _a})
    # model = ols("z ~ x + y", data).fit()
    # print(model.summary())

    # print "\n"

    # print 60*"-"
    # print "B fit"
    # print 60*"-"

    # data = pandas.DataFrame({'x': _chi, 'y': _eta, 'z': _b})
    # model = ols("z ~ x + y", data).fit()
    # print(model.summary())

    return None

visualize_OLVP_grid('output-2016_10_04_23_22_25')
