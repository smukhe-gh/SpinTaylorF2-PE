#===============================================================================
# SM: Code to visualize
#     SpinTaylorT2 and SpinTaylorF2
#===============================================================================

from pycbc import types, fft, waveform
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
import numpy as np
import matplotlib
from pycbc.filter import match
from pycbc.psd import aLIGOZeroDetHighPower
import STF2_waveform as STF2
from time import localtime, strftime

goldenratio = 2 / (1 + 5**.5)
matplotlib.rcParams.update({
        "font.size": 10.0,
        "axes.titlesize": 10.0,
        "axes.labelsize": 10.0,
        "xtick.labelsize": 10.0,
        "ytick.labelsize": 10.0,
        "legend.fontsize": 10.0,
        "figure.figsize": (15.3, 15.3*goldenratio),
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "text.usetex": True     # render all text with TeX
})

DATA = np.load("")

BANDS   = [None, 2, 1, 0, -1 ,-2]

if (1):
    for j in range(6):
        temp = 610 + j + 1
        ax = plt.subplot(temp)
        if j == 0:
            plt.plot(FD[0], FD[j+1], 'r-', alpha=0.8, linewidth=0.5, label="SpinTaylorF2 SB: %r"%BANDS[j])
        else:
            plt.plot(FD[0], FD[j+1], 'b-', alpha=0.8, linewidth=0.5, label="SpinTaylorF2 SB: %r"%BANDS[j])

        plt.subplots_adjust(hspace = .001)
        temp = tic.MaxNLocator(3)
        ax.yaxis.set_major_locator(temp)
        plt.legend(fontsize='x-small', loc='lower right', frameon=False)

    plt.savefig('./immediate/FD_waveforms.pdf')
    plt.close()

    TD = TD[:, 450000:]

    for j in range(6):
        temp = 610 + j + 1
        ax = plt.subplot(temp)
        if j == 0:
            plt.plot(TD[0], TD[j+1], 'r-', alpha=0.8, linewidth=0.5, label="SpinTaylorF2 SB: %r"%BANDS[j])
        else:
            plt.plot(TD[0], TD[j+1], 'b-', alpha=0.8, linewidth=0.5, label="SpinTaylorF2 SB: %r"%BANDS[j])

        plt.subplots_adjust(hspace = 0.25)
        temp = tic.MaxNLocator(3)
        ax.yaxis.set_major_locator(temp)
        plt.ylabel(r'$h_{(2,2)}$')
        plt.legend(fontsize='x-small', loc='lower right', frameon=False)

    plt.xlabel(r'$t$')
    plt.savefig('./immediate/TD_waveforms.pdf')
    plt.close()
