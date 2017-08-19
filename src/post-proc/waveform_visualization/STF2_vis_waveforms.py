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

#==============================================================================
options = {

    'ALPHA0' : 0.001,
    'KAPPA'  : -0.3,
    'CHI1'   : 0.5,
    'PSIJ'   : 0.001,
    'M1'     : 50.0,
    'PHI0'   : 0.001,
    'M2'     : 1.4,
    'THETAJ' : 0.7853981633974483,
    'APPROX' : 'SpinTaylorF2',
    'DEL_F'  : 1./256.,
    'F_MIN'  : 20.,
    'F_INJ'  : 20.,
    'F_MAX'  : 2000.,
    'BAND'   : 2,
    }

SAMPLES = int(options['F_MAX']/options['DEL_F']) + 1
BANDS   = [None, 2, 1, 0, -1 ,-2]
KAPPA   = np.linspace(-0.5, 1.0, 5)
CHI1    = np.linspace( 0.1, 1.0, 5)

string  = 'plot_param'

#==============================================================================

FD = np.zeros((7, SAMPLES))
TD = np.zeros((7, SAMPLES))

"""
Chaning things here to iterate over
chi for a fixed value of kappa.
"""
for index, chi1 in enumerate(CHI1):

    options['CHI1'] = chi1
    print "Generating SpinTaylorF2  \t chi1: %1.1f \t kappa: %1.1f" %(chi1, options["KAPPA"])

    FDW = STF2.generate_template(**options)
    TDW = FDW.to_timeseries()

    FD[index + 1] = np.real(FDW)
    TD[index + 1] = TDW[len(TDW)/2 - 1:]

FD[0] = FDW.sample_frequencies
TD[0] = TDW[len(TDW)/2 - 1:].sample_times

print "\nGenerating visualizations."

if string == "plot_SB":
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

    plt.savefig('./FD_waveforms.pdf')
    plt.close()

    TD = TD[:, 300000:]

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
    plt.savefig('./TD_waveforms.pdf')
    plt.close()

if string == "plot_param":
    TD = TD[:, 450000:]

    # TODO: Pull this range from the length of TD or FD.
    for j in range(5):
        temp = 610 + j + 1
        ax = plt.subplot(temp)
        if j % 2 == 0:
            plt.plot(TD[0], TD[j+1], 'r-', alpha=0.8, linewidth=0.5, label=r"$SB=%r, \chi_1=%1.2f, \kappa=%1.2f$"%(options["BAND"], CHI1[j], options["KAPPA"]))
        else:
            plt.plot(TD[0], TD[j+1], 'b-', alpha=0.8, linewidth=0.5, label=r"$SB=%r, \chi_1=%1.2f, \kappa=%1.2f$"%(options["BAND"], CHI1[j], options["KAPPA"]))

        plt.subplots_adjust(hspace = 1.0)
        temp = tic.MaxNLocator(3)
        ax.yaxis.set_major_locator(temp)
        plt.ylabel(r'$h_{(2,2)}$')
        plt.legend(fontsize='x-small', loc='lower left', frameon=False)

    # Change the file name here if you like.
    plt.xlabel(r'$t$')
    plt.savefig('./immediate/SB_%1.2f_kappa_%1.2f.pdf'%(options["BAND"], options["KAPPA"]))
    plt.close()
