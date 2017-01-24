#===============================================================================
# SM: Code to visualize SpinTaylorF2 
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
      # "subplots.left": 0.2,
      # "subplots.right": 0.75,
      # "subplots.bottom": 0.15,
      # "subplots.top": 0.75,
        "savefig.dpi": 600,
        "text.usetex": True     # render all text with TeX
})

options = {

    'ALPHA0' : 0.001,
    'KAPPA'  : 1.0,
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
    'BAND'   : None,

    'N'      : 50,
    'M'      : 50,

    'V_MASS1_RANGE' : [2.4, 50.0],
    'V_CHI1_RANGE'  : [0.20, 0.80],

    'V_KAPPA_RANGE'  : [-0.500, 0.999],
    'V_THETAJ_RANGE' : [0.001, 3.14],

    'OUTPUT_DIR'     : "output-%s" %strftime("%Y_%m_%d_%H_%M_%S", localtime()),
    'GENERATE_PLOTS' : 1,
    'THRESHOLD'      : 0.3}

SAMPLES = int(options['F_MAX']/options['DEL_F']) + 1
BANDS   = [None, 2, 1, 0, -1 ,-2]
KAPPA   = [0.0, 1.0]

FD = np.zeros((7, SAMPLES))
TD = np.zeros((7, SAMPLES))


for index, kappa in enumerate(KAPPA):

    options['KAPPA'] = kappa
    print "Generating SpinTaylorF2  \t kappa: %r" %kappa

    FDW = STF2.generate_template(**options)
    TDW = FDW.to_timeseries()

    FD[index + 1] = np.real(FDW)
    TD[index + 1] = TDW[len(TDW)/2 - 1:]

FD[0] = FDW.sample_frequencies
TD[0] = TDW[len(TDW)/2 - 1:].sample_times

# for j in range(6):
#     temp = 610 + j + 1
#     ax = plt.subplot(temp)
#     if j == 0:
#         plt.plot(FD[0], FD[j+1], 'r-', alpha=0.8, linewidth=0.5, label="SpinTaylorF2 SB: %r"%BANDS[j])
#     else:
#         plt.plot(FD[0], FD[j+1], 'b-', alpha=0.8, linewidth=0.5, label="SpinTaylorF2 SB: %r"%BANDS[j])

#     plt.subplots_adjust(hspace = .001)
#     temp = tic.MaxNLocator(3)
#     ax.yaxis.set_major_locator(temp)
#     plt.legend(fontsize='x-small', loc='lower right', frameon=False)

# plt.savefig('./FD_waveforms.pdf')
# plt.close()

# TD = TD[:, 300000:]

# for j in range(6):
#     temp = 610 + j + 1
#     ax = plt.subplot(temp)
#     if j == 0:
#         plt.plot(TD[0], TD[j+1], 'r-', alpha=0.8, linewidth=0.5, label="SpinTaylorF2 SB: %r"%BANDS[j])
#     else:
#         plt.plot(TD[0], TD[j+1], 'b-', alpha=0.8, linewidth=0.5, label="SpinTaylorF2 SB: %r"%BANDS[j])

#     plt.subplots_adjust(hspace = 0.25)
#     temp = tic.MaxNLocator(3)
#     ax.yaxis.set_major_locator(temp)
#     plt.ylabel(r'$h_{(2,2)}$')
#     plt.legend(fontsize='x-small', loc='lower right', frameon=False)

# plt.xlabel(r'$t$')
# plt.savefig('./TD_waveforms.pdf')
# plt.close()

TD = TD[:, 450000:]

for j in range(2):
    temp = 610 + j + 1
    ax = plt.subplot(temp)
    if j == 0:
        plt.plot(TD[0], TD[j+1], 'r-', alpha=0.8, linewidth=0.5, label=r"$\kappa=%r$"%KAPPA[j])
    else:
        plt.plot(TD[0], TD[j+1], 'b-', alpha=0.8, linewidth=0.5, label=r"$\kappa=%r$"%KAPPA[j])

    plt.subplots_adjust(hspace = 1.0)
    temp = tic.MaxNLocator(3)
    ax.yaxis.set_major_locator(temp)
    plt.ylabel(r'$h_{(2,2)}$')
    plt.legend(fontsize='x-small', loc='lower left', frameon=False)

plt.xlabel(r'$t$')
plt.savefig('./TD_waveforms_comparison.pdf')
plt.close()
