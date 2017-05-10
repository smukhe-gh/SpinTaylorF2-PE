#===============================================================================
# SM: Code to visualize
#     SpinTaylorT2 and SpinTaylorF2
#===============================================================================

import matplotlib
matplotlib.use('Agg')
from pycbc import types, fft, waveform
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
import numpy as np
import matplotlib
from pycbc.filter import match
from pycbc.psd import aLIGOZeroDetHighPower
import STF2_waveform as STF2
from time import localtime, strftime

if(0):
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

options = {

    'ALPHA0' : 0.001,
    'KAPPA'  : 0.201,
    'CHI1'   : 0.9001,
    'PSIJ'   : 0.001,
    'M1'     : 50.0,
    'PHI0'   : 0.001,
    'M2'     : 1.4,
    'THETAJ' : 1.58,
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

#========================================================================
string  = 'plot_SB' #'plot_param'
#========================================================================

FD = np.zeros((7, SAMPLES))
TD = np.zeros((7, SAMPLES))

FF = np.zeros((7, SAMPLES))
TF = np.zeros((7, SAMPLES))

for index, SB in enumerate(BANDS):

    options['BAND'] = SB
    print "Generating SpinTaylorF2  chi1: %1.1f \t kappa: %1.1f \t SB: %r" %(options["CHI1"], options["KAPPA"],SB)

    FDW, hc = STF2.generate_template(**options)

    #new additions
    TDW = FDW.to_timeseries()

    FD[index + 1] = np.real(FDW)
    TD[index + 1] = TDW[len(TDW)/2 - 1:]

FD[0] = FDW.sample_frequencies
TD[0] = TDW[len(TDW)/2 - 1:].sample_times

np.savez("./immediate/waveform_prec_RC_%1.2f.npz"%options["THETAJ"],
        frequency_domain = FD,
        time_domain      = TD,
        options          = options)
