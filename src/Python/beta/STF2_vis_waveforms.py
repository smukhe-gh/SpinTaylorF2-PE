#===============================================================================
# SM: Code to test waveform generation and computing the match capabilities for 
#     SpinTaylorT2 and SpinTaylorF2
#===============================================================================

from pycbc import types, fft, waveform
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
import numpy as np
import matplotlib
from pycbc.filter import match
from pycbc.psd import aLIGOZeroDetHighPower

goldenratio = 2 / (1 + 5**.5)
matplotlib.rcParams.update({
        "font.size": 8.0,
        "axes.titlesize": 8.0,
        "axes.labelsize": 8.0,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "legend.fontsize": 8.0,
        "figure.figsize": (14.3,14.3*goldenratio),
        "figure.dpi": 200,
      # "subplots.left": 0.2,
      # "subplots.right": 0.75,
      # "subplots.bottom": 0.15,
      # "subplots.top": 0.75,
        "savefig.dpi": 600,
        "text.usetex": True     # render all text with TeX
})

side_band = [None, 2, 1, 0, -1 ,-2]

#===============================================================================
# Generating FD waveforms for SpinTaylorF2

CHECK = []

for j, i in enumerate(side_band):

    print "Generating SpinTaylorF2 \t sideband: %r" %i

    sptilde, sctilde = waveform.get_fd_waveform(approximant="SpinTaylorF2",
                             mass1=14., mass2=1.4,
                             spin1x=0.0463165745854, spin1z=-0.000336979973602, spin1y=0.497850039031,
                             inclination=0., delta_f=1.0/256,
                             sideband=i,
                             f_lower=20.,
                             f_max=2000. )

    f = sptilde.get_sample_frequencies()
    
    CHECK.append(sptilde)

    temp = 610 + j + 1
    ax = plt.subplot(temp)
    if j == 0:
        plt.plot(f, sptilde, 'r-', alpha=0.8, linewidth=0.8, label="SpinTaylorF2 SB: %r"%i)
    else:
        plt.plot(f, sptilde, 'b-', alpha=0.8, linewidth=0.8, label="SpinTaylorF2 SB: %r"%i)

    plt.subplots_adjust(hspace = .001)
    temp = tic.MaxNLocator(3)
    ax.yaxis.set_major_locator(temp)
    plt.legend(fontsize='x-small', loc='upper left', frameon=False)

plt.tight_layout()
plt.xlabel("time")
plt.savefig("../output/plots/STF2_FD_harmonics.pdf")
plt.close()

WF_sum = CHECK[1] + CHECK[2] + CHECK[3] + CHECK[4] + CHECK[5]
WF_ful = CHECK[0]

plt.plot(f, CHECK[0])
plt.show()

print np.mean(np.array(WF_sum - WF_ful))