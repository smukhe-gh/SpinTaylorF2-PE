import matplotlib.pyplot as plt
from pycbc.waveform import get_td_waveform
import matplotlib

goldenratio = 2. / (1 + 5**.5)
matplotlib.rcParams.update({
        "font.size": 20.0,
        "axes.titlesize": 20.0,
        "axes.labelsize": 20.0,
        "xtick.labelsize": 20.0,
        "ytick.labelsize": 20.0,
        "legend.fontsize": 18.0,
        "figure.figsize": (30.0, 20.3*goldenratio),
        "figure.dpi": 300,
         # "subplots.left": 0.2,
         # "subplots.right": 0.75,
         # "subplots.bottom": 0.15,
         # "subplots.top": 0.75,
        "savefig.dpi": 300,
        "text.usetex": True
})


hp, hc = get_td_waveform(approximant='SpinTaylorT2',
                        mass1=14, mass2=1.4, #spin1z=0.8,
                        delta_t=1.0/1000.0, f_lower=30)

hp1, hc1 = get_td_waveform(approximant='SpinTaylorT2',
                        mass1=14, mass2=1.4, spin1z=0.8,
                        delta_t=1.0/1000.0, f_lower=30)

hp2, hc2 = get_td_waveform(approximant='SpinTaylorT2',
                           mass1=14, mass2=1.4, spin1z=-0.8,
                           delta_t=1.0/1000.0, f_lower=30)

hp3, hc3 = get_td_waveform(approximant='SpinTaylorT2',
                           mass1=14, mass2=1.4,
                           spin1x=-0.7, spin1y=-0.7,
                           delta_t=1.0/1000.0, f_lower=30)

lim1=(-11.0,0.1)

plt.figure(1)
plt.subplot(4,1,1)
plt.title(r'A')
plt.plot(hp.sample_times,hc,label='Non-spinning, $~|\mathbf{S}|=0.0$',c='g')
plt.ylabel(r'$h_{+}$')
plt.xlim(lim1)
plt.locator_params(axis='y',nbins=6)
plt.legend(loc='upper left',frameon=False)

plt.subplot(4,1,2)
plt.title(r'B')
plt.plot(hp1.sample_times,hc1,label='Aligned spin, $~S_{1z}=0.8$',c='r')
plt.ylabel(r'$h_{+}$')
plt.xlim(lim1)
plt.ylim(-1.0*1e-19, 1.0*1e-19)
plt.locator_params(axis='y',nbins=6)
plt.legend(loc='upper left',frameon=False)

plt.subplot(4,1,3)
plt.title(r'C')
plt.plot(hp2.sample_times,hc2, label='Anti-aligned spin, $~S_{1z}=-0.8$', c='b')
plt.xlim(lim1)
plt.ylabel(r'$h_{+}$')
plt.locator_params(axis='y',nbins=6) #to specify number of ticks on both or any single axes
plt.legend(loc='upper left',frameon=False)

plt.subplot(4,1,4)
plt.title(r'D')
plt.plot(hp3.sample_times,hc3, label='Precessing,$~S_{1x}=S_{1y}=-0.7$',c='black')
plt.xlim(lim1)
plt.ylabel(r'$h_{+}$')
plt.xlabel(r'$t$')
plt.locator_params(axis='y',nbins=6)
plt.legend(loc='upper left', frameon=False)

# plt.subplots_adjust(top=1.85)
plt.tight_layout()
# plt.suptitle(r'SpinTaylorT2 waveforms for $m_1=14~M_{\odot},~m_2=1.4~M_{\odot}$:~spinning and non-spinning cases')
plt.savefig("./immediate/spinning_nonspinning_waveforms.pdf")
