import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import numpy as np
from numpy import pi

goldenratio = 2. / (1 + 5**.5)
matplotlib.rcParams.update({
        "font.size": 24.0,
        "font.family": 'serif',
        "font.weight": 'normal',
        "axes.titlesize": 24.0,
        "axes.labelsize": 24.0,
        "xtick.labelsize": 24.0,
        "ytick.labelsize": 24.0,
        "legend.fontsize": 24.0,
        "figure.figsize": (20.3, 10.3*goldenratio),
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "text.usetex": True})

fig = plt.figure()
GS1 = gridspec.GridSpec(1, 2)
ax1 = fig.add_subplot(GS1[0])
ax2 = fig.add_subplot(GS1[1])

def m1m2_to_mchirpeta(m1, m2):
	eta = m1*m2/(m1+m2)/(m1+m2)
	mchirp = (m1+m2) * pow(eta, 3./5.)
	return (mchirp, eta)

chi1 = np.linspace(0, 1.0, 100)
eta = np.linspace(0.01, 0.13, 100)

def kappa(chi1, eta, flag):
	if flag==1:
	    return (145.83*chi1 - 155.92)*eta**2.0 - (1.10*chi1 + 0.16)*eta + 0.08*chi1 + 0.50
	elif flag==0:
            return (49.98*chi1 - 57.13)*eta**2.0 + (1.96*chi1 - 2.39)*eta - 0.09*chi1 + 0.82

print np.amax(kappa(chi1, 0.08, 1))
print np.amax(kappa(chi1, 0.08, 0))

if(1):
    ax1.plot(chi1, kappa(chi1, 0.08, 0), 'm-', linewidth=2.2, label=r'$\kappa_{max}\rightarrow (c)$')
    ax1.plot(chi1, kappa(chi1, 0.08, 1), 'g-', linewidth=2.2, label=r'$\kappa_{max}\rightarrow (a)$')
    ax1.set_ylim([-0.5, 1])
    ax1.set_xlim([0, 1])
    ax1.set_title(r"$\eta = 0.08$")
    ax1.legend(loc=4,frameon=False)
    ax1.set_xlabel(r'$\chi_1$')
    ax1.set_ylabel(r'$\kappa$')

    ax2.plot(eta, kappa(0.98, eta, 0), 'm-', linewidth=2.2, label=r'$\kappa_{max}\rightarrow (c)$')
    ax2.plot(eta, kappa(0.98, eta, 1), 'g-', linewidth=2.2, label=r'$\kappa_{max}\rightarrow (a)$')
    ax2.set_ylim([-0.5, 1])
    ax2.set_xlim([eta[0], eta[-1]])
    ax2.set_title(r"$\chi_1 = 0.8$")
    ax2.legend(loc=4,frameon=False)
    ax2.set_xlabel(r'$\eta$')
    ax2.set_ylabel(r'$\kappa$')

GS1.tight_layout(fig)
plt.savefig('./immediate/kappa_bounds.pdf', bbox_inches='tight')
