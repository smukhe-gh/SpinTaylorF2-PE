import matplotlib.pyplot as plt
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
        # "figure.figsize": (15.3, 10.3*goldenratio),
"figure.figsize": (13.3, 7.3),
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "text.usetex": True})

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

if(0):
    plt.plot(chi1, kappa(chi1, 0.08, 0), 'b--', label=r'$\kappa_{max}\rightarrow C$')
    plt.plot(chi1, kappa(chi1, 0.08, 1), 'k-', label=r'$\kappa_{max}\rightarrow A$')
    plt.fill_between(chi1, -0.5, kappa(chi1, 0.08, 0), color='0.65')
    plt.fill_between(chi1, kappa(chi1, 0.08, 0), kappa(chi1, 0.08, 1), color='b', alpha=0.4)
    plt.ylim([-0.5, 1])
    plt.xlim([0, 1])
    plt.legend(loc=0,frameon=False)
    plt.xlabel(r'$\chi_1$')
    plt.ylabel(r'$\kappa$')
    plt.savefig("./temp_kappa.pdf", bbox_inches='tight')

if(1):
    plt.plot(eta, kappa(0.98, eta, 0), 'b--', label=r'$\kappa_{max}\rightarrow C$')
    plt.plot(eta, kappa(0.98, eta, 1), 'k-', label=r'$\kappa_{max}\rightarrow A$')
    plt.fill_between(eta, -0.5, kappa(0.98, eta, 0), color='0.65')
    plt.fill_between(eta, kappa(0.98, eta, 0), kappa(0.98, eta, 1), color='b', alpha=0.4)
    plt.ylim([-0.5, 1])
    plt.xlim([eta[0], eta[-1]])
    plt.legend(loc=0,frameon=False)
    plt.xlabel(r'$\eta$')
    plt.ylabel(r'$\kappa$')
    plt.savefig("./temp_kappa_eta.pdf", bbox_inches='tight')