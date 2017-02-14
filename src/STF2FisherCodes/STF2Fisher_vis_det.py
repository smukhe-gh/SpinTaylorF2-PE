import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import glob

matplotlib.rcParams.update({
        "font.size": 20.0,
        "axes.titlesize": 20.0,
        "axes.labelsize": 20.0,
        "xtick.labelsize": 20.0,
        "ytick.labelsize": 20.0,
        "legend.fontsize": 20.0,
        "figure.figsize": (14.3, 12.3),
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "text.usetex": True
})


# files = glob.glob("../../output/datasets/output-2017_02_10_11_59_30/fisher*")   #m=0
files = glob.glob("../../output/datasets/output-2017_02_13_11_02_08/fisher*")   #m=None
# files = glob.glob("../../output/datasets/output-2017_02_10_20_20_58/fisher*")   #m=20

plt.cm = plt.get_cmap('viridis')
thetaJ = np.linspace(0.2, np.pi - 0.2, 50)
kappa = np.linspace(-0.5, 0.8, 50)
DAT = []

for i, file in enumerate(files):
    data = np.load(file)
    if i==0:
        print data["OPTIONS"]
    DAT.append(np.log(np.sqrt(np.abs(data['FDET']))))

    # plt.xlabel(r"$\kappa$")
    # plt.ylabel(r"$\theta_{J}$")
    # plt.title(r'$\chi_{1}=%1.2f$'%data['CHI1'] + r' $\eta=%1.2f$'%data['ETA'] + r" $m=0, \log\sqrt{|\rm det|}$")
    # plt.contourf(thetaJ, kappa, np.log(np.sqrt(np.abs(data['FDET']))),cmap='viridis')
    # plt.colorbar()
    # plt.savefig('./immediate/fisher_det_m0_plane_%r.pdf'%i, bbbox_inches='tight')
    # plt.close()


fig, axes = plt.subplots(nrows=2, ncols=2)
DAT = np.array(DAT)
DAT[np.isnan(DAT)] = np.mean(DAT[np.isfinite(DAT)])

cmax = np.amax(DAT)
cmin = np.amin(DAT)

print "det_max:", cmax
print "det_min:", cmin

for i, ax in enumerate(axes.flat):
    data = np.load(files[i])
    if i==0:
        print data["OPTIONS"]

    ax.set_xlabel(r"$\kappa$")
    ax.set_ylabel(r"$\theta_{J}$")

    thetaJ = np.linspace(0.2, np.pi - 0.2, 10)
    kappa = np.linspace(-0.5, 0.8, 10)
    ax.set_title(r'$\chi_{1}=%1.2f$'%data['CHI1'] + r' $\eta=%1.2f$'%data['ETA'])
    im = ax.contourf(thetaJ, kappa, np.log(np.sqrt(np.abs(data['FDET']))),cmap='viridis', vmin=cmin, vmax=cmax)
    plt.subplots_adjust(wspace=0.3, hspace=0.3)

fig.suptitle(r"$m=\rm None$")
cb = fig.colorbar(im, ax=axes.ravel().tolist(), spacing='proportional')
plt.savefig('./immediate/fisher_det_m_None_grid.pdf', bbbox_inches='tight')
plt.close()