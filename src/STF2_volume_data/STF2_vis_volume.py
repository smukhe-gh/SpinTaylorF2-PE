#=================================================================
# Volume slice visualization.
# SM 9/2016
# FIXME: Check if colorbar changes with the plot.
# FIMME: Change the plotting slice.
#=================================================================

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

"""
The code below should allow you to explore the [chi1, kappa, thetaJ] 3D space.
The funtion cube_show_slider takes the entire npz file as input and the axis along
which you would like to slice the volume. In our case:
    axis 0: chi1
    axis 1: theta_J
    axis 2: kappa
Therefore, you can set the axis appropriately to see the variation in that plane.
Also, you can set what to plot in the first line of cube_show_slider. I've commented
out the line where I plot the ratio between SNR_02, and SNR_00.

The following data is available in the npz file:
    DATE         = Date and time when the simulation was started
    SNR_0F       = Full SpinTaylorF2 SNR
    SNR_02       = SNR of m=2 sideband
    SNR_00       = SNR of m=0 sideband
    OLVP_0F_P2   = Overlap of Full SpinTaylorF2 with m=2
    OLVP_0F_P1   = Overlap of Full SpinTaylorF2 with m=1
    OLVP_0F_P0   = Overlap of Full SpinTaylorF2 with m=0
    OLVP_0F_M1   = Overlap of Full SpinTaylorF2 with m=-1
    OLVP_0F_M2   = Overlap of Full SpinTaylorF2 with m=-2
    OLVP_0F_P2P0 = Overlap of Full SpinTaylorF2 with m=2 + m=0
    THETAJ       = thetaJ axis
    KAPPA        = kappa axis
    CHI1         = chi1 axis
    ETA          = mass ratio
    OPTIONS      = the entire set of options used to create the dataset

Also, note that this data was created with 50 x 50 x 50 points, and therefore the
resolution is not that great.

Additionally, you can also choose the mass ratio you would want to look at: just
load the required file.
"""

def cube_show_slider(cube, axis=0, **kwargs):

    """
    Set what you want to plot here.
    """
    # cube = data['SNR_00']
    cube = np.divide(data['SNR_02'], data['SNR_00'])

    # add to mask regions
    for index, value in np.ndenumerate(cube):
    #    # if value < 1:
       if np.abs(value - 0.98) > 0.1:
           cube[index] = np.nan


    fig = plt.figure()
    ax = plt.subplot(111)
    fig.subplots_adjust(left=0.25, bottom=0.25)

    if axis == 0:
        im = cube[0, :, :]  #select the thetaJ, kappa plane
    elif axis == 1:
        im = cube[:, 0, :]  #select the chi1, kappa plane
    elif axis == 2:
        im = cube[:, :, 0]  #select the chi1, thetaJ plane


    x = np.array([0, 10, 20, 30, 40, 49])
    y = np.array([1, 10, 20, 30, 40, 50])
    ax.set_yticks(x)
    ax.set_xticks(y)

    if axis == 0:
        ax.set_xticklabels([r"$%.2f$"%data['KAPPA'][i] for i in x])
        ax.set_yticklabels([r"$%.2f$"%data['THETAJ'][50-i] for i in y])
    elif axis == 1:
        ax.set_xticklabels([r"$%.2f$"%data['CHI1'][i] for i in x])
        ax.set_yticklabels([r"$%.2f$"%data['THETAJ'][50-i] for i in y])
    elif axis == 2:
        ax.set_xticklabels([r"$%.2f$"%data['CHI1'][i] for i in x])
        ax.set_yticklabels([r"$%.2f$"%data['KAPPA'][50-i] for i in y])

    if axis == 0:
        ax.set_xlabel(r"$\kappa$")
        ax.set_ylabel(r"$\theta_{J}$")
    elif axis == 1:
        ax.set_xlabel(r"$\chi_{1}$")
        ax.set_ylabel(r"$\theta_{J}$")
    elif axis == 2:
        ax.set_xlabel(r"$\chi_{1}$")
        ax.set_ylabel(r"$\kappa$")


    # display image
    l = ax.imshow(np.flipud(im))
    fig.colorbar(l)

    # define slider
    axcolor = 'white'
    ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)

    if axis == 0:
        slider = Slider(ax, r"$\chi_{1}$", 0, 49,
                    valinit=0, valfmt='%i')
    elif axis == 1:
        slider = Slider(ax, r"$\kappa$", 0, 49,
                    valinit=0, valfmt='%i')
    elif axis == 2:
        slider = Slider(ax, r"$\theta_{J}$", 0, 49,
                    valinit=0, valfmt='%i')

    def update(val):
        ind = int(slider.val)
        if axis == 0:
            l.set_data(np.flipud(cube[ind, :, :]), **kwargs)
        elif axis == 1:
            l.set_data(np.flipud(cube[:, ind, :]), **kwargs)
        elif axis == 2:
            l.set_data(np.flipud(cube[:, :, ind]), **kwargs)

        fig.canvas.draw()

    slider.on_changed(update)
    plt.show()

data = np.load('../../output/datasets/output-2016_10_19_19_24_50/overlaps_eta_0.04_chi1_0.80_N_50.npz')
cube_show_slider(data, axis=1)
