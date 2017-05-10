#===============================================================================
# Generate SpinTaylorF2 waveform with sidebands.
# Soham M  9/2016
#===============================================================================

import numpy as np
import pycbc
from pycbc.types import TimeSeries, FrequencySeries, zeros
from pycbc.waveform import get_fd_waveform, get_td_waveform
from pycbc.waveform import td_approximants, fd_approximants
from lal import MTSUN_SI
import STF2_to_lal_coords as coords

def generate_template(**options):

    incl, psi0, spin1 = coords.to_lal_coords(options['M1'], options['M2'], \
                        options['CHI1'], options['KAPPA'], options['THETAJ'], \
                        options['PSIJ'], options['ALPHA0'], options['F_INJ'])

    # FIXME: Number of samples of the waveform.
    nsamples  = int((options['F_MAX'])/options['DEL_F']) + 1

    hpluss, hcross = get_fd_waveform(

        approximant = options['APPROX'], #Input

        mass1 = options['M1'], #Input
        mass2 = options['M2'], #Input

        delta_f = options['DEL_F'],    #Input
        f_lower = options['F_MIN'],    #Input
        f_final = options['F_MAX'],    #Input, can be turned off.

        distance    = 400.0,              #Default
        inclination = incl,               #to_lal_coords
        coa_phase   = options['PHI0'],    #Input

        spin1x = spin1[0],    #to_lal_coords
        spin1y = spin1[1],    #to_lal_coords
        spin1z = spin1[2],    #to_lal_coords
        spin2x = 0.0,         #Default
        spin2y = 0.0,         #Default
        spin2z = 0.0,         #Default

        phase_order     = 7,   #Default
        spin_order      = 6,   #Default
        amplitude_order = 0,   #Default

        sideband = options['BAND']   #Input
        )

    #--------------------------------------------------
    # NOTE: We're only plotting the hpluss strain here.
    #--------------------------------------------------

    # sin2Y, cos2Y = np.sin(2.*psi0), np.cos(2.*psi0)

    # hp   = pycbc.DYN_RANGE_FAC*hpluss
    # hc   = pycbc.DYN_RANGE_FAC*hcross

    # waveform =  hp*cos2Y + hc*sin2Y
    # waveform.resize(nsamples)

    return hpluss, hcross

