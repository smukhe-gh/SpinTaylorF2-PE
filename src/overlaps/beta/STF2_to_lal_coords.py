#===============================================================================
# SM 13/9/2016
# The following code is used to convert kappa, thetaJ to lal_coordinates.
# The code returns inclination angle, psi0, and spin components of M1.

# Note : This version of the code works on lal-version: 6.16.1.1. It may not work in earlier versions (you may have to remove or modify line 53) since the
# coordinate coversion implemented in LALSimInspiral.c is different in the earlier versions, namely: 6.14.0.1
# For a detailed review of the coordinates used in lal, see DCC LIGO-T1600045.
#===============================================================================

import numpy as np
import pycbc
from pycbc.types import TimeSeries, FrequencySeries, zeros
from pycbc.waveform import get_fd_waveform, get_td_waveform
from pycbc.waveform import td_approximants, fd_approximants
from lal import MTSUN_SI
import matplotlib.pyplot as plt

import lalsimulation
import lal
import numpy
from numpy import sqrt, double, complex128
from math import pow, log, cos, sin, acos, atan2

from pycbc.setuputils import pkg_config_header_strings
from pycbc.types import FrequencySeries, zeros
import pycbc.pnutils
from pycbc.waveform.utils import ceilpow2

import matplotlib.pyplot as plt

def m1m2_to_mchirpeta(m1, m2):
	eta    = m1*m2/(m1+m2)/(m1+m2)
	mchirp = (m1+m2) * np.pow(eta, 3./5.)
	return (mchirp, eta)

def mchirpeta_to_m1m2(mchirp, eta):
	mtotal = mchirp * np.pow(eta, -3./5)
	fac    = np.sqrt(1. - 4.*eta)
	return (mtot * (1. + fac) / 2., mtot * (1. - fac) / 2.)

def rotateY(lst, angle):
	sinR, cosR = np.sin(angle), np.cos(angle)
	return [lst[0]*cosR + lst[2]*sinR, lst[1], -lst[0]*sinR + lst[2]*cosR]

def rotateZ(lst, angle):
	sinR, cosR = np.sin(angle), np.cos(angle)
	return [lst[0]*cosR - lst[1]*sinR, lst[0]*sinR + lst[1]*cosR, lst[2]]

def to_lal_coords(m1, m2, chi1, kappa, thetaJ, psiJ, alpha0, f_inj):

    v0 = np.power(np.pi*MTSUN_SI*(m1+m2)*f_inj, 1./3.)
    gamma = m1*chi1*v0/m2
    denom = np.sqrt(1. + 2.*kappa*gamma + gamma*gamma)
    sinB, cosB = gamma*np.sqrt(1.-kappa*kappa)/denom, (1. + kappa*gamma)/denom
    sinA, cosA = np.sin(alpha0), np.cos(alpha0)

    lhat = [sinB*cosA, sinB*sinA, cosB]
    shat = [-sinB*cosA/gamma, -sinB*sinA/gamma, (kappa+gamma)/denom]

    # lhat = rotateZ(rotateY(lhat, thetaJ), psiJ)
    shat = rotateZ(rotateY(shat, thetaJ), psiJ)

    psi0 = np.arctan2(lhat[1], lhat[0])
    incl = np.arccos(lhat[2])
    shat = rotateZ(shat, -psi0)

    shat = rotateY(shat, -incl) # Additional rotation.

    spin = chi1*np.array(shat)  # Returns spin, not Shat.


    # print 'to_lal_coords: chi1\t', chi1
    # print 'to_lal_coords: kappa\t', kappa
    
    # print 'to_lal_coords: lhat\t', lhat
    # print 'to_lal_coords: spin\t', spin

    return (incl, psi0, spin, lhat)

def spintaylorf2_coords(**kwds):

    #SM: Changed nothing in this code. 

    f_lower = double(kwds['f_lower'])
    delta_f = double(kwds['delta_f'])
    distance = double(kwds['distance'])
    mass1 = double(kwds['mass1'])
    mass2 = double(kwds['mass2'])
    spin1x = double(kwds['spin1x'])
    spin1y = double(kwds['spin1y'])
    spin1z = double(kwds['spin1z'])
    phi0 = double(kwds['coa_phase'])               #Orbital Phase at coalescence
    phase_order = int(kwds['phase_order'])
    amplitude_order = int(kwds['amplitude_order'])
    inclination = double(kwds['inclination'])
    lnhatx = sin(inclination)
    lnhaty = 0.
    lnhatz = cos(inclination)
    psi = 0.

    tC= -1.0 / delta_f
    M = mass1 + mass2
    eta = mass1 * mass2 / (M * M)
    m_sec = M * lal.MTSUN_SI
    piM = lal.PI * m_sec

    vISCO = 1. / sqrt(6.)
    fISCO = vISCO * vISCO * vISCO / piM
    f_max = ceilpow2(fISCO)
    n = int(f_max / delta_f + 1)
    kmax = int(fISCO / delta_f)
    kmin = int(numpy.ceil(f_lower / delta_f))
    kmax = kmax if (kmax<n) else n

    #####Calculate the Orientation#####
    v0 = pow(piM *  kmin * delta_f,1./3)
    chi = sqrt(spin1x**2+spin1y**2+spin1z**2)
    kappa = (lnhatx*spin1x+lnhaty*spin1y+lnhatz*spin1z)/chi if (chi > 0.)  else 1.
    Jx0 = mass1*mass2*lnhatx/v0 + mass1*mass1*spin1x
    Jy0 = mass1*mass2*lnhaty/v0 + mass1*mass1*spin1y
    Jz0 = mass1*mass2*lnhatz/v0 + mass1*mass1*spin1z
    thetaJ = acos(Jz0 / sqrt(Jx0**2+Jy0**2+Jz0**2))
    psiJ = atan2(Jy0, -Jx0) # FIXME: check that Jy0 and Jx0 are not both 0

    return thetaJ, kappa, lnhatx, lnhaty, lnhatz

