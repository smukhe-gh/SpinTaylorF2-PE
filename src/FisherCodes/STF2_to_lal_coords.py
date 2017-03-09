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

    lhat = rotateZ(rotateY(lhat, thetaJ), psiJ)
    shat = rotateZ(rotateY(shat, thetaJ), psiJ)

    psi0 = np.arctan2(lhat[1], lhat[0])
    incl = np.arccos(lhat[2])

    # lhat = rotateZ(shat, -psi0) # FIXME: I don't think this is required.
    
    shat = rotateZ(shat, -psi0)
    shat = rotateY(shat, -incl)

    print "to_lal_coords shat", shat
    
    spin = chi1*np.array(shat)  # Returns spin, not Shat.

    return (incl, psi0, spin)



