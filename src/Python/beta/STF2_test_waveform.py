#===============================================================================
# Generate SpinTaylorF2 waveform with sidebands.
# Soham M  9/2016
# TODO: Add support for generating SpinTaylorT2
#===============================================================================

import numpy as np
import pycbc
from pycbc.types import TimeSeries, FrequencySeries, zeros
from pycbc.waveform import get_fd_waveform, get_td_waveform
from pycbc.waveform import td_approximants, fd_approximants
from lal import MTSUN_SI

#TODO: Check these functions; taken from STF2_Precessing.py
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

#TODO: Check how and why this works. See LIGO-T1500606-v2.
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
    
    shat = rotateZ(shat, -psi0)
    shat = rotateY(shat, -incl) #IMP: New rotation added.
    
    # SM: modified to return the spin vector instead of Shat 
    spin = chi1*np.array(shat)
    
    return (incl, psi0, spin)
  
def generate_template(**options):
                            
    incl, psi0, spin1 = to_lal_coords(options['M1'], options['M2'], \
                        options['CHI1'], options['KAPPA'], options['THETAJ'], \
                        options['PSIJ'], options['ALPHA0'], options['F_INJ'])
        
    print 'm1', options['M1']
    print 'm2', options['M2']

    print 'thetaJ', options['THETAJ']
    print 'inclination', incl
    print 'psi0', psi0
    print 'spin', spin1


    hpluss, hcross = get_fd_waveform(
    
        approximant = options['APPROX'], #Input
    
        mass1 = options['M1'], #Input
        mass2 = options['M2'], #Input
        
        delta_f = options['DEL_F'],    #Input
        f_lower = options['F_MIN'],    #Input
       
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
    
    sin2Y, cos2Y = np.sin(2.*psi0), np.cos(2.*psi0)

    hp   = pycbc.DYN_RANGE_FAC*pycbc.DYN_RANGE_FAC*hpluss
    hc   = pycbc.DYN_RANGE_FAC*pycbc.DYN_RANGE_FAC*hcross
    freq = hp.get_sample_frequencies().data

    return freq, hp, hc

#====================================================
# Test zone
#====================================================

# options = {
#     'APPROX' : 'SpinTaylorF2',
#     'M1'     : 10.0,
#     'M2'     : 1.4,
#     'DEL_F'  : 1./256.,
#     'F_MIN'  : 20,
#     'PHI0'   : 0.001,
#     'BAND'   : None,
#     'CHI1'   : 0.5,
#     'KAPPA'  : 0.5,
#     'THETAJ' : 0.7853981633974483,
#     'PSIJ'   : 0.001,
#     'ALPHA0' : 0.001,
#     'F_INJ'  : 20
# }
#
# hp, hx = generate_template(**options)
# freq   = hp.get_sample_frequencies().data
