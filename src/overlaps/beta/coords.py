#=======================================================================
# Compare coordinate conversion.
# SM 17/ 16
#=======================================================================

import numpy as np 
import STF2_to_lal_coords as tlc

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

#=======================================================================
# BEGIN CONTROL PANEL
#=======================================================================

options = {

    'ALPHA0' : 0.001,
    'KAPPA'  : 0.5,
    'CHI1'   : 0.5,
    'PSIJ'   : 0.001,
    'M1'     : 10.0,
    'PHI0'   : 0.001,
    'M2'     : 1.4,
    'THETAJ' : 0.7853981633974483,
    'APPROX' : 'SpinTaylorF2',
    'DEL_F'  : 1./256.,
    'F_MIN'  : 50.,
    'F_INJ'  : 50.,
    'F_MAX'  : 2000.,
    'BAND'   : None,

    'N'      : 50,
    'M'      : 3,

    'V_MASS1_RANGE' : [2.4, 50.0],
    'V_CHI1_RANGE'  : [0.20, 0.80],

    'V_KAPPA_RANGE'  : [-0.500, 0.999],
    'V_THETAJ_RANGE' : [0.001, 3.14]
    }

#=======================================================================
# END CONTROL PANEL
#=======================================================================


V_THETAJ = np.linspace(options['V_THETAJ_RANGE'][0], options['V_THETAJ_RANGE'][1], options['N'])
V_KAPPA  = np.linspace(options['V_KAPPA_RANGE'][0], options['V_KAPPA_RANGE'][1],  options['N'])

ZEROS = np.zeros((options['N'], options['N']))

for _thetaJ in xrange(options['N']):
    for _kappa in xrange(options['N']):

        options['THETAJ'] = V_THETAJ[_thetaJ]
        options['KAPPA']  = V_KAPPA[_kappa]

        # get inclination, psi0, and spin from to_lal_coords
        incl, psi0, spin, lhat, chi1, v0 = tlc.to_lal_coords(options['M1'], options['M2'], options['CHI1'], \
        	options['KAPPA'], options['THETAJ'], \
        	options['PSIJ'], options['ALPHA0'], options['F_INJ'])

        # set coords
        coords = {
        'f_lower' : options['F_MIN'],
        'delta_f' : options['DEL_F'],
        'distance' : 400.0,
        'mass1' : options['M1'],
        'mass2' : options['M2'],
        'spin1x' : spin[0],
        'spin1y' : spin[1],
        'spin1z' : spin[2],
        'coa_phase' : options['PHI0'],
        'phase_order' : 7,
        'amplitude_order' : 0,
        'inclination' : incl,
        'lnhat' : lhat,
        'v0' : v0
        }

        #pass these values to spintaylorF2 PyCBC code. 
        kappa, thetaJ, lnhat, spin_pycbc, chi_pycbc = tlc.spintaylorf2(**coords)

        # compute the difference:
        Q1 = kappa
        Q2 = V_KAPPA[_kappa]
        # Q1 = thetaJ
        # Q2 = V_THETAJ[_thetaJ]
        # Q1 = spin[2]
        # Q2 = spin_pycbc[2]
        # Q1 = lhat[0]
        # Q2 = lnhat[0]

        # Q1 = chi1
        # Q2 = chi_pycbc

        # if Q2 < 0 or Q1 < 0:
        #     print 'N'

        
        ZEROS[_thetaJ, _kappa] = np.abs(Q1 - Q2)
       

plt.figure(1)
plt.imshow(np.flipud(ZEROS))
plt.colorbar()
plt.show()
print ZEROS[ZEROS != 0.0]

# plt.plot(ZEROS[:,15])
# print np.argmax(ZEROS[:,15][5:45])
# plt.show()



