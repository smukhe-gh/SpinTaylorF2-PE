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
    'F_MIN'  : 20.,
    'F_INJ'  : 20.,
    'F_MAX'  : 2000.,
    'BAND'   : None,

    'N'      : 10,
    'M'      : 3,

    'V_MASS1_RANGE' : [2.4, 50.0],
    'V_CHI1_RANGE'  : [0.20, 0.80],

    'V_KAPPA_RANGE'  : [-0.500, 0.999],
    'V_THETAJ_RANGE' : [0.001, 3.14]
    }

#=======================================================================
# END CONTROL PANEL
#=======================================================================



#======================================================================================

V_THETAJ = np.linspace(options['V_THETAJ_RANGE'][0], \
options['V_THETAJ_RANGE'][1], options['N'])
V_KAPPA  = np.linspace(options['V_KAPPA_RANGE'][0],  \
options['V_KAPPA_RANGE'][1],  options['N'])

ZEROS_thetaJ = np.zeros((options['N'], options['N']))
ZEROS_kappa  = np.zeros((options['N'], options['N']))

ZEROS_lhatx  = np.zeros((options['N'], options['N']))
ZEROS_lhaty  = np.zeros((options['N'], options['N']))
ZEROS_lhatz  = np.zeros((options['N'], options['N']))

for _thetaJ in xrange(options['N']):
    for _kappa in xrange(options['N']):

        options['THETAJ'] = V_THETAJ[_thetaJ]
        options['KAPPA']  = V_KAPPA[_kappa]

        print '\n'
        
        incl, psi0, spin, lhat = tlc.to_lal_coords(options['M1'], options['M2'], options['CHI1'], \
        	options['KAPPA'], options['THETAJ'], \
        	options['PSIJ'], options['ALPHA0'], options['F_INJ'])

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
        'inclination' : incl
        }

        thetaJ, kappa, lnhatx, lnhaty, lnhatz = tlc.spintaylorf2_coords(**coords)

        ZEROS_thetaJ[_kappa, _thetaJ] = thetaJ - V_THETAJ[_thetaJ]
        ZEROS_kappa[_kappa, _thetaJ]  = kappa  - V_KAPPA[_kappa]

        ZEROS_lhatx[_kappa, _thetaJ] = int(lnhatx - lhat[0])
        ZEROS_lhaty[_kappa, _thetaJ] = int(lnhaty - lhat[1])
        ZEROS_lhatz[_kappa, _thetaJ] = int(lnhatz - lhat[2])


plt.figure(1)
plt.subplot(1,3,1)
plt.contourf(V_KAPPA, V_THETAJ, ZEROS_lhatx)
plt.colorbar()
plt.title('lhaty')
plt.xlabel('kappa')
plt.ylabel('thetaJ')

plt.subplot(1,3,2)
plt.contourf(V_KAPPA, V_THETAJ, ZEROS_lhaty)
plt.colorbar()
plt.title('lhaty')

plt.subplot(1,3,3)
plt.contourf(V_KAPPA, V_THETAJ, ZEROS_lhatz)
plt.colorbar()
plt.title('lhatz')

plt.show()



