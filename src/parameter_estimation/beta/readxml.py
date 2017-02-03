#===============================================================================
# Generates and modifies the XML injection file for multiple injections.
# Can modidy: mass1, mass2, eta, chirp, incl, spin1x, spin1y, spin1z, polarization.
# Parameters common to the injections are set in lalapps_inspinj.

#TODO: Any other parameters to replace?
# SM 2.10.2016
#===============================================================================

import xml.etree.ElementTree as ET
import numpy as np
import os

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

    'V_MASS1_RANGE' : [2.4, 50.0],
    'V_CHI1_RANGE'  : [0.20, 0.80],

    'V_KAPPA_RANGE'  : [-0.500, 0.999],
    'V_THETAJ_RANGE' : [0.001, 3.14],
	
	'FILE_NAME' 	 : 'injections.xml'
	}

#=======================================================================
# END CONTROL PANEL
#=======================================================================

#TODO: Fix parameteres to generate

def gen_xml(**options):
	os.system("lalapps_inspinj \
    --f-lower %r \
    --i-distr uniform \
    --m-distr componentMass \
    --min-mass1 %r \
	--max-mass1 %r \
    --min-mass2 1.4 \
    --max-mass2 1.4 \
    --coa-phase-distr uniform \
    --l-distr random \
    --waveform %r \
    --enable-spin \
    --min-spin1 %r \
	--max-spin1 %r \
    --min-spin2 0.0 \
    --max-spin2 0.0 \
    --d-distr uniform \
    --min-distance 600000 \
    --max-distance 600000 \
    --gps-start-time 1126258000 \
    --gps-end-time 1126259000 \
    --time-interval  0. \
    --time-step 500. \
    --seed 1990 \
    --amp-order 0 \
    --output ./%r" \
	
	%(
	options['F_MIN'], \
	options['V_MASS1_RANGE'][0], \
	options['V_MASS1_RANGE'][0], \
	options['APPROX'], \
	options['SPIN1_MIN'], \
	options['SPIN1_MAX'], \
	options['FILE_NAME'],	
	)
	
	return None

def parse_xml(injection_array, **options):
	tree = ET.parse('injections.xml')
	root = tree.getroot()

	headers = []
	for child in root:
		if child.attrib['Name'] == 'sim_inspiralgroup:sim_inspiral:table':
			for children in child:
				headers.append(children.attrib['Name'].split(':')[-1])

	for stream in root.iter('Stream'):
		if stream.attrib['Name'] == 'sim_inspiralgroup:sim_inspiral:table':
			process = stream.text.split('\n')
			for i, j in enumerate(process[1].split(',')):
					print "%s = %s"%(headers[i], j)

#write to the injection file.

#Condor submission follows.
