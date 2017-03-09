from numpy import load
import numpy as np
from pycbc import DYN_RANGE_FAC
from pycbc.psd import from_txt
from pycbc.types import FrequencySeries, zeros
import sys
import os.path

psd_loc = os.path.expanduser("PSDs")
#psd_loc = os.path.expanduser("/home/haris/Fisher/5VarRand/PSDs")

f_low = 10.

def load_psd(psd_choice, f_max, delta_f):
	""" Load PSD, options: early, HPZD """
	n_samples = int(f_max/delta_f)+1
        sample_rate = int(2*f_max)
	asd_txt_file = os.path.join(psd_loc, {'HPZD': "ZERO_DET_high_P.txt", 'early': "early_gaussian_asd.dat"}[psd_choice])
	psd_cache_file = os.path.join(psd_loc, "%s_%u_%u.npy" % (psd_choice, f_max, sample_rate))

	try: # For speed, use the numpy file
		temp = load(psd_cache_file)[:,1]
		psd = FrequencySeries(temp, delta_f=delta_f, dtype=np.double)
	except IOError:
		psd = from_txt(asd_txt_file, n_samples, delta_f, f_low)
		psd *= DYN_RANGE_FAC ** 2
		psd.save(psd_cache_file)

	return psd

