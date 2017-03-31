#======================================================================
# Soham M 3/2017
# Code to test and print output from STF2Fisher_SB_parallel.py
#======================================================================

import numpy as np

DATA = np.load("FISHER_SB_parallel_N2_test_SB.npz")
print "SHAPE: ", DATA["FISHER_DET"].shape

for slice in DATA["FISHER_DET"]:
	print slice