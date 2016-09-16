#=======================================================================
# Computes the OLVPs over the parameter space.
# Returns: Overlaps over thetaJ, Chi, eta and Kappa
# SM 9/ 16
# TODO: Add support for SpinTaylorF2
#=======================================================================

import numpy as np
from STF2_overlaps import compute_overlap
from numpy import sqrt
from time import localtime, strftime
import os

import STF2_vis_overlaps as vs
import STF2_vis_grid as vsg

def m1m2_to_mchirpeta(mass1, mass2):
        eta = mass1*mass2/(mass1+mass2)/(mass1+mass2)
        mc = (mass1+mass2) * pow(eta, 3./5.)
        return (mc, eta)

def mchirpeta_to_m1m2(mchirp, eta):
        mtot = mchirp * pow(eta, -3./5)
        fac = sqrt(1. - 4.*eta)
        return (mtot * (1. + fac) / 2., mtot * (1. - fac) / 2.)

#=======================================================================
# CONTROL PANEL: #TODO: Shift it outside the code.
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
    'BAND'   : None}

N = 10
M = 3

V_MASS1   = np.linspace(2.4, 50.0,  M)
V_CHI1    = np.linspace(0.5, 1.00,  M)

V_KAPPA   = np.linspace(-0.500, 0.999, N)
V_THETAJ  = np.linspace(0.001, 3.14, N)

generate_plots = 1

#=======================================================================
# Main loop for computations.
#=======================================================================

output_dir = "output-%s" %strftime("%Y_%m_%d_%H_%M_%S", localtime())

if not os.path.exists("../../output/datasets/%s" %output_dir):
    os.makedirs("../../output/datasets/%s" %output_dir)

OVLP = np.zeros((N, N, 9))

iter_O = 0

for _mass1 in xrange(M):
    for _chi1 in xrange(M): #generates a grid.

        iter_O = iter_O + 1

        options['M1']   = V_MASS1[_mass1]
        options['CHI1'] = V_CHI1[_chi1]

        iter_I   = 0
        MCHIRP, ETA = m1m2_to_mchirpeta(options['M1'], options['M2'])

        print 50*"-"
        print "eta: %r \t chi1: %r" %(ETA, V_CHI1[_chi1])
        print 50*"-" + "\n"

        filename = "overlaps_eta_%s_chi1_%s_N_%r.npz" \
        %('{:.2f}'.format(ETA),'{:.2f}'.format(V_CHI1[_chi1]) , N)

        for _thetaJ in xrange(N):
            for _kappa in xrange(N):    # generates each individual block.

                iter_I = iter_I + 1
                print "Iteration: [%r, %r] \t [%r, %r]" %(iter_I, N**2, iter_O, M**2)

                options['THETAJ'] = V_THETAJ[_thetaJ]
                options['KAPPA']  = V_KAPPA[_kappa]

                OVLP[_thetaJ][_kappa][:] = compute_overlap(**options)

        print "\n"

        #TODO: Do this better; this is a very crude way of doing this.

        np.savez("../../output/datasets/%s/%s" %(output_dir, filename),
                    DATE   = strftime("%Y-%m-%d %H:%M:%S", localtime()),
                    SNR_0F = OVLP[:, :, 0],
                    SNR_02 = OVLP[:, :, 1],
                    SNR_00 = OVLP[:, :, 2],
                    OLVP_0F_P2 = OVLP[:, :, 3],
                    OLVP_0F_P1 = OVLP[:, :, 4],
                    OLVP_0F_P0 = OVLP[:, :, 5],
                    OLVP_0F_M1 = OVLP[:, :, 6],
                    OLVP_0F_M2 = OVLP[:, :, 7],
                    OLVP_0F_P2P0 = OVLP[:, :, 8],
                    THETAJ = V_THETAJ,
                    KAPPA  = V_KAPPA,
                    CHI1   = V_CHI1[_chi1],
                    ETA    = ETA,
                    OPTIONS= options)

print "Finished creating datasets."

if generate_plots == 1:
    print "\nGenerating plots..."
    vs.visualize_OVLP(output_dir)
    vsg.visualize_OLVP_grid(output_dir)

