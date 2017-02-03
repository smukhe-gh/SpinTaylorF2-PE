#-----------------------------------------------------------------------
# Creating an injection using lalapps_inspinj 
#-----------------------------------------------------------------------

# Note: LALInference throws up errors with SpinTaylorF2. Can we fix this?

echo "============================================================================"
echo "Starting lalapps_inspinj"
echo "============================================================================"

lalapps_inspinj \
--f-lower 10 \
--i-distr uniform \
--m-distr fixMasses \
--fixed-mass1 50 \
--fixed-mass2 50 \
--coa-phase-distr uniform \
--l-distr random \
--waveform IMRPhenomPv2threePointFivePN \
--disable-spin \
--d-distr uniform \
--min-distance 600000 \
--max-distance 1000000 \
--gps-start-time 1126258000 \
--gps-end-time 1126258010 \
--time-interval  0.5 \
--time-step 5.0 \
--seed 1990 \
--amp-order 0 \
--output ../injections/injections.xml

#-----------------------------------------------------------------------
# Running lalinference_nest
#-----------------------------------------------------------------------

echo "============================================================================"
echo "Starting lalapps_chriplen" 
echo "============================================================================"

lalapps_chirplen --m1 50 --m2 50 --flow 10

echo "============================================================================"
echo "Starting lalinference_nest" 
echo "============================================================================"

lalinference_nest \
--L1-flow 10 \
--approx IMRPhenomPv2threePointFivePN \
--psdlength 63 \
--inj ../injections/injections.xml \
--V1-cache LALSimAdVirgo \
--nlive 1024 \
--V1-timeslide 0 \
--comp-max 80.0 \
--margphi  \
--srate 2048 \
--event 0 \
--seglen 8 \
--L1-channel L1:LDAS-STRAIN \
--H1-timeslide 0 \
--trigtime 1126258000 \
--tol 1.0 \
--psdstart 1126257920 \
--H1-cache LALSimAdLIGO \
--progress  \
--L1-timeslide 0 \
--H1-channel H1:LDAS-STRAIN \
--V1-channel V1:h_16384Hz \
--comp-min 30.0 \
--resume  \
--disable-spin  \
--V1-flow 10 \
--outfile ../nested_samples/lalinferencenest-0-V1H1L1-1126258000.0-0.hdf5 \
--L1-cache LALSimAdLIGO \
--randomseed 773032349 \
--dataseed 1234 \
--maxmcmc 2048 \
--H1-flow 10 --ifo V1  \
--ifo H1  \
--ifo L1  

#-----------------------------------------------------------------------
# Running lalapps_nest2pos
#-----------------------------------------------------------------------

echo "============================================================================"
echo "Starting lalapps_nest2pos" 
echo "============================================================================"

lalapps_nest2pos \
--pos ../posteriors/posterior_V1H1L1_1126258000-0.hdf5 ../nested_samples/lalinferencenest-0-V1H1L1-1126258000.0-0.hdf5





