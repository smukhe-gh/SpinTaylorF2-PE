# -----------------------------------------------------------------------
# Running cbcBayesPostProc.py
# -----------------------------------------------------------------------

#Throwing up errors -- Not a universal one. Other's don't seem to have
#reported such an error.

echo "============================================================================"
echo "Starting cbcBayesPostProc.py"
echo "============================================================================"

cbcBayesPostProc.py \
--inj ./injections.xml \
--snr ./lalinferencenest-0-V1H1L1-1126258000.0-0.hdf5_snr.txt \
--skyres 0.5 \
--outpath ./1126258000-0/V1H1L1 ./posterior_V1H1L1_1126258000-0.hdf5 \
--psdfiles ./lalinferencenest-0-V1H1L1-1126258000.0-0.hdf5V1-PSD.dat,./lalinferencenest-0-V1H1L1-1126258000.0-0.hdf5H1-PSD.dat,./lalinferencenest-0-V1H1L1-1126258000.0-0.hdf5L1-PSD.dat  \
--eventnum 0

