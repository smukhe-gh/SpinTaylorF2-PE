#Code to make injection files and get posterior samples.
import os

os.system('lalapps_inspinj \
--f-lower 20 \
--amp-order 0 \
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
--time-interval  0. \
--time-step 10. \
--seed 1990 \
--output ./test_injections.xml')

"""
lalapps_inspinj [options]
The following options are recognized.  Options not surrounded in []
are required. Defaults are shown in brackets
 [--help ]                 display this message
 [--verbose]               print progress information
 [--user-tag] usertag      set the usertag 
 [--output ] name          overwrite the standard file naming convention
 [--write-compress]        write a compressed xml file

Waveform details:
 [--seed] randomSeed       seed for random number generator (default : 1)
  --f-lower freq           lower cut-off frequency.
  --waveform wfm           set waveform type to wfm
  --amp-order              set PN order in amplitude

Time distribution information:
  --gps-start-time start   GPS start time for injections
  --gps-end-time end       GPS end time for injections
  --ipn-gps-time IPNtime   GPS end time for IPN trigger
  --t-distr timeDist       set the time step distribution of injections
                           fixed: fixed time step
                           uniform: uniform distribution
                           exponential: exponential distribution for Poisson process
  [--time-step] step       space injections by average of step seconds
                           (suggestion : 2630 / pi seconds)
  [--time-interval] int    distribute injections in an interval, int s
                           (default : 0 seconds)

Source distribution information:
  --l-distr  locDist       set the source location distribution,
                           locDist must be one of:
                           source: use locations from source-file
                           exttrig: use external trigger file
                           random: uses random locations
                           fixed: set fixed location
                           ipn: random locations from IPN skypoints
 [--longitude] longitude   read longitude if fixed value (degrees)
 [--latitude] latitude     read latitude if fixed value (degrees)
 [--d-distr] distDist      use a distribution over physical distance
                           source: take distance from galaxy source file
                           uniform: uniform distribution in distance
                           distancesquared: uniform distribution in distance^2
                           log10: uniform distribution in log10(d) 
                           volume: uniform distribution in volume
                           sfr: distribution derived from the SFR
 [--min-distance] DMIN     set the minimum (chirp) distance to DMIN kpc
 [--max-distance] DMAX     set the maximum (chirp) distance to DMAX kpc
                           min/max distance required if d-distr not 'source'
 [--source-file] sources   read source parameters from sources
                           requires enable/disable milkyway
 [--sourcecomplete] distance 
                           complete galaxy catalog out to distance (kPc)
 [--make-catalog]          create a text file of the completed galaxy catalog
 [--enable-milkyway] lum   enables MW injections, set MW luminosity
 [--disable-milkyway]      disables Milky Way injections
 [--dchirp-distr]          use a distribution over chirp distance
                           (normalized to a 1.4,1.4 Msun binary)
 [--z-distr]               use a distribution over redshift
                           currently only 'sfr' is supported
 [--local-rate] rho        set the local coalescence rate for --z-distr sfr
                           (suggestion: 1 per Mpc^3 per Myr)
 [--min-z]                 set the minimum redshift: at least 0.2 for sfr
 [--max-z]                 set the maximum redshift: at most 1.0 for sfr
 [--snr-distr]             use a distribution over expected (optimal) network SNR
                           uniform: uniform in SNR, log10: uniform in log10(SNR)
                           volume: uniform in 1/SNR^3
                           ( Setting max-snr == min-snr will allow you to choose a fixed SNR )
 [--ninja-snr]             use a NINJA waveform SNR calculation (if not set, use LALSimulation)
 [--min-snr] SMIN          set the minimum network snr
 [--max-snr] SMAX          set the maximum network snr
 [--min-coinc-snr] sm      Set the minimum SNR in two IFOs. Neglected if a single IFO is used
 [--ligo-psd] filename     Ascii, tab-separated file of frequency, value pairs to use for LIGO PSD in snr computation
 [--ligo-fake-psd] PSD     LALsimulation PSD fit to use instead of a file. Allowed values: LALLIGO, LALAdLIGO
 [--ligo-start-freq] freq  Frequency in Hz to use for LIGO snr computation
 [--virgo-psd] filename    Ascii, tab-separated file of frequency, value pairs to use for Virgo PSD in snr computation
 [--virgo-fake-psd] PSD    LALsimulation PSD fit to use instead of a file. Allowed values: LALVirgo, LALAdVirgo
 [--virgo-start-freq] freq Frequency in Hz to use for Virgo snr computation
 [--ifos] ifos             Comma-separated list of ifos to include in network SNR

  --i-distr INCDIST        set the inclination distribution, must be either
                           uniform: distribute uniformly over arccos(i)
                           gaussian: gaussian distributed in (i)
                           fixed: no distribution, fixed values of (i)
 [--polarization] psi      set the polarization angle for all injections (degrees)
 [--incl-std]  inclStd     std dev for gaussian inclination dist
 [--fixed-inc]  fixed_inc  value for the fixed inclination angle (in degrees) if '--i-distr fixed' is chosen.
 [--max-inc]  max_inc      value for the maximum inclination angle (in degrees) if '--i-distr uniform' is chosen. 
 [--coa-phase-distr] cDist set the coalescence phase distribution,
                           cDist must be one of:
                           uniform: use random, uniformly distributed coalescence phase [default]
                           fixed: set fixed coalescence phase
 [--fixed-coa-phase] phase set the coalescence phase (in degrees) for all injections if --coa-phase-distr=fixed
 [--ipn-file] ipnskypoints read IPN sky points from file
 [--exttrig-file] exttrig  XML file containing external trigger

Mass distribution information:
  --m-distr massDist       set the mass distribution of injections
                           must be one of:
                           source: using file containing list of mass pairs
                           nrwaves: using xml file with list of NR waveforms
                           (requires setting max/min total masses)
                           totalMass: uniform distribution in total mass
                           componentMass: uniform in m1 and m2
                           gaussian: gaussian mass distribution
                           log: log distribution in component mass
                           totalMassRatio: uniform distribution in total mass and
                           mass ratio m1 /m2
                           logTotalMassUniformMassRatio: log distribution in total mass
                           and uniform in mass ratio
                           totalMassFraction: uniform distribution in total mass and
                           in m1 /(m1+m2)
                           m1m2SquareGrid: component masses on a square grid
                           fixMasses: fix m1 and m2 to specific values
 [--ninja2-mass]           use the NINJA 2 mass-selection algorithm
 [--real8-ninja2]          when distributing by SNR for NINJA2, assume frames are REAL8
 [--mass-file] mFile       read population mass parameters from mFile
 [--nr-file] nrFile        read mass/spin parameters from xml nrFile
 [--min-mass1] m1min       set the minimum component mass to m1min
 [--max-mass1] m1max       set the maximum component mass to m1max
 [--min-mass2] m2min       set the min component mass2 to m2min
 [--max-mass2] m2max       set the max component mass2 to m2max
 [--min-mtotal] minTotal   sets the minimum total mass to minTotal
 [--max-mtotal] maxTotal   sets the maximum total mass to maxTotal
 [--fixed-mass1] fixMass1  set mass1 to fixMass1
 [--fixed-mass2] fixMass2  set mass2 to fixMass2
 [--fixed-mass2] fixMass2  set mass2 to fixMass2
 [--max-mtotal] maxTotal   sets the maximum total mass to maxTotal
 [--fixed-mass1] fixMass1  set mass1 to fixMass1
 [--fixed-mass2] fixMass2  set mass2 to fixMass2
 [--mean-mass1] m1mean     set the mean value for mass1
 [--stdev-mass1] m1std     set the standard deviation for mass1
 [--mean-mass2] m2mean     set the mean value for mass2
 [--stdev-mass2] m2std     set the standard deviation for mass2
 [--min-mratio] minr       set the minimum mass ratio
 [--max-mratio] maxr       set the maximum mass ratio
 [--mass1-points] m1pnt    set the number of grid points in the m1 direction if '--m-distr=m1m2SquareGrid'
 [--mass2-points] m2pnt    set the number of grid points in the m2 direction if '--m-distr=m1m2SquareGrid'

Spin distribution information:
  --disable-spin           disables spinning injections
  --enable-spin            enables spinning injections
                           One of these is required.
  [--spin-gaussian]        enable gaussian spin distribution
  --aligned                enforces the spins to be along the direction
                           of orbital angular momentum. Spin z-components are the only non-vanishing (unless '--axis-choice view' convention is chosen)
  [--axis-choice] choice   frame axis choice: 'angmomentum' (default) or 'view' to define convention for spin aligned case
  [--min-spin1] spin1min   Set the minimum spin1 to spin1min (0.0)
  [--max-spin1] spin1max   Set the maximum spin1 to spin1max (0.0)
  [--mean-spin1] spin1mean Set the mean for |spin1| distribution
  [--stdev-spin1] spin1std Set the standard deviation for |spin1|
  [--min-spin2] spin2min   Set the minimum spin2 to spin2min (0.0)
  [--max-spin2] spin2max   Set the maximum spin2 to spin2max (0.0)
  [--mean-spin2] spin2mean Set the mean for |spin2| distribution
  [--stdev-spin2] spin2std Set the standard deviation for |spin2|
  [--min-kappa1] kappa1min Set the minimum cos(S1.L_N) to kappa1min (-1.0)
  [--max-kappa1] kappa1max Set the maximum cos(S1.L_N) to kappa1max (1.0)
  [--min-abskappa1] abskappa1min 
                           Set the minimum absolute value of cos(S1.L_N)
                           to abskappa1min (0.0)
  [--max-abskappa1] abskappa1max 
                           Set the maximum absolute value of cos(S1.L_N) 
                           to abskappa1max (1.0)

Tapering the injection waveform:
  [--taper-injection] OPT  Taper the inspiral template using option OPT
                            (start|end|startend) 
  [--band-pass-injection]  sets the tapering method of the injected waveform

"""