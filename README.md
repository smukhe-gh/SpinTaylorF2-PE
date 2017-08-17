# Introduction
The code is designed as a wrapper to work with PyCBC and LALSuite functionalities to compute the SNR and overlaps of a waveform model, in particular, SpinTaylorF2. The code is written in Python, and uses process-level parallelization (using the joblib library) to compute things a bit faster. 
    
### Note
Many of the functions that are used in the code are now redundant after the recent updates to PyCBC which added code to handle the different sidebands of SpinTaylorF2 waveform. However, the code is still quite useful to run systematic overlap or SNR computations for a multi-dimensional parameter space.
 

