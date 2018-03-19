# SpinTaylorF2-PE
The code is designed as a wrapper to work with PyCBC [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1183449.svg)](https://doi.org/10.5281/zenodo.1183449) and LALSuite 
functionalities to compute the SNR and overlaps of a waveform model, 
in particular, SpinTaylorF2 [(1)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.89.044021). 
The code is written in Python, and uses 
process-level parallelization (using the joblib library) to compute 
things a bit faster. 
    
Many of the functions that are used in the code are now redundant after 
the recent updates to PyCBC which added code to handle the different 
sidebands of SpinTaylorF2 waveform. However, the code is still quite 
useful to run systematic overlap or SNR computations for a multi-dimensional parameter space.
 
---
(1) A. Lundgren and R. O’Shaughnessy, [https://journals.aps.org/prd/abstract/10.1103/PhysRevD.89.044021](Phys. Rev. D 89, 044021 (2014).)
