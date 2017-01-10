Readme file for 'SpinTaylorF2_PE'
Soham M, 18/2016

The code can:
    Compute SNR for SpinTaylorF2 : full waveform as well as sidebands.
    Compute overlaps between any two waveforms.
    Generate overlap plots over chi1 and eta. 

A few things to note:
    Run src/Python/STF2_compute.py to generate output.
    If you want to run the code in parallel (emabarassingly parallel for loop)
    you need the joblib Python module.

    Output is generated in the top directory, i.e. SpinTaylorF2_PE, with timestamped 
    folders. You'll find both the datasets and the plots here.

    To just plot existing data, pass the output directory directly to the fuction in 
    STF2_vis_grid.py and STF2_vis_overlaps.py
    
The code in beta:
    Is to compute the overlaps between SpinTaylorT2 and SpinTaylorF2 (sidebands).

In case you need to go back:
# Resets index to former commit; replace '56e05fced' with your commit code
git reset 56e05fced 

# Moves pointer back to previous HEAD
git reset --soft HEAD@{1}

git commit -m "Revert to 56e05fced"

# Updates working copy to reflect the new commit
git reset --hard
