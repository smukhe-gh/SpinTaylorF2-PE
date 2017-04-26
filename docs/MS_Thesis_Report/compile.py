#===========================================
# Code to compile pdflatex + bibtex
# files. To run: 
# python compile.py <filename without .tex>
# i.e. python compile.py draft
#===========================================

#!/usr/bin/python

import subprocess, sys

commands = [
    ['pdflatex', sys.argv[1] + '.tex'],
    ['bibtex', sys.argv[1] + '.aux'],
    ['pdflatex', sys.argv[1] + '.tex'],
    ['pdflatex', sys.argv[1] + '.tex']
]

for c in commands:
    subprocess.call(c)