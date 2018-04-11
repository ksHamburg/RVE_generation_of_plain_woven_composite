# RVE_generation_of_plain_woven_composite
Automatic generation and discretization of fully periodic RVEs of plain woven composites

These scripts are sublementary material for the publication:
N. Jendrysik, K. Schneider, S. Bargman: "Automatic generation and discretization of fully periodic RVEs of plain woven composites" in Journal of Composite Materials, SAGE

The provided scripts generate input files for the ABAQUS software package of plain woven composite RVEs for 
six one-dimensional strain loading cases and three one-dimensional stress loading
cases for different volume fractions.

You can run the "driver_script.py" via:

# python driver_script.py

required abaqus version:   6.14
required python version:   2.x
Required python libraries: os,sys,shutil,numpy,itertools,cPickle
Optional python libraries: scipy
