# RVE_generation_of_plain_woven_composite
Automatic generation and discretization of fully periodic RVEs of plain woven composites

These scripts are sublementary material for the publication:
N. Jendrysik, K. Schneider, S. Bargman: "Automatic generation and discretization of fully periodic RVEs of plain woven composites" in Journal of Composite Materials, SAGE

The provided scripts generate input files for the ABAQUS software package of plain woven composite RVEs for
six one-dimensional strain loading cases and three one-dimensional stress loading
cases for different volume fractions.

You can run the "driver_script.py" via:

```python driver_script.py```

## Requirements:
- required abaqus version:   6.14
- required python version:   2.x
- Required python libraries: os,sys,shutil,numpy,itertools,cPickle
- Optional python libraries: scipy

# Abaqus 2024
Since Abaqus 2024 the python interpreter of Abaqus is upgraded to version 3. Therefore the old scripts cannot be executed without any errors.
Therefore you can find the updated scripts in the folder ```Abaqus_2024```.
You can run the "driver_script.py" via:

```python3 driver_script.py```


## Requirements:
- required abaqus version:   2024
- required python version:   3.x
- Required python libraries: os,sys,shutil,numpy,itertools
- Optional python libraries: scipy,cPickle
