# MPID Protein Parameters for OpenMM

This repository contains MPID protein force field parameters mapped from the CHARMM Drude Force Field and example scripts for performing molecular dynamics simulations in OpenMM.

# Requirements

Python ≥ 3.9

OpenMM ≥ 7.5

MPIDPlugin (https://github.com/andysim/MPIDOpenMMPlugin)

# Usage
python examples/mpid_md.py 1ubq mpid 0 1

# Citation

Jing Huang, Andrew C. Simmonett, Frank C. Pickard, Alexander D. MacKerell, Bernard R. Brooks; Mapping the Drude polarizable force field onto a multipole and induced dipole model. J. Chem. Phys. 28 October 2017; 147 (16): 161702. https://doi.org/10.1063/1.4984113
