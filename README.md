# ITCF
This code takes in a path integral molecular dynamics simulation of a water monomer and calculates imaginary time correlation functions. The results can then be used to calculate the vibrational frequencies of the 3 normal modes (bend, symmetric stretch, and asymmetric stretch).

The imaginary time correlation functions calculated are:
- r(H1-O) + r(H2-O) : Used to calculate symmetric stretch mode
- r(H1-O) - r(H2-O) : Used to calculate asymmetric stretch mode
- theta (H1-O-H2) : Used to calculate the bending mode

This has been demonstrated in:
"Path integral molecular dynamic simulation of flexible molecular systems in their ground state: Application to the water dimer"
J. Chem. Phys. 148, 124116 (2018); https://doi.org/10.1063/1.5017532
Matthew Schmidt and Pierre-Nicholas Roy
