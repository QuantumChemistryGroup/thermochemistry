# thermochemistry
This repository **thermochemistry** is a collection of the Python 3 scripts to perform the gas-phase thermochemistry calculation based on molecular data within the following approximations:
1) translational partition function: *ideal gas*
2) rotational partition function: *rigid rotor* or *hindered rotor*
3) vibrational parition function: *hramonic oscillator* or *msRRHO approach of S. Grimme (DOI: 10.1002/chem.201200497)*
## Prerequisites
1) Linux operational system
2) Python 3
3) Numpy (if not installed, please install it with the following command: ```pip3 install numpy```)
## Running the code
### Generation of input files (*dat)
The easiest way to generate the *dat file with molecular parameters to perform the subsequent thermochemistry calculations id to run our script on the ORCA frequency job[^1]:
1) Frequencies scaling factor (1.011) and Harmonic Oscillator: ```python3 /path/to/script/td_input_orca.py m1_00_BC-PBE0-PBE.out 1.011 HO```
2) Frequencies scaling factor (1.011) and msRRHO (tau = 100, alpha = 4, enthapy is from HO) : ```python3 /path/to/script/td_input_orca.py m1_00_BC-PBE0-PBE.out 1.011 GR_100_4_0```
3) Frequencies scaling factor (1.011) and msRRHO (tau = 100, alpha = 4, enthapy is from msRRHO) : ```python3 /path/to/script/td_input_orca.py m1_00_BC-PBE0-PBE.out 1.011 GR_100_4_1```
4) For hindered rotor (*one* rotatable bond), see file ```examples/m1_01_NaH2O6_411-PBE0.dat``` (potential is given in kcal/mol), the corresponding output: ```examples/m1_01_NaH2O6_411-PBE0.dat```
### Thermochemistry calculation
Run the following command on the generated *dat file: ```python3 /path/to/script/thermochemistry_mmRRHO.py m1_00_BC-PBE0-PBE.dat```
The obtained *td file contains all necessary information on thermochemistry
[^1]: Please, see "examples" folder for m1_00_BC-PBE0-PBE.out (ORCA output), m1_00_BC-PBE0-PBE.dat (Harmonic Oscillator), m1_00_BC-PBE0-PBE.td (HO thermochemistry)
> [!IMPORTANT]
> **When using this code please cite the following publication:**
> 1) "Gas‐phase thermochemistry of noncovalent ligand–alkali metal ion clusters: An impact of low frequencies" A. A. Otlyotov, Y. Minenkov J. Comput. Chem. 2023, 44, 1807 – 1816 (DOI: 10.1002/jcc.27129)
> 2) "An Influence of Electronic Structure Theory Method, Thermodynamic and Implicit Solvation Corrections on the Organic Carbonates Conformational and Binding Energies"
A. S. Ryzhako, A. A. Tuma, A. A. Otlyotov, and Y. Minenkov, 2024, submitted.
> 
> **If using the hindered rotor approach please cite the following work:**
>
> 3) "Numerical Evaluation of Energy Levels and Wave Functions for Hindered Internal Rotation" G. Ercolani J. Chem. Educ. 2000, 77, 11, 1495 (DOI: 10.1021/ed077p1495)
>
> **If using the msRRHO approximation please cite the following publications:**
>
> 4) "Supramolecular Binding Thermodynamics by Dispersion-Corrected Density Functional Theory" S. Grimme, Chem. Eur. J 2012, 18, 9955 (DOI: 10.1002/chem.201200497)
> 5) "Calculation of absolute molecular entropies and heat capacities made simple" P. Pracht, S. Grimme, Chem. Sci., 2021,12, 6551-6568 (DOI: 10.1039/D1SC00621E)
> 
