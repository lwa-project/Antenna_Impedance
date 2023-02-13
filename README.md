Antenna_Impedance
===================
[![Paper]()]()

Description
-----------
Collection of scripts and files dealing with the LWA antenna impedance.

In the main directory are the general scripts:
* compute_IME.py: Reads in a single .s1p file containing antenna reflection coefficient data and computes the IME 
* plot_IME.py: Plots the output from computeIME.py along with the model IME from antenna simulations
* compute_IMF.py: Reads in a .npz file output by read_FEE_S_params.py and a collection of antenna .s2p files to compute the Impedance Mismatch Factor (IMF)
* plot_S_params.py: Takes in a collection of .s2p files and plots the full set of S-parameters for each antenna. Used to make the dipole-dipole measurement plots.
* plot_cross_coupling.py: Takes in a collection of .s2p files from an antenna-antenna cross coupling measurement and plots S12 and S21 for each of the four possible polarization permutations. 

Data Files
----------
In the Data/ directory there is a collection of various data files with their associated documentation.
Files include:
* BurnsZ.txt: Model IME from old antenna simulations
* A-NS_AntennaImpedance_062222.s1p: Antenna Impedance for A (N/S) polarization taken on 06/22/2022 at Whitham Reeve's house in Alaska
* B-EW_AntennaImpedance_062222.s1p: Antenna Impedance for B (E/W) polarization taken on 06/22/2022 at Whitham Reeve's house in Alaska
* A collection of dipole-dipole and antenna-antenna .s2p measurement files taken at the New Mexico stations (See "LWA_NM_Measurements_and_Procedures.pdf" for more information)

HX62A
-----
The Data/HX62A/ drectory contains:
* .s2p files containing S-parameter measurements for the HX62A present on the FEE test Fixture coupler PCB.
* average_HX62A_measurements.py: Reads in both .s2p files and returns the average S21 and S12 amplitudes in dB

FEE
---
The Data/FEE/ directory contains:
* A collection of .s1p files which contain S11 measurements that de-embed the HX62A.
* A collection of .s2p files which contains S-parameter measurements without the HX62A de-embedded.
* read_FEE_S_params.py: A script to read both the .s1p and .s2p files and write out proper S-pramaters (S11 with HX62A de-embdedd, correct S21 and S12, and S22).
* plot_FEE_Gain.py: A script to plot the FEE Gain files which are output by read_FEE_S_params.py

Example
-------
The following is a short example which shows how to generate the IMF correction for LWA data, which requires both LWA antenna S11 measurements and FEE S11 measurements.
The FEE S11 data with the HX62A de-embedded are located in `Data/FEE/S1P/`. The single antenna measurements for each station are located in the named directories in the `Data/` directory
under `Dipole-Dipole/` (e.g. the LWA-SV measurements are in `Data/LWA-SV/Dipole-Dipole/`.

First, the FEE S-parameters with the HX62A de-embedded have to be collected into a single .npz file. To do this, use `Data/FEE/read_FEE_S_params.py` for both polarizations:
```
python3 read_FEE_S_params.py -ns S1P/*A* S2p/*A*
```
and
```
python3 read_FEE_S_params.py -ns S1P/*B* S2p/*B*
```
This creates two files `FEE_S11_NS.npz` and `FEE_S11_EW.npz`, respectively. Finally, to compute the IMF correction using LWA1 antenna measurements and FEE S11 measurements use `compute_IMF.py`
for each polarization via
```
python3 compute_IMF.py -ns Data/FEE/FEE_S11_NS.npz Data/LWA1/Dipole-Dipole/*
```
and
```
python3 compute_IMF.py -ns Data/FEE/FEE_S11_EW.npz Data/LWA1/Dipole-Dipole/*
```
This creates two files named `IMF_NS.txt` and `IMF_EW.txt`, respectively, which contain the IMF as a function of frequency along with the $$16^{th}$$ and $$84^{th}$$ percentile uncertainty bounds.

Relevant Documentation and Papers
---------------------------------
[Hicks et al. (2012, PASP, 124, 1090)](https://iopscience.iop.org/article/10.1086/668121/pdf)
[Engineering Documentation for the LWA FEE (LWA Memo #190)](https://leo.phys.unm.edu/~lwa/memos/memo/lwa0190.pdf)
