Antenna_Impedance
===================
[![Paper]()]()

Description
-----------
Collection of scripts and files dealing with the LWA antenna impedance.

There are five main scripts in the repository:
* computeIME.py: Reads in a single .s1p file containing antenna reflection coefficient data and computes the IME 
* plotIME.py: Plots the output from computeIME.py along with the model IME from antenna simulations
* FEE_stats.py: Takes in a collection of .s2p files containing FEE S-parameter data and computes some simple statistics on S11. If an antenna .s2p file is also given, it will compute the IMF along with 1-sigma bounds.
* plot_S_params.py: Takes in a collection of .s2p files and plots the full set of S-parameters for each antenna. Used to make the dipole-dipole measurement plots.
* plot_cross_coupling.py: Takes in a collection of .s2p files from an antenna-antenna cross coupling measurement and plots S12 and S21 for each of the four possible polarization permutations. 
 
 Data Files
 ----------
 There is an assortment of data stored here for ease of access:
 * BurnsZ.txt: Model IME from old antenna simulations
 * A-NS_AntennaImpedance_062222.s1p: Antenna Impedance for A (N/S) polarization taken on 06/22/2022 at Whitham Reeve's house in Alaska
 * B-EW_AntennaImpedance_062222.s1p: Antenna Impedance for B (E/W) polarization taken on 06/22/2022 at Whitham Reeve's house in Alaska
 * A collection of .s2p files stored in Data/FEE_Bulk_Measurements/S2P/ (See "FEE_Bulk_Measurements.txt" for more information)
 * A collection of dipole-dipole and antenna-antenna .s2p measurement files taken at the New Mexico stations (See "LWA_NM_Measurements_and_Procedures.pdf" for more information)
 * Two measurement files of the S-parameters of the HX62A coupler which are needed to properly compute the FEE forward gain since the first FEE test fixture did not de-embed this. (See "HX62A_Insertion_Loss_Measurements.pdf")
