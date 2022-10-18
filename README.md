# Antenna_Impedance
Collection of scripts and files dealing with the LWA antenna impedance

Scripts
-------
There are three main scripts in the repository:
* computeIME.py: Reads in a single .s1p file containing antenna reflection coefficient data and computes the IME 
* plotIME.py: Plots the output from computeIME.py along with the model IME from antenna simulations
* FEE_stats.py: Takes in a collection of .s2p files containing FEE S-parameter data and computes some simple statistics on S11 and IMF
 
 Data Files
 ----------
 There is an assortment of data stored here for ease of access:
 * BurnsZ.txt: Model IME
 * A-NS_AntennaImpedance_062222.s1p: Antenna Impedance for A (N/S) polarization taken on 06/22/2022
 * B-EW_AntennaImpedance_062222.s1p: Antenna Impedance for B (E/W) polarization taken on 06/22/2022
 * A collection of .s2p files stored in Data/S2P/ (See "FEE_Bulk_Measurements.txt" for more information)
