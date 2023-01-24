#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#Simple script to averge the magnitude of S21 for the In-phase and Out-of-phase HX62A S-parameter measurements.

import numpy as np

#Read in the two files and get the data.
files = ['HX62Ax2_HybridCoupler_In-PhaseMeasurements_19Jan2023.s2p', 'HX62Ax2_HybridCoupler_Out-PhaseMeasurements_19Jan2023.s2p']

freqs, S21s, S12s = [], [], []
for file in files:
    f = open(file, 'r')

    freq = [] 
    S21_mag_db = []
    S12_mag_db = []
    
    while True:
        line = f.readline()
        
        #Skip the metadata header.
        if '!' in line or '#' in line:
            continue
        #Look for when we hit the end of the file.
        elif line == '':
            break
        else:
            l = line.split()

            freq.append(float(l[0]))
            S21_mag_db.append(float(l[3]))
            S12_mag_db.append(float(l[5]))

    f.close()

    #Convert from dB to linear so we can properly average
    freq = np.array(freq)   
    S12_mag_db = np.array(S12_mag_db)
    S21_mag_db = np.array(S21_mag_db)

    S21_mag = 10**(S21_mag_db / 20.0)
    S12_mag = 10**(S12_mag_db / 20.0)

    freqs.append(freq)
    S21s.append(S21_mag)
    S12s.append(S12_mag)

freqs = np.array(freqs) 

#Average and then convert back to dB.
#The extra factor of 1/2 is from the fact that these parameters were measured using two HX62As.
S21s = np.array(S21s)
S21 = np.mean(S21s, axis=0)
S21 = 20.0*np.log10(S21) /2.0

S12s = np.array(S12s)
S12 = np.mean(S12s, axis=0)
S12 = 20.0*np.log10(S12) / 2.0

#Write the average to file which can be read in as a "correction" to the FEE forward gain calculation in compute_IMF.py.
header = """HX62A Scattering Parameters
Freq [Hz]             |S21| [dB]             |S12| [dB]
"""
np.savetxt('HX62A_Insertion_Loss.txt', np.c_[freqs[0,:].view(float), S21.view(float), S12.view(float)], header=header)
