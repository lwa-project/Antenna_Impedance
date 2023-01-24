#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#Simple script to averge the magnitude of S21 for the In-phase and Out-of-phase HX62A S-parameter measurements.

import numpy as np

#Read in the two files and get the data.
files = ['HX62Ax2_HybridCoupler_In-PhaseMeasurements_19Jan2023.s2p', 'HX62Ax2_HybridCoupler_Out-PhaseMeasurements_19Jan2023.s2p']

freqs, S21s = [], []
for file in files:
    f = open(file, 'r')

    freq = [] 
    S21_mag_db, S21_phase_deg = [], []
    S12_mag_db, S12_phase_deg = [], []
    
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
            S21_phase_deg.append(float(l[4]))
            S12_mag_db.append(float(l[5]))
            S12_phase_deg.append(float(l[6]))

    f.close()

    freq = np.array(freq)
    freqs.append(freq)
   
    S12_mag_db = np.array(S12_mag_db)
    S12_phase_deg = np.array(S12_phase_deg)
    S21_mag_db = np.array(S21_mag_db)
    S21_phase_deg = np.array(S21_phase_deg)

    S21_mag = 10**(S21_mag_db / 20.0)
    S21s.append(S21_mag)

freqs = np.array(freqs) 
S21s = np.array(S21s)
S21 = np.mean(S21s, axis=0)
S21 = 20.0*np.log10(S21)

#Write the average to file which can be read in as a "correction" to the FEE forward gain calculation in compute_IMF.py.
header = """HX62A Scattering Parameters
Freq [Hz]             |S21| [dB]
"""
np.savetxt('HX62A_S21.txt', np.c_[freqs[0,:].view(float), S21.view(float)], header=header)
