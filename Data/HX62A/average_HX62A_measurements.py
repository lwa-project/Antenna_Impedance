#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#Simple script to averge the In-phase and Out-of-phase HX62A S-parameter measurements.

import numpy as np

#Read in the two files and get the data.
files = ['HX62Ax2_HybridCoupler_In-PhaseMeasurements_19Jan2023.s2p', 'HX62Ax2_HybridCoupler_Out-PhaseMeasurements_19Jan2023.s2p']

freqs, S11s, S21s, S12s, S22s = [], [], [], [], []
Sparams = {'11': [], '21': [], '12': [], '22': []}

for file in files:
    f = open(file, 'r')

    freq = [] 
    S11_mag_db, S11_phase_deg = [], []
    S21_mag_db, S21_phase_deg = [], []
    S12_mag_db, S12_phase_deg = [], []
    S22_mag_db, S22_phase_deg = [], []
    
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
            S11_mag_db.append(float(l[1]))
            S11_phase_deg.append(float(l[2]))
            S21_mag_db.append(float(l[3]))
            S21_phase_deg.append(float(l[4]))
            S12_mag_db.append(float(l[5]))
            S12_phase_deg.append(float(l[6]))
            S22_mag_db.append(float(l[7]))
            S22_phase_deg.append(float(l[8]))

    f.close()

    freq = np.array(freq)
    freqs.append(freq)
    
    S11_mag_db = np.array(S11_mag_db)
    S11_phase_deg = np.array(S11_phase_deg)
    S12_mag_db = np.array(S12_mag_db)
    S12_phase_deg = np.array(S12_phase_deg)
    S21_mag_db = np.array(S21_mag_db)
    S21_phase_deg = np.array(S21_phase_deg)
    S22_mag_db = np.array(S22_mag_db)
    S22_phase_deg = np.array(S22_phase_deg)

    mapper = {0: '11', 1: '21', 2: '12', 3: '22'}
    for i, (mag_db, phase_deg) in enumerate(zip([S11_mag_db,S21_mag_db,S12_mag_db,S22_mag_db], [S11_phase_deg,S21_phase_deg,S12_phase_deg,S22_phase_deg])):
        #Convert magnitude from dB to linear scale
        #and phase from degrees to radians.
        mag = 10**(mag_db / 20.0)
        phase = phase_deg * (np.pi/180.0)

        #Build the full complex S parameter and store it in the dictionary.
        S = mag*np.exp(1j*phase)
        Sparams[mapper[i]].append(S)

#Convert to arrays and do some simple statistics if more than 1 file was given.
freqs = np.array(freqs) 
S11s = np.array(Sparams['11'])
S21s = np.array(Sparams['21'])
S12s = np.array(Sparams['12'])
S22s = np.array(Sparams['22'])

S11 = np.mean(S11s, axis=0)
S21 = np.mean(S21s, axis=0)
S12 = np.mean(S12s, axis=0)
S22 = np.mean(S22s, axis=0)

#Write the average to file which can be read in as a "correction" to the FEE forward gain calculation in compute_IMF.py.
header = """HX62A Scattering Parameters
Freq [Hz]             Re(S11)                    Im(S11)                   Re(S21)                   Im(S21)                   Re(S12)                   Im(S12)                   Re(S22)                   Im(S22)
"""
np.savetxt('HX62A_S_Params.txt', np.c_[freqs[0,:].view(float), S11.view(float).reshape(S11.size, 2),S21.view(float).reshape(S21.size, 2),S12.view(float).reshape(S12.size, 2),S22.view(float).reshape(S22.size, 2)], header=header)
