#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

DATA_PATH = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(DATA_PATH, '..')

def main(args):

    #Read in the files and store all the data in lists.
    freqs, S11s, S21s, S12s, S22s = [], [], [], [], []

    Sparams = {'11': [], '21': [], '12': [], '22': []}
    for file in args.files:
        freq = [] 
        S11_mag_db, S11_phase_deg = [], []
        S21_mag_db, S21_phase_deg = [], []
        S12_mag_db, S12_phase_deg = [], []
        S22_mag_db, S22_phase_deg = [], []

        #If it is a S2P file, it's from the first run of measurements, so only grab S12, S21, and S22.
        if '.s2p' in file:
            f = open(file, 'r')
        
            #Read the data. 
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

                    #Get the frequencies and data.
                    freq.append(float(l[0]))

                    S21_mag_db.append(float(l[3]))
                    S21_phase_deg.append(float(l[4]))
                    S12_mag_db.append(float(l[5]))
                    S12_phase_deg.append(float(l[6]))
                    S22_mag_db.append(float(l[7]))
                    S22_phase_deg.append(float(l[8]))

            freq = np.array(freq)
            freqs.append(freq)
            
            S12_mag_db = np.array(S12_mag_db)
            S12_phase_deg = np.array(S12_phase_deg)
            S21_mag_db = np.array(S21_mag_db)
            S21_phase_deg = np.array(S21_phase_deg)
            S22_mag_db = np.array(S22_mag_db)
            S22_phase_deg = np.array(S22_phase_deg)


            #First, correct S21 and S12 for the HX62A insertion loss.
            hx = np.loadtxt(os.path.join(DATA_PATH, 'HX62A', 'HX62A_S_Params.txt'))
            hx21 = (hx[:,3] + 1j*hx[:,4]) / 2.0 #Factor of 2 since two HX62As were used to make the measurements.
            hx12 = (hx[:,5] + 1j*hx[:,6]) / 2.0

            hx21_db = 20.0*np.log10(np.abs(hx21))

            S21_mag_db += hx21_db

            S12 = ( 10**(S12_mag_db / 20.0) ) * np.exp(1j*(S12_phase_deg * np.pi/180.0))
            S12s.append(S12)

            S21 = ( 10**(S21_mag_db / 20.0) ) * np.exp(1j*(S21_phase_deg * np.pi/180.0))
            S21s.append(S21)

            S22 = ( 10**(S22_mag_db / 20.0) ) * np.exp(1j*(S22_phase_deg * np.pi/180.0))
            S22s.append(S22)

        #If it is a .s1p file, it's a new measurement with the hybrid coupler de-embedded and so only S11 was measured.
        elif '.s1p' in file:
            f = open(file, 'r')
 
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

                    #Get the frequencies and S11 info.
                    freq.append(float(l[0]))
                    S11_mag_db.append(float(l[1]))
                    S11_phase_deg.append(float(l[2]))

            freq = np.array(freq)
            freqs.append(freq)

            S11_mag_db = np.array(S11_mag_db)
            S11_phase_deg = np.array(S11_phase_deg)

            S11 = ( 10**(S11_mag_db / 20.0) ) * np.exp(1j*(S11_phase_deg * np.pi/180.0))
            S11s.append(S11)

        f.close()

    #Convert to arrays and do some simple statistics if more than 1 file was given.
    freqs = np.array(freqs) 
    S11s = np.array(S11s)
    S21s = np.array(S21s)
    S12s = np.array(S12s)
    S22s = np.array(S22s)

    #First, correct S21 and S12 for the HX62A insertion loss.
    #hx = np.loadtxt(os.path.join(DATA_PATH, 'HX62A', 'HX62A_S_Params.txt'))
    #hx21 = (hx[:,3] + 1j*hx[:,4]) / 2.0 #Factor of 2 since two HX62As were used to make the measurements.
    #hx12 = (hx[:,5] + 1j*hx[:,6]) / 2.0

    #S21s = np.abs(hx21)
    #S12s = hx12

    #Mean and standard deviation.
    S11 = np.mean(S11s, axis=0)
    S11_sigma = np.std(S11s, axis=0)
    S21 = np.mean(S21s, axis=0)
    S21_sigma = np.std(S21s, axis=0)
    S12 = np.mean(S12s, axis=0)
    S12_sigma = np.std(S12s, axis=0)
    S22 = np.mean(S22s, axis=0)
    S22_sigma = np.std(S22s, axis=0)

    #Compute the upper and lower 1 sigma bounds for each.
    S11u = np.abs(S11) + S11_sigma
    S11l = np.abs(S11) - S11_sigma
    S21u = np.abs(S21) + S21_sigma
    S21l = np.abs(S21) - S21_sigma
    S12u = np.abs(S12) + S12_sigma
    S12l = np.abs(S12) - S12_sigma
    S22u = np.abs(S22) + S22_sigma
    S22l = np.abs(S22) - S22_sigma
    
    #Plot 
    if not args.no_plot:
        fig, axes = plt.subplots(2, 2, num=1, sharex=True)
        axes = axes.flatten()
        fig.suptitle(f'LWA FEE S-Parameters (N={len(args.files)//2})', fontsize=14)
        for (ax, S, Su, Sl) in zip([axes[0], axes[1], axes[2], axes[3]], [S11, S21, S12, S22], [S11u, S21u, S12u, S22u], [S11l, S21l, S12l, S22l]):
            ax.plot(freqs[0,:]/1e6, 20.0*np.log10(np.abs(S))-20.0*np.log10(np.abs(hx21)) if ax == axes[1] else 20.0*np.log10(np.abs(S)), color='k', label='Mean')
            if len(args.files) > 1:
                ax.plot(freqs[0,:]/1e6, 20.0*np.log10(Su), 'r--', label=r'$1\sigma$ Bound')
                ax.plot(freqs[0,:]/1e6, 20.0*np.log10(Sl), 'r--')
                ax.fill_between(freqs[0,:]/1e6, 20.0*np.log10(Sl), 20.0*np.log10(Su), color='r', alpha=0.25)
            ax.legend(loc=0, fontsize=12)
            ax.tick_params(which='both', direction='in', length=6, labelsize=12)
            if ax == axes[2] or ax == axes[3]:
                ax.set_xlabel('Frequency [MHz]', fontsize=12)

        for (ax, param) in zip([axes[0], axes[1], axes[2], axes[3]], ['S11', 'S21', 'S12', 'S22']):
            ax.set_title(param, fontsize=12)
            ax.set_ylabel(param+' [dB]', fontsize=12)

        try:
            fig, ax = plt.subplots(1,1,num=2)
            ax.set_title('Impedance Mismatch Factor', fontsize=14)
            ax.plot(freqs[0,:]/1e6, IMF, 'k-', label='Mean')
            if len(args.files) > 1:
                ax.plot(freqs[0,:]/1e6, p16, 'r--', label=r'$16^{\rm{th}}$ and $83^{\rm{rd}}$ percentiles')
                ax.plot(freqs[0,:]/1e6, p83, 'r--')
                ax.fill_between(freqs[0,:]/1e6, p16, p83, color='r', alpha=0.25)

            if args.model:
                model = np.loadtxt('Data/BurnsZ.txt')
                imefreqs, re, im = model[:,0], model[:,1], model[:,2]
                Z = re + 1j*im
                ime = 1.0 - np.abs( (Z-100.0)/(Z+100.0) )**2
            
                intp = interp1d(imefreqs, ime, kind='cubic', bounds_error=False)
                ime = intp(freqs[0,:]/1e6)

                ax.plot(freqs[0,:]/1e6, ime, 'g-', label='Model IME')

            ax.legend(loc=0, fontsize=12)
            ax.set_xlabel('Frequency [MHz]', fontsize=12)
            ax.set_ylabel('IMF', fontsize=12)
            ax.tick_params(which='both', direction='in', length=6, labelsize=12)
        
        except NameError:
            plt.close(fig=2)
            pass
    
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Read in a collection of .s2p files and either compute the mean reflection coefficient or combine them with antenna data to compute the mean IMF.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('files', type=str, nargs='*',
            help='.s2p files to read in')
    parser.add_argument('-a', '--ants', type=str, nargs='*', default=None,
            help='.s2p file(s) containing antenna S-parameter data.')
    parser.add_argument('-n', '--no-plot', action='store_true',
            help='Do not plot the results')
    parser.add_argument('-m', '--model', action='store_true',
            help='Plot the IME model from Hicks et al 2012 with the IMF')
    parser.add_argument('-s', '--save', action='store_true',
            help='If antenna files are given, save the IMF results as a .txt file, else save FEE S11. \
                Also saves FEE forward gain as a separate .txt file.')
    
    args = parser.parse_args()
    main(args)
