#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

DATA_PATH = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(DATA_PATH, '..')

def main(args):

    #Read in the files and store all the data in lists.
    freqs, S11s, S21s, S12s, S22s = [], [], [], [], []
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
            hx = np.loadtxt(os.path.join(DATA_PATH, 'HX62A', 'HX62A_Insertion_Loss.txt'))
            hx21_db = hx[:,1] 
            hx12_db = hx[:,2] 

            S21_mag_db += np.abs(hx21_db)
            S12_mag_db += np.abs(hx12_db)

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

    #FEE Forward Gain.
    feeGain = np.abs(S21)
    
    #Plot
    if not args.no_plot: 
        fig, axes = plt.subplots(2, 2, num=1, sharex=True)
        axes = axes.flatten()
        fig.suptitle(f'LWA FEE S-Parameters (N={len(args.files)//2})', fontsize=14)
        for (ax, S, Su, Sl) in zip([axes[0], axes[1], axes[2], axes[3]], [S11, S21, S12, S22], [S11u, S21u, S12u, S22u], [S11l, S21l, S12l, S22l]):
            ax.plot(freqs[0,:]/1e6, 20.0*np.log10(np.abs(S)), color='k', label='Mean')
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
 
        plt.show()

    header1 = """LWA FEE S-Parameters
Originally, these were measured with a HX62A coupler on the coupler PCB in the FEE Test Fixture.
However, thanks to the new calibration scheme the following corrections have been achieved:
The S11 measurements de-embed the HX62A hyrbid coupler
The S21 and S12 measurements account for the HX62A insertion loss
Frequency [Hz]               Re(S11)                   Im(S11)               Re(S21)                   Im(S21)               Re(S12)                   Im(S12)               Re(S22)                   Im(S22)
    """

    data1 = np.c_[freqs[0,:].view(float), S11.view(float).reshape(S11.size, 2), S21.view(float).reshape(S21.size, 2), S12.view(float).reshape(S12.size, 2), S22.view(float).reshape(S22.size, 2)]

    header2 = """LWA FEE Forward Voltage Gain
Frequency [Hz]               Gain
    """

    data2 = np.c_[freqs[0,:].view(float), feeGain.view(float)]

    if args.save:
        if all('A' in os.path.basename(f) for f in args.files) or all('NS' in os.path.basename(f) for f in args.files):
            np.savetxt('FEE_S_Params_NS.txt', data1, header=header1)
            np.savetxt('FEE_Gain_NS.txt', data2, header=header2)
            np.savez('FEE_S11_NS.npz', freqs=freqs, S11s=S11s)

        elif all('B' in os.path.basename(f) for f in args.files) or all('EW' in os.path.basename(f) for f in args.files):
            np.savetxt('FEE_S_Params_EW.txt', data1, header=header1)
            np.savetxt('FEE_Gain_EW.txt', data2, header=header2)
            np.savez('FEE_S11_EW.npz', freqs=freqs, S11s=S11s)
    
        else:
            print('Unknown polarization of input files. Please check inputs and rerun.')
            sys.exit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Read in a collection of FEE S-parameter measurements from various runs to plot the results and combine them all into one nice file.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('files', type=str, nargs='+',
            help='.s1p and .s2p files to read in')
    parser.add_argument('-n', '--no-plot', action='store_true',
            help='Do not plot the FEE S-parameters')
    parser.add_argument('-s', '--save', action='store_true',
            help='Save the S-parameters and the forward gain')

    args = parser.parse_args()
    main(args)
