#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

DATA_PATH = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(DATA_PATH, 'Data')

def main(args):

    #Read in the files and store all the data in lists.
    freqs, S11s, S21s, S12s, S22s = [], [], [], [], []
    Sparams = {'11': [], '21': [], '12': [], '22': []}
    for file in args.files:
        f = open(file, 'r')
        
        #Read the data.
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
    
    #Compute IMF if one or more antenna files were given.
    if args.ants is not None:
        IMFs = []
        for file in args.ants:
            f = open(file, 'r')

            #Read the data.
            freq, mag_db, phase_deg = [], [], []
            while True:
                line = f.readline()
                if '!' in line or '#' in line:
                    continue
                elif line == '':
                    break
                else:
                    l = line.split()
                    freq.append(float(l[0]))
                    #Select S11 antenna data if polarization A (NS) FEE data is given.
                    #Else, choose S22 (which is antenna polarization B (EW).
                    if all('A' in os.path.basename(F) for F in args.files) or all('NS' in os.path.basename(F) for F in args.files):
                        mag_db.append(float(l[1]))
                        phase_deg.append(float(l[2]))
                    elif all('B' in os.path.basename(F) for F in args.files) or all('EW' in os.path.basename(F) for F in args.files):
                        mag_db.append(float(l[7]))
                        phase_deg.append(float(l[8]))
                    else:
                        print('Unknown polarization of FEE files. Unsure which antenna S parameter to use. Please recheck files and rerun.')
                        sys.exit()

            f.close()

            freq = np.array(freq)
            mag = 10**(np.array(mag_db)/20.0)
            phase = np.array(phase_deg)*(np.pi/180.0)

            antS11 = mag * np.exp(1j*phase)

            imfs = (1.0 - np.abs(S11s)**2) * (1.0 - np.abs(antS11)**2) / np.abs(1.0 - S11s*antS11)**2
            IMFs.append(imfs)

        IMFs = np.array(IMFs)
        
        IMF = np.mean(IMFs, axis=(0,1))
        p16, p83 = np.percentile(IMFs, q=[16, 83], axis=(0,1))

    #Convert S21 into a FEE forward voltage gain correction
    feeGain = np.abs(S21)

    #Load in the FEE Test Fixture HX62A forward gain and account for it.
    hx = np.loadtxt(os.path.join(DATA_PATH, 'HX62A', 'HX62A_S_Params.txt'))
    hxS21 = (hx[:,3] + 1j*hx[:,4]) / 2.0 #Factor of 2 since two HX62As were used to make the measurements.
    hxGain = np.abs(hxS21)

    feeGain -= hxGain

    #Plot the reflection coefficient of the FEE and the IMF, if computed.
    if not args.no_plot:
        fig, axes = plt.subplots(2, 2, num=1, sharex=True)
        axes = axes.flatten()
        fig.suptitle(f'LWA FEE S-Parameters (N={len(args.files)})', fontsize=14)
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

        try:
            fig, ax = plt.subplots(1,1,num=2)
            ax.set_title('Impedance Matching Factor', fontsize=14)
            ax.plot(freqs[0,:]/1e6, IMF, 'k-', label='Mean')
            if len(args.files) > 1:
                ax.plot(freqs[0,:]/1e6, p16, 'r--', label=r'$16^{\rm{th}}$ and $83^{\rm{rd}}$ percentiles')
                ax.plot(freqs[0,:]/1e6, p83, 'r--')
                ax.fill_between(freqs[0,:]/1e6, p16, p83, color='r', alpha=0.25)

            if args.model:
                model = np.loadtxt(os.path.join(DATA_PATH, 'BurnsZ.txt'))
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

    #Save, if requested.
    if args.save:
        try:
            header1 = f"""FEE S11 and IMF Data
Parameters Include:
* Measurement Frequencies in Hz
* Real and Imaginary Components of FEE S11
* Impedance Matching Factor Derived from FEE and Antenna S11 Measurements
* 16th and 83rd Percentiles of the IMF Distribution
Freq [Hz]              Re(S11)                   Im(S11)                   IMF                      IMF_P16                  IMF_P83
            """
            header3 = f"""FEE S21 and Forward Gain Data
Parameters Include:
* Measurement Frequencies in Hz
* Real and Imaginary Components of FEE S21
* FEE Forward Gain ( = |S21|^2 )
Freq [Hz]              Re(S21)                   Im(S21)                   Gain
            """
            if all('A' in os.path.basename(f) for f in args.files) or all('NS' in os.path.basename(f) for f in args.files):
                np.savetxt('IMF_NS.txt', np.c_[freqs[0,:].view(float), S11.view(float).reshape(S11.size, 2), IMF.view(float), p16, p83], header=header1)
                np.savetxt('FEE_Gain_NS.txt', np.c_[freqs[0,:].view(float), S21.view(float).reshape(S21.size, 2), (feeGain.view(float))**2], header=header3)
            elif all('B' in os.path.basename(f) for f in args.files) or all('EW' in os.path.basename(f) for f in args.files):
                np.savetxt('IMF_EW.txt', np.c_[freqs[0,:].view(float), S11.view(float).reshape(S11.size, 2), IMF.view(float), p16, p83], header=header1)
                np.savetxt('FEE_Gain_EW.txt', np.c_[freqs[0,:].view(float), S21.view(float).reshape(S21.size, 2), (feeGain.view(float))**2], header=header3)
            else:
                print('Unknown polarization of files. Please check inputs and rerun.')

        except NameError:
            header1 = f"""FEE S11 Data
Freq [Hz]              Re(S11)                   Im(S11)
            """
            header2 = f"""FEE S21 Data
Freq [Hz]              Re(S21)                   Im(S21)                   Gain
            """
            if all('A' in os.path.basename(f) for f in args.files) or all('NS' in os.path.basename(f) for f in args.files):
                np.savetxt('FEE_S11_NS.txt', np.c_[freqs[0,:].view(float), S11.view(float).reshape(S11.size, 2)], header=header1)
                np.savetxt('FEE_Gain_NS.txt', np.c_[freqs[0,:].view(float), S21.view(float).reshape(S21.size, 2), (feeGain.view(float))**2], header=header2)
            elif all('B' in os.path.basename(f) for f in args.files) or all('EW' in os.path.basename(f) for f in args.files):
                np.savetxt('FEE_S11_EW.txt', np.c_[freqs[0,:].view(float), S11.view(float).reshape(S11.size, 2)], header=header1)
                np.savetxt('FEE_Gain_EW.txt', np.c_[freqs[0,:].view(float), S21.view(float).reshape(S21.size, 2), (feeGain.view(float))**2], header=header2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Read in a collection of .s2p files and either compute the mean reflection coefficient or combine them with antenna data to compute the mean IMF.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('files', type=str, nargs='+',
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
