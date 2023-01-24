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

    #Read in the FEE S-parameters.
    fee = np.load(args.file[0])
    freqs, S11s = fee['freqs'], fee['S11s']

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
                if ('A' in os.path.basename(args.file[0])) or ('NS' in os.path.basename(args.file[0])):
                    mag_db.append(float(l[1]))
                    phase_deg.append(float(l[2]))
                elif ('B' in os.path.basename(args.file[0])) or ('EW' in os.path.basename(args.file[0])):
                    mag_db.append(float(l[7]))
                    phase_deg.append(float(l[8]))
                else:
                    print('Unknown polarization of FEE file. Unsure which antenna S parameter to use. Please recheck files and rerun.')
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

    if not args.no_plot:
        fig, ax = plt.subplots(1,1,num=1)
        ax.set_title('Impedance Mismatch Factor', fontsize=14)
        ax.plot(freqs[0,:]/1e6, IMF, 'k-', label='Mean')
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
        
        plt.show()

    #Save, if requested.
    if args.save:
        header = f"""FEE S11 and IMF Data
Parameters Include:
* Measurement Frequencies in Hz
* Real and Imaginary Components of FEE S11
* Impedance Matching Factor Derived from FEE and Antenna S11 Measurements
* 16th and 83rd Percentiles of the IMF Distribution
Freq [Hz]              IMF                      IMF_P16                  IMF_P83
            """
        if ('A' in os.path.basename(args.file[0])) or ('NS' in os.path.basename(args.file[0])):
            np.savetxt('IMF_NS.txt', np.c_[freqs[0,:].view(float), IMF.view(float), p16, p83], header=header)
        elif ('B' in os.path.basename(args.file[0])) or ('EW' in os.path.basename(args.file[0])):
            np.savetxt('IMF_EW.txt', np.c_[freqs[0,:].view(float), IMF.view(float), p16, p83], header=header)
        else:
            print('Unknown polarization of files. Please check inputs and rerun.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Read in a collection of .s2p files and either compute the mean reflection coefficient or combine them with antenna data to compute the mean IMF.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type=str, nargs=1,
            help='.npz file containing FEE S11 data')
    parser.add_argument('ants', type=str, nargs='+',
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
