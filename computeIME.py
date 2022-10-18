#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import argparse
import numpy as np
import matplotlib.pyplot as plt

def main(args):

    #Read in the files and skip the metadata headers.
    f = open(args.file, 'r')
    for i in range(12):
        f.readline()

    freqs, mag_db, phase_deg = [], [], []
    while True:
        line = f.readline()
        if line == '':
            break
        l = line.split()

        freqs.append(float(l[0]))
        mag_db.append(float(l[1]))
        phase_deg.append(float(l[2]))

    f.close()

    freqs = np.array(freqs)
    mag_db = np.array(mag_db)
    phase_deg = np.array(phase_deg)

    #Convert magnitude from dB to linear scale
    #and phase from degrees to radians.
    mag = 10**(mag_db / 20.0)
    phase = phase_deg * (np.pi/180.0)

    #Build the full complex reflection parameter.
    S11 = mag*np.exp(1j*phase)

    #Calculate IME. 
    IME = 1.0 - np.abs(S11)**2

    #Plot.
    if not args.no_plot:
        fig, ax = plt.subplots(1,1)
        ax.set_title('Impedance Mismatch Efficiency', fontsize=14)
        ax.plot(freqs/1e6, IME, '-')
        ax.set_xlabel('Frequency [MHz]', fontsize=12)
        ax.set_ylabel('Efficiency', fontsize=12)
        ax.tick_params(which='both', direction='in', length=6, labelsize=12)

        plt.show()

    #Save the results to a .txt file, if requested.
    if args.save:
        header = f"""IME computed from file: {args.file}
Freq [Hz]              IME
        """
        np.savetxt('IME.txt', np.c_[freqs, IME], header=header)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Read in .s1p from Keysight Technologies VNA to compute IME',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('file', type=str,
            help='.s1p files to read in')
    parser.add_argument('-n', '--no-plot', action='store_true',
            help='Do not plot the IME')
    parser.add_argument('-s','--save', action='store_true',
            help='Save the calculated IME as a text file')

    args = parser.parse_args()
    main(args)
