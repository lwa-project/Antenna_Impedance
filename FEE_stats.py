#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import argparse
import numpy as np
import matplotlib.pyplot as plt

def main(args):

    #Read in the files and store all the data in lists.
    freqs, S11s = [], [] 
    for file in args.files:
        f = open(file, 'r')
        
        #Skip the metadata header.
        for i in range(5):
            f.readline()

        #Read the data.
        freq, mag_db, phase_deg = [], [], []
        while True:
            line = f.readline()
            if line == '':
                break
            l = line.split()

            freq.append(float(l[0]))
            mag_db.append(float(l[1]))
            phase_deg.append(float(l[2]))

        f.close()

        freq = np.array(freq)
        mag_db = np.array(mag_db)
        phase_deg = np.array(phase_deg)

        #Convert magnitude from dB to linear scale
        #and phase from degrees to radians.
        mag = 10**(mag_db / 20.0)
        phase = phase_deg * (np.pi/180.0)

        #Build the full complex S11 parameter.
        S11 = mag*np.exp(1j*phase)

        freqs.append(freq)
        S11s.append(S11)

    #Convert to arrays and do some simple statistics.
    freqs = np.array(freqs) 
    S11s = np.array(S11s)

    #Mean, standard deviation, and min/max values.
    S11 = np.mean(S11s, axis=0)
    sigma = np.std(S11s, axis=0)
    minimum = np.min(20.0*np.log10(np.abs(S11s)), axis=0)
    maximum = np.max(20.0*np.log10(np.abs(S11s)), axis=0)

    #Compute the upper and lower 1 sigma bounds on S11.
    S11u = np.abs(S11) + sigma
    S11l = np.abs(S11) - sigma
    
    #Compute IMF if the antenna file was given.
    if args.ant is not None:
        f = open(args.ant[0], 'r')
        
        #Skip the metadata header.
        for i in range(12):
            f.readline()

        #Read the data.
        freq, mag_db, phase_deg = [], [], []
        while True:
            line = f.readline()
            if line == '':
                break
            l = line.split()

            freq.append(float(l[0]))
            mag_db.append(float(l[1]))
            phase_deg.append(float(l[2]))

        f.close()

        freq = np.array(freq)
        mag_db = np.array(mag_db)
        phase_deg = np.array(phase_deg)

        #Convert magnitude from dB to linear scale
        #and phase from degrees to radians.
        antmag = 10**(mag_db / 20.0)
        antphase = phase_deg * (np.pi/180.0)

        #Build the full complex S11 parameter for the antenna.
        antS11 = antmag*np.exp(1j*antphase)

        IMFs = (1.0 - np.abs(S11s)**2) * (1.0 - np.abs(antS11)**2) / np.abs(1.0 - S11s*antS11)**2
        
        IMF = np.mean(IMFs, axis=0)
        p16, p83 = np.percentile(IMFs, q=[16, 83], axis=0)

    #Plot the reflection coefficient of the FEE and the IMF, if computed.
    if not args.no_plot:
        fig, ax = plt.subplots(1,1,num=1)
        fig.suptitle(f'Mean Reflection Coefficient (N={len(args.files)})', fontsize=14)
        ax.plot(freqs[0,:]/1e6, 20.0*np.log10(np.abs(S11)), color='k', label='Average Reflection Coefficient')
        ax.plot(freqs[0,:]/1e6, 20*np.log10(S11u), 'r--', label=r'$1\sigma$ Bound')
        ax.plot(freqs[0,:]/1e6, 20*np.log10(S11l), 'r--')
        ax.fill_between(freqs[0,:]/1e6, 20*np.log10(np.abs(S11l)), 20*np.log10(np.abs(S11u)), color='r', alpha=0.25)
        #ax.plot(freqs[0,:]/1e6, 20*np.log10(minimum), 'k--')
        #ax.plot(freqs[0,:]/1e6, 20*np.log10(maximum), 'k--')
        ax.legend(loc=0, fontsize=12)
        ax.set_xlabel('Frequency [MHz]', fontsize=12)
        ax.set_ylabel(r'$\Gamma$ [dB]', fontsize=12)
        ax.tick_params(which='both', direction='in', length=6, labelsize=12)

        try:
            fig, ax = plt.subplots(1,1,num=2)
            fig.suptitle('Impedance Matching Factor', fontsize=14)
            ax.plot(freqs[0,:]/1e6, IMF, 'k-', label='Mean')
            ax.plot(freqs[0,:]/1e6, p16, 'r--', label=r'$16^{\rm{th}}$ and $83^{\rm{rd}}$ percentiles')
            ax.plot(freqs[0,:]/1e6, p83, 'r--')
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
            header = f"""FEE S11 Data
Freq [Hz]              S11              IMF
            """
            np.savetxt('FEE.txt', np.c_[freqs, S11, IMF], header=header)
            header = f"""FEE IMF Percentiles
Freq [Hz]              P16              P83
            """
            np.savetext('IMF_percentiles.txt', np.c_[freqs, p16, p83], header=header)

        except NameError:
            header = f"""FEE S11 Data
Freq [Hz]              S11
            """ 
            np.savetxt('FEE.txt', np.c_[freqs, S11], header=header)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Read in a collection of .s1p files and compute the mean reflection coefficient',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('files', type=str, nargs='*',
            help='.s1p files to read in')
    parser.add_argument('-a', '--ant', type=str, nargs=1, default=None,
            help='.s1p file containing antenna reflection coefficient data. \
            If this is not None, the output will be IMF, not IME')
    parser.add_argument('-n', '--no-plot', action='store_true',
            help='Do not plot the results')
    parser.add_argument('-s', '--save', action='store_true',
            help='Save the results as a .txt file')
    
    args = parser.parse_args()
    main(args)
