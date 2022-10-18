#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import argparse
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

def main(args):

    model = np.loadtxt('Data/BurnsZ.txt')
    freqs, re, im = model[:,0], model[:,1], model[:,2]
    Z = re + 1j*im
    model_IME = 1.0 - np.abs( (Z-100.0)/(Z+100.0) )**2

    intp = interp1d(freqs, model_IME, kind='cubic', bounds_error=False)

    freqs, IMEs = [], []
    for file in args.files:
        f = np.loadtxt(file)
        
        freq, IME = f[:,0], f[:,1]

        freqs.append(freq)
        IMEs.append(IME)

    fig, ax = plt.subplots(1,1)
    ax.set_title('Impedance Mismatch Efficiency', fontsize=16)
    for i in range(len(args.files)):
            ax.plot(freqs[i]/1e6, IMEs[i], '-', label='NS' if 'NS' in args.files[i] else 'EW')
    ax.plot(freqs[0]/1e6, intp(freqs[0]/1e6), '-', label='Model')
    ax.legend(loc=0, fontsize=12)
    ax.set_xlabel('Frequency [MHz]', fontsize=12)
    ax.set_ylabel('IME', fontsize=12)
    ax.tick_params(which='both', direction='in', length=6, labelsize=12)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot output from computeIME.py',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('files', type=str, nargs='*',
            help='.txt files output by computeIME.py')

    args = parser.parse_args()
    main(args)
