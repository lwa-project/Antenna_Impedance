#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import argparse
import numpy as np
import matplotlib.pyplot as plt

def main(args):

    #Read in the data.
    freqs, gains = [], []
    for file in args.files:
        f = np.loadtxt(file)
        freq, gain = f[:,0], f[:,-1]

        freqs.append(freq)
        gains.append(gain)

    freqs = np.array(freqs)
    gains = np.array(gains)

    if args.dB:
        gains = 20.0*np.log10(gains)
    else:
        gains = (gains.T / np.max(gains, axis=1)).T 

    mapper = {0: 'X pol', 1: 'Y pol'}

    #Plot.
    fig, ax = plt.subplots(1,1)
    ax.set_title('FEE Forward Gain', fontsize=14)
    for i in range(2):
        ax.plot(freqs[i]/1e6, gains[i], label=mapper[i])
    ax.set_xlabel('Frequency [MHz]', fontsize=12)
    ax.set_ylabel('Forward Gain' + ' [dB]' if args.dB else ' [lin.]', fontsize=12)
    ax.legend(loc=0, fontsize=12)
    ax.tick_params(which='both', direction='in', length=6, labelsize=12)

    plt.show()

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description='Plot the FEE forward gain output by plot_IMF.py',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('files', type=str, nargs=2,
        help='.txt file containing the FEE forward gain data')
    parser.add_argument('-d', '--dB', action='store_true',
        help='Plot the gain on a dB scale')

    args = parser.parse_args()
    main(args)

