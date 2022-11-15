#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import matplotlib.pyplot as plt

def main(args):
    #Ask the user if they want to fill in station/antenna info.
    #This will be used for plotting.
    station = input('What station was this data taken at? ')
    ant = input('Which antenna was this measure made at '+station+'? ')

    #Read in the data.
    f = open(args.file, 'r')
    freqs = []
    Sparams = {'S11': [], 'S21': [], 'S12': [], 'S22': []}
    mapper = {1: 'S11', 2: 'S21', 3: 'S12', 4: 'S22'}
    while True:
        line = f.readline()

        #Skip the metadata header and end when we finish the file.
        if '!' in line or '#' in line:
            continue
        elif line == '':
            break
        else:
            #Read the data and compute the S-parameters
            l = line.split()
            freqs.append(float(l[0]))

            for i in range(1,5):
                mag = 10**(float(l[2*i-1]) / 20.0)
                phase = float(l[2*i]) * (np.pi/180.0)
                S = mag * np.exp(1j * phase)
                Sparams[mapper[i]].append(S)
                
    freqs = np.array(freqs)

    #Build the plot.
    fig, axes = plt.subplots(2, 2, num=1, sharex=True)
    fig.suptitle(station+': Antenna '+ant, fontsize=16)

    axes = axes.flatten()
    for ax, (param, s) in zip(axes, Sparams.items()):
        ax.set_title(param, fontsize=14)
        ax.plot(freqs/1e6, 20.0*np.log10(np.abs(s)))
        ax.set_ylabel(param+' [dB]', fontsize=12)
        ax.tick_params(which='both', direction='in', length=6, labelsize=12)
        if ax == axes[2] or ax == axes[3]:
            ax.set_xlabel('Frequency [MHz]', fontsize=12)

    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot LWA Antenna S-parameters stored in a .s2p file from the FieldFox VNA',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type=str,
            help='.s2p file to read in')

    args = parser.parse_args()
    main(args)
