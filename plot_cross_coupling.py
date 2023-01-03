#!/usr/bin/env python3
# -*- conding: utf-8 -*-

import argparse
import numpy as np
import matplotlib.pyplot as plt


def main(args):
    #Ask the user which station this was done at.
    station = input('Which station were these data taken at? ')

    configs = {'ant1': [], 'ant2': [], 'pol1': [], 'pol2': []}
    S21s, S12s = [], []
    for file in args.files:
        print('Reading file: '+file)
        
        #Get some info from the user about the setup.
        ant1, ant2 = input('Which antennas? Please enter VNA Port 1 antenna first. [Ant1 Ant2] ').split(' ')
        pol1, pol2 = input('Which polarization combination? [Ant1_pol Ant2_pol] ').split(' ')

        #Add to the configurations dictionary.
        configs['ant1'].append(ant1)
        configs['ant2'].append(ant2)
        configs['pol1'].append(pol1)
        configs['pol2'].append(pol2)

        #Read the data.
        f = open(file, 'r')
        freqs, S12, S21 = [], [], []
        while True:
            line = f.readline()
            #Skip the metadata header and end when we hit the end of the file.
            if '!' in line or '#' in line:
                continue
            elif line == '':
                break
            else:
                l = line.split()
                freqs.append(float(l[0]))

                S21_mag = 10.0**(float(l[3]) / 20.0)
                S21_phase = float(l[4]) * (np.pi/180.0)
                S12_mag = 10.0**(float(l[5]) / 20.0)
                S12_phase = float(l[6]) * (np.pi/180.0)

                s21 = S21_mag * np.exp(1j * S21_phase)
                s12 = S12_mag * np.exp(1j * S12_phase)
                
                S21.append(s21)
                S12.append(s12)

        freqs = np.array(freqs)
        S21s.append(S21)
        S12s.append(S12)

    #Make the plots.
    fig, axes = plt.subplots(2, 2, num=1, sharex=True)
    fig.suptitle(station+' S21 Antenna Coupling Measurements', fontsize=14)
    axes = axes.flatten()
    for (ax, s21, a1, a2, p1, p2) in zip(axes, S21s, configs['ant1'], configs['ant2'], configs['pol1'], configs['pol2']):
        ax.plot(freqs/1e6, 20.0*np.log10(np.abs(s21)), '-')
        ax.set_ylim(-80.0, -25.0)
        ax.set_title('Antenna '+a1+' '+p1+' to Antenna '+a2+' '+p2, fontsize=12)
        ax.set_ylabel('S21 [dB]', fontsize=12)
        ax.tick_params(which='both', direction='in', length=6, labelsize=12)
        if ax == axes[2] or ax == axes[3]:
            ax.set_xlabel('Frequency [MHz]', fontsize=12)

    fig2, axes2 = plt.subplots(2, 2, num=2, sharex=True)
    fig2.suptitle(station+' S12 Antenna Coupling Measurements', fontsize=14)
    axes2 = axes2.flatten()
    for (ax, s12, a1, a2, p1, p2) in zip(axes2, S12s, configs['ant1'], configs['ant2'], configs['pol1'], configs['pol2']):
        ax.plot(freqs/1e6, 20.0*np.log10(np.abs(s12)), '-')
        ax.set_ylim(-80.0, -25.0)
        ax.set_title('Antenna '+a1+' '+p1+' to Antenna '+a2+' '+p2, fontsize=12)
        ax.set_ylabel('S12 [dB]', fontsize=12)
        ax.tick_params(which='both', direction='in', length=6, labelsize=12)
        if ax == axes2[2] or ax == axes2[3]:
            ax.set_xlabel('Frequency [MHz]', fontsize=12)

    plt.show()

if __name__ == '__main__':
    parser=argparse.ArgumentParser('Read in a collection of cross-antenna coupling measurements and plot them',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('files', type=str, nargs='*',
            help='.s2p files containing the cross-coupling measurements')

    args = parser.parse_args()
    main(args)
