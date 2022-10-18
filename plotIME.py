#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

model = np.loadtxt('BurnsZ.txt')
freqs, re, im = model[:,0], model[:,1], model[:,2]
Z = re + 1j*im
model_IME = 1.0 - np.abs( (Z-100.0)/(Z+100.0) )**2

intp = interp1d(freqs, model_IME, kind='cubic', bounds_error=False)

NS = np.loadtxt('IME_NS.txt')
freqs_NS, IME_NS = NS[:,0], NS[:,1]
EW = np.loadtxt('IME_EW.txt')
freqs_EW, IME_EW = EW[:,0], EW[:,1]

fig, ax = plt.subplots(1,1)
ax.set_title('Impedance Mismatch Efficiency', fontsize=16)
ax.plot(freqs_NS/1e6, IME_NS, '-', label='NS')
ax.plot(freqs_EW/1e6, IME_EW, '-', label='EW')
ax.plot(freqs_NS/1e6, intp(freqs_NS/1e6), '-', label='Model')
ax.legend(loc=0, fontsize=12)
ax.set_xlabel('Frequency [MHz]', fontsize=12)
ax.set_ylabel('IME', fontsize=12)
ax.tick_params(which='both', direction='in', length=6, labelsize=12)
plt.show()
