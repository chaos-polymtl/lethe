# SPDX-FileCopyrightText: Copyright (c) 2022-2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft, fftfreq

parser = argparse.ArgumentParser(
    description="Plot the drag and lift coefficients of the cylinder as a "
    "function of time and the FFT of the lift coefficient.")
parser.add_argument("folder", help="folder that contains the force.00.dat file")
parser.add_argument("--start-time", type=float, default=100.0,
                    help="time after which the vortex shedding is established "
                    "and the statistics are computed (default: 100)")
args = parser.parse_args()

data = np.loadtxt(args.folder + '/force.00.dat', skiprows=1)

time = data[:, 0]
CD = 2.0 * data[:, 1]
CL = 2.0 * data[:, 2]

# Drag and lift coefficients as a function of time
fig, (ax_CD, ax_CL) = plt.subplots(2, 1, sharex=True, figsize=(12, 6))
ax_CD.plot(time, CD, color='#1b9e77')
ax_CD.set_ylabel('$C_D$ [-]')
ax_CD.set_ylim(0.5, 2.0)  # Hide the very large drag of the impulsive start
ax_CD.grid()
ax_CL.plot(time, CL, color='#d95f02')
ax_CL.set_ylabel('$C_L$ [-]')
ax_CL.set_xlabel('Time [s]')
ax_CL.grid()
plt.savefig(args.folder + '/cylinderCDCL.png')

# Keep only the times where the vortex shedding is established
shedding = time >= args.start_time
time = time[shedding]
CD = CD[shedding]
CL = CL[shedding]

averagedCD = np.mean(CD)
averagedCL = np.mean(CL)

varCD = np.mean([abs(averagedCD - np.min(CD)), np.max(CD) - averagedCD])
varCL = np.mean([abs(averagedCL - np.min(CL)), np.max(CL) - averagedCL])

print("CD = " + str(averagedCD) + " +/- " + str(varCD))
print("CL = " + str(averagedCL) + " +/- " + str(varCL))

# The FFT needs a constant time step, which the force file does not always have
# (smaller first steps of BDF2, adaptive time stepping): interpolate CL on a
# uniform time grid
uniform_time = np.linspace(time[0], time[-1], len(time))
CL_uniform = np.interp(uniform_time, time, CL)

# One-sided spectrum of the lift coefficient
X = 2.0 * fft(CL_uniform) / len(CL_uniform)
freq = fftfreq(len(X), d=uniform_time[1] - uniform_time[0])
n_oneside = len(freq) // 2

# The shedding frequency is the peak of the spectrum, skipping the mean (f = 0)
# The Strouhal number is St = f D / U, with D = 1 and U = 1 in this example
shedding_frequency = freq[1 + np.argmax(np.abs(X[1:n_oneside]))]
print("St = " + str(shedding_frequency))

plt.figure(figsize=(12, 6))
plt.loglog(freq[:n_oneside], np.abs(X[:n_oneside]), color='#1b9e77')
plt.xlabel('Freq [Hz]')
plt.ylabel('FFT Amplitude [-]')
plt.grid(which='major')
plt.grid(which='minor', ls=':')
plt.axis([0.01, 10, 0.00001, 1])
plt.savefig(args.folder + '/cylinderFFT.png')
plt.show()
