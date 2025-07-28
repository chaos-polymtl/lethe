# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from cycler import cycler
from scipy.fft import rfft, rfftfreq
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d


# Set plot parameters
colors=['#1B9E77','#D95F02','#7570B3']
plt.rcParams['axes.prop_cycle'] = cycler(color = colors)
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize'] = (10,8)
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['lines.markersize'] = '11'
plt.rcParams['markers.fillstyle'] = "none"
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.columnspacing'] = 2
plt.rcParams['legend.handlelength'] = 3
plt.rcParams['legend.handletextpad'] = 0.2
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fancybox'] = False
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['font.size'] = '25'
plt.rcParams['font.family']='DejaVu Serif'
plt.rcParams['font.serif']='cm'
plt.rcParams['savefig.bbox']='tight'
plt.rcParams['legend.handlelength']=1

######################################################################

import sys
sys.path.append("$LETHE_PATH/contrib/postprocessing/")
from lethe_pyvista_tools import *

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the Pseudo-2d gas solid fluidized bed example')
parser.add_argument("-f", "--folder", type=str, help="Folder path. This folder is the folder which contains the .prm file.", required=True)
args, leftovers=parser.parse_known_args()

# Simulation folder
folder=args.folder

# Starting-ending vtu id
start = 200
end = 1200

# Load lethe data
pvd_particles = 'out_particles.pvd'
pvd_fluid     = 'out.pvd'
prm_file = 'gas-solid-fluidized-bed.prm'
particles = lethe_pyvista_tools(folder, prm_file, pvd_particles)
fluid = lethe_pyvista_tools(folder, prm_file, pvd_fluid)
time = np.array(particles.time_list)

# Get properties on the plane at 45mm above the floating wall
pressure = np.zeros(end-start)
void_fraction = np.zeros(end-start)
bed_height = np.zeros(end-start)
y = 0.004

for i in range(start, end):

    print(f"Time: {time[i]:.2f} s")
    df_fluid = fluid.get_df(i)
    df_particles = particles.get_df(i)
    df_loc = pd.DataFrame(np.copy(df_particles.points), columns=['x', 'y','z'])

    bed_height[i-start] = df_loc['z'].max()
    sampled_data = df_fluid.sample_over_line([0.042, y, 0.045], [0.048, y, 0.045])
    void_fraction[i - start] = np.mean(sampled_data['void_fraction'])
    pressure[i - start] = np.mean(sampled_data['pressure'])


# Plot mean pressure over time
reference_pressure = pd.read_csv('reference/relative_pressure.csv')
plt.figure()
plt.plot(reference_pressure['t'],reference_pressure['p'], label='Reference')
plt.plot(time[start:end],pressure-np.mean(pressure),label='Lethe')
plt.legend()
plt.xlabel('Time [s]')
plt.ylabel('Relative pressure [Pa]')
plt.xlim(3,6)
plt.ylim(-300,300)
plt.grid()
plt.subplots_adjust(left=0.2)
plt.savefig('pressure-fluctuations')
plt.show()

# Plot void fraction over time
reference_voidage = pd.read_csv('reference/voidage.csv')
plt.figure()
plt.plot(reference_voidage['t'],reference_voidage['voidage'], label='Reference')
plt.plot(time[start:end],void_fraction, label='Lethe')
plt.legend()
plt.xlabel('Time [s]')
plt.ylabel('Void fraction [-]')
plt.xlim(2,4)
plt.ylim(0.1,1.0)
plt.grid()
plt.savefig('void-fraction-fluctuations')
plt.show()

# Plot bed height
reference_height = pd.read_csv('reference/height.csv')
plt.figure()
plt.plot(reference_height['t'],reference_height['h'], label='Reference')
plt.plot(time[start:end],bed_height, label='Lethe')
plt.legend()
plt.xlabel('Time [s]')
plt.ylabel('Bed height [m]')
plt.xlim(2,4)
plt.ylim(0.08,0.18)
plt.grid()
plt.subplots_adjust(left=0.2)
plt.savefig('bed-height')
plt.show()

# Calculate Power Spectral Density of pressure
plt.rcParams['lines.linewidth'] = 3
# Interpolate signal to have dt=1e-3 instead of 1e-2 
# and so fs~1000 (like experiment) instead of fs~100
t=time[start:end]
p=pressure - np.mean(pressure)
f_interp = interp1d(t, p, kind='linear')
n_steps = len(t) * 10
x_new = np.linspace(t.min(), t.max(), n_steps)
y_new = f_interp(x_new)
dt = (t.max() - t.min()) / n_steps
fs = 1/dt

# Plot interpolation
plt.show()
plt.xlabel('Time [s]')
plt.ylabel('Relative pressure [Pa]')
plt.plot(t, p, 'b.', label='Original data')
plt.plot(x_new, y_new, 'r-', label='Cubic interpolation')
plt.legend()
plt.title("Cubic interpolation of pressure signal")
plt.grid()
plt.tight_layout()
plt.show()

# Calculate PSD with fft
Y = rfft(y_new)
frequencies = rfftfreq(len(y_new), dt)
psd = (dt/len(y_new) ) * np.abs(Y)**2
psd[1:len(y_new)//2] *= 2 
amplitudes = np.abs(Y)
amplitudes[0] = 0  # Ignore DC component
freq_dominante = frequencies[np.argmax(amplitudes)]
f = freq_dominante

# Apply filter to PSD
psd_smoothed = gaussian_filter1d(psd, sigma=1.5)

# Plot
reference_psd = pd.read_csv('reference/psd.csv')
plt.figure()
plt.loglog(reference_psd['f'],reference_psd['PSD'], label='Reference')
plt.loglog(frequencies, psd, label='Lethe')
plt.loglog(frequencies, psd_smoothed,label='Lethe-Gauss')
plt.xlim(0.5, 100)
plt.xlabel("Frequency [Hz]")
plt.ylabel(f'PSD $[Pa^2/Hz]$')
plt.grid()
plt.subplots_adjust(left=0.2, bottom=0.15)
plt.tight_layout()
plt.legend()
plt.savefig('pressure-psd')
plt.show() 
