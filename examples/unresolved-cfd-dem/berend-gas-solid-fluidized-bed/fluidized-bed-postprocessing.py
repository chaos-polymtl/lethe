# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from scipy import signal
from cycler import cycler

# Set plot parameters
colors=['#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02']
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

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the 2d-sandpile DEM example')
parser.add_argument("-f", "--folder", type=str, help="Folder path. This folder is the folder which contains the .prm file.", required=True)
args, leftovers=parser.parse_known_args()

# Simulation folder
folder=args.folder

# Starting vtu id (~2s at least)
start = 200
end = 1000

# Load lethe data
pvd_particles = 'out_particles.pvd'
pvd_fluid     = 'out.pvd'
prm_file = 'gas-solid-fluidized-bed.prm'
particles = lethe_pyvista_tools(folder, prm_file, pvd_particles)
fluid = lethe_pyvista_tools(folder, prm_file, pvd_fluid)
time = np.array(particles.time_list)

# Get mean pressures and void fractions on the plane at 45mm above the floating wall
pressure = np.zeros(end-start)
void_fraction = np.zeros(end-start)
bed_height = np.zeros(end-start)
y_values = [0.001, 0.003, 0.005, 0.007]

for i in range(start, end):

    df_fluid = fluid.get_df(i)
    sampled_pressures = []
    sampled_void_fractions = []

    for y in y_values:
        sampled_data = df_fluid.sample_over_line([0, y, 0.045], [0.09, y, 0.045])
        sampled_pressures.append(np.mean(sampled_data['pressure']))
        sampled_void_fractions.append(np.mean(sampled_data['void_fraction']))

    pressure[i - start] = np.mean(sampled_pressures)
    void_fraction[i - start] = np.mean(sampled_void_fractions)

    df_particles = particles.get_df(i)
    df_loc = pd.DataFrame(np.copy(df_particles.points), columns=['x', 'y','z'])
    bed_height[i-start] = df_loc['z'].max()



# Plot mean pressure and void fraction over time
reference_pressure = pd.read_csv('reference/relative_pressure.csv')
plt.figure()
plt.plot(time[start:end],pressure-np.mean(pressure),label='Lethe')
plt.plot(reference_pressure['t'],reference_pressure['p'], label='Ref')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Relative pressure (Pa)')
plt.xlim(3,6)
plt.ylim(-300,300)
plt.grid()
plt.subplots_adjust(left=0.2)
plt.savefig('pressure-fluctuations')
plt.show()

reference_voidage = pd.read_csv('reference/voidage.csv')
plt.figure()
plt.plot(time[start:end],void_fraction, label='Lethe')
plt.plot(reference_voidage['t'],reference_voidage['voidage'], label='Ref')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Void fraction')
plt.yticks(np.arange(0.1, 0.7, 0.1))
plt.xticks(np.arange(2, 4, 0.25))
plt.xlim(2,4)
plt.ylim(0.1,0.7)
plt.grid()
plt.savefig('void-fraction-fluctuations')
plt.show()


# Calculate Power Spectral Density of pressure
reference_psd = pd.read_csv('reference/psd.csv')
plt.figure()
dt = time[1] - time[0]
fs = 1 / dt
print(fs)
f, P = signal.periodogram(pressure-np.mean(pressure),fs)
plt.loglog(f, P, label='Lethe')
plt.loglog(reference_psd['f'],reference_psd['PSD'], label='Ref')
plt.legend()
plt.xlim(0.5, 100)
plt.ylim(1e-1,1e7)
plt.xlabel('frequency [Hz]')
plt.ylabel(f'PSD $[Pa^2/Hz]$')
plt.grid()
plt.subplots_adjust(left=0.2)
plt.savefig('pressure-psd')
plt.show()

# N = len(pressure)
# fft_vals = np.fft.fft(pressure-np.mean(pressure))
# fft_freqs = np.fft.fftfreq(N, 1/fs)
# psd = (1 / (fs * N)) * np.abs(fft_vals)**2
# psd[1:N//2] *= 2 
# plt.loglog(fft_freqs, psd, label='Lethe')
# plt.loglog(reference_psd['f'],reference_psd['PSD'], label='Ref')
# plt.legend()
# plt.xlim(0.5, 100)
# plt.ylim(1e-1,1e7)
# plt.xlabel('frequency [Hz]')
# plt.ylabel(f'PSD $[Pa^2/Hz]$')
# plt.grid()
# plt.subplots_adjust(left=0.2)
# plt.show()

# reference_height = pd.read_csv('reference/height.csv')
# plt.figure()
# plt.plot(time[start:end],bed_height, label='Lethe')
# plt.plot(reference_height['t'],reference_height['h'], label='Ref')
# plt.legend()
# plt.xlabel('Time (s)')
# plt.ylabel('Bed height (m)')
# plt.yticks(np.arange(0.05, 0.25, 0.05))
# plt.xticks(np.arange(2, 4, 0.25))
# plt.xlim(2,4)
# plt.ylim(0.05,0.25)
# plt.grid()
# plt.subplots_adjust(left=0.2)
# plt.savefig('bed-height')
# plt.show()