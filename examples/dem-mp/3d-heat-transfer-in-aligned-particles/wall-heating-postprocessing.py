# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from cycler import cycler

# Set plot parameters
colors=['#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02']
plt.rcParams['axes.prop_cycle'] = cycler(color = colors)
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize'] = (10,8)
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.columnspacing'] = 2
plt.rcParams['legend.handlelength'] = 1
plt.rcParams['legend.handletextpad'] = 0.2
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fancybox'] = False
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['font.size'] = '20'
plt.rcParams['font.family']='DejaVu Serif'
plt.rcParams['font.serif']='cm'
plt.rcParams['savefig.bbox']='tight'

######################################################################

import sys
sys.path.append("$LETHE_PATH/contrib/postprocessing/")
from lethe_pyvista_tools import *

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the 3d-thermal-linear case')
parser.add_argument("-f", "--folder", type=str, help="Folder path. This folder is the folder which contains the .prm file.", required=True)
args, leftovers=parser.parse_known_args()

# Simulation folder
folder=args.folder

# Load lethe data
pvd_name = 'out.pvd'
prm_file = 'wall-heating.prm'
ignore_data = ['type', 'diameter', 'volumetric contribution', 'torque', 'fem_torque','fem_force']
particle = lethe_pyvista_tools(folder, prm_file, pvd_name, ignore_data=ignore_data)
time = np.array(particle.time_list)

# Get position and temperature of particles at last frame
df_load = particle.get_df(len(time)-1)
df = pd.DataFrame(np.copy(df_load.points[:, 0]), columns=['x'])
df['temperature'] = df_load['temperature']
x = df['x'].to_numpy()
T_sim = df['temperature'].to_numpy()

# Calculate analytical solution
r = 0.005 # Particle radius
T_left = 0 # Temperature of left wall
T_right = 50 # Temperature of right wall
x_left = -0.045 # Position of left wall after it has stopped
x_right = 0.045 # Position of right wall after it has stopped
A = (T_left - T_right)/(x_left-x_right)
B = T_left - A*x_left
T_analytical = A*x + B

# Fit simulation results
a, b = np.polyfit(x, T_sim, deg=1)
T_fit = a*x + b

# Print some results
overlaps = 2*r -abs(x[1:] - x[:-1])
mean_overlap = np.mean(overlaps)
print('\n')
print('=' * 80)
print(f"Mean overlap : {mean_overlap:.6f} m")
print(f"Analytical solution : {A:.2f} x + {B:.2f}")
print(f"Simulation : {a:.2f} x + {b:.2f}")
print('=' * 80)
print('\n')

# Plot the temperatures at last frame
plt.figure()
plt.plot(x, T_sim, 'kx', label= "Lethe")
# plt.plot(x, T_fit, '--', alpha = 0.7, label= "Linear regression")
plt.plot(x, T_analytical, alpha=0.5, label= "Analytical solution")
plt.ylim(T_left,T_right)
plt.legend()
plt.grid()
plt.xlabel('x (m)')
plt.ylabel('Temperature (°C)')
plt.savefig('temperature-linear')
plt.show()