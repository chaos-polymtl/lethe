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

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the thermal equilibrium example')
parser.add_argument("-f", "--folder", type=str, help="Folder path. This folder is the folder which contains the .prm file.", required=True)
args, leftovers=parser.parse_known_args()

# Simulation folder
folder=args.folder

# Load lethe data
pvd_name = 'out.pvd'
prm_file = 'equilibrium.prm'
ignore_data = ['type', 'mass','omega','velocity']
particle = lethe_pyvista_tools(folder, prm_file, pvd_name, ignore_data=ignore_data)
time = np.array(particle.time_list)

# Particle radius
r = 0.5 * particle.get_df(0)['diameter'][0]

# Store the mean temperatures
mean_temperature_left_x = np.zeros(len(time))
mean_temperature_right_x = np.zeros(len(time))

# Values where the overlaps are stored
mean_overlap = np.zeros(len(time))

# Calculate the mean temperatures (left: x<0, right: x>0)
for i in range(len(time)):

    df_load = particle.get_df(i)
    df = pd.DataFrame(np.copy(df_load.points[:, 0]), columns=['x'])
    df['temperature'] = df_load['temperature']

    for j in range(len(df)):

        x = df['x'][j]
        T = df['temperature'][j]

        if x<0:
            mean_temperature_left_x[i] += T
        elif x>0:
            mean_temperature_right_x[i] += T

    mean_temperature_right_x[i] *= 2/len(df)
    mean_temperature_left_x[i] *= 2/len(df)

# Print some results
diff = np.abs(mean_temperature_left_x - mean_temperature_right_x)
indices = np.where(diff < 0.01)[0]
ids=["0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9"]
overlaps = [2*r -abs(df['x'][i]-df['x'][i+1]) for i in range(0,9)]
mean_overlap = np.mean(overlaps)
print('\n')
print('=' * 80)
if len(indices) > 0:
    time_about_equal = time[indices[0]]
    print(f"After {time_about_equal:.2f} s, the two mean temperatures have a difference of less than 1%.")
else:
    print("The two mean temperatures never reach a difference below 1% during the simulation. Simulation time should be increased")
print('=' * 80)
print("| Ids | Overlaps (m) at the end")
for i in range(9):
    print(f"| {ids[i]} | {overlaps[i]:.5e}")
print(f"Mean overlap: {mean_overlap:.5e} m.")
print('=' * 80)
print('\n')

# Plot the evolution of the mean temperatures
plt.figure()
plt.plot(time, mean_temperature_left_x, label= "Mean temperature x<0")
plt.plot(time, mean_temperature_right_x, label= "Mean temperature x>0")
plt.legend()
plt.grid()
plt.title("Evolution of the mean temperatures", pad=25)
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')
plt.savefig('mean-temperatures')
plt.show()






