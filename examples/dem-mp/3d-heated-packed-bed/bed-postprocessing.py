# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import matplotlib.lines as mlines
from cycler import cycler

# Set plot parameters
colors=['#1B9E77','#D95F02','#7570B3']
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

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the 3d-heated-packed-bed example')
parser.add_argument("-f", "--folder", type=str, help="Folder path. This folder is the folder which contains the .prm file.", required=True)
args, leftovers=parser.parse_known_args()

# Simulation folder
folder=args.folder

# Starting vtu id (heating starts)
start = 139 # 7s

# Height (z) of the heating wall
h_top = 0.2

# Same fixed parameters as reference
particle_radius=0.0032
heights = [0.040 , 0.060, 0.073] # Heights for mean temperatures
particle_number = 8849

# Load lethe data
prm_file = 'heat-packed-bed.prm'
pvd_name = 'out.pvd'
ignore_data = ['type', 'omega','specific_heat']
particle = lethe_pyvista_tools(folder, prm_file, pvd_name, ignore_data=ignore_data)
time = np.array(particle.time_list)

# Get kinetic knergy at the end of the first two stages 
K = []
for i in [118,138]:
    df_load = particle.get_df(i)
    vx = df_load['velocity'][:, 0]
    vy = df_load['velocity'][:, 1]
    vz = df_load['velocity'][:, 2]
    v_squared = vx*vx + vy*vy +vz*vz
    m = df_load['mass']
    K.append(0.5*np.sum(m*v_squared))

# Print the kinetic energy and the packed bed porosity
particle_volume = 4.0/3.0 * np.pi * (particle_radius)**3
total_volume = 0.105 * 0.1 * h_top
porosity = 1 - (particle_number * particle_volume) / total_volume
print('\n')
print('=' * 80)
print(f"At t = {time[118]} s, kinetic energy: {K[0]} J")
print(f"At t = {time[138]} s, kinetic energy: {K[1]} J")
print(f"Total porosity: {porosity:.4f} %")
print('=' * 80)
print('\n')

# Sampling point for the mean temperatures 
# (chosen in the middle of the packed bed to represent the thermocouples in the experiment)
x_sample = 0.05
y_sample = 0.05

# Array where the mean temperatures of the packed bed at each height are stored
time -= time[start] # remove time for loading stages
temperature = np.zeros((len(time)-start,3))

# Distance to define which particles are considered around the sampling points 
# (radius of the particles)
D = particle_radius 

for i_h, h in enumerate(heights):
    for i in range(start, len(time)):

        df_load = particle.get_df(i)
        df = pd.DataFrame(np.copy(df_load.points), columns=['x', 'y','z'])
        df['temperature'] = df_load['temperature']

        df_filtered = df.copy()

        # Compute the distance between each particle and the sampling point
        df_filtered['dist_x'] = ((df_filtered['x']-x_sample)** 2) ** 0.5
        df_filtered['dist_y'] = ((df_filtered['y']-y_sample)** 2) ** 0.5
        df_filtered['dist_z'] = ((df_filtered['z']-(h_top-h))** 2) ** 0.5
        
        # Keep the particles close to the sampling point
        df_to_sample = df_filtered[(df_filtered['dist_x'] < 2*D) & (df_filtered['dist_y'] < 2*D) & (df_filtered['dist_z'] < D)]

        # Get the mean temperature around the sampling point
        temperature[i-start][i_h] = df_to_sample['temperature'].mean()

# Plot the evolution of the temperature of the packed bed at three different heights.
plt.figure()
for i_h, h in enumerate(heights):
    # Read temperature data from paper
    reference_data_sim = pd.read_csv('reference/steel_Beaulieu_' + f"{h:.3f}" + '_rough.csv')
    reference_data_exp = pd.read_csv('reference/steel_Beaulieu_' + f"{h:.3f}" + '_experiment.csv')

    # Read and draw error
    reference_data_error_1 = pd.read_csv('reference/steel_Beaulieu_' + f"{h:.3f}" + '_error1.csv')
    reference_data_error_2 = pd.read_csv('reference/steel_Beaulieu_' + f"{h:.3f}" + '_error2.csv')
    t1, T1 = reference_data_error_1['t'], reference_data_error_1['T']
    t2, T2 = reference_data_error_2['t'], reference_data_error_2['T']
    T2_interp = np.interp(t1, t2, T2)
    plt.fill_between(t1, T1, T2_interp, color=colors[i_h], alpha=0.2, linewidth=0.7)

    # Plot
    line1, = plt.plot(time[start:],temperature[:, i_h], color=colors[i_h],label= f"Lethe")
    line2, = plt.plot(reference_data_sim['t'],reference_data_sim['T'], '--', color=colors[i_h], label= f"Beaulieu simulation")
    line3, = plt.plot(reference_data_exp['t'],reference_data_exp['T'], '.', color=colors[i_h], label= f"Beaulieu experiment")

# Add color legend
color_one   = mlines.Line2D([], [], color=colors[0], marker='s', linestyle='None', label=f"$h_1 = {heights[0]:.3f}$ m")
color_two   = mlines.Line2D([], [], color=colors[1], marker='s', linestyle='None', label=f"$h_2 = {heights[1]:.3f}$ m")
color_three = mlines.Line2D([], [], color=colors[2], marker='s', linestyle='None', label=f"$h_3 = {heights[2]:.3f}$ m")
legend = plt.legend(handles=[line1, line2, line3, color_one, color_two, color_three],loc='upper left',fontsize=18)
legend.legend_handles[0].set_color('k')
legend.legend_handles[1].set_color('k')
legend.legend_handles[2].set_color('k')

plt.grid()
plt.ylim(19,34)
plt.xlabel('Time (s) after packing')
plt.ylabel('Temperature (°C)')
plt.savefig('mean-temperatures')
plt.show()