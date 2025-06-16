# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import numpy as np
import pandas as pd
import argparse

import sys
sys.path.append("$LETHE_PATH/contrib/postprocessing/")
from lethe_pyvista_tools import *

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the 3d-thermal-packed-bed example')
parser.add_argument("-f", "--folder", type=str, help="Folder path. This folder is the folder which contains the .prm file.", required=True)
parser.add_argument("--start", type=int, help="Starting vtu. Last vtu of the loading part of the simulation.", required=True)
parser.add_argument("--htop", type=float, help="Height of the top wall. Height of the top wall after it has stopped moving.", required=True)

args, leftovers=parser.parse_known_args()

# Simulation folder
folder=args.folder

# Starting vtu id (heating starts)
start = args.start

# Height of the top wall once still
h_top = args.htop

# Same fixed parameters as reference
particle_radius=0.0032
prm_file = 'heat-packed-bed.prm'
# Heights for mean temperatures
heights = [0.040 , 0.060, 0.073]
# Number of particles
particle_number = 8849
ylim=(19,34)
    
# Load lethe data
pvd_name = 'out.pvd'
ignore_data = ['type', 'mass', 'velocity' , 'omega']
particle = lethe_pyvista_tools(folder, prm_file, pvd_name, ignore_data=ignore_data)
time = np.array(particle.time_list)
time -= time[start]

# Sampling point for the mean temperatures 
# (chosen in the middle of the packed bed to represent the thermocouple in the experiment)
x_sample = 0.05
y_sample = 0.05

# Values where the mean temperature of the packed bed at each height is stored
temperature = np.zeros((len(time)-start,3))

# Distance to define which particles are considered around the sampling points 
# (diameter of the particles)
D = 2*particle_radius 

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
        # Keep the particles close to the sampled point
        df_to_sample = df_filtered[(df_filtered['dist_x'] < 2*D) & (df_filtered['dist_y'] < 2*D) & (df_filtered['dist_z'] < D)]
        # Get the mean temperature around the sampling point
        temperature[i-start][i_h] = df_to_sample['temperature'].mean()

# Print packed bed porosity
particle_volume = 4.0/3.0 * np.pi * (particle_radius)**3
total_volume = 0.102 * 0.1 * h_top
porosity = 1 - (particle_number * particle_volume) / total_volume
print('\n')
print('=' * 80)
print(f"Total porosity (%) : {porosity:.4f}")
print('=' * 80)
print('\n')
# Write the porosity on a file
with open(folder + "/porosity.txt", "w") as file:
    file.write("Porosity: \n")
    file.write(f"{porosity:.4f} %\n")

# Export simulation data to csv
temperature_data = pd.DataFrame({'time': time[start:]})
for i_h, h in enumerate(heights):
    temperature_data[f'T_{h:.3f}'] = temperature[:, i_h]
temperature_data.to_csv(folder +'/temperatures.csv')