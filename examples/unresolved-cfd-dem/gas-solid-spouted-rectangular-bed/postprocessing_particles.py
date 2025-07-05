# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import argparse
from cycler import cycler
from collections import defaultdict
from functions_particles import *

# Set plot parameters
plt.rcParams['lines.markersize'] = 11
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.fancybox'] = False
plt.rcParams['font.size'] = 20
plt.rcParams['legend.handlelength'] = 2
plt.rcParams['lines.linewidth'] = 3
colors = ['#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02']
plt.rcParams['axes.prop_cycle'] = cycler(color = colors)

######################################################################

import sys
sys.path.append("$LETHE_PATH/contrib/postprocessing/")
from lethe_pyvista_tools import *

parser = argparse.ArgumentParser(description = 'Arguments for the post-processing of the gas-solid-spouted-rectangular-bed unresolved CFD-DEM example')
parser.add_argument("-f", "--folder", type = str, help = "Folder path. This folder is the folder which contains the .prm file.", required = True)
args, leftovers = parser.parse_known_args()

# Simulation folder
folder = args.folder

# Starting time step for post-processing
start = 200

# Load Lethe data
pvd_name = 'out_particles.pvd'
prm_file = 'gas-solid-spouted-rectangular-bed.prm'
ignore_data = ['ID', 'diameter', 'fem_force', 'fem_torque', 'mass', 'omega', 'type', 'volumetric_contribution']
particle = lethe_pyvista_tools(folder, prm_file, pvd_name, ignore_data = ignore_data, first = start)
time_list = particle.time_list

############################################################

# Heights
h0 = 0.2  # Static bed height in meters
ratio_h = np.array([0.15, 0.30, 0.45, 0.60, 0.75, 0.80])
h_list = ratio_h * h0

# Get min and max x for grid spacing
first_df = pd.DataFrame(np.copy(particle.get_df(0).points), columns = ['x', 'y', 'z'])
x_min, x_max = first_df['x'].min(), first_df['x'].max()
x_centers = np.linspace(x_min, x_max, 200)

# Define search windows
dx = 0.004
dy = 0.004

# Initialize stats
stats = defaultdict(lambda: {'count': 0, 'mean': 0.0, 'sq_diff_accumulator': 0.0})

############################################################
# Main loop over time steps

for t in range(len(time_list)):
    # Get positions and velocity array
    df_positions = pd.DataFrame(np.copy(particle.get_df(t).points), columns = ['x', 'y', 'z'])
    velocities = np.copy(particle.get_df(t)['velocity'])
    # Get the magnitude of the velocities
    velocity_mag = np.linalg.norm(velocities, axis = 1)

    for H in h_list:
        # Group by x at this height and apply Welford's algorithm
        grouped_by_x = group_velocities_by_x(H, df_positions, velocity_mag, x_centers, dx = dx, dy = dy)
        welford(grouped_by_x, stats)

############################################################
# Organize results in a dictionary

simulation_data = {} # Dictionary to store simulation data

for H in h_list:
    entries = []
    for (h, x), profile_stats in stats.items():
        if np.isclose(h, H):
            # Calculate the standard deviation
            std_dev = standard_deviation(profile_stats)
            entries.append([x, profile_stats['mean'], std_dev, profile_stats['count']])

    entries.sort(key = lambda e: e[0])
    
    # Store the results: x position, mean velocity magnitude, standard deviation, and count
    simulation_data["simulation_data_" + str(int(H / h0 * 100)).zfill(3)] = np.array(entries) if entries else np.empty((0, 4))

############################################################
# Experimental data

paper_data = {} # Dictionary to store experimental data

# Read velocity magnitude data from paper
for val in ratio_h:
    key = "paper_data_" + str(int(val * 100)).zfill(3)
    file_path = "data/spouted_bed_data_" + str(int(val * 100)).zfill(3) + ".csv"
    paper_data[key] = pd.read_csv(file_path)

############################################################
# Plot results

keys = [["paper_data_015", "paper_data_045", "paper_data_075"],
        ["paper_data_030", "paper_data_060", "paper_data_080"]]

color_map = dict(zip((keys[0] + keys[1]), colors))

fig, axs = plt.subplots(2, 1, figsize=(12, 18))

custom_lines = [
    Line2D([0], [0], color = 'black', lw = 3, linestyle = '-', label = 'Simulation'),
    Line2D([0], [0], marker = 'o', color = 'w', markerfacecolor = 'black', markersize = 11, label = 'Yue et al.')
]

# Plot group 1
plot_data(keys, 0, simulation_data, paper_data, axs, color_map, custom_lines)

# Clear extra labels for second plot
custom_lines = custom_lines[:2]

# Plot group 2
plot_data(keys, 1, simulation_data, paper_data, axs, color_map, custom_lines)

plt.tight_layout()
plt.show()