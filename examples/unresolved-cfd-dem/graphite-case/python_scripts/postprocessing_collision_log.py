# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import pandas as pd
from cycler import cycler
from functions_collision_log import *

# Set plot parameters
plt.rcParams['lines.markersize'] = 11
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.fancybox'] = False
plt.rcParams['font.size'] = 20
plt.rcParams['legend.handlelength'] = 2
plt.rcParams['lines.linewidth'] = 3
colors = ['#7570B3','#E7298A','#66A61E','#E6AB02']
plt.rcParams['axes.prop_cycle'] = cycler(color = colors)


# Set simulations parameters for compariason between different collisions logs from different simulations
# Define the variables for how many different collision logs you want to analyze
inlet_velocity_1 = 5 # m/s
distance_between_nozzle_and_wall_1 = 0.01 # m

######################################################################

# Read the CSV file containing collision statistics
# Load how many different collision logs you want to analyze
csv_file_1 = pd.read_csv('collision_statistics.csv')

# Function to calculate collision frequency by particle diameter
df_collision_frequency_1 = collision_frequency_by_diameter(csv_file_1)

# Function to plot collision frequency by particle diameter
plot_collision_frequency(df_collision_frequency_1, output_png = 'collision_frequency_histogram.png')

# Function to calculate the mean speed difference between the start and end of a collision, classified by particle diameter
mean_delta_speed_1 = mean_speed_difference(csv_file_1)
