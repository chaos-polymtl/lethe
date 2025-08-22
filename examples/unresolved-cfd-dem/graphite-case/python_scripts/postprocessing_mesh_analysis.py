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
from functions_mesh_analysis import *

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

parser = argparse.ArgumentParser(description = 'Arguments for the post-processing of the mesh analysis of the graphite case')
parser.add_argument("-f_1", "--folder_1", type = str, help = "Folder path for the mesh 1. It contains the .prm file.", required = True)
parser.add_argument("-f_2", "--folder_2", type = str, help = "Folder path for the mesh 2. It contains the .prm file.", required = True)
#parser.add_argument("-f_3", "--folder_3", type = str, help = "Folder path for the mesh 3. It contains the .prm file.", required = True)

args, leftovers = parser.parse_known_args()

# Simulation folders
folder_1 = args.folder_1
folder_2 = args.folder_2
#folder_3 = args.folder_3

# Starting time step for post-processing
start_1 = 13
start_2 = 37
#start_3 = 0

# Load Lethe data
pvd_name = 'out.pvd'
prm_file = 'pure_flow.prm'
ignore_data = ['average_pressure', 'dummy_rss', 'pressure', 'q_criterion', 'reynolds_normal_stress', 'reynolds_shear_stress_uv', 'reynolds_shear_stress_uw', 'reynolds_shear_stress_vw', 'subdomain', 'turbulent_kinetic_energy', 'velocity', 'velocity_divergence', 'vorticity']
flow_mesh_1 = lethe_pyvista_tools(folder_1, prm_file, pvd_name, ignore_data = ignore_data, first = start_1)
flow_mesh_2 = lethe_pyvista_tools(folder_2, prm_file, pvd_name, ignore_data = ignore_data, first = start_2)
#flow_mesh_3 = lethe_pyvista_tools(folder_3, prm_file, pvd_name, ignore_data = ignore_data, first = start)

# Extract positions and velocities for each mesh
x_positions_mesh_1, avg_velocity_mesh_1 = get_positions_and_velocities(flow_mesh_1)
x_positions_mesh_2, avg_velocity_mesh_2 = get_positions_and_velocities(flow_mesh_2)
#x_positions_mesh_3, avg_velocity_mesh_3 = get_positions_and_velocities(flow_mesh_3)

# Plotting the results
plt.figure(figsize = (10, 6))
plt.plot(x_positions_mesh_1, avg_velocity_mesh_1, label = 'Mesh 1', color = colors[0], marker = 'o', markersize = 3)
plt.plot(x_positions_mesh_2, avg_velocity_mesh_2, label = 'Mesh 2', color = colors[1], marker = 'o', markersize = 3)
#plt.plot(x_positions_mesh_3, avg_velocity_mesh_3, label = 'Mesh 3', color = colors[2], marker = 'o', markersize = 3)

plt.xlabel('x position (m)')
plt.ylabel('Average Velocity (m/s)')
plt.title('Average Velocity Magnitude at y=0 and z=0.015 for Different Meshes')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('mesh_analysis.png')
plt.show()