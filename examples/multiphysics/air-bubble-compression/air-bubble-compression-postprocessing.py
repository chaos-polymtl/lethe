# SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#############################################################################
'''Importing Libraries'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

import sys
sys.path.append("$LETHE_PATH/contrib/postprocessing/")
from lethe_pyvista_tools import *

#############################################################################

#############################################################################
# Run script : python3 path_to_air-bubble-compression-postprocessing.py path_to_case prm_file_name

#Take case path as argument
parser = argparse.ArgumentParser(description='Arguments for the post-processing of the air bubble compression example')
parser.add_argument("simulation_path", type=str, help="Path to the folder which contains the .prm file")
parser.add_argument("prm_file_name", type=str, help="Name of the .prm file")
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
args, leftovers=parser.parse_known_args()

simulation_path = args.simulation_path
prm_file_name = args.prm_file_name

# Name the save path
save_path = simulation_path

# Create the fluids object
fluids = lethe_pyvista_tools(simulation_path, prm_file_name, 'air-bubble-compression.pvd')

# Get active times
time_list = fluids.time_list

# Create list to collect air density and pressure values
density = []
pressure = []

# Sampling points
sample_point_a = [-0.5, 0, 0]
sample_point_b = [0.5, 0, 0]

# Read PVTU files
for i in range(len(fluids.list_vtu)):
    # Store results in 'df'
    df = fluids.get_df(i)

    # Extract density and pressure at the center of the bubble
    sampled_data = df.sample_over_line(sample_point_a,sample_point_b,resolution=2)
    density_over_line = pd.DataFrame(sampled_data["density_00"])
    pressure_over_line = pd.DataFrame(sampled_data["pressure"])
    density.append(density_over_line.values[1])
    pressure.append(pressure_over_line.values[1])

# Analytical solutions
density_ref = fluids.prm_dict["density_ref"]
volumetric_flow_rate = 0.0025*0.1*4
bubble_initial_radius = 0.03
R = fluids.prm_dict["R"]
T = fluids.prm_dict["T"]
# Time
tf = time_list[len(time_list)-1]
analytical_time = np.linspace(0,tf,1000)
# Density
analytical_density=[density_ref/(1-(volumetric_flow_rate*t/(np.pi*bubble_initial_radius*bubble_initial_radius))) for t in analytical_time]
# Pressure
analytical_pressure=[(rho-density_ref)*R*T for rho in analytical_density]

# Fontsize for plots
plt.rcParams['font.size'] = '16'

# Figure name
figure_name = fluids.prm_dict["output name"]

# Plot density
fig0 = plt.figure(figsize=(9.5, 6))
ax0 = fig0.add_subplot(111)
plt.plot(time_list, density, "sk", mfc="none", linewidth=2, label=f'Density – Lethe')
plt.plot(analytical_time, analytical_density, "--r", linewidth=2, label=f'Density – Analytical')
plt.xlabel('Time [s]')
plt.ylabel(r'Density [kg/m$^3$]')
plt.legend(loc="best")
plt.tight_layout()
if (not args.validate):
    plt.savefig(f'{simulation_path}/figure-' + figure_name + '-density.png', dpi=300)
else:
    plt.savefig(f'{simulation_path}/figure-' + figure_name + '-density.pdf')

# Plot pressure
fig1 = plt.figure(figsize=(9.5, 6))
ax1 = fig1.add_subplot(111)
plt.plot(time_list, pressure, "sk", mfc="none", linewidth=2, label=f'Pressure – Lethe')
plt.plot(analytical_time, analytical_pressure, "--r", linewidth=2, label=f'Pressure – Analytical')
plt.xlabel('Time [s]')
plt.ylabel('Pressure [Pa]')
plt.legend(loc="best")
plt.tight_layout()
if (not args.validate):
    plt.savefig(f'{simulation_path}/figure-' + figure_name + '-pressure.png', dpi=300)
    plt.show()
else:
    plt.savefig(f'{simulation_path}/figure-' + figure_name + '-pressure.pdf')

    # Store the plotted values, with the analytical solution at the same
    # times, so that they can be compared between validation runs
    np.savetxt(f'{simulation_path}/solution-density-pressure.dat',
               np.column_stack((time_list,
                                np.ravel(density),
                                np.interp(time_list, analytical_time, analytical_density),
                                np.ravel(pressure),
                                np.interp(time_list, analytical_time, analytical_pressure))),
               header="t density density_analytical pressure pressure_analytical")
