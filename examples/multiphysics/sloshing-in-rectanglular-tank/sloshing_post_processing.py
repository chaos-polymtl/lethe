#############################################################################
"""
Postprocessing code for the sloshing in a fixed rectangular tank example.
This code extracts and compares the evolution of the relative amplitude of
the wave.
"""
#############################################################################

#############################################################################
'''Importing Libraries'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys
sys.path.append("../../../contrib/postprocessing/")
from lethe_pyvista_tools import *

#############################################################################

#############################################################################
# Run script : python3 path_to_sloshing_post_processing.py path_to_case prm_file_name

#Take case path as argument
simulation_path = sys.argv[1]
prm_file_name = sys.argv[2]

# Name the save path
save_path = simulation_path

# Create the fluids object
fluids = lethe_pyvista_tools(simulation_path, prm_file_name)

# Get the pvd name
pvd_name = fluids.prm_dict["output name"]

# Read data
fluids.read_lethe_to_pyvista(f'{pvd_name}.pvd')

# Set phase_limit to search for height values
phase_limit = 0.5

# Get active times
time_list = fluids.time_list

# Create list to hold wave heights
H = []

# Create beginning and end points for height lines
H_a = [0,-0.02, 0]
H_b = [0, 0.02, 0]

# Read VTU files
for i in range(len(fluids.list_vtu)):
    # Store results in 'df'
    df = fluids.df[i]

    # Extract phase values and points over the center line
    sampled_data = df.sample_over_line(H_a, H_b, resolution=5000)
    phase_over_line = pd.DataFrame(sampled_data["phase"])
    points_over_line = pd.DataFrame(sampled_data.points)

    # Find highest point of fluid 1
    fluid1_points = points_over_line[phase_over_line[0] > phase_limit].values
    H_max = fluid1_points[0][1]
    for point in fluid1_points:
        if point[1] > H_max:
            H_max = point[1]
    if i == 0:
        amplitude_0 = H_max

    H.append(H_max)

# Calculate relative
relative_amplitude = [h/amplitude_0 for h in H]

# Validation data
nu = fluids.prm_dict["kinematic viscosity"]
g = 1
Re = np.sqrt(g)/nu
validation_file_name = f"{simulation_path}/analytical_solution_Re{Re:.0f}.csv"
analytical_values = pd.read_csv(validation_file_name)
analytical_time = np.array(analytical_values['x'])
analytical_solution = np.array(-analytical_values[f"Re{Re:.0f}"])

# Figures
plt.rcParams['font.size'] = '16'
fs1 = 16
fs2 = 14

# Figure name
output_path = fluids.prm_dict["output path"]
figure_name = output_path[9:-1]

# Plot heights
fig0 = plt.figure(figsize=(19, 12))
ax0 = fig0.add_subplot(111)
plt.plot(time_list, relative_amplitude, "k", linewidth=2, label=f'Re = {Re:.0f} - Lethe')
plt.plot(analytical_time, analytical_solution, "--r", linewidth=2, label=f'Re = {Re:.0f} - Analytical')
plt.xlabel('Time')
plt.ylabel('Relative amplitude')
plt.legend(loc="best")
plt.savefig(f'{simulation_path}/figure_' + figure_name + '.png')
plt.show()