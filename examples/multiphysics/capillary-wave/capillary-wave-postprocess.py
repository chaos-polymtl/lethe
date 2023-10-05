#############################################################################
"""
Postprocessing code for the capillary wave example.
This code extracts and compares the evolution of the relative amplitude of
the wave with the analytical solution from Prosperetti [1].

*** Before running this script, make sure to generate a file with analytical
solution using capillary-wave-prosperetti-solution.py ***
"""
#############################################################################

#############################################################################
'''Importing Libraries'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("$LETHE_PATH/contrib/postprocessing/")
from lethe_pyvista_tools import *
#############################################################################

#############################################################################
# Run script: python3 path_to_capillary-wave-postprocess.py path_to_case prm_file_name path_to_analytical_solution_csv_file

# Check the number of input arguments
if len(sys.argv)!=4:
    print("********************************************************************************\n"
          "Incorrect number of arguments\n"
          "Run script with: \n"
          "\t path_to_capillary-wave-postprocess.py path_to_case prm_file_name path_to_analytical_solution_csv_file\n"
          "********************************************************************************\n")
    exit(1)

#############################################################################
#----------------------------------
# Read, extract and compute
# quantities of interest
#----------------------------------
# Parse arguments
simulation_path = sys.argv[1]
prm_file_name = sys.argv[2]
analytical_path = sys.argv[3]

# Name of the pvd file
time_step_multiplier = prm_file_name[19:-4]
pvd_file_name = prm_file_name[:18] + f"-{time_step_multiplier}.pvd"
time_step_multiplier=time_step_multiplier.replace('_', '.')

# Create the fluids object
fluids = lethe_pyvista_tools(simulation_path, prm_file_name, pvd_file_name)

# Set phase_limit to search for height values
phase_limit = 0.5

# Get active times
time_list = fluids.time_list

# Create list to hold wave heights
H = []

# Create beginning and end points for height lines
H_a = [0,-2e-6, 0]
H_b = [0, 2e-6, 0]

# Read VTU files
for i in range(len(fluids.list_vtu)):
    # Store results in 'df'
    df = fluids.get_df(i)

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

# Calculate relative amplitude
relative_amplitude = [h/amplitude_0 for h in H]

# Read verification data from Prosperetti [1]
validation_file_name = analytical_path
analytical_values = pd.read_csv(validation_file_name)
analytical_time = np.array(analytical_values['t'])
analytical_solution = np.array(analytical_values['a/a0'])

#############################################################################
#----------------------------------
# Plot
#----------------------------------
# Figure name
output_path = fluids.prm_dict["output path"]
figure_name = output_path[9:].replace('/', '')
figure_name = figure_name.replace('.', '')

# Plot
plt.rcParams['font.size'] = '18'
fig0 = plt.figure(figsize=(12, 8))
ax0 = fig0.add_subplot(111)
plt.plot(time_list, relative_amplitude, "s", mfc='none', label=f'Lethe — $\Delta t= {time_step_multiplier} \Delta t_\sigma$')
plt.plot(analytical_time, analytical_solution, "--k", linewidth=2, label=f'Prosperetti (1981)')
plt.xlabel('Time')
plt.ylabel('Relative amplitude')
plt.ticklabel_format(axis='both', style='sci', scilimits=(4,-9))
plt.legend(loc="best")
plt.tight_layout()
plt.savefig(f'{simulation_path}/figure-' + figure_name + '.png',dpi=300)

#----------------------------------
# Write csv file with values
#----------------------------------
lethe_df = pd.DataFrame({'t': time_list, 'a/a0': relative_amplitude})
lethe_df.to_csv(f'{simulation_path}/lethe-' + figure_name + '.csv', index=False)

#############################################################################
#----------------------------------
# References
#----------------------------------

# [1] A. Prosperetti, “Motion of two superposed viscous fluids,” Phys. Fluids, vol. 24, no. 7, pp. 1217–1223, Jul. 1981, doi: 10.1063/1.863522.