#############################################################################
"""
Postprocessing code for the water injection in a closed cell example.
This code extracts and compares the evolution of the density and pressure
in the air
"""
#############################################################################

#############################################################################
'''Importing Libraries'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys
sys.path.append("../../../../contrib/postprocessing/")
from lethe_pyvista_tools import *

#############################################################################

#############################################################################
# Run script : python3 path_to_water-injection-in-a-closed-cell-postprocessing.py path_to_case prm_file_name

#Take case path as argument
simulation_path = sys.argv[1]
prm_file_name = sys.argv[2]

# Name the save path
save_path = simulation_path

# Create the fluids object
fluids = lethe_pyvista_tools(simulation_path, prm_file_name, 'water-injection-in-a-closed-cell.pvd')

# Get active times
time_list = fluids.time_list

# Create list to collect air density and pressure values
density = []
pressure = []

# Sampling points
sample_point_a = [-0.5, 0, 0]
sample_point_b = [0.5, 0, 0]

# Read VTU files
for i in range(len(fluids.list_vtu)):
    # Store results in 'df'
    df = fluids.get_df(i)

    # Extract density and pressure at the center of the cavity
    sampled_data = df.sample_over_line(sample_point_a,sample_point_b,resolution=2)
    density_over_line = pd.DataFrame(sampled_data["density_00"])
    pressure_over_line = pd.DataFrame(sampled_data["pressure"])
    density.append(density_over_line.values[1])
    pressure.append(pressure_over_line.values[1])

# Analytical solutions
density_ref = fluids.prm_dict["density_ref"]
volumetric_flow_rate = 0.1
initial_air_height = 0.1
R = fluids.prm_dict["R"]
T = fluids.prm_dict["T"]
# Time
tf = time_list[len(time_list)-1]
analytical_time = np.linspace(0,tf,1000)
# Density
analytical_density=[density_ref/(1-(volumetric_flow_rate*t/initial_air_height)) for t in analytical_time]
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
plt.savefig(f'{simulation_path}/figure-' + figure_name + '-density.png', dpi=300)

# Plot pressure
fig1 = plt.figure(figsize=(9.5, 6))
ax1 = fig1.add_subplot(111)
plt.plot(time_list, pressure, "sk", mfc="none", linewidth=2, label=f'Pressure – Lethe')
plt.plot(analytical_time, analytical_pressure, "--r", linewidth=2, label=f'Pressure – Analytical')
plt.xlabel('Time [s]')
plt.ylabel('Pressure [Pa]')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig(f'{simulation_path}/figure-' + figure_name + '-pressure.png', dpi=300)

plt.show()