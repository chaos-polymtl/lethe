#############################################################################
"""
Postprocessing automation tool.

By: Victor Oliveira Ferreira, Audrey Collard-Daigneault
Date: May 5th, 2022
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
'''Simulation properties'''

#Take case path as argument
simulation_path = sys.argv[1]
prm_file_name = sys.argv[2]

# Name the save path
save_path = simulation_path

# Create the particle object
pvd_name = 'hopper.pvd'
ignore_data = ['type', 'diameter', 'volumetric contribution']
particle = lethe_pyvista_tools(simulation_path, prm_file_name, pvd_name, ignore_data=ignore_data)

#############################################################################
# Beginning of flow (after loading particles)
start = int(particle.prm_dict['end time'] /
            (particle.prm_dict['output frequency'] * particle.prm_dict['time step'])) - 1

# Position and normal vector of the outlet
hopper_outlet = 0
if int(particle.prm_dict['nx']):
    normal_vect = 0
elif int(particle.prm_dict['ny']):
    normal_vect = 1
elif int(particle.prm_dict['nz']):
    normal_vect = 2
else:
    sys.exit('No normal vector for Hopper outlet.')

# Create a list to store the number of particles below the outlet
number_of_particles_below = []
mass_discharge = []
volume = 4./3. * np.pi * (particle.prm_dict['diameter']/2.)**3

# Correction regards data in paper (nb of particle)
# Thickness of hopper domain was increased and number of particle was multiplied since the simulated case has no
# periodic boundaries unlike simulations in the paper.
# The correction factor is corresponding to the multiplication of the width/the number of particle
n_part_paper = 6790
correction_factor = particle.prm_dict['number']/n_part_paper

# Total mass = nρV/factor
n_particle = particle.prm_dict["number"]
density = particle.prm_dict["density particles"]
print(f'Total mass in hopper : {n_particle * density * volume / correction_factor * 1000:.2f} g.')

# Create a list to store the "flow rate" of particles
rate = []

# Loop through all results
for i in range(len(particle.list_vtu)):
    # Store results in 'df'
    df = particle.get_df(i)

    # Select the data (if particle is completely under hopper outlet)
    vertical_position = df.points[:, normal_vect]
    vertical_position_below = vertical_position[vertical_position < (hopper_outlet - particle.prm_dict['diameter']/2.)]

    # Number of particles below and mass discharge
    number_of_particles_below.append(len(vertical_position_below))
    mass_discharge.append(len(vertical_position_below) * density * volume)

    # Calculate the rate of discharge
    if i == 0:
        rate.append(0.0)
    else:
        rate.append((number_of_particles_below[i] - number_of_particles_below[i - 1])/particle.prm_dict['time step'])

# Export data to csv
data = pd.DataFrame({'time': particle.time_list, 'rate': rate,
                    'number_of_particles': number_of_particles_below, 'mass_discharge': mass_discharge})

# Export data to csv
data.to_csv(save_path + '/results_' + pvd_name + '.csv')

# Read data from paper
paper_data = pd.read_csv('paper_data.csv')

# Find range to calculate rate (this part is kind of headcoded)
p0 = start + int(0.25/(particle.prm_dict['output frequency'] * particle.prm_dict['time step']))
p1 = p0 +    int(0.5 /(particle.prm_dict['output frequency'] * particle.prm_dict['time step']))

# Calculate mass flow rate
p = np.polyfit([value - particle.time_list[p0] for value in particle.time_list[p0:p1]],
               [value * 1000 / correction_factor  for value in mass_discharge[p0:p1]], 1)
print(f'Mass flow rate is : {p[0]:.2f} g/s.')

# Plot results
plt.plot(data['time'][start:] - data['time'][start], data['mass_discharge'][start:] * 1000 / correction_factor, label="Lethe DEM")
plt.plot(paper_data['time'], paper_data['discharge'], '.k', label="Anand et al.")
plt.xlabel('Time (s)')
plt.ylabel('Mass discharged from the hopper (g)')
plt.legend()
plt.grid()
plt.savefig('figure_' + pvd_name + '.png')
plt.show()


