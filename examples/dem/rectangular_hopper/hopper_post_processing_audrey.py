#############################################################################
"""
Postprocessing automation tool.

By: Victor Oliveira Ferreira
Date: May 5th, 2022
"""
#############################################################################

#############################################################################
'''Importing Libraries'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Lethe_pyvista_tools import *

import sys
#############################################################################

#############################################################################
'''Simulation properties'''

#Take case path as argument
simulation_path = sys.argv[1]

# Name the save path
save_path = simulation_path

# Create the particle object
particle = Lethe_pyvista_tools(simulation_path, prm_file_name='hopper_wiebke.prm')

# Read data
particle.read_lethe_to_pyvista('sim3.pvd')


#############################################################################
# Beginning of flow
start = int(particle.prm_dict['end time'] / (particle.prm_dict['output frequency'] * particle.prm_dict['time step']))

# Position of the outlet
hopper_outlet = 0

# Create a list to store the number of particles below the outlet
number_of_particles_below = []
mass_discharge = []
volume = 4./3. * np.pi * (particle.prm_dict['diameter']/2.)**3
print(f'Total mass {volume * particle.prm_dict["density particles"]*particle.prm_dict["number"]}')

# Create a list to store the "flow rate" of particles
rate = []

# Loop through all results
for i in range(len(particle.list_vtu)):

    # Store results in 'df'
    exec(f'df = particle.df_{i}')

    vertical_position =\
        df.points[:, 2]

    # Select the data
    vertical_position_below = vertical_position[vertical_position < hopper_outlet]

    # Number of particles above
    number_of_particles_below.append(len(vertical_position_below))
    mass_discharge.append(number_of_particles_below[-1] * particle.prm_dict['density particles'] * volume)

    if i > 0:
        rate.append((number_of_particles_below[i] - number_of_particles_below[i - 1])/particle.prm_dict['time step'])
    
tab = pd.concat([pd.Series(particle.time_list[1:]), pd.Series(rate), pd.Series(number_of_particles_below),
                 pd.Series(mass_discharge)], axis = 1)
tab.columns = ['time', 'rate', 'number_of_particles', 'mass_discharge']
tab.to_csv(save_path + '/results_.csv')


paper_data = pd.read_csv('paper_data.csv', ',')

plt.plot([value - particle.time_list[start] for value in particle.time_list[start:]],
         [value * 1000 for value in mass_discharge[start:]], label="Lethe DEM results")
plt.plot(paper_data['time'], paper_data['discharge'], label="Anand et al. results")


p = np.polyfit([value - particle.time_list[start+int(0.5/(particle.prm_dict['output frequency'] * particle.prm_dict['time step']))]
                for value in particle.time_list[start+int(0.5/(particle.prm_dict['output frequency'] * particle.prm_dict['time step'])):
                                                start+int(0.75/(particle.prm_dict['output frequency'] * particle.prm_dict['time step']))]],
               [value * 1000 for value in mass_discharge[start+int(0.25/(particle.prm_dict['output frequency'] * particle.prm_dict['time step'])):
                                                         start+int(0.5/(particle.prm_dict['output frequency'] * particle.prm_dict['time step']))]], 1)
print(p)
plt.xlabel('Time (s)')
plt.ylabel('Mass discharged from the hopper (g)')
plt.legend()
plt.grid()
plt.show()

