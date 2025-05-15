# SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

from lethe_pyvista_tools import *
import matplotlib.pyplot as plt
from pyvista_utilities import *
from log_utilities import *

#############################
# Path and filenames
#############################
simulation_path = "./"
prm_filename = "pneumatic-conveying.prm"
pvd_name = "cfd_dem"
log_filename = "pneumatic-log.out"

#############################
# Read data and process
#############################
# In this post-processing, we need the ID of the particles, the position, and the velocity, the other data is not
# necessary and takes a lot of memory
ignore_data = ['omega', 'type', 'diameter', 'mass', 'fem_force', 'fem_torque']

# Create the particle object
particles = lethe_pyvista_tools(f'{simulation_path}', prm_filename, f'{pvd_name}_particles.pvd',
                                ignore_data=ignore_data, step=5)
prm_dict = particles.prm_dict
dataframe = []

for i in range(len(particles.list_vtu)):
    df = particles.get_df(i)
    data = pd.DataFrame()
    # Get the ID of the particles
    data['ID'] = df['ID'].astype(int)

    # Get the position of the particles
    data['x'] = df.points[:, 0]
    data['y'] = df.points[:, 1]
    data['z'] = df.points[:, 2]

    # Get the velocity of the particles
    data['u'] = df['velocity'][:, 0]
    data['v'] = df['velocity'][:, 1]
    data['w'] = df['velocity'][:, 2]

    dataframe.append(data)

#############################
# Post-processing
#############################
# We specify the volume of the triangulation since the mesh is not a perfect cylinder
# This value is at the beginning of the output of the simulation
volume_triangulation = 0.00540043
L = 1.0

# Create the post-processing and log object
postprocessing = PneumaticConveyingUtilities(dataframe, prm_dict, particles.time_list, L, volume_triangulation)
log = LogUtilities(f"{simulation_path}/{log_filename}")

# Get the slug length, the solid mass flow rate, the average particle velocity, and the slug velocity
postprocessing.calculate_slug_length(void_fraction_threshold=0.55)
postprocessing.calculate_solid_mass_flow_rate()
postprocessing.calculate_average_particle_velocity()
postprocessing.calculate_slug_velocity(sample_step=3)
postprocessing_dataframe = postprocessing.get_data()

# Get the fluid velocity and the correctional volumetric force (beta)
log.flow_monitoring()
log_dataframe = log.get_log_data()

#############################
# Print results
#############################
# Calculate the time-averaged values
start = 4.0 # ~quasi-steady state
quasisteasy_postprocessing_dataframe = postprocessing_dataframe[postprocessing_dataframe['time'] > start]
quasisteasy_log_dataframe = log_dataframe[log_dataframe['time'] > start]

print('\nTime-averaged in quasi-steady state values:')
print(f"Fluid velocity:             {quasisteasy_log_dataframe['fluid_velocity'].mean():2.3f} ± "
      f"{quasisteasy_log_dataframe['fluid_velocity'].std():2.3f} m/s")
print(f"Beta force:                 {quasisteasy_log_dataframe['beta'].mean():2.0f}  ± "
        f"{quasisteasy_log_dataframe['beta'].std():2.0f}    m/s²")
print(f"Slug velocity:              {quasisteasy_postprocessing_dataframe['slug_velocity'].mean():2.3f} ± "
        f"{quasisteasy_postprocessing_dataframe['slug_velocity'].std():2.3f} m/s")
print(f"Particles in slug velocity: {quasisteasy_postprocessing_dataframe['average_velocity'].mean():2.3f} ± "
        f"{quasisteasy_postprocessing_dataframe['average_velocity'].std():2.3f} m/s")
print(f"Solid mass flow rate:       {quasisteasy_postprocessing_dataframe['mass_flow_rate'].mean():2.2f}          kg/s")
print(f"Slug length:                {quasisteasy_postprocessing_dataframe['distance'].mean():2.2f}  ± "
        f"{quasisteasy_postprocessing_dataframe['distance'].std():2.2f}  m")

print('\nResult from u_s = 0.967 u_p + 0.5(gD)**0.5:')
u_p = quasisteasy_postprocessing_dataframe["average_velocity"].mean()
u_s = quasisteasy_postprocessing_dataframe["slug_velocity"].mean()
calc_u_s = 0.967*u_p+ 0.5 * (9.81 * 0.084) ** 0.5
print(f'u_s = {calc_u_s:2.3} m/s')
print(f'Difference: {abs(calc_u_s - u_s)/u_s:.3f} %')

#############################
# Plot results
#############################
# All plots in one figure
fig = plt.figure(figsize=(10, 9))

# Average fluid velocity and beta force
flow_plot = fig.add_subplot(311)
beta_plot = flow_plot.twinx()

flow_color = 'tab:orange'
flow_plot.plot(log_dataframe['time'], log_dataframe['fluid_velocity'], color=flow_color)
flow_plot.plot(log_dataframe['time'], 3 * np.ones(len(log_dataframe)),
               color=flow_color, linestyle='--', label='Target: 3 m/s')

beta_color = 'tab:purple'
beta_plot.plot(log_dataframe['time'], log_dataframe['beta'], color=beta_color)

# Slug velocity and particles in slug
slug_plot = fig.add_subplot(312)
particles_plot = slug_plot.twinx()

slug_color = 'k'
slug_plot.plot(postprocessing_dataframe['time'], postprocessing_dataframe['slug_velocity'], color=slug_color)

particle_color = 'tab:blue'
particles_plot.plot(postprocessing_dataframe['time'], postprocessing_dataframe['average_velocity'],
                    color=particle_color)

# Solid mass flow rate and average particle velocity
mass_plot = fig.add_subplot(313)
length_plot = mass_plot.twinx()

mass_color = 'tab:red'
mass_plot.plot(postprocessing_dataframe['time'], postprocessing_dataframe['mass_flow_rate'], color=mass_color)

length_color = 'tab:green'
length_plot.plot(postprocessing_dataframe['time'], postprocessing_dataframe['distance'], color=length_color)

# Formatting shenanigans
flow_plot.set(ylim=[1.5, 4.0], xlim=[0, 5])
flow_plot.set_ylabel('Average fluid velocity (m/s)', color=flow_color)
flow_plot.tick_params(axis='y', labelcolor=flow_color)
flow_plot.grid(color='k', linestyle='-', linewidth=0.4)
flow_plot.xaxis.set_major_locator(plt.MaxNLocator(11))
flow_plot.axvspan(0, start, color='gray', alpha=0.15, lw=0)
flow_plot.get_legend()

beta_plot.set(ylim=[0, 5000])
beta_plot.set_ylabel('Beta force (m/s²)', color=beta_color)
beta_plot.tick_params(axis='y', labelcolor=beta_color)
beta_plot.axes.xaxis.set_ticklabels([])

slug_plot.set(ylim=[0, 2.5], xlim=[0, 5])
slug_plot.set_ylabel('Slug velocity (m/s)', color=slug_color)
slug_plot.tick_params(axis='y', labelcolor=slug_color)
slug_plot.grid(color='k', linestyle='-', linewidth=0.4)
slug_plot.xaxis.set_major_locator(plt.MaxNLocator(11))
slug_plot.axvspan(0, start, color='gray', alpha=0.15, lw=0)

particles_plot.set(ylim=[0, 2.5], xlim=[0, 5])
particles_plot.set_ylabel('Particle velocity (m/s)', color=particle_color)
particles_plot.tick_params(axis='y', labelcolor=particle_color)
particles_plot.axes.xaxis.set_ticklabels([])

mass_plot.set(ylim=[0, 2.5], xlim=[0, 5])
mass_plot.set_xlabel('Time (s)')
mass_plot.set_ylabel('Solid mass flow rate (kg/s)', color=mass_color)
mass_plot.tick_params(axis='y', labelcolor=mass_color)
mass_plot.grid(color='k', linestyle='-', linewidth=0.4)
mass_plot.xaxis.set_major_locator(plt.MaxNLocator(11))
mass_plot.axvspan(0, start, color='gray', alpha=0.15, lw=0)

length_plot.set(ylim=[0.1, 0.6], xlim=[0, 5])
length_plot.set_ylabel('Slug length (m)', color=length_color)
length_plot.tick_params(axis='y', labelcolor=length_color)
length_plot.xaxis.set_major_locator(plt.MaxNLocator(11))
length_plot.set_xlabel('Time (s)')

fig.savefig(f'./pneumatic-conveying-data.png', bbox_inches='tight', dpi=300)
plt.tight_layout()
plt.show()
