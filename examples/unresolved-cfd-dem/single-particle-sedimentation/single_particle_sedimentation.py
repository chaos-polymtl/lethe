# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


# This script is used to compare the simulation results of a single particle sedimentation to the analytical solution.
# It will plot the results of the particle position and velocity as a function of time and compare it to the analytical solution given the particle and fluid pvd files paths.

import numpy as np
import pandas as pd
import pyvista as pv
from argparse import ArgumentParser

# Define argument for pvd file path of both particle and fluid
parser = ArgumentParser()
parser.add_argument("--particle_pvd", type=str, required=True)
parser.add_argument("--fluid_pvd", type=str, required=True)
args = parser.parse_args()

import matplotlib.pyplot as plt

import sys
sys.path.append("./")

from single_particle_sedimentation_funcs import *
define_plot_style()

# Parameters for the analytical solution
convert_m_to_cm = 100
# Data:
dp = 0.002663
Vp = 4/3 * np.pi * (0.5 * dp)**3
rhop = 1029.
Mp = rhop * Vp

# Constants:
nu = 0.0000008379
rhof = 996.7775
mu = nu * rhof
g = -9.81

# Parameters
dt = 0.05
t_end = 2
initial_position = 0.08
initial_urel = 0.

# Calculate the analytical solution
analytical_time = np.arange(0, t_end, dt)
analytical_velocity = np.zeros_like(analytical_time)
analytical_position = np.zeros_like(analytical_time)
analytical_position[0] = initial_position

i = 1
current_t = dt
while current_t <= t_end:
    # Calculate forces
    Fg = fg(Mp, g)
    Fb = fb(Vp, rhof, g)
    Fd = dallavalle(dp, rhof, mu, analytical_velocity[i-1])

    # Calculate acceleration
    a = (Fb - Fd + Fg) / Mp

    # Calculate velocity
    analytical_velocity[i] = analytical_velocity[i-1] + a * dt

    # Calculate position
    analytical_position[i] = analytical_position[i-1] + analytical_velocity[i] * dt

    # Calculate time
    current_t += dt
    analytical_time[i] = current_t

    i += 1

# Get simulated and analytical particle position
time, position = get_time_position_from_simuation(args.particle_pvd)

# Plot position as a function of time and compare to analytical solution
plt.plot(time, position*convert_m_to_cm, label="Simulation", lw=0, marker='o')
plt.plot(analytical_time, analytical_position*convert_m_to_cm, label="Analytical", lw=0, marker='s')
plt.xlabel("Time [s]")
plt.ylabel("Particle position [cm]")
plt.legend()
plt.show()

# Get simulated and analytical particle velocity
time, velocity = get_time_and_velocity_from_simulation(args.particle_pvd)

# Compare results to analytical from RMSE
rmse = compare_to_analytical(analytical_velocity, velocity)

# Plot velocity as a function of time and compare to analytical solution
plt.plot(time, -velocity*convert_m_to_cm, label=f"Simulation - RMSE = {convert_m_to_cm*rmse:.2f} cm/s", lw=0, marker='o')
plt.plot(analytical_time, -analytical_velocity*convert_m_to_cm, label="Analytical", lw=0, marker='s')
plt.xlabel("Time [s]")
plt.ylabel("Particle velocity [cm/s]")
plt.legend()
plt.show()

# Define points to find sampled void fraction
time_points_to_sample_voidfraction = [1, len(time)//2, -1]

# For time points
for i in time_points_to_sample_voidfraction:
    # Get center void fraction from simulation
    time, center_void_fraction = get_center_void_fraction_from_simulation(args.fluid_pvd, time_point=i)

    # Plot center void fraction as a function of position in y
    plt.plot(center_void_fraction.points[:, 1]*convert_m_to_cm, center_void_fraction['void_fraction'], label=f"{time:.2f} s")
plt.xlabel("Position in y [cm]")
plt.ylabel("Center void fraction")
plt.legend(title="Time")
plt.show()

