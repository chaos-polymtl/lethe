# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np
import pyvista as pv
import pandas as pd

from matplotlib import pyplot as plt

# Set the matplotlib style

def define_plot_style():
    from cycler import cycler

    colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']

    plt.rcParams['axes.prop_cycle'] = cycler(color = colors)
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['figure.figsize'] = (10,8)
    plt.rcParams['figure.autolayout'] = True
    plt.rcParams['lines.linewidth'] = 5
    plt.rcParams['lines.markersize'] = '11'
    plt.rcParams['markers.fillstyle'] = "none"
    plt.rcParams['lines.markeredgewidth'] = 2
    plt.rcParams['legend.columnspacing'] = 2
    plt.rcParams['legend.handlelength'] = 3
    plt.rcParams['legend.handletextpad'] = 0.2
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.fancybox'] = False
    plt.rcParams['xtick.major.width'] = 2
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 2
    plt.rcParams['font.size'] = '25'
    plt.rcParams['font.family']='DejaVu Serif'
    plt.rcParams['font.serif']='cm'
    plt.rcParams['legend.fontsize']='20'
    plt.rcParams['savefig.bbox']='tight'
    plt.rcParams['legend.handlelength']=1

# Definition of functions
def fg(Mp, g):
    return Mp * g

def fb(Vp, rhof, g):
    return -Vp * rhof * g

def dallavalle(dp, rhof, mu, urel):
    Rep = np.max([0.0000001, np.abs(urel) * rhof * dp / mu])
    beta = np.pi * dp**2 / 8 * rhof * (0.63 + 4.8/(Rep**(1/2)))**2 * np.abs(urel)
    fd = beta * urel
    return fd

def get_time_and_velocity_from_simulation(path_to_pvd):
    reader = pv.PVDReader(path_to_pvd)

    time = reader.time_values

    velocity = np.zeros(len(reader.time_values))
    for i in range(len(reader.time_values)):
        reader.set_active_time_point(i)
        data = reader.read()
        velocity[i] = data[0]['velocity'][:, 1][0]

    return time, velocity

def get_time_position_from_simulation(path_to_pvd):
    reader = pv.PVDReader(path_to_pvd)

    time = reader.time_values

    position = np.zeros(len(reader.time_values))
    for i in range(len(reader.time_values)):
        reader.set_active_time_point(i)
        data = reader.read()
        position[i] = data[0].points[:, 1]

    return time, position

def compare_to_analytical(analytical_velocity, simulated_velocity):
    # Interpolate analytical velocity to the same time points as the simulation
    # First, check what is the ratio between the number of elements in the analytical and simulated velocity
    ratio = len(analytical_velocity) / len(simulated_velocity)
    
    # Get analytical velocity at the same time points as the simulation
    analytical_velocity_interpolated = np.zeros(len(simulated_velocity))
    for i in range(len(simulated_velocity)):
        analytical_velocity_interpolated[i] = analytical_velocity[int(i * ratio)]
    
    # Calculate the RMSE
    rmse = np.sqrt(np.mean((simulated_velocity - analytical_velocity_interpolated)**2))

    return rmse

def get_center_void_fraction_from_simulation(path_to_pvd, time_point =-1):
    reader = pv.PVDReader(path_to_pvd)

    time = reader.time_values

    reader.set_active_time_value(time[time_point])
    data = reader.read()
    line_lower_bound = [0, data.bounds[2], 0]
    line_upper_bound = [0, data.bounds[3], 0]
    sample = data[0].sample_over_line(line_lower_bound, line_upper_bound, resolution=1000)
    sample.points = sample.points

    return time[time_point], sample
