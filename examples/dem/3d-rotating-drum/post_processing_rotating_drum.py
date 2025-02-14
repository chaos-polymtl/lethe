# SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from lethe_pyvista_tools import *
from numpy.ma.core import arctan
from numpy.ma.extras import average

parser = argparse.ArgumentParser(description='Arguments for calculation of the velocity profile in the rotating drum')
parser.add_argument("-f", "--folder", type=str, help="Folder path", required=True)
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
args, leftovers=parser.parse_known_args()

folder=args.folder

# Experimental date for depth as a function of velocity
depth_exp =  np.loadtxt(folder + "/" +"depth_velocity.dat", skiprows=1)

# Experimental date extracted using Web plot digitizer
free_surface_exp = np.loadtxt(folder + "/" +"free_surface_velocity.dat", skiprows=1)

# Starting vtu id. Here, we are interested in the last 100 vtu
start = 800

# Load lethe data
pvd_name = 'out.pvd'
ignore_data = ['type', 'volumetric contribution', 'torque', 'fem_torque',
               'fem_force']
particle = lethe_pyvista_tools("./", "rotating-drum.prm",
                               pvd_name, ignore_data=ignore_data)
time = particle.time_list

# We fix the angle at 27 deg like is was mentioned in the experimental article
angle = 27. * 2. * np.pi / 360.

# H is the normal distance between the free surface and the center of the cylinder
H = -0.03
R=0.12

# Points where the local averages are made in the depth of the cylinder
y_depth = np.linspace(-R, H, 100) # -0.12 is the cylinder radius
depth_vel_value_x = np.zeros_like(y_depth, float)
depth_vel_value_y = np.zeros_like(y_depth, float)
x_depth = 0.
D = 1. * 0.003 # D is the distance in which particle are considered for the local average


# Point where the local averages are made on the free surface of the cylinder
x_free_surface = np.linspace (-0.95*R, 0.95*R,100)
vel_value_x_free_surface = np.zeros_like(x_free_surface, float)
vel_value_y_free_surface = np.zeros_like(x_free_surface, float)
y_free_surface = H


# Change frame of reference. (loc is for local)
# The X direction is parallel to the free surface and pointing down.
x_loc_unit = np.array([-np.cos(angle), -np.sin(angle)])
y_loc_unit = np.array([-np.sin(angle), np.cos(angle)])

for i in range(start, len(time)):
    df_load = particle.get_df(i)

    # We do a scalar product to find the velocity in the new frame of reference.
    df = pd.DataFrame(np.copy(df_load.points), columns=['x', 'y', 'z'])
    df['v_x_loc'] = pd.DataFrame(df_load['velocity'][:, 1]) * x_loc_unit[
        0] + pd.DataFrame(df_load['velocity'][:, 2]) * x_loc_unit[1]

    # Same here, scalar product
    df['v_y_loc'] = pd.DataFrame(df_load['velocity'][:, 1]) * y_loc_unit[
        0] + pd.DataFrame(df_load['velocity'][:, 2]) * y_loc_unit[1]


    # The cylinder is 0.36 m long from -0.18 to 0.18
    cut_off = 0.175
    condition = (df['x'] < cut_off) & (df['x'] > -cut_off)
    df = df[condition]

    df['x_loc'] = (df['y'] ) * x_loc_unit[0] + (df['z']) * \
                  x_loc_unit[1]
    df['y_loc'] = (df['y'] ) * y_loc_unit[0] + (df['z']) * \
                  y_loc_unit[1]

    df_filtered = df.copy()

    for index, j in enumerate(y_depth):
        # Compute the distance between each particle and the local average calculation point.
        df_filtered['dist'] = (df_filtered['x_loc'] ** 2. + (df_filtered['y_loc'] - j) ** 2.) ** 0.5
        condition_2 = (df_filtered['dist'] < D)
        df_filtered_filtered = df_filtered[condition_2]

        # Number of particle in the local average zone
        n = len(df_filtered_filtered)

        # Average velocity in the local average
        tot_vel_x = df_filtered_filtered['v_x_loc'].sum()
        depth_vel_value_x[index] += tot_vel_x / n
        tot_vel_y = df_filtered_filtered['v_y_loc'].sum()
        depth_vel_value_y[index] += tot_vel_y / n

    for index, xx in enumerate(x_free_surface):
        # Compute the distance between each particle and the local average calculation point.
        df_filtered['dist_free_surface'] = ((df_filtered['x_loc']-xx) ** 2. + (df_filtered['y_loc']-y_free_surface) ** 2.) ** 0.5
        condition_2 = (df_filtered['dist_free_surface'] < D)
        df_filtered_filtered = df_filtered[condition_2]

        # Number of particle in the local average zone
        n = len(df_filtered_filtered)


        # Average velocity in the local average
        tot_vel_x = df_filtered_filtered['v_x_loc'].sum()
        vel_value_x_free_surface[index] += tot_vel_x / n
        tot_vel_y = df_filtered_filtered['v_y_loc'].sum()
        vel_value_y_free_surface[index] += tot_vel_y / n

# Average over time of the averages
depth_vel_value_x = depth_vel_value_x / (len(time) - start)
depth_vel_value_y = depth_vel_value_y / (len(time) - start)
vel_value_x_free_surface = vel_value_x_free_surface / (len(time) - start)
vel_value_y_free_surface = vel_value_y_free_surface / (len(time) - start)

plt.plot(depth_vel_value_x, y_depth, label="Lethe")
plt.plot(-depth_exp[:,0],
         depth_exp[:,1], "xk", label="RPT")
plt.xlabel("u (m/s)")
plt.ylabel("y (m)")
plt.plot(y_depth * 1.214749159, y_depth, "-.", label=r"$\omega r$")
plt.legend()
if (args.validate):     
    plt.savefig("lethe-rotating-drum-comparison-depth.pdf")
else:
    plt.savefig("lethe-rotating-drum-comparison-depth.png",dpi=300)
    plt.show()

plt.plot(x_free_surface/R, vel_value_x_free_surface/np.max(vel_value_x_free_surface),  label="Lethe")
plt.plot(free_surface_exp[:,0],
         free_surface_exp[:,1], "xk", label="RPT")
plt.xlabel("${x}/{R}$")
plt.ylabel("${v_x}/{v_{max}}$ (m)")
plt.ylim(-0.01,1.1)
plt.xlim(-1.01,1.01)
plt.legend()
if (args.validate):     
    plt.savefig("lethe-rotating-drum-comparison-free-surface.pdf")
else:
    plt.savefig("lethe-rotating-drum-comparison-free-surface.png",dpi=300)
    plt.show()




