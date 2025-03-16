# SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from sklearn.metrics import r2_score

from lethe_pyvista_tools import *

parser = argparse.ArgumentParser(description='Arguments for calculation of the velocity profile in the rotating drum')
parser.add_argument("-f", "--folder", type=str, help="Folder path", required=True)
parser.add_argument("--prm", type=str, help="prm file", required=True)
parser.add_argument("--automatic", action="store_true", default=False, required=False)
args, leftovers=parser.parse_known_args()

folder=args.folder

# Starting vtu id. Here, we are interested in the last 100 vtu
start = 400

# Number of sampled particles
n_particle_sample=10

# Distance to wall, used to extract the particles
cut_off = 0.07

# Radius of the cylinder
R=0.069

# Sample point in x
n_y_sample = 8

# Load lethe data
pvd_name = 'out.pvd'
ignore_data = ['type', 'volumetric contribution', 'torque', 'fem_torque',
               'fem_force']
particle = lethe_pyvista_tools(folder, args.prm,
                               pvd_name, ignore_data=ignore_data)
time = particle.time_list


# Points where the local averages are made in the depth of the cylinder
y_sample = np.linspace(-0.5*R,0.5*R,n_y_sample)

# Values where the height are stored
z_sample = np.zeros(n_y_sample)

# Number of sampled files
n_files = 0

D = 1. * 0.003 # D is the distance in which particle are considered for the local average

for i in range(start, len(time)):
    df_load = particle.get_df(i)

    # We do a scalar product to find the velocity in the new frame of reference.
    df = pd.DataFrame(np.copy(df_load.points), columns=['x', 'y', 'z'])

    # We take the particles close to the walls
    condition = (df['x'] > cut_off) | (df['x'] < -cut_off)
    df = df[condition]
    df_filtered = df.copy()

    for index, y_loc in enumerate(y_sample):
        # Compute the distance between each particle and the location of the x sampling point
        df_filtered['dist'] = ((df_filtered['y']-y_loc)** 2.) ** 0.5

        df_to_sample = df_filtered[df_filtered['dist'] < D]

        # Take the n_sample highest particles
        df_sampled = df_to_sample.nlargest(n_particle_sample, 'z')

        z_sample[index] += df_sampled['z'].mean()

    n_files += 1



z_sample = z_sample / n_files
p = np.polyfit(y_sample, z_sample,1)



R2 = r2_score(z_sample,np.polyval(p,y_sample))
print("R2: ", R2)
print("Slope: ",p[0])
print("Angle: ",np.arctan(p[0])*180/np.pi)

with open("post-process.txt", "w") as file:
    file.write("Angle\n")
    file.write(f"{np.arctan(p[0])*180/np.pi}\n")

if (not args.automatic):
    plt.plot(y_sample, z_sample,'s',label='Sampled points')
    plt.plot(y_sample, np.polyval(p,y_sample),'--',label='Linear fit')
    plt.ylabel("Height [m]")
    plt.xlabel("x [m]")
    plt.legend()
    plt.savefig("angle_calculation.png",dpi=300)
    plt.show()







