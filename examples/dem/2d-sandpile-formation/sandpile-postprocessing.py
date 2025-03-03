# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from sklearn.metrics import r2_score

# Set plot parameters
plt.rcParams['figure.figsize'] = (10,8)
plt.rcParams['lines.markersize'] = '11'
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.fancybox'] = False
plt.rcParams['font.size'] = 20
plt.rcParams['font.family']='DejaVu Serif'
plt.rcParams['legend.handlelength']=1.5
plt.rcParams['lines.linewidth'] = 2

######################################################################


import sys
sys.path.append("$LETHE_PATH/contrib/postprocessing/")
from lethe_pyvista_tools import *

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the 2d-sandpile DEM example')
parser.add_argument("-f", "--folder", type=str, help="Folder path. This folder is the folder which contains the .prm file.", required=True)
parser.add_argument("--prm", type=str, help="prm file", required=True)
parser.add_argument("--regression", action="store_true", default=False, help="Plot the least squares regression",required=False)
parser.add_argument("--rollingmethod", type=str, choices=["constant", "viscous", "epsd"], help="Rolling resistance method. Must be one of: constant, viscous, or epsd.", required=True)
args, leftovers=parser.parse_known_args()

# Simulation folder
folder=args.folder

# Rolling resistance method
rollingmethod = args.rollingmethod

# Starting vtu id for the height of the pile
start = 500

# Number of sampled particles for the angle
n_angle_sample = 8

# Height of the bottom part of the mesh
h = 0.57

# Load lethe data
pvd_name = 'out.pvd'
ignore_data = ['type', 'volumetric contribution', 'torque', 'fem_torque',
               'fem_force']
particle = lethe_pyvista_tools(folder, args.prm,pvd_name, ignore_data=ignore_data)
time = particle.time_list


# Sampling points for the angle of repose
x_sample = np.linspace(-0.5,-0.1,n_angle_sample)

# Values where the heights are stored for the angle of repose
y_sample = np.zeros(n_angle_sample)

# Values where the height of the pile is stored
height = np.zeros(len(time)-start)

D = 0.0089 # D is the distance in which particle are considered around the sampling points (diameter of the biggest particle)

for i in range(start, len(time)):

    df_load = particle.get_df(i)
    df = pd.DataFrame(np.copy(df_load.points), columns=['x', 'y','z'])

    height[i-start] += df['y'].max() + h

    if i==len(time)-1 :

        for index, x_loc in enumerate(x_sample):

            df_filtered = df.copy()

            # Compute the distance between each particle and the sampling point
            df_filtered['dist'] = ((df_filtered['x']-x_loc)** 2) ** 0.5
            
            # Keep the particles close to the sampled point
            df_to_sample = df_filtered[df_filtered['dist'] < D]

            # Take the highest particle around the sampled one
            df_sampled = df_to_sample.nlargest(1, 'y')
            y_sample[index] += df_sampled['y']


# Results of the regression
p = np.polyfit(x_sample, y_sample,1)

R2 = r2_score(y_sample,np.polyval(p,x_sample))
print("R2: ", R2)
print("Slope: ",p[0])
print("Angle: ",np.arctan(p[0])*180/np.pi)

if (args.regression):
    # Plot the least squares regression used to calculate the angle
    plt.figure()
    plt.plot(x_sample, y_sample,'s',label='Sampled points')
    plt.plot(x_sample, np.polyval(p,x_sample),'--',label='Linear fit')
    plt.legend()
    plt.show()


# Write the angle on a file
with open(folder+"/data/angle_" + rollingmethod + ".txt", "w") as file:
    file.write("Angle\n")
    file.write(f"{np.arctan(p[0])*180/np.pi}\n")


# Read height data from paper
paper_data = pd.read_csv('data/reference/extraction_model_' + rollingmethod + '.csv')


# Plot the evolution of the height of the pile
plt.figure()
plt.plot(time[start:],height, label= "Lethe-DEM " + rollingmethod)
plt.plot(paper_data['x'],paper_data['Curve1'], '--', label= "Ai2010 " + rollingmethod)
plt.legend()
plt.grid()
plt.title("Evolution of the height of the pile with model " + rollingmethod, pad=25)
plt.xlabel('Time (s)')
plt.ylabel('Height of the pile (m)')
plt.yticks(np.arange(0.1, 0.32, 0.02))
plt.xticks(np.arange(0, 60, 10))
plt.savefig(folder + '/data/figure-height-comparison-' + rollingmethod)
plt.show()


# Export simulation data to csv
data_height = pd.DataFrame({'time': time[start:], 'height':height})
data_height.to_csv(folder + '/data/height_' + rollingmethod + '.csv')







