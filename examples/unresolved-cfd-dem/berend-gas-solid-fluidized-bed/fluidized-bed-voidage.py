# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from scipy import signal
from cycler import cycler

# Set plot parameters
colors=['#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02']
plt.rcParams['axes.prop_cycle'] = cycler(color = colors)
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize'] = (10,8)
plt.rcParams['lines.linewidth'] = 4
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
plt.rcParams['savefig.bbox']='tight'
plt.rcParams['legend.handlelength']=1

######################################################################

import sys
sys.path.append("$LETHE_PATH/contrib/postprocessing/")
from lethe_pyvista_tools import *

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the 2d-sandpile DEM example')
parser.add_argument("-f", "--folder", type=str, help="Folder path. This folder is the folder which contains the .prm file.", required=True)
args, leftovers=parser.parse_known_args()

# Simulation folder
folder=args.folder

# Starting vtu id
start = 200
end   = 650

# Load lethe data
pvd_particles = 'out_particles.pvd'
prm_file = 'gas-solid-fluidized-bed.prm'
particles = lethe_pyvista_tools(folder, prm_file, pvd_particles)
time = np.array(particles.time_list)

# Get void fractions on the plane at 45mm above the floating wall
voidage = np.zeros(end-start)
n_sample = 200
x_sample = np.linspace(0,0.09,n_sample)
h = 0.045
particle_diameter = 0.001545
D = 1.0/4.0 * particle_diameter

for i in range(start, end):
    n_detected = 0
    df_load = particles.get_df(i)
    df = pd.DataFrame(np.copy(df_load.points), columns=['x', 'y','z'])

    for index, x_loc in enumerate(x_sample):

        df_filtered = df.copy()

        # Compute the distance between each particle and the sampling point
        df_filtered['dist_x'] = ((df_filtered['x']-x_loc)** 2) ** 0.5
        df_filtered['dist_z'] = ((df_filtered['z']-h)** 2) ** 0.5
        
        # Keep the particles on the the sampling point
        df_sampled = df_filtered[(df_filtered['dist_x'] < D) & (df_filtered['dist_z'] < D)]

        # If there are particles on the sampling point
        if (len(df_sampled)>0):
            n_detected +=1
    voidage[i-start] = 1 - n_detected / n_sample


# Plot voidage over time
reference_voidage = pd.read_csv('reference/voidage.csv')
plt.figure()
plt.plot(time[start:end],voidage, label='Lethe')
plt.plot(reference_voidage['t'],reference_voidage['voidage'], label='Ref')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Voidage')
plt.xlim(2,4)
plt.ylim(0.1,0.7)
plt.grid()
plt.savefig('voidage-fluctuations')
plt.show()