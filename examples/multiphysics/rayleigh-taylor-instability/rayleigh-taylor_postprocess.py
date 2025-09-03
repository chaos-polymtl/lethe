# SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for rayleigh-taylor-instability example
This code extracts the y position of the bubble and the spike and compares it
to the results of He et al (1999)
"""

#-------------------------------------------
# Modules
#-------------------------------------------
import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import pyvista as pv
import tqdm

import os
import sys
sys.path.append("$LETHE_PATH/contrib/postprocessing/")
from lethe_pyvista_tools import *


#--------------------------------------------
# Plot parameters
#--------------------------------------------
from cycler import cycler

colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']

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

#--------------------------------------------
# Main
#--------------------------------------------

#To make it work, type "python3 rayleigh-taylor_postprocess.py -f . -p rayleigh-taylor-instability.prm" into the terminal.
parser = argparse.ArgumentParser(description='Arguments for the validation of the 2D rising bubble benchmark')
parser.add_argument("-f", "--folder", type=str, help="Path to the folder in which the simulation is run. This is the folder that contains the prm file.", required=True)
parser.add_argument("-p", "--prm", type=str, help="Name of the prm file", required=True)

args, leftovers=parser.parse_known_args()


#Load reference data from He et al (1999)
ref_data_file = "ref_He_et_al_data.txt"
ref_data = np.loadtxt(ref_data_file,skiprows=1)
ref_time = ref_data[:,0]
ref_y_bubble = ref_data[:,1]
ref_y_spike = ref_data[:,2]

#Constants
H = 0.25
g = 9.81

#Set phase_limit to search for y values
phase_limit = 0.5

#Take case path as argument and store it
simulation_path = args.folder
prm_file_name=args.prm
pvd_file_name="out.pvd"

#Define list of VTU files
sim = lethe_pyvista_tools(simulation_path, prm_file_name, pvd_file_name)

# Get active times
time_list = sim.time_list

time_list = [x * ((g/H)**0.5) for x in time_list]

#Create a list to fill with maximum y in which phase < phase_limit
y_bubble_list = []

#Create a list to fill with minimum y in which phase > phase_limit
y_spike_list = []

#Create beginning and end points for spike line
a_spike = [0.125, 0, 0]
b_spike = [0.125, 1, 0]

#Create beginning and end points for bubble line
a_bubble = [0, 0, 0]
b_bubble = [0, 1, 0]

#Read VTU files
pbar = tqdm(total = len(sim.list_vtu), desc="Calculate Spike and Bubble")
for i in tqdm(range(len(sim.list_vtu))):

    step = sim.get_df(i)

    # Extract the spike and the bubble
    sampled_data_spike = step.sample_over_line(a_spike, b_spike, resolution=1000)
    phase_spike = pd.DataFrame((sampled_data_spike["filtered_phase"].copy()))
    points_spike = pd.DataFrame(sampled_data_spike.points.copy())

    sampled_data_bubble = step.sample_over_line(a_bubble, b_bubble, resolution=1000)
    phase_bubble = pd.DataFrame(sampled_data_bubble["filtered_phase"].copy())
    points_bubble = pd.DataFrame(sampled_data_bubble.points.copy())

    #Find min 'y' in phase > phase_limit (SPIKE)
    fluid1_points = points_spike[phase_spike[0] > phase_limit].values
    y_min = fluid1_points[0][1]
    for points in fluid1_points:
        if points[1] < y_min:
            y_min = points[1]
    y_spike_list.append(y_min)

    #Find max 'y' in phase < phase_limit (BUBBLE)
    fluid0_points = points_bubble[phase_bubble[0] < phase_limit].values
    y_max = fluid0_points[0][1]
    for points in fluid0_points:
        if points[1] > y_max:
            y_max = points[1]
    y_bubble_list.append(y_max)



#Figure
fig0 = plt.figure()
ax0 = fig0.add_subplot(111)

ax0.plot(time_list, y_bubble_list,  label="Bubble - Lethe")
line, = ax0.plot(ref_time, ref_y_bubble,linestyle="dotted", label="Bubble - He et al (1999)")
line.set_path_effects([pe.Stroke(linewidth=5.5, foreground='black'),
                       pe.Normal()])

ax0.plot(time_list, y_spike_list,  label="Spike - Lethe")
line, = ax0.plot(ref_time, ref_y_spike, linestyle="dotted",   label=r"Spike - He \textit{et al.} (1999)")
line.set_path_effects([pe.Stroke(linewidth=5.5, foreground='black'),
                       pe.Normal()])

ax0.set_ylabel(r'$y$')
ax0.set_xlabel(r'$t^* = t \sqrt{g/H}$')
ax0.set_xlim([0, 4.5])
ax0.set_ylim([0., 1.2])
ax0.legend()
plt.tight_layout()
fig0.savefig('./spike_and_bubble_evolution.png')
plt.show()


#Plot the mass of fluid throughout the simulation

t,area,geo_area = np.loadtxt(sim.path_output+"mass_conservation_information.dat", usecols=(0,7,8),skiprows=1,unpack=True)
plt.plot(t,area/area[0],label="Algebraic measure")
plt.plot(t,geo_area/geo_area[0],label="Geometric measure")
plt.legend()

plt.xlabel(r"Time $(t)$ [s]")
plt.ylabel(r"Relative area $\left(\frac{A}{A_0}\right)$")
plt.tight_layout()
plt.savefig('./mass_of_fluid_1.png')
plt.show()
