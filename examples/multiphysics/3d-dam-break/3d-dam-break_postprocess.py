"""
Postprocessing code for the 3D dam-break with obstacle example
This code extracts the height of fluid1 and compares them to
the experimental results of MAritime Research Institute
Netherlands (MARIN)
"""

#-------------------------------------------
# Modules
#-------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv

import os
import sys
import csv

#--------------------------------------------
# Main
#--------------------------------------------

#Load reference data from MARIN
exp_data_file = "experimental_data.txt"
exp_data = np.loadtxt(exp_data_file, skiprows=1)
exp_time = exp_data[:,0]
exp_H3 = exp_data[:,1]
exp_H2 = exp_data[:,2]
exp_H1 = exp_data[:,3]
exp_H4 = exp_data[:,4]

#Set phase_limit to search for y values
phase_limit = 0.5

#Take case path as argument and store it
output_path = sys.argv[1]

#Define list of VTU files
list_vtu = os.listdir(output_path)
list_vtu = [x for x in list_vtu if ("vtu" in x and "pvtu" not in x and "boundaries" not in x)]

# Sort VTU files to ensure they are in the same order as the time step
list_vtu = sorted(list_vtu)

# Read the pvd file to extract the times
reader = pv.get_reader(f"{output_path}/3d-dam-break.pvd")

# Get active times
time_list = reader.time_values

#Create lists to fill with fluid 1 height values
H1 = []
H2 = []
H3 = []
H4 = []

#Create beginning and end points for height lines
H1_a = [0.496, 0, 0]
H1_b = [0.496, 0, 1]
H2_a = [0.992, 0, 0]
H2_b = [0.992, 0, 1]
H3_a = [1.488, 0, 0]
H3_b = [1.488, 0, 1]
H4_a = [2.638, 0, 0]
H4_b = [2.638, 0, 1]

#Read VTU files
for vtu_file in list_vtu:
    sim = pv.read(f"{output_path}/{vtu_file}")

    #Get heights
    r = 1000
    sampled_H1_data = sim.sample_over_line(H1_a, H1_b, resolution=r)
    phase_H1 = pd.DataFrame(sampled_H1_data["phase"])
    points_H1 = pd.DataFrame(sampled_H1_data.points)
    sampled_H2_data = sim.sample_over_line(H2_a, H2_b, resolution=r)
    phase_H2 = pd.DataFrame(sampled_H2_data["phase"])
    points_H2 = pd.DataFrame(sampled_H2_data.points)
    sampled_H3_data = sim.sample_over_line(H3_a, H3_b, resolution=r)
    phase_H3 = pd.DataFrame(sampled_H3_data["phase"])
    points_H3 = pd.DataFrame(sampled_H3_data.points)
    sampled_H4_data = sim.sample_over_line(H4_a, H4_b, resolution=r)
    phase_H4 = pd.DataFrame(sampled_H4_data["phase"])
    points_H4 = pd.DataFrame(sampled_H4_data.points)

    #Find highest fluid1 points
    H1_fluid1_points = points_H1[phase_H1[0] > phase_limit].values
    if len(H1_fluid1_points) > 0:
        H1_max = H1_fluid1_points[0][2]
        for point in H1_fluid1_points:
            if point[2] > H1_max:
                H1_max = point[2]
    else:
        H1_max = 0
    H1.append(H1_max)

    H2_fluid1_points = points_H2[phase_H2[0] > phase_limit].values
    if len(H2_fluid1_points) > 0:
        H2_max = H2_fluid1_points[0][2]
        for point in H2_fluid1_points:
            if point[2] > H2_max:
                H2_max = point[2]
    else:
        H2_max = 0
    H2.append(H2_max)

    H3_fluid1_points = points_H3[phase_H3[0] > phase_limit].values
    if len(H3_fluid1_points) > 0:
        H3_max = H3_fluid1_points[0][2]
        for point in H3_fluid1_points:
            if point[2] > H3_max:
                H3_max = point[2]
    else:
        H3_max = 0
    H3.append(H3_max)

    H4_fluid1_points = points_H4[phase_H4[0] > phase_limit].values
    if len(H4_fluid1_points) > 0:
        H4_max = H4_fluid1_points[0][2]
        for point in H4_fluid1_points:
            if point[2] > H4_max:
                H4_max = point[2]
    else:
        H4_max = 0
    H4.append(H4_max)

#Figures
tf = time_list[len(time_list)-1]
lw = 1.5
lw_dashed = 2
plt.rcParams['font.size'] = '16'
fs1 = 16
fs2 = 14

#Heights
fig0 = plt.figure(figsize=(24.5, 15.5))
ax0 = fig0.add_subplot(221)
ax0.plot(time_list, H1, '-k', linewidth=lw, label="H1 - Lethe")
ax0.plot(exp_time, exp_H1, 'r', linestyle="dashed",  linewidth=lw_dashed,  label="H1 - MARIN")
ax0.set_ylabel(r'H1 [m]', fontsize=fs1)
ax0.set_xlabel(r't [s]', fontsize=fs1)
ax0.set_xlim([0, tf])
ax0.legend(loc="best")
for label in (ax0.get_xticklabels() + ax0.get_yticklabels()):
	label.set_fontsize(fs2)

ax1 = fig0.add_subplot(222)
ax1.plot(time_list, H2, '-k', linewidth=lw, label="H2 - Lethe")
ax1.plot(exp_time, exp_H2, 'r', linestyle="dashed",  linewidth=lw_dashed,  label="H2 - MARIN")
ax1.set_ylabel(r'H2 [m]', fontsize=fs1)
ax1.set_xlabel(r't [s]', fontsize=fs1)
ax1.set_xlim([0, tf])
ax1.legend(loc="best")
for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(fs2)

ax2 = fig0.add_subplot(223)
ax2.plot(time_list, H3, '-k', linewidth=lw, label="H3- Lethe")
ax2.plot(exp_time, exp_H3, 'r', linestyle="dashed",  linewidth=lw_dashed,  label="H3 - MARIN")
ax2.set_ylabel(r'H3 [m]', fontsize=fs1)
ax2.set_xlabel(r't [s]', fontsize=fs1)
ax2.set_xlim([0, tf])
ax2.legend(loc="best")
for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
	label.set_fontsize(fs2)

ax3 = fig0.add_subplot(224)
ax3.plot(time_list, H4, '-k', linewidth=lw, label="H4 - Lethe")
ax3.plot(exp_time, exp_H4, 'r', linestyle="dashed",  linewidth=lw_dashed,  label="H4 - MARIN")
ax3.set_ylabel(r'H4 [m]', fontsize=fs1)
ax3.set_xlabel(r't [s]', fontsize=fs1)
ax3.set_xlim([0, tf])
ax3.legend(loc="best")
for label in (ax3.get_xticklabels() + ax3.get_yticklabels()):
	label.set_fontsize(fs2)

fig0.savefig('./H1_to_H4_evolution.png')
plt.show()
