"""
Postprocessing code for the 3D dam-break with obstacle example
This code extracts the height of fluid1 and compares them to
the experimental results of Maritime Research Institute
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

#Set phase_limit to search for height values
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
for i in range(1,5):
    exec(f'H{i} = []')

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

    for i in range(1,5):
        #Get heights
        exec(f'sampled_H{i}_data = sim.sample_over_line(H{i}_a, H{i}_b, resolution=1000)')
        exec(f'phase_H{i} = pd.DataFrame(sampled_H{i}_data["phase"])')
        exec(f'points_H{i} = pd.DataFrame(sampled_H{i}_data.points)')

        #Find highest fluid1 points
        exec(f'H{i}_fluid1_points = points_H{i}[phase_H{i}[0] > phase_limit].values')
        exec(f'if len(H{i}_fluid1_points) > 0: \n'
             f'\t H{i}_max = H{i}_fluid1_points[0][2]\n'
             f'\t for point in H{i}_fluid1_points: \n'
             f'\t \t if point[2] > H{i}_max:\n'
             f'\t \t \t H{i}_max = point[2]\n'
             f'else:\n'
             f'\t H{i}_max = 0\n'
             f'H{i}.append(H{i}_max)')

#Figures
tf = time_list[len(time_list)-1]
plt.rcParams['font.size'] = '16'
fs1 = 16
fs2 = 14

#Heights
fig0 = plt.figure(figsize=(24.5, 15.5))
for i in range(1,5):
    exec(f'ax{i-1} = fig0.add_subplot(22{i})')
    exec(f'ax{i-1}.plot(time_list, H{i}, "-k", linewidth=1.5, label="H{i} - Lethe")')
    exec(f'ax{i-1}.plot(exp_time, exp_H{i}, "--r",  linewidth=2,  label="H{i} - MARIN")')
    exec(f'ax{i-1}.set_ylabel(r"H{i} [m]", fontsize=fs1)')
    exec(f'ax{i-1}.set_xlabel(r"t [s]", fontsize=fs1)')
    exec(f'ax{i-1}.set_xlim([0, tf])')
    exec(f'ax{i-1}.legend(loc="best")')
    exec(f'for label in (ax{i-1}.get_xticklabels() + ax{i-1}.get_yticklabels()):\n'
         f'\t label.set_fontsize(fs2)')
fig0.savefig('./H1_to_H4_evolution_2.png')
plt.show()