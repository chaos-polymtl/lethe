"""
Postprocessing code for rayleigh-taylor-instability example
This code extracts the y position of the bubble and the spike and compares it
to the results of He et al (1999)
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

#To make it work, type "python3 rayleigh-taylor_postprocess.py ./output/adaptive/" or
#"python3 rayleigh-taylor_postprocess.py ./output/constant /" into the terminal.


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
output_path = sys.argv[1]



#Define list of VTU files
list_vtu = os.listdir(output_path)
list_vtu = [x for x in list_vtu if ("vtu" in x and "pvtu" not in x)]

# Sort VTU files to ensure they are in the same order as the time step
list_vtu = sorted(list_vtu)

# Read the pvd file to extract the times
if "constant" in list_vtu[0]:
    reader = pv.get_reader("output/constant/rayleigh-taylor-constant-sharpening.pvd")
else:
    reader = pv.get_reader("output/adaptive/rayleigh-taylor-adaptive-sharpening.pvd")

# Get active times
time_list = reader.time_values
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
for vtu_file in list_vtu:
    sim = pv.read(f"{output_path}/{vtu_file}")

    sampled_data_spike = sim.sample_over_line(a_spike, b_spike, resolution=1000)
    phase_spike = pd.DataFrame(sampled_data_spike["filtered_phase"])
    points_spike = pd.DataFrame(sampled_data_spike.points)

    sampled_data_bubble = sim.sample_over_line(a_bubble, b_bubble, resolution=1000)
    phase_bubble = pd.DataFrame(sampled_data_bubble["filtered_phase"])
    points_bubble = pd.DataFrame(sampled_data_bubble.points)

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
ax0.plot(time_list, y_spike_list, '-k', linewidth=2, label="Spike - Lethe")
ax0.plot(ref_time, ref_y_spike, 'r', linestyle="dashed",  linewidth=2.5,  label="Spike - He et al (1999)")
ax0.plot(time_list, y_bubble_list, color='#2e62ff', linewidth=2, label="Bubble - Lethe")
ax0.plot(ref_time, ref_y_bubble, color='#fc8c03',  linewidth=2.5, linestyle="dashed", label="Bubble - He et al (1999)")
ax0.set_ylabel(r'$y$')
ax0.set_xlabel(r'$t^* = t \sqrt{g/H}$')
ax0.set_xlim([0, 4.5])
ax0.set_ylim([0.1, 0.8])
ax0.legend(loc="upper left")
plt.title("Spike and bubble evolution with {} interface sharpening".format(output_path[9:-1]))
fig0.savefig('./spike_and_bubble_evolution_{}.png'.format(output_path[9:-1]))
plt.show()


#Plot the mass of fluid along the simulation

#Functions to turn .dat data in numpy array
#Credits to Lucka Barbeau for the two functions below !
def is_number(s):
    try:
        float(s)
    except ValueError:
        return False
    return True
def read_my_data(results_path):
    force_file = open(results_path, 'r')
    list_of_list_of_vars_name = [[]];
    list_of_list_of_vars = [[]];
    fix_vars_name = False

    nb_set_of_vars = 0;
    for line in force_file.readlines():
        i = 0;
        for word in line.split():
            if is_number(word):
                fix_vars_name = True
                list_of_list_of_vars[nb_set_of_vars][i] = np.append(list_of_list_of_vars[nb_set_of_vars][i],
                                                                    float(word))
                i += 1
            else:
                if word != 'particle':
                    if fix_vars_name:
                        nb_set_of_vars += 1
                        list_of_list_of_vars_name.append([])
                        list_of_list_of_vars.append([])
                        fix_vars_name = False
                    list_of_list_of_vars_name[nb_set_of_vars].append(word)
                    list_of_list_of_vars[nb_set_of_vars].append(np.array([]))
    return list_of_list_of_vars_name, list_of_list_of_vars

list_of_list_of_vars_name,list_of_list_of_vars=read_my_data(output_path + "VOF_monitoring_fluid_1.dat")


plt.plot(list_of_list_of_vars[0][0],list_of_list_of_vars[0][2])
plt.title("Evolution of the mass of fluid 1 with {} interface sharpening".format(output_path[9:-1]))
plt.xlabel("Time (s)")
plt.ylabel("Mass of fluid 1 (kg)")
plt.ylim(37.4,37.6)
plt.savefig('./mass_of_fluid_1_{}.png'.format(output_path[9:-1]))
plt.show()


