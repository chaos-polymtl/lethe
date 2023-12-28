#############################################################################
"""
Postprocessing code for melting cavity example

"""
#############################################################################

'''Parameters'''

#############################################################################
'''Importing Libraries'''
from math import pi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import cm
import pyvista as pv

import glob

import os
import sys

#############################################################################

#############################################################################
'''Functions'''

#############################################################################

#Constants
H=0.714
alpha =  0.040516842071415184 # Used for the time scaling
St = 0.041 # Used for the time scaling
time_scaling = St * alpha 
y_top	= 0.999 * H
y_mid	= 0.5 * H
y_bot	= 0.001 * H

#Set temperature to search for maximum x
temperature_limit = 29.75

#Take case path as argument and store it
output_path = sys.argv[1]

# Read the pvd file to extract the times
reader = pv.get_reader(f"{output_path}/melting.pvd")
# Get active times
time_list = reader.time_values

#Define list of VTU files
list_vtu = os.listdir(output_path)
list_vtu = [x for x in list_vtu if  ("vtu" in x and "pvtu" not in x) ]

# Sort VTU files to ensure they are in the same order as the time step
list_vtu = sorted(list_vtu)

#Create a list to fill with maximum x in which phase > phase_limit
x_list_top = []
x_list_mid = []
x_list_bot = []

#Create beginning and end points for spike line
a_top = [0, y_top, 0]
b_top = [1, y_top, 0]

a_mid = [0,y_mid,0]
b_mid = [1,y_mid,0]

a_bot = [0, y_bot, 0]
b_bot = [1, y_bot, 0]


#Read vtu data
for vtu_file in (list_vtu):
    sim = pv.read(f"{output_path}/{vtu_file}")
    print("Processing file: ", vtu_file)

    sampled_data_top = sim.sample_over_line(a_top, b_top, resolution=1000)
    temperature_top = pd.DataFrame(sampled_data_top["temperature"])
    points_top = pd.DataFrame(sampled_data_top.points)

    sampled_data_bot = sim.sample_over_line(a_bot, b_bot, resolution=1000)
    temperature_bot = pd.DataFrame(sampled_data_bot["temperature"])
    points_bot = pd.DataFrame(sampled_data_bot.points)

    sampled_data_mid = sim.sample_over_line(a_mid, b_mid, resolution=1000)
    temperature_mid = pd.DataFrame(sampled_data_mid["temperature"])
    points_mid = pd.DataFrame(sampled_data_mid.points)
    
    #Find min 'x' in temperature > temperature_limit (top)
    fluid1_points = points_top[temperature_top[0] < temperature_limit].values
    x_min = fluid1_points[0][0]
    for points in fluid1_points:
        if points[0] < x_min:
            x_min = points[0]
    x_list_top.append(x_min)

    #Find min 'x' in temperature > temperature_limit (mid)
    fluid1_points = points_mid[temperature_mid[0] < temperature_limit].values
    x_min = fluid1_points[0][0]
    for points in fluid1_points:
        if points[0] < x_min:
            x_min = points[0]
    x_list_mid.append(x_min)

    #Find min 'x' in temperature > temperature_limit (bot)
    fluid1_points = points_bot[temperature_bot[0] < temperature_limit].values
    x_min = fluid1_points[0][0]
    for points in fluid1_points:
        if points[0] < x_min:
            x_min = points[0]
    x_list_bot.append(x_min)


#Experimental data from Gau and Viskanta 1986
top_x_exp = [0.00853, 0.01279, 0.01931, 0.02575, 0.03428, 0.04298, 0.0536, 0.06438]
top_y_exp = [0.11508, 0.18212, 0.28603, 0.31508, 0.38994, 0.50391, 0.61117, 0.72514]

center_x_exp = [0.00836, 0.01279, 0.01931, 0.02584, 0.0342, 0.04298, 0.0536, 0.06455]
center_y_exp = [0.11508, 0.16313, 0.20894, 0.23911, 0.29609, 0.36089, 0.45475, 0.52961]

bottom_x_exp = [0.00853, 0.01271, 0.01931, 0.02575, 0.0342, 0.04306, 0.0536, 0.06455]
bottom_y_exp = [0.1095, 0.12849, 0.14078, 0.16425, 0.21117, 0.22346, 0.26592, 0.26592]

melted_x_exp = [0.00946, 0.01433, 0.02113, 0.02877, 0.03774, 0.04739, 0.05956, 0.07164, 0.08104, 0.09069]
melted_y_exp = [0.11972, 0.16112, 0.20476, 0.22713, 0.2786, 0.34014, 0.43077, 0.5035, 0.56392, 0.61664]

#Simulation data by Blais and Ilinca 2018
melted_x_blais=[0.0009408610966497919,0.0026209691082628825,0.004569890217872624 ,0.0068548364163561895,0.011155915715326668,	
                0.018010754310777356,0.024731181999040716,0.034005377700162306,0.03891129082781425,0.043682793048278865,0.04811827889624767 ,0.0532258062055861,0.02963709512669266 ,0.05759408877905576 , 		
                0.06451613064900559,0.07432795690430948,0.08266128497149776,0.08931451374345097,0.09596774251540417,0.09899193344975654 ]	
melted_y_blais= [0.0372093373985142, 0.0669768073173256, 0.0855814156911383, 0.1079070181302468, 0.1423255798170668, 0.1897674397560892, 0.2344186446343063, 0.2874418750203834, 0.3162790662601486, 0.3441860391463121,0.37209301203247547,0.39999998491863886,0.2613953388414236,	0.42418608439039535,0.4651162952033166,	0.5153488705285886,	0.5600000150813612,	0.5934883885773017,	0.6279069804268439,	0.6437209337398515]	



#Make the time_list and x_list dimensionless
time_list = np.array(time_list) * time_scaling

#Plot outputs for solid front position
fig0 = plt.figure()
ax0 = fig0.add_subplot(111)
ax0.plot(time_list, x_list_top, 'g-.', label="Simulation (top)")
ax0.plot(top_x_exp, top_y_exp, 'gs', label="Experiment (top) - G&V (1986)")
ax0.plot(time_list, x_list_mid, 'r', label="Simulation (center)")
ax0.plot(center_x_exp, center_y_exp, 'ro', label="Experiment (center) - G&V (1986)")
ax0.plot(time_list, x_list_bot, 'b--', label="Simulation (bottom)")
ax0.plot(bottom_x_exp, bottom_y_exp, 'b^', label="Experiment (bottom) - G&V (1986)")
ax0.set_ylabel(r'$\delta x$')
ax0.set_xlabel(r'$\tau$')
ax0.set_xlim([0, 0.07])
ax0.legend(loc="upper left")
fig0.savefig(f'./xmax-t.png')
plt.show()
plt.close()

#Load melted volume
t_list,melt_frac=np.loadtxt(output_path+"/liquid_fraction.dat",unpack=True,skiprows=1)

#Plot melted volume
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(t_list * time_scaling , melt_frac, '-k', label="Simulation")
ax1.plot(melted_x_exp, melted_y_exp, 'ro', label="Experiment - G&V (1986)")
ax1.plot(melted_x_blais, melted_y_blais, 'bs', label="Simulation - B&I (2018)")
ax1.set_ylabel(r'melted volume fraction')
ax1.set_xlabel(r'$\tau$')
ax1.set_xlim([0, 0.07])
ax1.legend(loc="upper left")
fig1.savefig(f'./melted-volume-fraction.png')
plt.show()
plt.close()
