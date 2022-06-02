#############################################################################
"""
Postprocessing code for melting_cavity example

"""
#############################################################################

'''Parameters'''

#############################################################################
'''Importing Libraries'''
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv
import os
import sys
#############################################################################

#############################################################################
'''Functions'''

#############################################################################

#Constants
tolerance 	= 0.001
max_y		= 0.714
y_top		= 0.85 * max_y
y_center	= 0.6 * max_y
y_bottom	= 0.18 * max_y
n_files 	= 400
t_to_tau_factor	= 0.00000164

#Take case path as argument and store at currentPath
currentPath = sys.argv[1]

#Define list of VTK files and time list:
list_vtu = os.listdir(currentPath)
list_vtu = [x for x in list_vtu if "pvd" not in x ]
list_vtu = [x for x in list_vtu if "pvtu" not in x ]

#Create a list of time_steps
data = pd.read_csv("output/melting.pvd",sep='"',header=5, usecols=[1],skiprows=[38,39]) 
data.columns = ["a"] 
time_list = data['a']
print(data['a'])

#Set temperature_threshold to search for liquid-solid interface
T_threshold = 105.05

#Sort list
time_list, list_vtu = (list(t) for t in zip(*sorted(zip(time_list, list_vtu))))

#Create a list to fill with the interface location
x_list_top = []
x_list_center = []
x_list_bottom = []
melted_volume_percentage = []

#Read vtu data
for i in range(0, len(list_vtu)):
    #Read DF from VTK files
    exec(f'df_{i} = pv.read(f\'{currentPath}/{list_vtu[i]}\')')

    #Select a data to apply the slice   
    exec(f'df = df_{i}')

    #Find max 'x' in which T > T_threshold
    top_filter = abs(df.points[:,1] - y_top) < tolerance
    filtered_temperatures = pd.Series(pd.DataFrame(df['temperature'])[top_filter][0]).reset_index(drop = True)
    top_points_temperatures = pd.concat([pd.DataFrame(df.points[top_filter][: , 0]), filtered_temperatures] , axis = 1, ignore_index = True)
    
    center_filter = abs(df.points[:,1] - y_center) < tolerance
    filtered_temperatures = pd.Series(pd.DataFrame(df['temperature'])[center_filter][0]).reset_index(drop = True)
    center_points_temperatures = pd.concat([pd.DataFrame(df.points[center_filter][: , 0]), filtered_temperatures] , axis = 1, ignore_index = True)
    
    bottom_filter = abs(df.points[:,1] - y_bottom) < tolerance
    filtered_temperatures = pd.Series(pd.DataFrame(df['temperature'])[bottom_filter][0]).reset_index(drop = True)
    bottom_points_temperatures = pd.concat([pd.DataFrame(df.points[bottom_filter][: , 0]), filtered_temperatures] , axis = 1, ignore_index = True)

    x_list_top.append(max(top_points_temperatures[top_points_temperatures[1] > T_threshold][0], default=0))
    x_list_center.append(max(center_points_temperatures[center_points_temperatures[1] > T_threshold][0], default=0))
    x_list_bottom.append(max(bottom_points_temperatures[bottom_points_temperatures[1] > T_threshold][0], default=0))
    
    melted_volume_percentage.append(sum(df['temperature'] > T_threshold) / df['temperature'].size)
    
    #Clean variables
    del top_points_temperatures, center_points_temperatures, bottom_points_temperatures
    del top_filter, center_filter, bottom_filter
    del filtered_temperatures
    
    print("Finished reading ", round(i / n_files * 100, 2), "% of the files")
    
#Sort x_lists, since the time_list is extracted from the .pvd file, the x_lists and time_list are not consistent
x_list_top.sort()
x_list_center.sort()
x_list_bottom.sort()
melted_volume_percentage.sort()

#Experimental data from Gau and Viskanta 1986
top_x_exp = [0.00853, 0.01279, 0.01931, 0.02575, 0.03428, 0.04298, 0.0536, 0.06438]
top_y_exp = [0.11508, 0.18212, 0.28603, 0.31508, 0.38994, 0.50391, 0.61117, 0.72514]

center_x_exp = [0.00836, 0.01279, 0.01931, 0.02584, 0.0342, 0.04298, 0.0536, 0.06455]
center_y_exp = [0.11508, 0.16313, 0.20894, 0.23911, 0.29609, 0.36089, 0.45475, 0.52961]

bottom_x_exp = [0.00853, 0.01271, 0.01931, 0.02575, 0.0342, 0.04306, 0.0536, 0.06455]
bottom_y_exp = [0.1095, 0.12849, 0.14078, 0.16425, 0.21117, 0.22346, 0.26592, 0.26592]

melted_x_exp = [0.00946, 0.01433, 0.02113, 0.02877, 0.03774, 0.04739, 0.05956, 0.07164, 0.08104, 0.09069]
melted_y_exp = [0.11972, 0.16112, 0.20476, 0.22713, 0.2786, 0.34014, 0.43077, 0.5035, 0.56392, 0.61664]

#Make the time_list and x_list dimensionless
time_list = [x * t_to_tau_factor for x in time_list]

#Plot outputs
fig0 = plt.figure()
ax0 = fig0.add_subplot(111)
ax0.plot(time_list, x_list_top, 'g-.', label="Simulation (top)")
ax0.plot(top_x_exp, top_y_exp, 'gs', label="Experiment (top) - G&V (1986)")
ax0.plot(time_list, x_list_center, 'r', label="Simulation (center)")
ax0.plot(center_x_exp, center_y_exp, 'ro', label="Experiment (center) - G&V (1986)")
ax0.plot(time_list, x_list_bottom, 'b--', label="Simulation (bottom)")
ax0.plot(bottom_x_exp, bottom_y_exp, 'b^', label="Experiment (bottom) - G&V (1986)")
ax0.set_ylabel(r'$\delta x$')
ax0.set_xlabel(r'$\tau$')
ax0.set_xlim([0, 0.07])
ax0.legend(loc="upper left")
fig0.savefig(f'./xmax_t.png')
plt.show()
plt.close()

#Plot melted volume
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(time_list, melted_volume_percentage, '-k', label="Simulation")
ax1.plot(melted_x_exp, melted_y_exp, 'ro', label="Experiment - G&V (1986)")
ax1.set_ylabel(r'melted volume fraction')
ax1.set_xlabel(r'$\tau$')
ax1.set_xlim([0, 0.07])
ax1.legend(loc="upper left")
fig1.savefig(f'./melted_volume_fraction.png')
plt.show()
plt.close()
