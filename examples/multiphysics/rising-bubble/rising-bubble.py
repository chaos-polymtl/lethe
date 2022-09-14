#############################################################################
"""
Postprocessing code for rising-bubble example

"""
#############################################################################

'''Importing Libraries'''
from math import pi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv

import glob

import os
import sys
#############################################################################

#Take case path as argument and store it
output_path = sys.argv[1]

# Read the pvd file to extract the times
reader = pv.get_reader("output/rising-bubble.pvd")
# Get active times
time_list = reader.time_values



#Define list of VTU files
list_vtu = os.listdir(output_path)
list_vtu = [x for x in list_vtu if  ("vtu" in x and "pvtu" not in x) ]

# Sort VTU files to ensure they are in the same order as the time step
list_vtu = sorted(list_vtu)

#Set phase_limit to search for maximum x
phase_limit = 0.5

#Create lists to fill with average y, dy, dt and bubble rise velocity
y_list = []
y_diff_list = []
t_diff_list = []
bubble_rise_vel = []

#Read vtu data
for i in range(0, len(list_vtu)):
    #Read DF from VTK files
    exec(f'df_{i} = pv.read(f\'{output_path}/{list_vtu[i]}\')')

    #Select a data to apply the slice   
    exec(f'df = df_{i}')

    #find average 'y' by averaging minimum and maximum y of the bubble
    points_y = pd.DataFrame(df.points[:, 1])
    phase  = pd.DataFrame(df['phase'])
    
    y_max = max(points_y[phase[0] > phase_limit].values)[0]
    y_min = min(points_y[phase[0] > phase_limit].values)[0]
    y_mean = (y_min + y_max) * 0.5
    y_list.append(y_mean)
       

# This array could be smoothed to obtain a better velocity field
smoothed_y= np.array(y_list)

#Calculate bubble rise velocity
length = len(smoothed_y)
bubble_rise_vel = [None] * length
for i in range(length):
    if i == 0 :
        bubble_rise_vel[i] = (smoothed_y[1] - smoothed_y[0]) / (time_list[1]-time_list[0])
    elif i == (length-1) :
        bubble_rise_vel[i] = (smoothed_y[-1] -  smoothed_y[-2]) / (time_list[-1]-time_list[-2])
    else :
        print (i)
        dt_0 = time_list[i]-time_list[i-1]
        dt_1 = time_list[i+1]-time_list[i]
        bubble_rise_vel[i] = (smoothed_y[i+1] - (dt_0/dt_1)*smoothed_y[i-1] - (1-dt_0/dt_1)*smoothed_y[i]) / (dt_1*(1+dt_0/dt_1))


#Data from Zahedi, Kronbichler and Kreiss (2012)
x_ref = [0.047, 0.139, 0.232 ,0.324 ,0.416 ,0.508 ,0.601 ,0.693 ,0.785 ,0.878 ,0.97 ,1.062 ,1.154 ,1.247 ,1.339 ,1.431 ,1.524 ,1.616 ,1.708 ,1.8 ,1.893 ,1.985 ,2.077 ,2.17 ,2.262 ,2.354 ,2.446 ,2.539 ,2.631 ,2.723 ,2.816 ,2.908 ,2.979]
y_ref = [0.501 ,0.505 ,0.513 ,0.525 ,0.54 ,0.557 ,0.576 ,0.597 ,0.618 ,0.641 ,0.663 ,0.685 ,0.707 ,0.729 ,0.75 ,0.77 ,0.791 ,0.811 ,0.83 ,0.849 ,0.868 ,0.886 ,0.904 ,0.922 ,0.94 ,0.958 ,0.976 ,0.993 ,1.011 ,1.029 ,1.047 ,1.064 ,1.078]
x_vel = [0.001 ,0.03 ,0.056 ,0.081 ,0.106 ,0.136 ,0.17 ,0.203 ,0.237 ,0.271 ,0.304 ,0.338 ,0.372 ,0.41 ,0.452 ,0.494 ,0.544 ,0.607 ,0.687 ,0.779 ,0.872 ,0.964 ,1.056 ,1.149 ,1.241 ,1.333 ,1.426 ,1.518 ,1.61 ,1.702 ,1.795 ,1.887 ,1.979 ,2.072 ,2.164 ,2.256 ,2.348 ,2.441 ,2.533 ,2.625 ,2.718 ,2.81 ,2.902 ,2.978]
y_vel = [0.004 ,0.017 ,0.029 ,0.041 ,0.053 ,0.066 ,0.081 ,0.095 ,0.109 ,0.122 ,0.135 ,0.147 ,0.158 ,0.17 ,0.182 ,0.194 ,0.205 ,0.218 ,0.229 ,0.237 ,0.24 ,0.241 ,0.239 ,0.236 ,0.231 ,0.227 ,0.222 ,0.217 ,0.213 ,0.208 ,0.205 ,0.201 ,0.198 ,0.196 ,0.194 ,0.192 ,0.192 ,0.191 ,0.191 ,0.192 ,0.192 ,0.193 ,0.193 ,0.194]

fig0 = plt.figure()
ax0 = fig0.add_subplot(111)
ax0.plot(time_list, smoothed_y, '-ok', label="Simulation")
ax0.plot(x_ref, y_ref, 'ro',label="Reference - Zahedi, Kronbichler and Kreiss (2012)")
ax0.set_ylabel(r'Bubble center height')
ax0.set_xlabel(r'$t$')
ax0.legend(loc="upper left")
fig0.savefig(f'./ymean-t.png')
plt.show()

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(time_list, bubble_rise_vel , '-ok', label="Simulation")
ax1.plot(x_vel, y_vel, 'ro',label="Reference - Zahedi, Kronbichler and Kreiss (2012)")
ax1.set_ylabel(r'Rise velocity')
ax1.set_xlabel(r'$t$')
ax1.legend(loc="upper left")
ax1.legend(loc=4)
fig1.savefig(f'./bubble-rise-velocity.png')
plt.show()
