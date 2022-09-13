#############################################################################
"""
Postprocessing code for gls_VOF_dam-break_Martin_and_Moyce example

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
from pyvista.plotting.renderer import CameraPosition

from scipy.linalg.special_matrices import dft
#############################################################################

#############################################################################
'''Functions'''

#############################################################################

#Constants
L1 = 3.5
g  = 1.0
time_correction = 0.175

#Take case path as argument and store at currentPath
currentPath = sys.argv[1]

#Define list of VTK files and time list:
list_vtu = os.listdir(currentPath)
list_vtu = [x for x in list_vtu if "pvd" not in x ]
list_vtu = [x for x in list_vtu if "pvtu" not in x ]
list_vtu = [x for x in list_vtu if "vtu" in x ]


#Create a list of time_steps
data = pd.read_csv("output/dam-break_VOF.pvd",sep='"',header=5, usecols=[1],skiprows=[38,39]) 
data.columns = ["a"] 
time_list = data['a']
print(data['a'])

#Set phase_limit to search for maximum x
phase_limit = 0.1

#Sort list
time_list, list_vtu = (list(t) for t in zip(*sorted(zip(time_list, list_vtu))))

#Create a list to fill with maximum x in which phase > phase_limit
x_list = []

#Read vtu data
for i in range(0, len(list_vtu)):
    #Read DF from VTK files
    exec(f'df_{i} = pv.read(f\'{currentPath}/{list_vtu[i]}\')')

    #Select a data to apply the slice   
    exec(f'df = df_{i}')

    #find max 'x' in which phase > 0
    points = pd.DataFrame(df.points[:, 0])
    phase  = pd.DataFrame(df['phase'])

    x_max = max(points[phase[0] > phase_limit].values)[0]
    x_list.append(x_max)

print('List of maximum values for x:')
print(x_list)

#Sort x_list, since the time_list is extracted from the .pvd file, the x_list and time_list are not consistent
x_list.sort()

#Experimental data from Martin & Moyce 1952
x_exp = [0, 0.41, 0.84, 1.19, 1.43, 1.63, 1.82, 1.97, 2.2, 2.32, 2.5, 2.64, 2.82, 2.96]
y_exp = [1, 1.11, 1.23, 1.44, 1.67, 1.89, 2.11, 2.33, 2.56, 2.78, 3, 3.22, 3.44, 3.67]

#Make the time_list and x_list dimensionless
time_list = [x * ((2 * g / L1) ** 0.5) for x in time_list]
x_list = [x / L1 for x in x_list]

#Time-correction
time_list = [x + time_correction for x in time_list]

fig0 = plt.figure()
ax0 = fig0.add_subplot(111)
ax0.plot(time_list, x_list, '-ok', label="Simulation")
ax0.plot(x_exp, y_exp, 'ro',label="Experiment - Martin and Moyce (1952)")
ax0.set_ylabel(r'$\delta$')
ax0.set_xlabel(r'$\tau$')
ax0.set_xlim([0, 3.5])
ax0.set_ylim([1, 4])
ax0.legend(loc="upper left")
fig0.savefig(f'./xmax-t.png')
plt.show()
