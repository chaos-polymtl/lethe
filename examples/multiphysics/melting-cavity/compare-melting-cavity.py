#############################################################################
"""
Postprocessing code to compare cases in the melting cavity example

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

# List of name for the cases used in the comparison
case_list=["Lethe - Viscous penality", "Lethe - Darcy penality", "Lethe - Hybrid penality"]






#Experimental data from Gau and Viskanta 1986
melted_x_exp = [0.00946, 0.01433, 0.02113, 0.02877, 0.03774, 0.04739, 0.05956, 0.07164, 0.08104, 0.09069]
melted_y_exp = [0.11972, 0.16112, 0.20476, 0.22713, 0.2786, 0.34014, 0.43077, 0.5035, 0.56392, 0.61664]

#Simulation data by Blais and Ilinca 2018
melted_x_blais=[0.0009408610966497919,0.0026209691082628825,0.004569890217872624 ,0.0068548364163561895,0.011155915715326668,	
                0.018010754310777356,0.024731181999040716,0.034005377700162306,0.03891129082781425,0.043682793048278865,0.04811827889624767 ,0.0532258062055861,0.02963709512669266 ,0.05759408877905576 , 		
                0.06451613064900559,0.07432795690430948,0.08266128497149776,0.08931451374345097,0.09596774251540417,0.09899193344975654 ]	
melted_y_blais= [0.0372093373985142, 0.0669768073173256, 0.0855814156911383, 0.1079070181302468, 0.1423255798170668, 0.1897674397560892, 0.2344186446343063, 0.2874418750203834, 0.3162790662601486, 0.3441860391463121,0.37209301203247547,0.39999998491863886,0.2613953388414236,	0.42418608439039535,0.4651162952033166,	0.5153488705285886,	0.5600000150813612,	0.5934883885773017,	0.6279069804268439,	0.6437209337398515]	


#Plot melted volume
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(melted_x_exp, melted_y_exp, 'o', label="Experiment - G&V (1986)")
ax1.plot(melted_x_blais, melted_y_blais, 's', label="Simulation - B&I (2018)")
ax1.set_ylabel(r'melted volume fraction')
ax1.set_xlabel(r'$\tau$')
ax1.set_xlim([0, 0.07])

#Load melted volume
for c in range(1,len(sys.argv)):
    t_list,melt_frac=np.loadtxt(sys.argv[c]+"/liquid_fraction.dat",unpack=True,skiprows=1)
    ax1.plot(t_list * time_scaling , melt_frac, '-', label=case_list[c-1])


ax1.legend(loc="upper left")
fig1.savefig(f'./comparison-melted-volume-fraction.png')
plt.show()
plt.close()
