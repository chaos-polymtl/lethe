"""
Postprocessing code for Jurin's law example
This code extracts the difference in height between the meniscus
and the fluid on the side.
"""
# -------------------------------------------
# Modules
# -------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv
from natsort import os_sorted
from math import *

import os
import sys

# --------------------------------------------
# Main
# --------------------------------------------

# To make it work, type "python3 rayleigh-taylor_postprocess.py ./output/adaptive/" or
# "python3 rayleigh-taylor_postprocess.py ./output/constant/" into the terminal.

def get_deltaH(output_path,prm):
    print(output_path)
    
    phase_limit = prm.phase_limit
    g = prm.g
    # Take case path as argument and store it

    # Define list of VTU files
    list_vtu = os.listdir(output_path)
    pvd_file = [x for x in list_vtu if (x.endswith('.pvd'))]
    list_vtu = [x for x in list_vtu if ("pvtu" in x)]
    list_vtu = set(list_vtu)

    # Sort VTU files to ensure they are in the same order as the time step
    list_vtu = os_sorted(list_vtu)
    # Read the pvd file to extract the times

    reader = pv.get_reader(output_path + "/" + pvd_file[0])

    # Get active times
    time_list = reader.time_values
    #  Create lists to fill the height of different points of the meniscus (defined by y such that phase(y) < phase_limit)
    meniscus_height = []

    # Create lists to fill the height of different points of the right side (defined by y such that phase(y) < phase_limit)
    side_height=[]
    side_heights=[[]]

    # Create beginning and end points for meniscus lines and side lines
    a_meniscus = [1, 0, 0]
    b_meniscus = [1, 8, 0]

    a_side_left = [2.1, 0, 0] # eventuellement décaler légèrement pour que ce soit dans le domaine
    b_side_left = [2.1, 8, 0]

    side_lines = [[[2, 0, 0], [2, 8, 0]]]

    deltaH = []

    # Read VTU files
    for j,vtu_file in enumerate(list_vtu):
        sim = pv.read(f"{output_path}/{vtu_file}")
    
        sampled_data_meniscus = sim.sample_over_line(a_meniscus, b_meniscus, resolution=1000)
        phase_meniscus = pd.DataFrame(sampled_data_meniscus["phase_order"])
        points_meniscus = pd.DataFrame(sampled_data_meniscus.points)
    
        # Find min 'y' in phase > phase_limit (MENISCUS)
        fluid1_points_meniscus = points_meniscus[phase_meniscus[0] < phase_limit].values
        y_min_meniscus = fluid1_points_meniscus[0][1]
        points_heights=[points[1] for points in fluid1_points_meniscus]
        meniscus_height.append(np.max(points_heights))
        
        for i in range(0, 1):
            sampled_data_side = sim.sample_over_line(side_lines[i][0], side_lines[i][1],     resolution=1000)
            phase_side = pd.DataFrame(sampled_data_side["phase_order"])
            points_side = pd.DataFrame(sampled_data_side.points)
            
            # Find min 'y' in phase > phase_limit (SIDE)
            fluid1_points_side = points_side[phase_side[0] < phase_limit].values
            y_min_side = fluid1_points_side[0][1]
            for points in fluid1_points_side:
                if points[1] > y_min_side:
                    y_min_side = points[1]
            side_heights[i].append(y_min_side)
            
        side_height.append(np.mean([side_heights[i][-1] for i in range(len(side_heights))]))
        deltaH.append(-side_height[-1] + meniscus_height[-1])
        
    return deltaH, time_list
        
def analytical_solution(prm,angle):
    g = prm.g
    mu = prm.mu
    sigma = prm.sigma
    rho_l=prm.rho_l
    r = prm.r
    
    angle = np.pi*angle/180 #given in degrees in the prm file
    
    angle_correction = (r/2*cos(angle))*(2-sin(angle)-asin(cos(angle))/cos(angle)) 
    
    analytical_deltaH = (sigma*cos(angle))/(rho_l*g*r) - angle_correction
    print(analytical_deltaH)
    
    return analytical_deltaH*1000
      
#deltaH,time_list = get_deltaH('verifications/150',prm)
# Figure
#fig0 = plt.figure()
#ax0 = fig0.add_subplot(111)
#ax0.plot(time_list, deltaH, '-k', linewidth=2, label="Delta H")
#ax0.plot(time_list, meniscus_height, '-r', linewidth=2, label="Meniscus height")
#ax0.plot(time_list, side_height, '-b', linewidth=2, label="Side height")
#ax0.set_ylabel(r'$\Delta h$')
#ax0.set_xlabel(r'$t$')
# ax0.set_xlim([0, 4.5])
# ax0.set_ylim([0.0, 0.008])
#ax0.legend(loc="upper left")
#plt.title("DeltaH evolution in time")
#fig0.savefig('./delta_H.png')
#plt.show()



