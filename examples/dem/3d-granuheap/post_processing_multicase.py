# SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

from PIL import Image
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from matplotlib.colors import LinearSegmentedColormap
from post_processing_functions import *

# Path of the granuheap experimental result
exp_path = 'experimental_result.png' 
# Name of directory for each simulation (see CASE_PREFIX from the launch_lethe_locally.py file used)
num_name = 'wetsand'
# Name of simulation output (see OUTPUT NAME set in the simulation subsection of the parameter file)
num_output = 'granuheap'
# Output path (see OUTPUT PATH set in the simulation subsection of the parameter file)
out_path = 'output'
# Definition of variable parameters
parameter1_name = 'Surface Energy'
parameter1 = [0.0010, 0.0040, 0.0070, 0.0100]
parameter2_name = 'Rolling Friction'
parameter2 = [0.70, 0.50, 0.30]
# Number of pixels in height and width of your experimental support (to adjust if you change experimental result)
height_exp = 60
width_exp = 85


############################## SCALING ######################################
print('\n ******************************** SCALING **************************** \n')

# Scaling experiment & numeric
scale = 0.07 # Never use 0 (it will produce a white image)
scaling, height, width, to_print = scaling_paraview(width_exp, f'{num_name}_{parameter1[0]:.4f}_{parameter2[0]:.2f}/', num_output, out_path, scale)

# If scaling is a succes, map processing for each simulation
scale_list = [scale]
width_list =[width]
i = 0
while scaling is False and i < 25: # This loop is necessary to scale. The scaling use a maximum of 25 iterations.
    if i < 3:
        if width_list[i] > width_exp:
            scale += 0.01
        elif width_list[i] < width_exp:
            scale -= 0.01
        scaling, height, width, to_print = scaling_paraview(width_exp, f'{num_name}_{parameter1[0]:.4f}_{parameter2[0]:.2f}/', num_output, out_path, scale)
    else:
        sorted_indices = sorted(range(len(width_list)), key=lambda k: width_list[k])
        width_list_sorted = [width_list[i] for i in sorted_indices]
        scale_list_sorted = [scale_list[i] for i in sorted_indices]
        spline = UnivariateSpline(width_list_sorted, scale_list_sorted, k=3)
        x = width_exp
        scale = float(spline(x).item())
        scaling, height, width, to_print = scaling_paraview(width_exp, f'{num_name}_{parameter1[0]:.4f}_{parameter2[0]:.2f}/', num_output, out_path, scale)
    scale_list.append(scale)
    width_list.append(width)
    i += 1
print(to_print)

############################## NUMERIC ######################################
print('\n ******************************** NUMERIC **************************** \n')

if scaling:
    for n in parameter1:
        for u in parameter2:
            path = f'{num_name}_{n:.4f}_{u:.2f}/'
            # Take pictures in paraview
            picture_paraview(path, num_output, out_path, scale)
            # Create a list of file paths 
            list_files  = find_files(f'{num_name}_{n:.4f}_{u:.2f}/map')
            # Create a map with all files in the list
            map_processing(list_files, f'map/map_{n:.4f}_{u:.2f}.png', height, height_exp)

############################## DIFFERENCE ######################################
    print('\n ******************************** DIFFERENCE **************************** \n')
    error_matrix = np.empty((len(parameter2), len(parameter1)))
    image_matrice = np.empty((len(parameter2), len(parameter1)), dtype=object)
            
    for i, n in enumerate(parameter1):
        for j, u in enumerate(parameter2):
            difference_images(exp_path, f'map/map_{n:.4f}_{u:.2f}.png', f'im/image_difference_{n:.4f}_{u:.2f}.png')
            error_matrix[j, i] = RMSE(exp_path, f'map/map_{n:.4f}_{u:.2f}.png')
            image_matrice[j, i] = f'im/image_difference_{n:.4f}_{u:.2f}.png'


############################## POST-PROCESSING ######################################
    print('\n ******************************** POST-PROCESSING **************************** \n') 
    # HeatMap for RMSE 
    sns.set_theme() 
    plt.figure(figsize=(10, 8))  
    sns.heatmap(error_matrix, annot=True, cmap='YlOrBr', fmt='.2f', linewidths=.5, xticklabels=parameter1, yticklabels=parameter2 )
    plt.title('Root Mean Squared Error (RMSE) Heatmap ')
    plt.xlabel(parameter1_name)
    plt.ylabel(parameter2_name)
    plt.savefig('im/error_values_heatmap.png')
    # plt.show()

    # Shape profile error (Experimental-Numeric)
    fig, axs = plt.subplots(len(parameter2), len(parameter1), figsize=(10, 8), constrained_layout=True)
    colors = [(1, 0, 0), (1, 1, 1), (0, 0, 1)] 
    cmap_name = 'red_white_blue'
    cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=256)
    for i, n in enumerate(parameter1):
        for j, u in enumerate(parameter2):
            image_path = image_matrice[j, i]  
            image = plt.imread(image_path)  
            im = axs[j, i].imshow(image, cmap=cm, vmin=-1, vmax=1)  
            axs[j, i].set_xticks([])
            axs[j, i].set_yticks([])
            if i == 0:
                axs[j, i].set_ylabel(parameter2[j])
            if j == len(parameter2) -1:
                axs[j, i].set_xlabel(parameter1[i])
    cbar = fig.colorbar(im, ax=axs.ravel().tolist(), orientation='vertical', fraction=0.05, pad=0.04)  
    cbar.set_label('Signed Error')  
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])  
    plt.suptitle('Profile shape error (Experimental-Numeric)')

    fig.supxlabel(parameter1_name, fontsize=12, ha='center')
    fig.supylabel(parameter2_name, fontsize=12, va='center', rotation='vertical')
    plt.savefig('im/profile_shape_error.png')


    plt.show()
