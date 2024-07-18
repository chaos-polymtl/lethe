from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from paraview.simple import *
from scipy.interpolate import UnivariateSpline
from post_processing_functions import *

# Path to the granuheap experimental result
exp_path = 'experimental_result.png' 
# Name of simulation output (see OUTPUT NAME set in the simulation subsection of the parameter file)
num_output = 'granuheap'
# Output path (see OUTPUT PATH set in the simulation subsection of the parameter file)
out_path = 'output'
# Number of pixels in height and width of your experimental support (to adjust if you change experimental result)
height_exp = 60
width_exp = 85

############################## SCALING ######################################
print('\n ******************************** SCALING **************************** \n')

# Scaling experiment & numeric
scale = 0.07 # Never use 0 (it will produce a white image)
scaling, height, width, to_print = scaling_paraview(width_exp, f'', num_output, out_path, scale)

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
        scaling, height, width, to_print = scaling_paraview(width_exp, f'', num_output, out_path, scale)
    else:
        sorted_indices = sorted(range(len(width_list)), key=lambda k: width_list[k])
        width_list_sorted = [width_list[i] for i in sorted_indices]
        scale_list_sorted = [scale_list[i] for i in sorted_indices]
        spline = UnivariateSpline(width_list_sorted, scale_list_sorted, k=3)
        x = width_exp
        scale = float(spline(x).item())
        scaling, height, width, to_print = scaling_paraview(width_exp, f'', num_output, out_path, scale)
    scale_list.append(scale)
    width_list.append(width)
    i += 1
print(to_print)


############################## NUMERIC ######################################
print('\n ******************************** NUMERIC **************************** \n')

# If scaling is True, continue with the map processing
if scaling:
    # Take pictures in paraview
    picture_paraview(f'', num_output, out_path, scale)
    # Create a list of file path
    list_files  = find_files(f'map')
    # Create a map with all files in the list
    map_processing(list_files, f'map.png', height, height_exp)

############################## DIFFERENCE ######################################
    print('\n ******************************** DIFFERENCE **************************** \n')
    difference_images(exp_path, f'map.png', f'image_difference.png')
    error = RMSE(exp_path, f'map.png')

