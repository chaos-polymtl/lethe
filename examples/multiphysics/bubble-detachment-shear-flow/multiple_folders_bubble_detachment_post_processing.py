"""
Quick user guide : 

This python script is used to plot the contour of the bubble at detachment time.

How to use this script?

1- Create a first directory (we will call it outputs but feel free to name it as you like).
2- For each simulation, create a new directory in outputs : output1, output2, output3,...
3- Run the simulations and store the results in their correct directories
4- Execute the script with the path of outputs as argument :
   python3 multiple_folders_new_bubble_detachment_post_processing /PATH/TO/OUTPUTS/DIR /PATH/TO/STORE/RESULTS
5- Enjoy your plots. 

"""

#-------------------------------------------
# Modules
#-------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import pyvista as pv
import os
import sys
from natsort import os_sorted

import scienceplots
plt.style.use(['science','ieee'])


from functions_bubble_detachment_post_processing import bubble_volume, get_numerical_detachment_time, get_contour_at_fixed_time, get_volume_at_time

#For controlling font sized globally
SMALL_SIZE = 5
MEDIUM_SIZE = 8
BIGGER_SIZE = 15

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', labelsize=8)              # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Get first and second argument
rootdir = sys.argv[1]

folder_name_list = []
root, dirs, files = next(os.walk(rootdir, topdown=True))

# The following lines allow for a user to type outputs/ or outputs as argument of the script.
for dir in dirs:
    if (root[-1] == "/"):
        folder_name_list.append(str(root + dir))
    else:
        folder_name_list.append(str(root + "/" + dir))

folder_name_list = os_sorted(folder_name_list)
dirs = os_sorted(dirs)

# Create legend entries 
pparams=dict(xlabel=r'$x\text{ [mm]}$', ylabel=r'$y\text{ [mm]}$')

with plt.style.context(['science', 'ieee']): # Syntax to use a specific plot style with the scienceplots module
    fig = plt.figure()
    ax = fig.add_subplot(111)


    for i in range(len(folder_name_list)): # Loop on the different directories in the output directory
         print(f'Reading directory : {folder_name_list[i]}')
         
         # The three following lines gather the detachment time, detachment volume and the contour of the bubble at detachment.
         time, detachment_time,x,y = get_numerical_detachment_time(folder_name_list[i])
         absolute_volume, time_integrated_flux_volume, analytical_volume_array, time = bubble_volume(
         folder_name_list[i])
         detachment_volume = get_volume_at_time(absolute_volume,time,detachment_time)
         
         # Displays the detachment time and detachment volume for each directory
         print(f'For the folder : {folder_name_list[i]} \n' +
               f'Detachment time : t_det = {detachment_time} s \n' +
               f'Detachment volume : V_det = {detachment_volume} m³')
         label_loop =r' ($t_{det}=$'+f'{detachment_time:2f} s)'
         
         # Plot and save the contour of the bubble at detachment
         ax.scatter(x, y, s=0.5, marker=".", label=label_loop)


    # Perform some modifications of the legend for aesthetics.
    ax.set_xlim(left=-1)
    ax.set(**pparams)
    handles, labels = plt.gca().get_legend_handles_labels()
    fig.legend(loc='outside center right',frameon = True,edgecolor='k',prop={'size': MEDIUM_SIZE},ncol=1, fancybox=False, bbox_to_anchor=(0.80, -0.15))
    fig.savefig('./' + savedir + '_bubble_contour.pdf', format="pdf", dpi=500)
    fig.savefig('./' + savedir + '_bubble_contour.png', format="png", dpi=500)


