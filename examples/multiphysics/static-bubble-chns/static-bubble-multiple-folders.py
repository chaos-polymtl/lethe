# SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Quick user guide: 

This python script is used to plot a slice of the pressure field and the difference between the averaged pressure outside of the bubble and inside the bubble.

How to use this script?

1- Create a first directory (we will call it outputs but feel free to name it as you like).
2- For each simulation, create a new directory in outputs : output1, output2, output3,...
   For this step, it may be easier to use the "generate_cases_locally.py" script to generate
   all the directories. Then launch the simulations with "launch_lethe_locally.py" and then 
   use "static-bubble-multiple-folders.py"
3- Run the simulations and store the results in their correct directories
4- Execute the script with the path of outputs as argument :
   python3 multiple_folders_new_bubble_detachment_post_processing /PATH/TO/OUTPUTS/DIR 
5- Enjoy your plots. 

This should be called as follows : python3 static-bubble-multiple-folders.py /PATH/TO/OUTPUTS/DIR
"""
# -------------------------------------------
# Modules
# -------------------------------------------

from postprocessing_static_bubble import get_pressure_difference, analytical_solution,get_pressure_slice
import numpy as np
import os
import sys
import scienceplots
from natsort import os_sorted
import matplotlib.pyplot as plt

# For plotting nice Latex-style plots and controlling the sizes of the texts on
# the plot
plt.style.use(['science', 'ieee'])

SMALL_SIZE = 10
MEDIUM_SIZE = 10
BIGGER_SIZE = 15
plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=8)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=8)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Class containing relevant physical properties of the problem
class parametres():
    sigma = 1
    radius = 0.15

    def set_radius(self, r):
        self.radius = r


prm = parametres()

if (len(sys.argv)==0):
    raise Exception(f"No outputs folder was specified as an argument of the python script.")
rootdir = sys.argv[1]

folder_name_list = []
root, dirs, files = next(os.walk(rootdir, topdown=True))

for dir in dirs:
    folder_name_list.append(str(root + "/" + dir))

folder_name_list = os_sorted(folder_name_list)
dirs = os_sorted(dirs)

# Syntax specific to the SciencePlots module 
with plt.style.context(['science', 'ieee']):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pparams = dict(xlabel=r'$\text{Rayon} (R) \text{[m]}$',
                   ylabel=r'$\Delta p \text{[Pa]}$')

plots = []

for i in range(len(folder_name_list)):
    # Get the radius of the case from the directory name. (Works if the subdirectories were generated using the generate_cases_locally.py script)
    radius = float(dirs[i][2:])
    prm.set_radius(radius)
    # Get the numerical pressure difference from each directory
    pressure_diff = get_pressure_difference(folder_name_list[i], prm)
    # Store the plot handle with label "Lethe"
    plot_handle, = ax.plot(radius, pressure_diff, 'k*', label="Lethe")
    if i == 0:  # Store the handle only once
        plots.append(plot_handle)

# Store the analytical solution plot handle
radius_array, pressure_diff_analytical = analytical_solution(prm)
analytical_handle, = ax.plot(radius_array, pressure_diff_analytical, label=r'$\Delta p$ analytique (Young-Laplace)')
plots.append(analytical_handle)

# Set the legend with the stored handles
ax.legend(handles=plots, loc='upper right', frameon=True, edgecolor='k', prop={'size': 8}, ncol=1)

ax.set(**pparams)
ax.set_xlim([0, 0.5])
ax.set_ylim([2, 10])
# Save the figures
fig.savefig('pressure-difference.pdf', format="pdf", dpi=500)
fig.savefig('pressure-difference.png', format="png", dpi=500)

# Syntax specific to the SciencePlots module     
with plt.style.context(['science', 'ieee']):
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    pparams1 = dict(xlabel=r'$x \text{[m]}$',
                    ylabel=r'$\mathrm{Pression} \text{[Pa]}$')
    # Get the pressure profile for each subdirectory of the main output directory
    for i in range(len(folder_name_list)):
        p, x = get_pressure_slice(folder_name_list[i])
        label_loop = r'$R = $ ' + dirs[i][2:] + r'$ \text{[m]}$'
        ax2.plot(x, p,label=label_loop )
    ax2.legend(loc='upper left', frameon=True, edgecolor='k',
               prop={'size': 8.5}, ncol=1)
    ax2.set(**pparams1)
    ax2.grid(which='major', color='lightgrey', linestyle='--',alpha=0.8)
    fig2.savefig('pressure-profile.pdf', format="pdf", dpi=500)
    fig2.savefig('pressure-profile.png', format="png", dpi=500)
