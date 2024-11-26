# SPDX-FileCopyrightText: Copyright (c) 2022 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for 2D-backward-facing-step example
Computes velocity distributions at outlet and
compares it with analytical solution (Poiseuille)
"""

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import re
import os
import argparse

########################################
########################################

# EXAMPLE TO COMPUTE VELOCITY DISTRIBUTION AT OUTLET
# NOTE : THIS FILE MUST BE IN THE "2d-backward-facing-step" DIRECTORY TO 
#        WORK PROPERLY

# Parse Reynolds number
parser = argparse.ArgumentParser(description='Arguments to compute the velocity distribution at outlet')
parser.add_argument("-Re", "--Reynolds", type=int, help="Reynolds number (only 100 and 1000 supported)", required=True)

args, leftovers = parser.parse_known_args()
Re = args.Reynolds

#Define folder path according to Re
folder = "./Reynolds" + str(Re) +'/'

# VARIABLES
L_out = 50                      # Outlet length
L_in = 15                       # Inlet length
h_in = 0.5                      # Half height of inlet
h_out = 1                       # Half height of outlet
u_moy_in = 1                    # Inflow mean velocity
u_moy_out = 0.5                 # Outflow mean velocity
plt.rc('axes', labelsize=14)
L_out = L_in + L_out

# DATA EXTRACTION
#  Iterate through files in directory, we find the maximum number of files and the last file
max_number = -1
pattern = re.compile(r'backward_facing_step_output\.(\d+)\.00000.vtu')

for filename in os.listdir(folder):
    match = pattern.match(filename)
    if match:
        number = int(match.group(1))
        if number > max_number:
            max_number = number

n_files_to_read = 0
if (Re == 100):
    n_files_to_read = max_number # All of the files
elif (Re == 1000):
    n_files_to_read = 1  # Only the last file

# EXACT SOLUTION
k_out = -3/2*u_moy_out/h_out**2
y_out_an = np.linspace(-h_out, h_out, 1001)
u_out_an = k_out*(y_out_an**2-h_out**2)

# LOOP OVER VTU FILES AND PLOT RELEVANT SOLUTIONS
a = np.array([L_out, 0, 0])
b = np.array([L_out, 2, 0])
plt.figure()

for i in range(0, n_files_to_read):
    if (Re == 100):
        file = (folder + '/backward_facing_step_output.' 
                   + f'{i+1:05d}' + '.00000.vtu')
        plt.title("Re = 100")
    elif (Re == 1000):
        file = (folder + '/backward_facing_step_output.' 
                   + f'{max_number:05d}' + '.00000.vtu')
        plt.title("Re = 1000")

    data = pv.read(file)
    data .set_active_vectors("velocity")

    # Data
    data_outlet = data.sample_over_line(a, b, resolution=1000)
    y_outlet = data_outlet["Distance"]
    u_outlet = data_outlet["velocity"][:,0]
    # Graph
    if (Re == 100):
        plt.plot(u_outlet, y_outlet, label="Mesh " + str(i+1))
    elif (Re == 1000):
        plt.plot(u_outlet, y_outlet, label="Mesh " + str(max_number))

# Plot reference solution
plt.plot(u_out_an, y_outlet, label="Analytical Poiseuille Flow", linestyle='--',
             linewidth=2.0, color='k')

# Graph setup
plt.xlabel(r'$u_{out}$')
plt.ylabel("y/h")
plt.legend(fontsize=8)
plt.axis([0,1.5,0,2])

# Save plot
if (Re == 100):
    plt.savefig('Reynolds100-poiseuille.png')
elif (Re == 1000):
    plt.savefig('Reynolds1000-poiseuille.png')