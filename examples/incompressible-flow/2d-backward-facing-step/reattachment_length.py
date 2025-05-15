# SPDX-FileCopyrightText: Copyright (c) 2022, 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for 2D-backward-facing-step example
Computes the reattachment length (x_r) for Re = 100 and Re = 1000
for several meshes using a bisection algorithm
"""

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import re
import os
import argparse
import sys

########################################
########################################

# EXAMPLE TO COMPUTE x_r (REATTACHMENT LENGTH) FOR Re = 100 and Re = 1000
# NOTE : THIS FILE MUST BE IN THE "2d-backward-facing-step" DIRECTORY TO 
#        WORK PROPERLY

# Parse Reynolds number
parser = argparse.ArgumentParser(description='Arguments to compute the velocity distribution at outlet')
parser.add_argument("-Re", "--Reynolds", type=int, help="Reynolds number (only 100 and 1000 supported)", required=True)

args, leftovers = parser.parse_known_args()
Re = args.Reynolds

if (Re != 100 and Re != 1000):
    sys.exit("This Reynolds number is not supported. Only Re = 100 and Re = 1000.")

#Define folder path according to Re
folder = "./Reynolds" + str(Re)

# VARIABLES
L_in = 15           # Inlet length


# Initial guesses depending on the Reynolds number and reference reattachment values
# by Erturk (2008) (see attached .txt file for other Reynolds x_r values)
if (Re == 100):
    x_inf = 1.5
    x_sup = 3
    x_r_ref = 2.922    
elif (Re == 1000):
    x_inf = 12
    x_sup = 14
    x_r_ref = 13.121 

tol = 1e-12                  # Bisection stop criterion

# DATA EXTRACTION
max_number = -1
pattern = re.compile(r'backward_facing_step_output\.(\d+)\.00000.vtu')

for filename in os.listdir(folder):
    match = pattern.match(filename)
    if match:
        number = int(match.group(1))
        if number > max_number:
            max_number = number

n_files_to_read = max_number
x_r = []

for i in range(0, n_files_to_read):
    # Read all files except the 0 file
    file = (folder + '/backward_facing_step_output.' 
                   + f'{i+1:05d}' + '.00000.vtu')
    data = pv.read(file)
    data.set_active_vectors("velocity")

    # Initial guesses
    x1 = L_in + x_inf
    x2 = L_in + x_sup
    j = 0
    # Bisection loop
    while abs(x2-x1)>tol and j<=20:
        xm = (x1+x2)/2

        # Profiles extraction
        a = np.array([x1, 0, 0])
        b = np.array([x1, 0.01, 0])
        profil1 = data.sample_over_line(a, b, resolution=1000)
        u1 = profil1["velocity"][:,0]
        y = profil1["Distance"]

        a = np.array([xm, 0, 0])
        b = np.array([xm, 0.01, 0])
        profilm = data.sample_over_line(a, b, resolution=1000)
        um = profilm["velocity"][:,0]

        a = np.array([x2, 0, 0])
        b = np.array([x2, 0.01, 0])
        profil2 = data.sample_over_line(a, b, resolution=1000)
        u2 = profil2["velocity"][:,0]


        # Derivatives
        dudy1 = (-u1[2] + 4*u1[1] - 3*u1[0]) / (y[2] - y[0])
        dudym = (-um[2] + 4*um[1] - 3*um[0]) / (y[2] - y[0])
        dudy2 = (-u2[2] + 4*u2[1] - 3*u2[0]) / (y[2] - y[0])

        # Signs check
        if (abs(dudym) < 0.1*tol):
            break
        if (dudy1*dudym < 0):
            x2 = xm
        else:
            x1 = xm

        j+=1

    x_r.append(xm - L_in)   # Add distance from the step
    print("Output file " + str(i+1) + ", reattachment point x_r: " , xm - L_in)

# Calculate relative error using last reattachment point and ref values
e = abs(x_r_ref - x_r[-1])/x_r_ref
print("Reference value Erturk (2008): ", x_r_ref)
print("Relative error abs(x_r_ref - x_r)/x_r_ref = ", e)
