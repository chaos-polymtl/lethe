# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for 2D-sudden-expansion-flow example
Computes velocity distributions at different cross sections and
compares it with literature results and analytical solution (Poiseuille)
"""

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import re
import os
import argparse

########################################
########################################

# NOTE: This file must be in the "2d-sudden-expansion-flow" directory to work properly

# Parse Reynolds number
parser = argparse.ArgumentParser(description='Arguments to compute the velocity distribution at outlet')
parser.add_argument("-Re", "--Reynolds", type=int, help="Reynolds number (only 70 and 610 supported)", required=True)

args, leftovers = parser.parse_known_args()
Re = args.Reynolds

if (Re != 70 and Re != 610):
    sys.exit("This Reynolds number is not supported. Only Re = 70 and Re = 610.")

#Define folder path according to Re
folder = "./Reynolds" + str(Re) +'/'

# Variables
L_out = 1.0                     # Outlet length
L_in = 0.08                     # Inlet length
L_out = L_in + L_out            # Total length
h = 0.01                        # constriction diameter

# Plotting text size
plt.rc('axes', labelsize=14)
plt.rc('axes', titlesize=14)
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('legend', fontsize=14)

# Load reference data
if (Re == 70):
    # Durst et al, Re=70, numerical
    ref_Re70_data = [f"./Reynolds70/Reference-Data/Durst-Re70-{x}mm.txt"
                     for x in [10, 35, 50, 140]]
    ref_Re70 = [np.loadtxt(file, delimiter=",", skiprows=1) for file in ref_Re70_data]
    # Durst et al, Re=70, experimental
    refexp_Re70_data = [f"./Reynolds70/Reference-Data/Durst-Re70exp-{x}mm.txt"
                    for x in [10, 35, 50, 140]]
    refexp_Re70 = [np.loadtxt(file, delimiter=",", skiprows=1) for file in refexp_Re70_data]
    # Kanna et al, Re=70, numerical
    ref2_Re70_data = [f"./Reynolds70/Reference-Data/Kanna-Re70-{x}mm.txt"
                        for x in [35, 50, 140]]
    ref2_Re70 = [np.loadtxt(file, delimiter=",", skiprows=1) for file in ref2_Re70_data]

if (Re == 610):
    # Durst et al, Re=610, numerical
    ref_Re610_data = [f"./Reynolds610/Reference-Data/Durst-Re610-{x}mm.txt"
                    for x in [10, 35, 50, 140, 200, 290, 440, 600]]
    ref_Re610 = [np.loadtxt(file, delimiter=",", skiprows=1) for file in ref_Re610_data]
    # Durst et al, Re=610, experimental
    refexp_Re610_data = [f"./Reynolds610/Reference-Data/Durst-Re610exp-{x}mm.txt"
                    for x in [10, 35, 50, 140, 200, 290, 440, 600]]
    refexp_Re610 = [np.loadtxt(file, delimiter=",", skiprows=1) for file in refexp_Re610_data]

# Data extraction
# Iterate through files in directory, we find the maximum number of files and the last file
max_number = -1
pattern = re.compile(r'sudden-expansion-flow-output\.(\d+)\.00000.vtu')

for filename in os.listdir(folder):
    match = pattern.match(filename)
    if match:
        number = int(match.group(1))
        if number > max_number:
            max_number = number

u_moy_in = 1            # Inflow mean velocity
u_moy_out = 0.5         # Outflow mean velocity
u_max_Re70 = 1.48635    # Maximum velocity at Re=70
u_max_Re610 = 1.32786   # Maximum velocity at Re=610

# Exact solution
k_out = -3/2*u_moy_out/h**2
y_out_an = np.linspace(-h, h, 1001)
u_out_an = k_out*(y_out_an**2-h**2)

# Velocity profiles along the domain
file = (folder + '/sudden-expansion-flow-output.' 
        + f'{max_number:05d}' + '.00000.vtu')
if (Re == 70):
    a = np.array([[L_in-h,        0,    0],
                  [L_in+3.5*h,    0,    0],
                  [L_in+5*h,      0,    0],
                  [L_in+14*h,     0,    0]])
    b = np.array([[L_in-h,        2*h,  0],
                  [L_in+3.5*h,    2*h,  0],
                  [L_in+5*h,      2*h,  0],
                  [L_in+14*h,     2*h,  0]])
elif (Re == 610):
    a = np.array([[L_in-h,        0,    0],
                  [L_in+3.5*h,    0,    0],
                  [L_in+5*h,      0,    0],
                  [L_in+14*h,     0,    0],
                  [L_in+20*h,     0,    0],
                  [L_in+29*h,     0,    0],
                  [L_in+44*h,     0,    0],
                  [L_in+60*h,     0,    0]])    
    b = np.array([[L_in-h,        2*h,  0],
                  [L_in+3.5*h,    2*h,  0],
                  [L_in+5*h,      2*h,  0],
                  [L_in+14*h,     2*h,  0],
                  [L_in+20*h,     2*h,  0],
                  [L_in+29*h,     2*h,  0],
                  [L_in+44*h,     2*h,  0],
                  [L_in+60*h,     2*h,  0]])

data = pv.read(file)
data .set_active_vectors("velocity")

for row in range(a.shape[0]):
    plt.figure()
    # Data
    data_cross_section = data.sample_over_line(a[int(row),:], b[int(row),:], resolution=1000)
    y = data_cross_section["Distance"]
    u = data_cross_section["velocity"][:,0]
    
    if(Re == 70):
        # Graph
        plt.plot(y, u/u_max_Re70, label="Lethe", linewidth=2.0, color='c')
        # Plot reference solutions
        plt.plot(ref_Re70[int(row)][:,0], ref_Re70[int(row)][:,1], label="Durst (1993) num.", 
                 marker='s', markerfacecolor='none', linestyle='', linewidth=2.0, color='k')
        plt.plot(refexp_Re70[int(row)][:,0], refexp_Re70[int(row)][:,1], label="Durst (1993) exp.", 
                 marker='o', markerfacecolor='none', linestyle='', linewidth=2.0, color='r')
        if(int(row) != 0):
            plt.plot(ref2_Re70[int(row)-1][:,0], ref2_Re70[int(row)-1][:,1], label="Kanna (2005) num.", 
                 marker='^', markerfacecolor='none', linestyle='', linewidth=2.0, color='b')

        plt.xticks(rotation=90)
        plt.xlabel("y [m]")
        plt.ylabel(r'$u/u_{max}$')
        plt.legend()
        plt.title(f'Re = 70, x = {int(a[int(row),0]*1e3)}mm')
        plt.axis([0,2*h,0,1])
        plt.savefig(f'Reynolds70-{row}.png', bbox_inches='tight')
    elif(Re == 610):
        # Graph
        plt.plot(y, u/u_max_Re610, label="Lethe", linewidth=2.0, color='c')
        # Plot reference solutions
        plt.plot(ref_Re610[int(row)][:,0], ref_Re610[int(row)][:,1], label="Durst (1993) num.", 
                 marker='s', markerfacecolor='none', linestyle='', linewidth=2.0, color='k')
        plt.plot(refexp_Re610[int(row)][:,0], refexp_Re610[int(row)][:,1], label="Durst (1993) exp.", 
                 marker='o', markerfacecolor='none', linestyle='', linewidth=2.0, color='r')
        plt.xticks(rotation=90)
        plt.xlabel("y [m]")
        plt.ylabel(r'$u/u_{max}$')
        plt.legend()
        plt.title(f'Re = 610, x = {int(a[int(row),0]*1e3)}mm')
        plt.axis([0,2*h,-0.25,1])
        plt.savefig(f'Reynolds610-{row}.png', bbox_inches='tight')

# Velocity profile at outlet
# Comparison with Poiseuille analytical solution
a_outlet = np.array([L_out, 0,  0])
b_outlet = np.array([L_out, 2*h, 0])
plt.figure()

# Data
data_cross_section = data.sample_over_line(a_outlet, b_outlet, resolution=1000)
y_outlet = data_cross_section["Distance"]
u_outlet = data_cross_section["velocity"][:,0]
# Graph
plt.plot(y_outlet, u_outlet, label="Lethe", linewidth=2.0, color='c')
# Plot reference solution
plt.plot(y_out_an+h, u_out_an, label="Analytical Poiseuille Flow", linestyle='--',
         linewidth=2.0, color='k')
# Graph setup
plt.xticks(rotation=90)
plt.xlabel("y [m]")
plt.ylabel(r'$u$')
plt.legend()
plt.axis([0,2*h,0,1])


# Save plot
if (Re == 70):
    plt.savefig('Reynolds70-poiseuille.png', bbox_inches='tight')
elif (Re == 610):
    plt.savefig('Reynolds610-poiseuille.png', bbox_inches='tight')
