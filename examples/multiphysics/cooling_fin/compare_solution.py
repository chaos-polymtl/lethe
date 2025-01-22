# SPDX-FileCopyrightText: Copyright (c) 2022, 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for cylindrical fin
"""

# Modules
#-------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import argparse


#--------------------------------------------
# Main
#--------------------------------------------

parser = argparse.ArgumentParser(description='Arguments for the validation of the fin')
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
args, leftovers=parser.parse_known_args()

#--------------------------------------------
# Analytical solution
#--------------------------------------------
h = 10
k = 100
L = 0.2
R = 0.01
P = 2*np.pi*R
A = np.pi*R**2
T_inf = 20
T_b = 100
x_analytical = np.linspace(0,L,100)
m = np.sqrt(P*h/k/A)
T_analytical = T_inf + (T_b-T_inf) * np.cosh(m*(L-x_analytical))/np.cosh(m*L)




#--------------------------------------------
# Simulation Data
# Load VTU file
vtu_file="output/out.00004.00000.vtu"
sim = pv.read(vtu_file)
sim.set_active_scalars("temperature")

# Create begin and end point of line
a = [-0.1, 0, 0]
b = [0.1, 0, 0]

# Extract all field over the line using pyvista
sampled_data=sim.sample_over_line(a, b, resolution=50)

# Get u component of the velocity from sampled data
x = sampled_data["Distance"]
T = sampled_data["temperature"][:]

plt.plot(x,T,'o',mfc='None',label="Lethe")
plt.plot(x_analytical,T_analytical,color="black",label="Analytical solution")

plt.xlabel("x")
plt.ylabel("T")
plt.legend()
plt.show()