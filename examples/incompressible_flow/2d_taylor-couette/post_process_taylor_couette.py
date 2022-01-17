"""
Postprocessing code for 2D lid-driven cavity
This code extracts the u velocity along the y axis at x=0.5 and compares
it to the results of Ghia et al.

This comparison is similar to what is done in the step-57 of deal.II

Author : Bruno Blais
"""

# Modules
#-------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import cm
import pyvista as pv

import glob

import os
import sys

#--------------------------------------------
# Main
#--------------------------------------------

# Load VTU file
vtu_file="couette.0001.0000.vtu"
sim = pv.read(vtu_file)
sim.set_active_vectors("velocity")

# Create begin and end point of line
a = [0.25, 0, 0]
b = [1., 0, 0]

# Extract all field over the line using pyvista
sampled_data=sim.sample_over_line(a, b, resolution=500)

# Get u component of the velocity from sampled data
r = sampled_data["Distance"]+0.25
u = sampled_data["velocity"][:,1]

# Calculate analytical solution
omega = 1
kappa = 0.25
r_a=np.linspace(0.25,1,10)
R = 1 
mu = 1

u_a = omega * R * kappa * (R / r_a - r_a /R) / (1/kappa - kappa)

torque = 4*np.pi*mu*omega*R**2*(kappa**2/(1-kappa**2))

print("The torque is : ", torque)

plt.plot(r,u,'-',label="Lethe")
plt.plot(r_a,u_a,'s',ms=8,mfc="None",markeredgecolor="black",label="Analytical solution")
plt.ylabel("$u_{\\theta}$")
plt.xlabel("R")
plt.legend()
plt.tight_layout()
plt.savefig("lethe_analytical_taylor_couette_comparison.png", dpi=300)
plt.show()
