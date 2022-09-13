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
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import cm
import pyvista as pv

import glob

import os
import sys
from pyvista.plotting.renderer import CameraPosition

#--------------------------------------------
# Main
#--------------------------------------------

# Load reference data
# Note, this reference data is taken from step-57 of deal.II
ref_data_ghia = "../ref-2d-ghia-u.txt"
raw_ghia=np.loadtxt(ref_data_ghia,skiprows=2)
y_ghia=raw_ghia[:,0]
u_ghia=raw_ghia[:,6]

ref_data_erturk = "../ref-2d-erturk-u.txt"
raw_erturk=np.loadtxt(ref_data_erturk,skiprows=2)
y_erturk=raw_erturk[:,0]
u_erturk=raw_erturk[:,4]

# Load VTU file
vtu_file="out.0320.0000.vtu"
sim = pv.read(vtu_file)
sim.set_active_vectors("velocity")

# Create begin and end point of line
a = [0.5, 0, 0]
b = [0.5, 1, 0]

# Extract all field over the line using pyvista
sampled_data=sim.sample_over_line(a, b, resolution=1000)

# Get u component of the velocity from sampled data
y = sampled_data["Distance"]
u = sampled_data["velocity"][:,0]

plt.plot(u,y,label="Lethe")
plt.plot(u_ghia,y_ghia,'s',ms=8,mfc="None",markeredgecolor="black",label="Ghia 1982")
plt.plot(u_erturk,y_erturk,'o',ms=8,mfc="None",markeredgecolor="black",label="Erturk 2005")

plt.xlabel("u")
plt.ylabel("y")
plt.legend()
plt.savefig("lethe-ghia-re-7500-comparison.png", dpi=300)
plt.show()


plt.plot(u,y,label="Lethe")
plt.plot(u_ghia,y_ghia,'s',ms=8,mfc="None",markeredgecolor="black",label="Ghia 1982")
plt.plot(u_erturk,y_erturk,'o',ms=8,mfc="None",markeredgecolor="black",label="Erturk 2005")

plt.xlabel("u")
plt.ylabel("y")
plt.xlim([0.2,1])
plt.ylim([0.9,1])
plt.legend()
plt.savefig("lethe-ghia-re-7500-comparison-zoom.png", dpi=300)
plt.show()
