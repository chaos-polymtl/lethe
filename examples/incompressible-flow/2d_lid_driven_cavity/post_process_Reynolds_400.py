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
import pyvista as pv

#--------------------------------------------
# Main
#--------------------------------------------

# Load reference data
# Note, this reference data is taken from step-57 of deal.II
ref_data = "ref-2d-ghia-u.txt"
raw_ghia=np.loadtxt(ref_data,skiprows=2)
y_ghia=raw_ghia[:,0]
u_ghia=raw_ghia[:,2]

# Load VTU file
vtu_file="out.0001.0000.vtu"
sim = pv.read(vtu_file)
sim.set_active_vectors("velocity")

# Create begin and end point of line
a = [0.5, 0, 0]
b = [0.5, 1, 0]

# Extract all field over the line using pyvista
sampled_data=sim.sample_over_line(a, b, resolution=100)

# Get u component of the velocity from sampled data
y = sampled_data["Distance"]
u = sampled_data["velocity"][:,0]

plt.plot(u,y,label="Lethe")
plt.plot(u_ghia,y_ghia,'s',ms=8,mfc="None",markeredgecolor="black",label="Ghia 1982")
plt.xlabel("u")
plt.ylabel("y")
plt.legend()
plt.savefig("lethe-ghia-re-400-comparison.png", dpi=300)
plt.show()
