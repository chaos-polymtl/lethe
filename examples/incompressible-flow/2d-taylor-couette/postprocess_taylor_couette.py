"""
Postprocessing code for the 2D taylor-couette example
This code extracts the data from the vtu file and plots it against the analytical solution

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
plt.savefig("lethe-analytical-taylor-couette-comparison.png", dpi=300)
plt.show()
