# SPDX-FileCopyrightText: Copyright (c) 2022 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for the 2D taylor-couette example
This code extracts the data from the vtu file and plots it against the analytical solution
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

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the 2d Taylor-Couette flow')
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
args, leftovers=parser.parse_known_args()

# Load VTU file
vtu_file="couette.00001.00000.vtu"
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


plt.plot(r,u,'-',label="Lethe")
plt.plot(r_a,u_a,'s',ms=8,mfc="None",markeredgecolor="black",label="Analytical solution")
plt.ylabel("$u_{\\theta}$")
plt.xlabel("R")
plt.legend()
plt.tight_layout()

if (not args.validate):
  print("The torque is : ", torque)
  plt.savefig("lethe-analytical-taylor-couette-comparison.png", dpi=300)
  plt.show()

else:
  # Combine the vectors into an array (columns)
  data = np.column_stack((r,u))

  # Output file name
  output_file = "solution.dat"

  # Save the data to a .dat file with space as the separator
  np.savetxt(output_file, data, fmt="%.6f",header="r u", delimiter=" ") 
  plt.savefig("lethe-analytical-taylor-couette-comparison.pdf")

