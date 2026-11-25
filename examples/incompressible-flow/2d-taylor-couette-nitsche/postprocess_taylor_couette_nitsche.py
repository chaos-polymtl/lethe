# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for the 2D Taylor-Couette flow with a Nitsche immersed boundary
This code reads the torques on the inner and outer cylinders and the L2 error of the
velocity for every mesh, and compares the torque against the analytical solution
"""

# Modules
#-------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import argparse

#--------------------------------------------
# Main
#--------------------------------------------

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the 2d Taylor-Couette flow with a Nitsche immersed boundary')
parser.add_argument("-f", "--folder", type=str, help="Folder which contains the output files of the simulation", default="./")
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
args, leftovers=parser.parse_known_args()

folder = args.folder

# Torque on the inner cylinder (Nitsche solid 0) and on the outer cylinder (boundary 0)
cells, T_z_solid = np.loadtxt(folder + "/torque_solid_00.dat", skiprows=1, usecols=(0,3), unpack=True)
T_z_wall = np.loadtxt(folder + "/torque.00.dat", skiprows=1, usecols=(3))

# L2 error of the velocity with respect to the analytical solution
error_velocity = np.loadtxt(folder + "/L2Error.dat", skiprows=1, usecols=(1))

# Analytical torque on the inner cylinder
omega = 1
kappa = 0.25
R = 1
mu = 1
torque = -4*np.pi*mu*omega*R**2*(kappa**2/(1-kappa**2))
relative_error = np.abs((T_z_solid-torque)/torque)

print("The analytical torque on the inner cylinder is : ", torque)
for n, T in zip(cells, T_z_solid):
  print(f"{int(n):8d} cells: T_z = {T:.6f}")

# Plot the torques
fig0 = plt.figure()
plt.semilogx(cells, -T_z_solid, 'o-', label="Lethe, inner cylinder")
plt.semilogx(cells, T_z_wall, 's-', label="Lethe, outer cylinder")
plt.axhline(-torque, color='k', ls='--', label="Analytical solution")
plt.xlabel("Number of cells")
plt.ylabel("$|T_z|$")
plt.legend()
plt.tight_layout()

# Plot the L2 error of the velocity, with a first-order slope (h is proportional to cells^-1/2)
fig1 = plt.figure()
plt.loglog(cells, error_velocity, 'o-', label="Lethe")
plt.loglog(cells, error_velocity[0]*np.sqrt(cells[0]/cells), 'k--', label="First order")
plt.xlabel("Number of cells")
plt.ylabel("L2 error of the velocity")
plt.legend()
plt.tight_layout()

if (not args.validate):
  fig0.savefig("lethe-taylor-couette-nitsche-torque.png", dpi=300)
  fig1.savefig("lethe-taylor-couette-nitsche-error.png", dpi=300)
  plt.show()

else:
  # Store the plotted data so that they can be compared between validation runs
  np.savetxt("solution-torque-error.dat",
             np.column_stack((cells, T_z_solid, T_z_wall, relative_error, error_velocity)),
             header="cells T_z_solid T_z_wall T_z_solid_relative_error error_velocity")
  fig0.savefig("lethe-taylor-couette-nitsche-torque.pdf")
  fig1.savefig("lethe-taylor-couette-nitsche-error.pdf")
