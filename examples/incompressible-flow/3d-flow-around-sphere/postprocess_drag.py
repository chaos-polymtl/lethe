# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for the 3D flow around a sphere
This code reads the force on the sphere for every simulation, computes the drag
coefficient and compares it with the reference values at Re = 0.1 and Re = 150
"""

# Modules
#-------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

#--------------------------------------------
# Main
#--------------------------------------------

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the 3D flow around a sphere')
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
args, leftovers = parser.parse_known_args()

# Flow conditions of the parameter files
U = 1                           # Inlet velocity
D = 1                           # Sphere diameter
A = np.pi * D**2 / 4            # Frontal area of the sphere

# Reference drag coefficients: Stokes law (24/Re) at Re = 0.1 and the correlation
# quoted in the documentation of this example at Re = 150
reynolds = [0.1, 150]
C_D_reference = {0.1: 240., 150: 0.889}

# Output folder of every simulation (see the output path of the parameter files)
cases = [("matrix-based", ""), ("matrix-free", "-matrix-free")]

Re_list, solver_list, cells_list, C_D_list, error_list = [], [], [], [], []
for Re in reynolds:
  for solver, suffix in cases:
    folder = f"sphere-{Re:g}{suffix}"
    filename = os.path.join(folder, "force.00.dat")
    if not os.path.isfile(filename):
      print(f"Warning: {filename} does not exist, the {solver} case at Re = {Re:g} is skipped")
      continue
    # Force on the sphere (boundary id 0) at the end of the simulation.
    # Lethe solves with a unit density, so the force is already divided by the density.
    cells, f_x = np.loadtxt(filename, skiprows=1, usecols=(0, 1), ndmin=2)[-1]
    C_D = 2 * f_x / (U**2 * A)
    error = (C_D - C_D_reference[Re]) / C_D_reference[Re]
    print(f"Re = {Re:g}, {solver}: {int(cells)} cells, C_D = {C_D:.6f}, "
          f"reference = {C_D_reference[Re]:g}, relative difference = {error:+.4f}")
    Re_list.append(Re); solver_list.append(solver); cells_list.append(cells)
    C_D_list.append(C_D); error_list.append(error)

# Plot the drag coefficient of every simulation, normalized by the reference value
fig0 = plt.figure()
markers = {"matrix-based": 'o', "matrix-free": 's'}
for solver, _ in cases:
  indices = [i for i, s in enumerate(solver_list) if s == solver]
  plt.semilogx([Re_list[i] for i in indices],
               [C_D_list[i] / C_D_reference[Re_list[i]] for i in indices],
               markers[solver], mfc='None', ms=10, label=f"Lethe, {solver}")
plt.axhline(1, color='k', ls='--', label="Reference")
plt.xlabel("$Re$")
plt.ylabel(r"$C_D / C_{D,\mathrm{ref}}$")
plt.legend()
plt.tight_layout()

if (not args.validate):
  fig0.savefig("lethe-sphere-drag-coefficient.png", dpi=300)
  plt.show()

else:
  # Store the plotted data so that they can be compared between validation runs.
  # The solver is stored as 0 for matrix-based and 1 for matrix-free.
  solver_index = [0 if s == "matrix-based" else 1 for s in solver_list]
  np.savetxt("solution-drag-coefficient.dat",
             np.column_stack((Re_list, solver_index, cells_list, C_D_list,
                              [C_D_reference[Re] for Re in Re_list], error_list)),
             header="Re solver(0:matrix-based,1:matrix-free) cells C_D C_D_reference relative_difference")
  fig0.savefig("lethe-sphere-drag-coefficient.pdf")
