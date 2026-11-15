# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for the 3D flow around a sphere
This code reads the force on the sphere for every simulation and computes the drag
coefficient. It compares it with the reference values at Re = 0.1 and Re = 150 and,
at Re = 0.1, extrapolates the drag coefficient to an infinite domain from the
simulations on the domains enlarged twice and four times.
"""

# Modules
#-------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

# Plot font and colors
font = {'weight' : 'normal',
        'size'   : 13}
plt.rc('font', **font)
colors = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']

#--------------------------------------------
# Functions
#--------------------------------------------

def drag_coefficient(folder):
  """
  Return the number of cells and the drag coefficient at the end of the simulation
  whose results are in folder, or None if the simulation did not produce them.
  The force on the sphere (boundary id 0) is divided by the density, since Lethe
  solves with a unit density.
  """
  filename = os.path.join(folder, "force.00.dat")
  if not os.path.isfile(filename):
    print(f"Warning: {filename} does not exist, this simulation is skipped")
    return None
  cells, f_x = np.loadtxt(filename, skiprows=1, usecols=(0, 1), ndmin=2)[-1]
  return cells, 2 * f_x / (U**2 * A)

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

# Reference drag coefficients: Stokes law (24/Re) and Oseen's correction
# (24/Re (1 + 3/16 Re)) at Re = 0.1, and the correlation quoted in the
# documentation of this example at Re = 150
C_D_stokes = 240.
C_D_oseen = 24 / 0.1 * (1 + 3 / 16 * 0.1)
C_D_reference = {0.1: C_D_stokes, 150: 0.889}

# Simulations: (Reynolds number, solver, domain scale, output folder).
# The solver is 0 for matrix-based and 1 for matrix-free.
simulations = [(0.1, 0, 1, "sphere-0.1"),
               (0.1, 1, 1, "sphere-0.1-matrix-free"),
               (0.1, 1, 2, "sphere-0.1-matrix-free-x2"),
               (0.1, 1, 4, "sphere-0.1-matrix-free-x4"),
               (150, 0, 1, "sphere-150"),
               (150, 1, 1, "sphere-150-matrix-free")]
solver_names = ["matrix-based", "matrix-free"]

results = []
for Re, solver, scale, folder in simulations:
  result = drag_coefficient(folder)
  if result is None:
    continue
  cells, C_D = result
  error = (C_D - C_D_reference[Re]) / C_D_reference[Re]
  print(f"Re = {Re:g}, {solver_names[solver]}, domain scale {scale}: {int(cells)} cells, "
        f"C_D = {C_D:.6f}, reference = {C_D_reference[Re]:g}, relative difference = {error:+.4f}")
  results.append((Re, solver, scale, cells, C_D, C_D_reference[Re], error))
results = np.array(results)

# Domain size study at Re = 0.1 with the matrix-free solver. The disturbance of a
# Stokes flow decays like the inverse of the distance to the sphere, so the error
# due to the confinement of the domain decreases like 1/domain_scale. The drag
# coefficient is extrapolated linearly in 1/domain_scale to an infinite domain from
# the two largest domains.
study = results[(results[:, 0] == 0.1) & (results[:, 1] == 1)] if results.size else np.empty((0, 7))
study = study[np.argsort(study[:, 2])]
inverse_scale = 1 / study[:, 2]
C_D_study = study[:, 4]
extrapolation = None
if np.all(np.isin([2, 4], study[:, 2])):
  C_D_x2 = C_D_study[study[:, 2] == 2][0]
  C_D_x4 = C_D_study[study[:, 2] == 4][0]
  C_D_extrapolated = 2 * C_D_x4 - C_D_x2
  extrapolation = (C_D_extrapolated, (C_D_extrapolated - C_D_stokes) / C_D_stokes)
  print(f"Re = 0.1, C_D extrapolated to an infinite domain = {C_D_extrapolated:.6f}, "
        f"relative difference with Stokes law = {extrapolation[1]:+.4f}")

# Figure 1: drag coefficient as a function of the size of the domain at Re = 0.1
fig0 = plt.figure()
x_line = np.linspace(0, 1.05, 2)
plt.plot(x_line, C_D_stokes * np.ones(2), lw=2, ls='--', color=colors[0], label="Stokes law")
plt.plot(x_line, C_D_oseen * np.ones(2), lw=2, ls='--', color=colors[1], label="Oseen correction")
if extrapolation is not None:
  # Line through the two largest domains, down to an infinite domain (1/scale = 0)
  x_fit = np.array([0, 0.5])
  plt.plot(x_fit, C_D_extrapolated + (C_D_x2 - C_D_extrapolated) / 0.5 * x_fit, lw=2, ls=':', color=colors[2])
  plt.plot(0, C_D_extrapolated, '*', ms=14, color=colors[2], label="Extrapolation")
plt.plot(inverse_scale, C_D_study, 'ko-', lw=2, ms=7, label="Lethe")
plt.xlim(-0.05, 1.05)
plt.xlabel(r"$1 / \mathrm{domain\ scale}$ [-]")
plt.ylabel(r"$C_D$ [-]")
plt.title(r"$\mathrm{Re}=0.1$, matrix-free Q2-Q2")
plt.legend(loc='upper left')
plt.tight_layout()

# Figure 2: drag coefficient of every simulation on the base domain, normalized
# by the reference value
fig1 = plt.figure()
markers = ['o', 's']
for solver in (0, 1):
  selected = results[(results[:, 1] == solver) & (results[:, 2] == 1)] if results.size else np.empty((0, 7))
  plt.semilogx(selected[:, 0], selected[:, 4] / selected[:, 5], markers[solver], mfc='None', mew=2, ms=10,
               color=colors[solver], label=f"Lethe, {solver_names[solver]}")
if extrapolation is not None:
  plt.semilogx(0.1, C_D_extrapolated / C_D_stokes, '*', ms=14, color=colors[2], label="Lethe, extrapolated")
plt.axhline(1, color='k', ls='--', lw=2, label="Reference")
plt.xlabel(r"$\mathrm{Re}$ [-]")
plt.ylabel(r"$C_D / C_{D,\mathrm{ref}}$ [-]")
plt.legend()
plt.tight_layout()

if (not args.validate):
  fig0.savefig("lethe-sphere-drag-domain-size.png", dpi=300)
  fig1.savefig("lethe-sphere-drag-coefficient.png", dpi=300)
  plt.show()

else:
  # Store the plotted data so that they can be compared between validation runs
  np.savetxt("solution-drag-coefficient.dat", results,
             header="Re solver(0:matrix-based,1:matrix-free) domain_scale cells C_D C_D_reference relative_difference")
  if extrapolation is not None:
    np.savetxt("solution-drag-coefficient-extrapolated.dat", np.array([extrapolation]),
               header="C_D_extrapolated relative_difference_with_stokes")
  fig0.savefig("lethe-sphere-drag-domain-size.pdf")
  fig1.savefig("lethe-sphere-drag-coefficient.pdf")
