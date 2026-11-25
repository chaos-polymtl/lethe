# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for the 3D ribbon mixer simulated with a single rotating frame (SRF)
This code reads the torque on the rotating vessel wall at Re = 1, computes the power
number and compares it with the experimental data
"""

# Modules
#-------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import argparse


def read_experimental_data(filename):
  """
  Read the experimental power curves. The file is tab-separated, uses decimal commas
  and contains three (Re, Np) series, some of which have missing values.
  Returns a list of (Re, Np) arrays, one per series.
  """
  with open(filename, encoding="utf-8") as f:
    lines = f.read().splitlines()[1:]

  n_series = 3
  series = [([], []) for _ in range(n_series)]
  for line in lines:
    columns = line.replace(",", ".").split("\t")
    for i in range(n_series):
      if 2*i+1 < len(columns) and columns[2*i].strip() and columns[2*i+1].strip():
        series[i][0].append(float(columns[2*i]))
        series[i][1].append(float(columns[2*i+1]))
  return [(np.array(Re), np.array(Np)) for Re, Np in series]


#--------------------------------------------
# Main
#--------------------------------------------

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the 3D ribbon mixer simulated with a single rotating frame')
parser.add_argument("-f", "--folder", type=str, help="Folder which contains the output files of the simulation", default="./output")
parser.add_argument("-e", "--experimental", type=str, help="File which contains the experimental data", default="../Np_vs_Re/experimental-data.dat")
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
args, leftovers = parser.parse_known_args()

# Operating conditions of ribbon-gls-Re1.prm
D = 0.27                        # Impeller diameter
N = 10 / 2 / np.pi              # Rotation speed in rotations per second (omega_z = 10 rad/s)
nu = 0.11602395351399172        # Kinematic viscosity
Re = N * D**2 / nu

# Torque on the rotating vessel wall (boundary id 1) at the end of the steady simulation.
# Lethe solves with a unit density, so the torque is already divided by the density.
T_z = np.loadtxt(args.folder + "/torque.01.dat", skiprows=1, usecols=(3), ndmin=1)[-1]
Np = 2 * np.pi * np.abs(T_z) / (N**2 * D**5)

# Experimental power number at the same Reynolds number, interpolated in log-log
experimental = read_experimental_data(args.experimental)
Np_experimental = np.array([np.exp(np.interp(np.log(Re), np.log(Re_exp), np.log(Np_exp)))
                            for Re_exp, Np_exp in experimental])
relative_difference = (Np - Np_experimental) / Np_experimental

print(f"Re = {Re:.6f}, T_z = {T_z:.10e}, Np = {Np:.6f}, Np*Re = {Np*Re:.6f}")
for i, (Np_exp, difference) in enumerate(zip(Np_experimental, relative_difference)):
  print(f"Experiment {i+1}: Np = {Np_exp:.6f}, relative difference = {difference:+.4f}")

# Plot the power curves around the simulated Reynolds number
fig0 = plt.figure()
markers = ['o', '*', '^']
for i, (Re_exp, Np_exp) in enumerate(experimental):
  plt.loglog(Re_exp, Np_exp, markers[i], mfc='None', label=f"Experiment {i+1}")
plt.loglog(Re, Np, 's', color="black", label="Lethe (SRF)")
plt.xlim(1e-2, 1e2)
plt.xlabel("$Re$")
plt.ylabel("$N_p$")
plt.legend()
plt.tight_layout()

if (not args.validate):
  fig0.savefig("ribbon-mixer-srf-power-number.png", dpi=300)
  plt.show()

else:
  # Store the plotted data so that they can be compared between validation runs
  np.savetxt("solution-power-number.dat",
             np.column_stack(([Re], [T_z], [Np], [Np*Re], [Np_experimental], [relative_difference])).reshape(1, -1),
             header="Re T_z Np Np*Re Np_experiment_1 Np_experiment_2 Np_experiment_3 relative_difference_1 relative_difference_2 relative_difference_3")
  fig0.savefig("ribbon-mixer-srf-power-number.pdf")
