# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for the sedimentation of one particle example.
This code compares the velocity of the particle with the one measured by
ten Cate et al. [1] for the E4 experiment, as plotted by Barbeau et al. [2].
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Plot colors
colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the sedimentation of one particle')
parser.add_argument("-f", "--folder", type=str, help="Folder which contains the results of the simulation", default="out")
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
args, leftovers = parser.parse_known_args()

# Velocity and position of the particle written by Lethe
particle = pd.read_csv(f"{args.folder}/ib_force.00.dat", sep=r"\s+")
t = particle["time"].to_numpy()
v_y = particle["v_y"].to_numpy()
p_y = particle["p_y"].to_numpy()

# Experimental data of ten Cate et al. [1], as plotted by Barbeau et al. [2]
t_exp, v_y_exp = np.loadtxt("ten-cate-E4.dat", unpack=True)

plt.rcParams['font.size'] = '14'
fig0 = plt.figure(figsize=(7, 8))
ax0 = fig0.add_subplot(111)
ax0.plot(t, v_y, color=colors[0], lw=2, label="Lethe")
ax0.plot(t_exp, v_y_exp, 'o', mfc='none', color='k', label="Experiment - ten Cate et al. (2002)")
ax0.set_xlabel("Time [s]")
ax0.set_ylabel(r"$v_y$ [cm/s]")
ax0.legend(loc="upper right")
plt.tight_layout()

if (not args.validate):
    plt.savefig("velocity-comparison.png", dpi=300)
    plt.show()
else:
    plt.savefig("sedimentation-velocity.pdf")
    # Store the velocity and the position of the particle, as well as its
    # velocity at the times of the experiment, so that they can be compared
    # between validation runs
    np.savetxt("solution-particle.dat",
               np.column_stack((t, v_y, p_y)),
               header="t v_y p_y")
    np.savetxt("solution-velocity-experiment.dat",
               np.column_stack((t_exp, v_y_exp, np.interp(t_exp, t, v_y))),
               header="t v_y_experiment v_y")

# [1] A. ten Cate, C. H. Nieuwstad, J. J. Derksen and H. E. A. Van den Akker,
# "Particle imaging velocimetry experiments and lattice-Boltzmann simulations
# on a single sphere settling under gravity", Phys. Fluids 14 (2002) 4012-4025,
# doi:10.1063/1.1512918.
# [2] L. Barbeau, S. Golshan, J. Deng, S. Etienne, C. Beguin and B. Blais,
# "High-order moving immersed boundary and its application to a resolved
# CFD-DEM model", Comput. Fluids 268 (2024) 106094,
# doi:10.1016/j.compfluid.2023.106094.
