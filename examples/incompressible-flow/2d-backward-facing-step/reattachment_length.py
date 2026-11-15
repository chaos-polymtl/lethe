# SPDX-FileCopyrightText: Copyright (c) 2022-2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for 2D-backward-facing-step example
Computes the reattachment length (x_r) for Re = 100 and Re = 1000
for every mesh adaptation step from the velocity gradient on the bottom wall
"""

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import re
import os
import argparse
import sys

########################################
########################################

# EXAMPLE TO COMPUTE x_r (REATTACHMENT LENGTH) FOR Re = 100 and Re = 1000
# NOTE : THIS FILE MUST BE IN THE "2d-backward-facing-step" DIRECTORY TO
#        WORK PROPERLY. THE SIMULATIONS MUST OUTPUT THE VELOCITY GRADIENT
#        (set output velocity gradient = true in subsection post-processing)

# Parse arguments
parser = argparse.ArgumentParser(description='Arguments to compute the reattachment length')
parser.add_argument("-Re", "--Reynolds", type=int, help="Reynolds number (only 100 and 1000 supported)", required=True)
parser.add_argument("--validate", action="store_true", help="Launches the script in validation mode. This will log the content of the graph and prevent the display of figures", default=False)
args, leftovers = parser.parse_known_args()
Re = args.Reynolds

if (Re != 100 and Re != 1000):
    sys.exit("This Reynolds number is not supported. Only Re = 100 and Re = 1000.")

# Define folder path according to Re
folder = "./Reynolds" + str(Re)

# VARIABLES
L_in = 15            # Inlet length (position of the step)
L_out = 50           # Outlet length (length of the bottom wall downstream of the step)

# Reference reattachment values by Erturk (2008)
# (see attached benchmark-data.txt file for other Reynolds x_r values)
if (Re == 100):
    x_r_ref = 2.922
elif (Re == 1000):
    x_r_ref = 13.121

# The wall shear is sampled on a line slightly above the bottom wall (y = 0), so
# that every sampling point lies inside the domain
offset = 1e-8
n_samples = 100000

# DATA EXTRACTION
max_number = -1
pattern = re.compile(r'backward_facing_step_output\.(\d+)\.00000.vtu')

for filename in os.listdir(folder):
    match = pattern.match(filename)
    if match:
        number = int(match.group(1))
        if number > max_number:
            max_number = number

if max_number < 1:
    sys.exit("No output file found in " + folder + ".")

n_files_to_read = max_number
cells = []
x_r = []

for i in range(0, n_files_to_read):
    # Read all files except the 0 file
    file = (folder + '/backward_facing_step_output.'
                   + f'{i+1:05d}' + '.00000.vtu')
    data = pv.read(file)
    if "velocity_gradient" not in data.array_names:
        sys.exit("The velocity gradient is not in " + file + ". Set "
                 "'output velocity gradient = true' in the post-processing subsection.")

    # Wall shear du/dy along the bottom wall, downstream of the step. The
    # component 1 of the velocity gradient is du/dy (xy component).
    a = np.array([L_in + offset, offset, 0])
    b = np.array([L_in + L_out - offset, offset, 0])
    profile = data.sample_over_line(a, b, resolution=n_samples)
    valid = profile["vtkValidPointMask"].astype(bool)
    x = profile.points[valid, 0] - L_in
    dudy = profile["velocity_gradient"][valid, 1]

    # The wall shear is negative under the recirculation zone and the flow
    # reattaches where it becomes positive again. The reattachment point is the
    # end of the longest negative stretch (the primary recirculation zone), so
    # that a corner eddy at the foot of the step or a smaller secondary
    # recirculation zone further downstream is not mistaken for it.
    negative = (dudy < 0).astype(int)
    starts = np.where(np.diff(np.concatenate(([0], negative))) == 1)[0]
    ends = np.where(np.diff(np.concatenate((negative, [0]))) == -1)[0]
    reattaching = ends < len(dudy) - 1
    if not np.any(reattaching):
        sys.exit("No reattachment point found in " + file + ".")
    starts, ends = starts[reattaching], ends[reattaching]
    k = ends[np.argmax(x[ends] - x[starts])]
    x_reattachment = x[k] - dudy[k] * (x[k+1] - x[k]) / (dudy[k+1] - dudy[k])

    cells.append(data.n_cells)
    x_r.append(x_reattachment)
    print("Output file " + str(i+1) + " (" + str(data.n_cells) + " cells), reattachment point x_r: ", x_reattachment)

cells = np.array(cells)
x_r = np.array(x_r)

# Calculate relative error using last reattachment point and ref values
relative_error = np.abs(x_r_ref - x_r) / x_r_ref
print("Reference value Erturk (2008): ", x_r_ref)
print("Relative error abs(x_r_ref - x_r)/x_r_ref = ", relative_error[-1])

# Plot the reattachment length for every mesh adaptation step
fig0 = plt.figure()
plt.semilogx(cells, x_r, 'o-', label="Lethe")
plt.axhline(x_r_ref, color='k', ls='--', label="Erturk (2008)")
plt.xlabel("Number of cells")
plt.ylabel("$x_r$")
plt.title(f"Re = {Re}")
plt.legend()
plt.tight_layout()

if (not args.validate):
    plt.show()

else:
    # Store the plotted data so that they can be compared between validation runs
    np.savetxt(f"solution-reattachment-length-Re{Re}.dat",
               np.column_stack((cells, x_r, relative_error)),
               header="cells x_r relative_error")
    fig0.savefig(f"lethe-reattachment-length-Re{Re}.pdf")
