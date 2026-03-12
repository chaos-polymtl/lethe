#############################################################################
"""
SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
This script extracts particles properties and writes them to a particles.input
file for "file" insertion method.
"""
#############################################################################
'''Importing Libraries'''
import argparse
import sys
import numpy as np
from lethe_pyvista_tools import *
#############################################################################

# Create the particle object
# Argument parser
parser = argparse.ArgumentParser(description="Extract particle data from Lethe output")
parser.add_argument("prm_file_name", help="Name of the .prm file")
args = parser.parse_args()

prm_file_name = args.prm_file_name

pvd_name = 'out.pvd'
output_file_name = "particles.input"
particle = lethe_pyvista_tools("./", prm_file_name, pvd_name)
#############################################################################

df = particle.get_df(-1)

# Positions
p_x = df.points[:, 0]
p_y = df.points[:, 1]
p_z = df.points[:, 2]

# Velocity
v_x = df["velocity"][:, 0]
v_y = df["velocity"][:, 1]
v_z = df["velocity"][:, 2]

# Omega
w_x = df["omega"][:, 0]
w_y = df["omega"][:, 1]
w_z = df["omega"][:, 2]

# Diameter
diameters = df["diameter"][:]

with open(output_file_name, 'w') as file:
    # Write content to the file
    file.write(
        "p_x; p_y; p_z; v_x  ; v_y; v_z; w_x; w_y; w_z; diameters; \n")

    for px, py, pz, vx, vy, vz, wx, wy, wz, d in zip(p_x, p_y, p_z, v_x, v_y, v_z, w_x,
                                                                                   w_y, w_z, diameters):
        file.write(
            f"{px}; {py}; {pz}; {vx}; {vy}; {vz}; {wx}; {wy}; {wz}; {d}; \n")

print("Job is finished, the particle properties are saved in the file: " + output_file_name)