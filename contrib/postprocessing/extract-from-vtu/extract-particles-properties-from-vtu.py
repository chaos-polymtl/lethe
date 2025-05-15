#############################################################################
"""
Extraction particles properties and write them an particles.input file for
"file" insertion method.
"""
#############################################################################
'''Importing Libraries'''
import sys
import numpy as np
from lethe_pyvista_tools import *
#############################################################################

# Create the particle object
prm_file_name = "01-01-00.prm"
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

print("Job is finish")