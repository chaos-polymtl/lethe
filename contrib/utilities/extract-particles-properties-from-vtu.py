#############################################################################
"""
Extraction particles properties and write them an particles.input file for
"from_file" insertion method.
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

# FEM force
fem_f_x = df["fem_force"][:, 0]
fem_f_y = df["fem_force"][:, 1]
fem_f_z = df["fem_force"][:, 2]

# FEM force
fem_t_x = df["fem_torque"][:, 0]
fem_t_y = df["fem_torque"][:, 1]
fem_t_z = df["fem_torque"][:, 2]

with open('example.txt', 'w') as file:
    # Write content to the file
    file.write('Hello, this is a text file.\n')
    file.write('This is another line of text.')

with open(output_file_name, 'w') as file:
    # Write content to the file
    file.write(
        "p_x; p_y; p_z; v_x  ; v_y; v_z; w_x; w_y; w_z; diameters; fem_force_x; fem_force_y; fem_force_z; "
        "fem_torque_x; fem_torque_y; fem_torque_z \n")

    for px, py, pz, vx, vy, vz, wx, wy, wz, d, ffx, ffy, ffz, ftx, fty, ftz in zip(p_x, p_y, p_z, v_x, v_y, v_z, w_x,
                                                                                   w_y, w_z, diameters, fem_f_x,
                                                                                   fem_f_y, fem_f_z, fem_t_x,
                                                                                   fem_t_y, fem_t_z):
        file.write(
            f"{px}; {py}; {pz}; {vx}; {vy}; {vz}; {wx}; {wy}; {wz}; {d}; {ffx}; {ffy}; {ffz}; {ftx}; {fty}; {ftz} \n")

print("Job is finish")