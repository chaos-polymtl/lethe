# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import argparse
import os
import sys
sys.path.append("$LETHE_PATH/contrib/postprocessing/")
from lethe_pyvista_tools import *

######################################################################
# Define function to write flow disturbance data on .dat file

def write_flow_disturbance_data(fluid,particle,folder):

    df_fluid = fluid.get_df(1)

    # The particle may be moved between 0 and h for x,y,z to assess how its position in the cell affects the results
    x, y, z = particle.prm_dict['list x'], particle.prm_dict['list y'], particle.prm_dict['list z']
    dp = particle.prm_dict['diameter']

    # Get grid cell length from grid arguments (this works for a subdivided_hyper_cube)
    initial_refinement = fluid.prm_dict['initial refinement']
    grid_args = fluid.prm_dict['grid arguments'].split(': ')
    subdivisions =  int(grid_args[0])
    corner_1 = float(grid_args[1])
    corner_2 = float(grid_args[2])
    h = (corner_2 - corner_1) / subdivisions / 2**(initial_refinement)

    # Get fluid velocity from closest point to the particle
    point_id = df_fluid.find_closest_point([x,y,z])
    u_h = df_fluid.point_data["velocity"][point_id][0]

    # Write dp/h, Rep, |u_h - u_inf | / u_inf and particle position on file
    nu_f = fluid.prm_dict['kinematic viscosity']
    u_inf = 1 # should not change

    file_name = folder + 'velocity_disturbance.dat'
    header = 'dp/h,Rep,|u_h-u_inf|/u_inf\n'
    line = f'{dp/h:.6f},{u_inf*dp/nu_f:.4f},{abs(u_h-u_inf)/u_inf:.10f}\n'

    file_exists = os.path.exists(file_name)
    with open(file_name, 'a') as file:
        if not file_exists:
            file.write(header)
        file.write(line)


######################################################################


'''
This post-processing code can be used to get data from multiple single-particle-flow cases of dp/h, Rep = (u_inf*dp)/nu_f, |u_h-u_inf|/u_inf when varying the cell size, the particle diameter, the fluid kinematic viscosity and the position of the particle (preferably within the same cell, where (0,0,0) is).
'''

parser = argparse.ArgumentParser(description='Arguments for the post-processing of the simulation')
parser.add_argument("-f", "--folder", type=str, help="Folder path. This folder is the folder where all the cases are generated.", required=True)
args, leftovers=parser.parse_known_args()

# Example folder
ex_folder=args.folder

# File names
particle_prm_file = 'initial-particle.prm'
pvd_file     = 'out.pvd'
coupling_prm_file = 'single-particle-flow.prm'

# Loop through all entries in the example folder
for entry in os.listdir(ex_folder):

    # Check if the entry matches the prefix
    if entry.startswith('single_particle_flow_case_'):

        case_folder = os.path.join(ex_folder, entry)
        print(f"Processing folder {case_folder}.")

        if (not os.path.exists(f"{case_folder}/output/{pvd_file}")) or (not os.path.exists(f"{case_folder}/output_dem/{pvd_file}")):
            print(f"Skipping {case_folder} — {pvd_file} not found for one of the simulations.\n")
            continue

        fluid = lethe_pyvista_tools(ex_folder, f"{case_folder}/{coupling_prm_file}", pvd_file)
        if len(fluid.time_list) < 2 :
            print(f"Skipping this folder. Not enough iterations.\n")
            continue
        particle = lethe_pyvista_tools(ex_folder, f"{case_folder}/{particle_prm_file}", pvd_file)
        write_flow_disturbance_data(fluid,particle,ex_folder)

        print(f"Finished processing folder {case_folder}.\n")