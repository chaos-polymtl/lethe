# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Script to launch different cases (as folders) of the same problem 
with different levels of refinement and polynomial degrees for the velocity and the pressure.
"""
import os
import argparse
import re

PATH = os.getcwd()
PRM_FILE = 'mms_2d_steady.prm'

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--num_cores", type=int, default=1, help="Number of cores to use (default: 1)")
args = parser.parse_args()
print(f"Using {args.num_cores} cores.")

LETHE_EXEC = f'mpirun -np {args.num_cores} lethe-fluid'

print(LETHE_EXEC)

for root, directories, files in os.walk(PATH):
    pattern = r'ref_(\d+)_degu_(\d+)_degp_(\d+)'
    match = re.search(pattern, os.path.basename(root))  # Extract folder name details
    if match:
       os.chdir(root)
       refinement_level = int(match.group(1))
       velocity_approx_degree = int(match.group(2))
       pressure_approx_degree = int(match.group(3))
       # Calculate the number of cells
       num_cells = 2**refinement_level * 2**refinement_level * 8         
       # Print the message with the simulation parameters
       print(f'Running simulation with velocity approximation degree {velocity_approx_degree}, '
             f'pressure approximation degree {pressure_approx_degree}, '
             f'and refinement level {refinement_level}, corresponding to {num_cells} cells.')         
       os.system(f'{LETHE_EXEC} {PRM_FILE}')