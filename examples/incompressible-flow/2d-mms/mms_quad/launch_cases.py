# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Script to launch different cases (as folders) of the same problem 
with different polynomial degrees for the velocity and the pressure (mesh is adapted automatically).
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

for root, directories, files in os.walk(PATH):
    pattern = r'degu_(\d+)_degp_(\d+)'
    match = re.search(pattern, os.path.basename(root))  # Extract folder name details
    if match:
       os.chdir(root)
       velocity_approx_degree = int(match.group(1))
       pressure_approx_degree = int(match.group(2))
       # Print the message with the simulation parameters
       print(f'Running simulation with velocity approximation degree {velocity_approx_degree}, '
             f'pressure approximation degree {pressure_approx_degree},  at increasing refinement levels')         
       os.system(f'{LETHE_EXEC} {PRM_FILE}')