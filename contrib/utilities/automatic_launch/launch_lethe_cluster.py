"""
Summary: Script to launch all the Lethe simulations on a cluster, one job for every cases.
"""

import os

PATH = os.getcwd()

# User input
CASE_PREFIX = 'cylinder_u_'
PRM_FILE = 'cylinder.prm'
SHELL_FILE = 'launch_lethe.sh'

for root, directories, files in os.walk(PATH):

    if CASE_PREFIX in root:

        os.chdir(root)

        case_name = root.split('/')[-1]
        os.system(f'sbatch -J {case_name} {SHELL_FILE}')