"""
Summary: Script to launch all the Lethe simulations on a cluster, one job for every cases.
"""

import os

PATH = os.getcwd()

# File to look for in each case folder
SHELL_FILE = 'launch_lethe.sh'
case_count = 0

# Walk through the directory structure
for root, directories, files in os.walk(PATH):
    if SHELL_FILE in files:
        case_count+=1
        os.chdir(root)
        os.system(f'sbatch {SHELL_FILE}')
        os.chdir(PATH)

print(f'Number of cases launched = {case_count}')


