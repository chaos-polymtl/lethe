# SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Summary: Script to generate different cases (as folders) of the same problem 
with one parameters being changed.

This script generates multiples folders of the example of the flow around a cylinder.
In each folder, we change the input velocity (u) in the parameter file in order to 
have a different case of the same problem.
"""

import jinja2
import os
import numpy as np
import shutil

PATH = os.getcwd()

# User input
CASE_PREFIX = 'wetsand_'
PRM_FILE = 'granuheap_multicase.prm'
MESH_FILE_1 = 'support.msh'
MESH_FILE_2 = 'cylinder.msh'
number_of_cases = 4

# Generation of data points
energy_first = 0.0010
energy_last = 0.0100
energy = np.linspace(energy_first, energy_last, number_of_cases)

rolling_friction_first = 0.3
rolling_friction_last = 0.7
rolling_friction = np.linspace(rolling_friction_first, rolling_friction_last, number_of_cases-1)

# Create Jinja template
templateLoader = jinja2.FileSystemLoader(searchpath=PATH)
templateEnv = jinja2.Environment(loader=templateLoader)
template = templateEnv.get_template(PRM_FILE)

# Generation of different cases
for v in rolling_friction:
    for u in energy:

        case_folder_name = f'{CASE_PREFIX}{u:.4f}_{v:.2f}'

        if os.path.exists(case_folder_name) and os.path.isdir(case_folder_name):
            shutil.rmtree(case_folder_name)
        
        # Insert the velocity in the prm template with Jinja2 and render it
        parameters = template.render(energy = u, rolling_friction = v)

        # Create the folder of the case and put the prm template in it
        case_path = f'{PATH}/{case_folder_name}'
        os.mkdir(case_path)
        shutil.copy(f'{PATH}/{PRM_FILE}', f'{case_path}/{PRM_FILE}')
        
        # Copy the mesh file (in order to launch Lethe in seperate folders)
        shutil.copy(f'{PATH}/{MESH_FILE_1}', f'{case_path}/{MESH_FILE_1}')
        shutil.copy(f'{PATH}/{MESH_FILE_2}', f'{case_path}/{MESH_FILE_2}')

        # Write a unique prm file with the prm template being updated
        with open(f'{case_path}/{PRM_FILE}', 'w') as f:
            f.write(parameters)

        print(f'{case_folder_name} generated successfully!')
