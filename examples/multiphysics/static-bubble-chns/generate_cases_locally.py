# SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Summary: Script to generate different cases (as folders) of the same problem 
with one parameters being changed.

This script generates multiples folders of the example of the static bubble using the CHNS model.
In each folder, we change the input radius (r) and adapt the mobility in the parameter file in order to have a different case of the same problem.
"""

import jinja2
import os
import numpy as np
import shutil

PATH = os.getcwd()

# User input
CASE_PREFIX = 'R_'
PRM_FILE = 'static-bubble-chns.prm'

# Generation of data points

radii = [0.15,0.20,0.25,0.30,0.35,0.40]

sigma = 1
t_sim = 2
epsilon = 0.01381

# Create Jinja template
templateLoader = jinja2.FileSystemLoader(searchpath=PATH)
templateEnv = jinja2.Environment(loader=templateLoader)
template = templateEnv.get_template(PRM_FILE)

# Generation of different cases
for r in radii:
    mobility = r*r*epsilon/(10*sigma*t_sim)
    case_folder_name = f'{CASE_PREFIX}{r:.2f}'

    if os.path.exists(case_folder_name) and os.path.isdir(case_folder_name):
        shutil.rmtree(case_folder_name)
    
    # Insert the radii in the prm template with Jinja2 and render it
    parameters = template.render(R = r, mobility=mobility)

    # Create the folder of the case and put the prm template in it
    case_path = f'{PATH}/{case_folder_name}'
    os.mkdir(case_path)
    shutil.copy(f'{PATH}/{PRM_FILE}', f'{case_path}/{PRM_FILE}')
    
    # Write a unique prm file with the prm template being updated
    with open(f'{case_path}/{PRM_FILE}', 'w') as f:
        f.write(parameters)

    print(f'{case_folder_name} generated successfully!')
