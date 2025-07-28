# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Script to generate different cases (as folders) of the same single-particle-flow 
example with some parameters changed.

"""

######################################################################
# Import Libraries

import jinja2
import os
import numpy as np
import shutil
import argparse

######################################################################

PATH = os.getcwd()

parser = argparse.ArgumentParser(description='Arguments for generating the cases')
parser.add_argument("-n", "--number", type=int, help="Number of cases to generate.", required=True)
args, leftovers=parser.parse_known_args()

# Number of cases
number_of_cases=args.number

# User input
CASE_PREFIX = 'single_particle_flow_case_'
PARTICLE_PRM_FILE = 'initial-particle-template.tpl'
COUPLING_PRM_FILE = 'single-particle-flow-template.tpl'

# Generate parameter values
diameter_range     = (1e-6,1e-1)
viscosity_range    = (1e-9,1e-1)
# grad_div_range     = (,)
divisions_range = (1,10)
refinement_range   = (1,4)

diameter = np.linspace(*diameter_range, number_of_cases)
viscosity = np.random.uniform(*viscosity_range, size=number_of_cases)
l2_smoothing_length = 2*diameter
divisions = np.sort(np.random.randint(*divisions_range, size=number_of_cases))
refinement = np.sort(np.random.randint(*refinement_range, size=number_of_cases))
# particle location, with (0,0,0), the particle is usually in the corner of a cell
x,y,z = 0.0, 0.0, 0.0
grad_div = 0.1
velocity_order = 1
pressure_order = 1

# Create Jinja template
templateLoader = jinja2.FileSystemLoader(searchpath=PATH)
templateEnv = jinja2.Environment(loader=templateLoader)
particle_template = templateEnv.get_template(PARTICLE_PRM_FILE)
coupling_template = templateEnv.get_template(COUPLING_PRM_FILE)

# Generation of different cases
for i in range(number_of_cases):

		case_folder_name = f'{CASE_PREFIX}{i}'

		if os.path.exists(case_folder_name) and os.path.isdir(case_folder_name):
				shutil.rmtree(case_folder_name)
		
		# Insert the parameters in the prm template with Jinja2 and render it
		particle_parameters = particle_template.render(diameter = diameter[i], divisions = divisions[i], refinement = refinement[i], x=x, y=y, z=z, case_folder = case_folder_name)
		coupling_parameters = coupling_template.render(diameter = diameter[i], divisions = divisions[i], refinement = refinement[i], l2_smoothing_length = l2_smoothing_length[i], kinematic_viscosity = viscosity[i], grad_div = grad_div, velocity_order = velocity_order, pressure_order = pressure_order, case_folder = case_folder_name)

		# Create the folder of the case
		case_path = f'{PATH}/{case_folder_name}'
		os.mkdir(case_path)

		# Write a prm file from the template
		with open(f'{case_path}/initial-particle.prm', 'w') as f:
				f.write(particle_parameters)
		with open(f'{case_path}/single-particle-flow.prm', 'w') as f:
				f.write(coupling_parameters)