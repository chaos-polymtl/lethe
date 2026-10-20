# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Script to generate the cases (as folders) of the parametric study of the single
particle sedimentation example. The mesh size and the L2 smoothing length of the
void fraction are varied. Each folder contains the parameter file of the particle
insertion (lethe-particles) and of the CFD-DEM simulation (lethe-fluid-particles).
"""

import jinja2
import os
import shutil

PATH = os.getcwd()

# User input
CASE_PREFIX = 'case_'
TEMPLATE_FILES = ['initial-particles.tpl', 'single-particle-sedimentation.tpl']
PRM_FILES = ['initial-particles.prm', 'single-particle-sedimentation.prm']

# Particle diameter
dp = 0.002663

# Number of cells in the x and z directions. The domain is 20 particle
# diameters wide, hence the cell size over the particle diameter is 20/n. An odd
# number of cells is used for all meshes so that the particle always falls along
# a line of cell centers.
n_cells = [15, 11, 9, 7]

# L2 smoothing length of the void fraction over the particle diameter
l2_over_dp = [2, 5, 10]

# Create Jinja templates
templateLoader = jinja2.FileSystemLoader(searchpath=PATH)
templateEnv = jinja2.Environment(loader=templateLoader, keep_trailing_newline=True)
templates = [templateEnv.get_template(f) for f in TEMPLATE_FILES]

# Generation of the different cases
for n in n_cells:
    for l in l2_over_dp:
        case_folder_name = f'{CASE_PREFIX}nx_{n}_Ldp_{l}'

        if os.path.exists(case_folder_name) and os.path.isdir(case_folder_name):
            shutil.rmtree(case_folder_name)

        # Create the folder of the case
        case_path = f'{PATH}/{case_folder_name}'
        os.mkdir(case_path)

        # The domain is twice as high as it is wide
        for template, prm_file in zip(templates, PRM_FILES):
            parameters = template.render(nx=n, ny=2 * n, nz=n,
                                         l2_smoothing_length=f'{l * dp:g}')

            # Write a unique prm file with the prm template being updated
            with open(f'{case_path}/{prm_file}', 'w') as f:
                f.write(parameters)

        print(f'{case_folder_name} generated successfully!')
