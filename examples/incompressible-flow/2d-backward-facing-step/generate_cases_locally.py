# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Summary: Script to generate the cases (as folders) of the Reynolds number sweep of the
backward-facing step.

Every case is a folder reynolds_sweep/Re_XXXX that contains the parameter file
generated from the template backward-facing-step.tpl, in which the kinematic viscosity
is set to obtain the Reynolds number of the case (Re = 2/nu). By default, the Reynolds
numbers are those of the benchmark of Erturk (2008), given in benchmark-data.txt.
"""

import argparse
import os
import shutil

import jinja2

PATH = os.path.dirname(os.path.abspath(__file__))

# User input
SWEEP_FOLDER = 'reynolds_sweep'
CASE_PREFIX = 'Re_'
TEMPLATE_FILE = 'backward-facing-step.tpl'
PRM_FILE = 'backward-facing-step.prm'

parser = argparse.ArgumentParser(description='Generate the cases of the Reynolds number sweep of the backward-facing step')
parser.add_argument("-Re", "--reynolds", type=int, nargs='+', help="Reynolds numbers of the cases",
                    default=[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
args, leftovers = parser.parse_known_args()

# Create Jinja template
templateLoader = jinja2.FileSystemLoader(searchpath=PATH)
templateEnv = jinja2.Environment(loader=templateLoader, keep_trailing_newline=True)
template = templateEnv.get_template(TEMPLATE_FILE)

# Generation of the cases
for reynolds in args.reynolds:
    case_folder_name = f'{CASE_PREFIX}{reynolds:04d}'
    case_path = os.path.join(PATH, SWEEP_FOLDER, case_folder_name)

    if os.path.isdir(case_path):
        shutil.rmtree(case_path)
    os.makedirs(case_path)

    # The Reynolds number is based on twice the step height (Re = 2/nu)
    parameters = template.render(kinematic_viscosity=f'{2 / reynolds:.10g}', reynolds=reynolds)

    with open(os.path.join(case_path, PRM_FILE), 'w') as f:
        f.write(parameters)

    print(f'{SWEEP_FOLDER}/{case_folder_name} generated successfully!')
