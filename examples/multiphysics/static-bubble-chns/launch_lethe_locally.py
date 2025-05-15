# SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Summary: Script to launch all the Lethe simulations, one at a time locally.
"""

import os

PATH = os.getcwd()

# User input
CASE_PREFIX = 'R_'
PRM_FILE = 'static-bubble-chns.prm'
LETHE_EXEC = 'lethe-fluid'

for root, directories, files in os.walk(PATH):

    if CASE_PREFIX in root:

        os.chdir(root)

        os.system(f'mpirun -np 8 {LETHE_EXEC} {PRM_FILE}')
