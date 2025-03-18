# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Script to launch different cases (as folders) of the same problem 
with different polynomial degrees for the velocity and the pressure (mesh is adapted automatically).
"""
import os

PATH = os.getcwd()
PRM_FILE = 'mms_2d_steady.prm'
LETHE_EXEC = 'mpirun -np 6 lethe-fluid'

for root, directories, files in os.walk(PATH):
    if PRM_FILE in files and root != PATH:
                os.chdir(root)
                os.system(f'{LETHE_EXEC} {PRM_FILE}' + '> output.txt')