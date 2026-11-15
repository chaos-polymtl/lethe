# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Summary: Script to launch the Lethe simulations of the Reynolds number sweep of the
backward-facing step locally, one at a time. The cases must have been generated with
generate_cases_locally.py.
"""

import argparse
import os
import subprocess

PATH = os.path.dirname(os.path.abspath(__file__))

# User input
SWEEP_FOLDER = 'reynolds_sweep'
CASE_PREFIX = 'Re_'
PRM_FILE = 'backward-facing-step.prm'
LETHE_EXEC = 'lethe-fluid'

parser = argparse.ArgumentParser(description='Launch the cases of the Reynolds number sweep of the backward-facing step')
parser.add_argument("-np", "--nproc", type=int, help="Number of processes of every simulation", default=4)
args, leftovers = parser.parse_known_args()

sweep_path = os.path.join(PATH, SWEEP_FOLDER)
cases = sorted(case for case in os.listdir(sweep_path) if case.startswith(CASE_PREFIX))

for case in cases:
    case_path = os.path.join(sweep_path, case)
    print(f'Launching {case} on {args.nproc} processes')
    # The output of the simulation is written in the log file of the case
    with open(os.path.join(case_path, 'log'), 'w') as log:
        subprocess.run(['mpirun', '-np', str(args.nproc), LETHE_EXEC, PRM_FILE],
                       cwd=case_path, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
