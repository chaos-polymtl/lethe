#!/bin/bash
# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the contact between two
# sedimenting particles (drafting, kissing and tumbling, section 4.3 of
# Barbeau et al., Computers & Fluids 268 (2024) 106094) in Lethe.
# The case is simulated with Q1Q1 elements on a mesh of h = 1/128 cm.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Discretization that is simulated. q2q1-h128 is also supported, but it is
# about 40 times more expensive than q1q1-h128 and is not part of the
# validation for now.
discretization="q1q1-h128"

# Store filenames of all plots in a variable (space-separated)
plots="two-particles-velocity-$discretization.pdf two-particles-trajectories-$discretization.pdf"

# Store filenames of all data files in a variable (space-separated)
data="solution-particle-0-$discretization.dat solution-particle-1-$discretization.dat solution-metrics-$discretization.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=128

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/sedimentation-2-particles"
action="mpirun -np $n_proc lethe-fluid-sharp sedimentation-2-particles-$discretization.prm"

recreate_folder "$folder"

# Remove the outputs of a previous run, so that they are not reported again
# if this run fails. This must be done before the simulation, which writes
# some of these files.
rm -f $plots $data
rm -rf "out-$discretization"

{ time $action ; } &> "$folder/log-$discretization"

# Process the simulation
python3 ./post-process-sedimentation-2-particles.py -d "$discretization" --folder "out-$discretization" --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Add the metrics and the plots to the report
add_to_report "$folder" "solution-metrics-$discretization.dat" $plots
