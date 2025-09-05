# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the cuboid sedimentation case

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="cuboid-sedimentation-velocity.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="cuboid-sedimentation-velocity.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/sedimentation-1-cuboid"
action_sharp="mpirun -np $n_proc lethe-fluid-sharp sedimentation-1-cuboid.prm"

recreate_folder "$folder"

{ time $action_sharp ; } &> "$folder/log-sharp"


# Process the simulation
python3 ./post-process-sedimentation-1-cuboid.py --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf
