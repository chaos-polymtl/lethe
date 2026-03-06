# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the Fichera oven example in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="fichera_oven_convergence.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-fichera-oven.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/fichera-oven"
action_fichera="mpirun -np $n_proc lethe-fluid fichera_oven.prm" 

recreate_folder "$folder"

{ time $action_fichera ; } &> "$folder/log-fichera"

# Process the simulation
python3 ./fichera_oven.py --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf


