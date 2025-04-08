# SPDX-FileCopyrightText: Copyright (c) 2020-2023 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the 2D lid-driven-cavity benchmark
# in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="lethe-ghia-re-400-comparison.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution.dat"

# Default path
output_root="./"

# Default number of cores
n_proc=1

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/2d-lid-driven-cavity"
action="mpirun -np $n_proc lethe-fluid cavity.prm" 
recreate_folder "$folder"

{ time $action ; } &> "$folder/log"

python3 post_process_Reynolds_400.py --validate

        # Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf



