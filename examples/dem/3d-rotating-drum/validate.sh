# SPDX-FileCopyrightText: Copyright (c) 2020-2023 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the 2D Taylor-Couette benchmark
# in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="lethe-rotating-drum-comparison-depth.pdf lethe-rotating-drum-comparison-free-surface.pdf"

# Default path
output_root="./"

# Default number of cores
n_proc=1

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/3d-rotating-drum"
action="mpirun -np $n_proc lethe-particles load-rotating-drum.prm"
recreate_folder "$folder"

{ time $action ; } &> "$folder/log_loading" & pid=$!

show_animation $pid "---> Loading rotating drum"
wait $pid
# Loading done

action="mpirun -np $n_proc lethe-particles rotating-drum.prm"

{ time $action ; } &> "$folder/log_rotating" & pid=$!
show_animation $pid "---> Running simulation"
wait $pid

python3 post_processing_rotating_drum.py  -f ./ --validate

# Copy the information to the log folder
cp $plots $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf



