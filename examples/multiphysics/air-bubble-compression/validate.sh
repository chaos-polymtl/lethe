# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the air bubble compression
# example in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="figure-air-bubble-compression-density.pdf figure-air-bubble-compression-pressure.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-density-pressure.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/air-bubble-compression"
action="mpirun -np $n_proc lethe-fluid air-bubble-compression.prm"

recreate_folder "$folder"

# Remove the outputs of a previous run, so that they are not reported again
# if this run fails. This must be done before the simulation, which writes
# some of these files.
rm -f $plots $data
rm -f output/air-bubble-compression*

{ time $action ; } &> "$folder/log"

# Process the simulation
python3 ./air-bubble-compression-postprocessing.py . air-bubble-compression.prm --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf
