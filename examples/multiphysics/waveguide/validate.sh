# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the waveguide example in
# Lethe. The errors with respect to the analytical solution are read from the
# log of the simulation.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="waveguide-convergence.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-errors.dat solution-convergence-rates.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/waveguide"
action="mpirun -np $n_proc lethe-fluid waveguide.prm"

recreate_folder "$folder"

# Remove the outputs of a previous run, so that they are not reported again
# if this run fails. This must be done before the simulation, which writes
# some of these files.
rm -f $plots $data
rm -f degree_1/*

{ time $action ; } &> "$folder/log"

# Process the simulation
python3 ./waveguide-postprocessing.py --log "$folder/log" --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf
