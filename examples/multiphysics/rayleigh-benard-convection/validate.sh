# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the Rayleigh-Bénard
# convection example at Ra=10^4 in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="solution-rayleigh-uy.pdf solution-rayleigh-xv.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-rayleigh_10k-uy.dat solution-rayleigh_10k-xv.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/rayleigh-benard-convection"
action="mpirun -np $n_proc lethe-fluid rayleigh-benard-convection-Ra10k.prm"

recreate_folder "$folder"

# Remove the outputs of a previous run, so that they are not reported again
# if this run fails. This must be done before the simulation, which writes
# some of these files.
rm -f $plots $data
rm -f output_10k/rayleigh-benard_convection*

{ time $action ; } &> "$folder/log"

# Process the simulation
python3 ./rayleigh-benard-convection.py -f ./output_10k -Ra 10k --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf
