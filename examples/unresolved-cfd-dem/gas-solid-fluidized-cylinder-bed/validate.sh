# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the cylindrical gas-solid fluidized
# bed benchmark in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="pressure-drop-Re-mf-qcm-si.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-pressure-drop-Re.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/gas-solid-fluidized-bed-cylindrical"
action_packing="mpirun -np $n_proc lethe-particles packing-particles.prm" 
action_cfddem="mpirun -np $n_proc lethe-fluid-particles-matrix-free mf-fluidized-bed-modelA-semi-implicit.prm" 

recreate_folder "$folder"

{ time $action_packing ; } &> "$folder/log-packing"
{ time $action_cfddem ; } &> "$folder/log-cfddem"


# Process the simulation
python3 ./plot-pressure.py --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf


