# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the pseudo-2D gas-solid fluidized
# bed benchmark in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="pressure-fluctuations.pdf void-fraction-fluctuations.pdf bed-height.pdf pressure-psd.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-relative-pressure.dat solution-void-fraction.dat solution-bed-height.dat solution-pressure-power-spectral-density.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/pseudo-2d-gas-solid-fluidized-bed"
action_packing="mpirun -np $n_proc lethe-particles packing-particles.prm" 
action_cfddem="mpirun -np $n_proc lethe-fluid-particles gas-solid-fluidized-bed.prm" 

recreate_folder "$folder"

{ time $action_packing ; } &> "$folder/log-packing"
{ time $action_cfddem ; } &> "$folder/log-cfddem"


# Process the simulation
python3 ./fluidized-bed-postprocessing.py -f ./ --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf


