# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the cylindrical packed bed
# example in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="pressure_drop_void_fraction.pdf velocity_void_fraction.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-profiles.dat solution-pressure-drop.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/cylindrical-packed-bed"
action_packing="mpirun -np $n_proc lethe-particles packing-particles.prm"
action_vans="mpirun -np $n_proc lethe-fluid-vans flow-in-bed.prm"

recreate_folder "$folder"

# Remove the outputs and checkpoints of a previous run, so that they are not
# reported again if this run fails. This must be done before the simulation,
# which writes some of these files.
rm -f $plots $data
rm -f output_dem/out* output/out* dem.* dem_*

{ time $action_packing ; } &> "$folder/log-packing"
{ time $action_vans ; } &> "$folder/log-vans"

# Process the simulation
python3 ./post_process_packed_bed.py --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf
