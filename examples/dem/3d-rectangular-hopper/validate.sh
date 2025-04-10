# SPDX-FileCopyrightText: Copyright (c) 2020-2023 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the rectangular hopper
# in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="hopper-flow-rate.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="mass-and-discharge-rate.txt solution.dat"

# Default path
output_root="./"

# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/3d-rectangular-hopper"
action="mpirun -np $n_proc lethe-particles hopper.prm" 
recreate_folder "$folder"

# Recreate the mesh
gmsh -3 hopper_structured.geo > "$folder/log-mesh"

{ time $action ; } &> "$folder/log" & pid=$!
show_animation $pid "---> Running regular simulation"
wait $pid

python3 hopper_post_processing.py -f . --prm hopper.prm --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300 -pointsize 12 text:mass-and-discharge-rate.txt mass-and-discharge-rate.pdf

magick -density 300  $output_root/report.pdf mass-and-discharge-rate.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf


folder="$output_root/3d-rectangular-hopper-periodic"
action="mpirun -np $n_proc lethe-particles hopper_periodic.prm" 
recreate_folder "$folder"

# Recreate the mesh
gmsh -3 hopper_structured_periodic.geo > "$folder/log-mesh"

{ time $action ; } &> "$folder/log" & pid=$!
show_animation $pid "---> Running periodic BC simulation" 
wait $pid

python3 hopper_post_processing.py -f . --prm hopper_periodic.prm --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300 -pointsize 12 text:mass-and-discharge-rate.txt mass-and-discharge-rate.pdf

magick -density 300  $output_root/report.pdf mass-and-discharge-rate.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf


