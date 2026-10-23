# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the heated packed bed
# example in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="mean-temperatures.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-mean-temperatures.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/3d-heated-packed-bed"
action_load="mpirun -np $n_proc lethe-particles load-packed-bed.prm"
action_flip="mpirun -np $n_proc lethe-particles flip-packed-bed.prm"
action_heat="mpirun -np $n_proc lethe-particles heat-packed-bed.prm"

recreate_folder "$folder"

# Remove the outputs and checkpoints of a previous run, so that they are not
# reported again if this run fails. This must be done before the simulation,
# which writes some of these files.
rm -f $plots $data
# The output folder is removed as a whole, since the number of files it
# contains on many cores exceeds the maximal argument list length of rm.
rm -rf output dem_mp*

# Create the meshes of the solid surfaces
for surface in square-top square-rake square-side; do
  gmsh -2 $surface.geo -o $surface.msh >> "$folder/log-mesh"
done

{ time $action_load ; } &> "$folder/log-load"
{ time $action_flip ; } &> "$folder/log-flip"
{ time $action_heat ; } &> "$folder/log-heat"

# Process the simulation
python3 ./bed-postprocessing.py -f ./ --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf
