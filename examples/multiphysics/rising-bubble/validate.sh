# SPDX-FileCopyrightText: Copyright (c) 2020-2023 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the 2D lid-driven cavity benchmark
# in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="bubble-barycenter-case1.pdf bubble-velocity-case1.pdf bubble-contour-case1.pdf geo-mass-conservation-case1.pdf global-mass-conservation-case1.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-contour-case1-Geometric.dat solution-contour-case1-PDE-based.dat solution-contour-case1-Projection.dat solution-barycenter-case1-Geometric.dat solution-barycenter-case1-PDE-based.dat solution-barycenter-case1-Projection.dat solution-velocity-case1-Geometric.dat solution-velocity-case1-PDE-based.dat solution-velocity-case1-Projection.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/rising-bubble"
action_proj="mpirun -np $n_proc lethe-fluid rising-bubble-proj.prm" 
action_geo="mpirun -np $n_proc lethe-fluid rising-bubble-geo.prm" 
action_alge="mpirun -np $n_proc lethe-fluid rising-bubble-alge.prm" 
recreate_folder "$folder"

{ time $action_proj ; } &> "$folder/log-proj"
{ time $action_geo ; } &> "$folder/log-geo"
{ time $action_alge ; } &> "$folder/log-alge"


# Process the simulation
python3 ./rising-bubble.py -p rising-bubble-proj -g rising-bubble-geo -a rising-bubble-alge -c 1 --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf


