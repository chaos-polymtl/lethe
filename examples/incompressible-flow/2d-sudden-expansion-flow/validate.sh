# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the 2D sudden-expansion
# flow example in Lethe. Both the steady (Re = 70) and the transient
# (Re = 610) cases are launched.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="sudden-expansion-Re70.pdf sudden-expansion-Re610.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-Re70-profiles.dat solution-Re70-outlet.dat solution-Re610-profiles.dat solution-Re610-outlet.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/2d-sudden-expansion-flow"

recreate_folder "$folder"

# Remove the outputs of a previous run, so that they are not reported again
# if this run fails. This must be done before the simulation, which writes
# some of these files.
rm -f $plots $data
rm -f Reynolds70/sudden-expansion-flow-output* Reynolds610/sudden-expansion-flow-output*

# Create the mesh
gmsh -2 two-dimensional-jet.geo -o two-dimensional-jet.msh > "$folder/log-mesh"

# The parameter files must be launched from their folder, since they refer to
# the mesh as ../two-dimensional-jet.msh
for Re in 70 610; do
  action="mpirun -np $n_proc lethe-fluid Reynolds$Re-jet.prm"
  { time (cd Reynolds$Re && $action) ; } &> "$folder/log-Re$Re"
done

# Process the simulations
for Re in 70 610; do
  python3 ./velocity_profiles.py -Re $Re --validate
done

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf
