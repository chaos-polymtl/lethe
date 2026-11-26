# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the 3D ribbon mixer
# simulated with a single rotating frame (SRF) in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="ribbon-mixer-srf-power-number.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-power-number.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/3d-ribbon-mixer-srf"
action="mpirun -np $n_proc lethe-fluid ribbon-gls-Re1.prm"

recreate_folder "$folder"

# The simulation and its post-processing are run from the Re1 folder
cd Re1

# Remove the outputs of a previous run, so that they are not reported again
# if this run fails. This must be done before the simulation, which writes
# some of these files.
rm -f $plots $data
rm -f output/mixer* output/torque*

# Create the mesh
gmsh -3 ../diff-step-mesh.geo -o diff-step-mesh.msh > "$folder/log-mesh"

{ time $action ; } &> "$folder/log"

# Process the simulation
python3 ./postprocess_power_number.py --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

cd ..

# Add the plots to the report
add_to_report "$folder" $plots
