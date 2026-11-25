# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the single particle
# sedimentation example in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="particle-position.pdf particle-velocity.pdf center-void-fraction.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-particle.dat solution-rmse.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/single-particle-sedimentation"
action_dem="mpirun -np $n_proc lethe-particles initial-particles.prm"
action_cfddem="mpirun -np $n_proc lethe-fluid-particles single-particle-sedimentation.prm"

recreate_folder "$folder"

# Remove the outputs and checkpoints of a previous run, so that they are not
# reported again if this run fails. This must be done before the simulation,
# which writes some of these files.
rm -f $plots $data
rm -f output_dem/out* output/result_* dem.* dem_*

{ time $action_dem ; } &> "$folder/log-dem"
{ time $action_cfddem ; } &> "$folder/log-cfddem"

# Process the simulation
python3 ./single_particle_sedimentation.py --particle_pvd output/result__particles.pvd --fluid_pvd output/result_.pvd --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf
