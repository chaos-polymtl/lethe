# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the capillary wave
# example in Lethe. Only the case with a time step of 15 times the capillary
# time-step constraint is launched. It is much faster than the 0.95 case and,
# since its stability relies on the GLS stabilization of the Navier-Stokes
# equations, it is more sensitive to the finite element formulation.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="capillary-wave-amplitude.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-amplitude.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/capillary-wave"
action="mpirun -np $n_proc lethe-fluid capillary-wave-TSM-15.prm"

recreate_folder "$folder"

# Remove the outputs of a previous run, so that they are not reported again
# if this run fails. This must be done before the simulation, which writes
# some of these files.
rm -f $plots $data
rm -f output-TSM-15/capillary-wave-TSM-15*

{ time $action ; } &> "$folder/log"

# Process the simulation
python3 ./capillary-wave-postprocess.py . capillary-wave-TSM-15.prm ./capillaryWaveData_rho_1_nu_5e-6_prosperetti.csv --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf
