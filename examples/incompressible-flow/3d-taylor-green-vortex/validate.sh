# SPDX-FileCopyrightText: Copyright (c) 2020-2023 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the 2D lid-driven cavity benchmark
# in Lethe.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh


# Store filenames of all plots in a variable (space-seperated)
plots=" dissipation_comparison.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution_enstrophy.dat solution_kinetic_energy.dat"

# Default path
output_root="./"

# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/3d-taylor-green-vortex"
action="mpirun -np $n_proc lethe-fluid-matrix-free tgv-matrix-free-validation.prm"
recreate_folder "$folder"

{ time $action ; } &> "$folder/log"
python3 calculate_dissipation_rate_constant_cfl.py -f output/  -i kinetic_energy.dat
python3 plot_dissipation_rate.py -ke output/ke_rate.dat -ens output/enstrophy.dat -v 0.000625 --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf


