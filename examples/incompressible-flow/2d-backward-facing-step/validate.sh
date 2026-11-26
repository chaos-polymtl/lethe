# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the 2D backward-facing step
# in Lethe. The reattachment length is computed for ten Reynolds numbers
# (Re = 100, 200, ..., 1000) and compared with the benchmark of Erturk (2008).

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="lethe-reattachment-length-sweep.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-reattachment-length-sweep.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

# Reynolds numbers of the sweep
reynolds_numbers="100 200 300 400 500 600 700 800 900 1000"

folder="$output_root/2d-backward-facing-step"

recreate_folder "$folder"

# Remove the outputs of a previous run, so that they are not reported again
# if this run fails. This must be done before the simulations, which write
# some of these files.
rm -f $plots $data
rm -rf reynolds_sweep

# Generate the cases of the sweep
python3 ./generate_cases_locally.py -Re $reynolds_numbers > "$folder/log-generate"

# Launch the simulations one at a time
for Re in $reynolds_numbers; do
  case_name=$(printf "Re_%04d" $Re)
  (cd reynolds_sweep/$case_name && { time mpirun -np $n_proc lethe-fluid backward-facing-step.prm ; } &> "$folder/log-$case_name")
done

# Process the simulations
python3 ./postprocess_reynolds_sweep.py --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Add the plots to the report
add_to_report "$folder" $plots
