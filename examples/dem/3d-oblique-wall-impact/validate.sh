# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the oblique wall impact
# example in Lethe. The 34 serial simulations are launched n_proc at a time.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="rebound.pdf omega.pdf coeff_restitution.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-rebound.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=34

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/3d-oblique-wall-impact"

recreate_folder "$folder"

# Remove the outputs and the parameter files of a previous run, so that they
# are not reported again if this run fails. This must be done before the
# simulation, which writes some of these files. Every simulation writes its
# results in a folder named after its angle.
rm -f $plots $data
for prm in run_oblique_impact_*.prm; do
  [ -f "$prm" ] || continue
  angle=${prm#run_oblique_impact_}
  rm -rf "./${angle%.prm}"
  rm -f "$prm"
done

# Generate the parameter files
python3 ./oblique_wall_impact_case_generator.py

# Launch the simulations. The output of all simulations is written in the same
# log, so that the wall time is that of the whole set of simulations.
run_cases() {
  ls run_oblique_impact_*.prm | xargs -P $n_proc -n 1 lethe-particles
}

{ time run_cases ; } &> "$folder/log"

# Process the simulations
python3 ./oblique_wall_impact_postprocessing.py --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf
