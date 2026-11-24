#!/bin/bash
# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the contact between two
# sedimenting particles (drafting, kissing and tumbling, section 4.3 of
# Barbeau et al., Computers & Fluids 268 (2024) 106094) in Lethe.
# The case is simulated with Q1Q1 elements on a mesh of h = 1/128 cm.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Discretizations that are simulated. q2q1-h128 is also supported, but it is
# about ten times more expensive than q1q1-h128 and is not part of the
# validation for now.
discretizations="q1q1-h128"

# Store filenames of all plots and data files in variables (space-separated)
plots=""
data=""
for d in $discretizations; do
  plots="$plots two-particles-velocity-$d.pdf two-particles-trajectories-$d.pdf"
  data="$data solution-particle-0-$d.dat solution-particle-1-$d.dat solution-metrics-$d.dat"
done

# Default path
output_root="."
# Default number of cores
n_proc=128

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/sedimentation-2-particles"

recreate_folder "$folder"

# Remove the outputs of a previous run, so that they are not reported again
# if this run fails. This must be done before the simulations, which write
# some of these files.
rm -f $plots $data
for d in $discretizations; do
  rm -rf "out-$d"
done

status=0
for d in $discretizations; do
  { time mpirun -np $n_proc lethe-fluid-sharp "sedimentation-2-particles-$d.prm" ; } &> "$folder/log-$d"
  simulation_status=$?
  if [ $simulation_status -ne 0 ]; then
    echo "Error: the $d simulation failed, see $folder/log-$d." >&2
    status=$simulation_status
  fi

  # Process the simulation, even if it stopped early, to show how far it went
  if ! python3 ./post-process-sedimentation-2-particles.py -d "$d" --folder "out-$d" --validate; then
    echo "Error: the post-processing of the $d simulation failed, see $folder/log-$d." >&2
    status=1
  fi
done

# Copy the information to the log folder. A file is missing if its
# post-processing failed.
for file in $plots $data; do
  if [ -f "$file" ]; then
    cp "$file" "$folder"
  fi
done

# Append the plots to the report, which only exists when this script is
# launched by validate_lethe.sh
existing_plots=""
for plot in $plots; do
  if [ -f "$plot" ]; then
    existing_plots="$existing_plots $plot"
  fi
done
report="$output_root/report.pdf"
previous_report=""
if [ -f "$report" ]; then
  previous_report="$report"
fi
if [ -n "$existing_plots" ]; then
  magick -density 300 $previous_report $existing_plots -quality 100 "$output_root/temporary.pdf"
  mv "$output_root/temporary.pdf" "$report"
fi

exit $status
