# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the 2D sandpile formation
# example in Lethe. The simulations with the three rolling resistance models
# are launched at the same time, each of them on a third of the n_proc cores.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Rolling resistance models
models="epsd constant viscous"

# Store filenames of all plots in a variable (space-seperated)
plots=""
# Store filenames of all data files in a variable (space-seperated)
data=""
for model in $models; do
  plots="$plots sandpile-height-$model.pdf"
  data="$data solution-height-$model.dat solution-angle-$model.dat"
done

# Default path
output_root="./"
# Default number of cores, shared between the three simulations
n_proc=12

# Parse command-line arguments
parse_command_line "$@"

# Number of cores of each simulation
n_proc_model=$((n_proc / 3))
if [ $n_proc_model -lt 1 ]; then
  n_proc_model=1
fi

folder="$output_root/2d-sandpile-formation"

recreate_folder "$folder"

# Remove the outputs and checkpoints of a previous run, so that they are not
# reported again if this run fails. This must be done before the simulation,
# which writes some of these files.
rm -f $plots $data
for model in $models; do
  rm -f output_${model}/sandpile_${model}* sandpile_${model}.* sandpile_${model}_*
done

# Create the mesh
gmsh -2 sandpile.geo -o sandpile.msh > "$folder/log-mesh"

# Launch the three simulations at the same time. Only the total wall time is
# measured, since the comparison of the validation runs sums the wall times
# of all the logs.
run_models() {
  for model in $models; do
    mpirun -np $n_proc_model lethe-particles sandpile-$model.prm &> "$folder/log-$model" &
  done
  wait
}

{ time run_models ; } &> "$folder/log-time"

# Process the simulations
for model in $models; do
  python3 ./sandpile-postprocessing.py --folder ./ --rollingmethod $model --validate
done

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf
