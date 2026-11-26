# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the 3D flow around a sphere
# at Re = 0.1 and Re = 150, with the matrix-based and the matrix-free solvers.
# At Re = 0.1, the matrix-free simulation is also carried out on domains
# enlarged twice and four times to extrapolate the drag coefficient to an
# infinite domain.

# Load the validation functions
. ../../../contrib/validation/validation_functions.sh

# Store filenames of all plots in a variable (space-seperated)
plots="lethe-sphere-drag-domain-size.pdf lethe-sphere-drag-coefficient.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="solution-drag-coefficient.dat solution-drag-coefficient-extrapolated.dat"

# Default path
output_root="./"
# Default number of cores
n_proc=16

# Parse command-line arguments
parse_command_line "$@"

folder="$output_root/3d-flow-around-sphere"

recreate_folder "$folder"

# Remove the outputs of a previous run, so that they are not reported again
# if this run fails. This must be done before the simulations, which write
# some of these files.
rm -f $plots $data
rm -rf sphere-0.1 sphere-150 sphere-0.1-matrix-free sphere-0.1-matrix-free-x2 \
       sphere-0.1-matrix-free-x4 sphere-150-matrix-free
rm -f sphere.msh sphere-x2.msh sphere-x4.msh

# Create the mesh of the original domain and the meshes of the enlarged domains
gmsh -3 sphere.geo -o sphere.msh > "$folder/log-mesh"
for scale in 2 4; do
  gmsh -setnumber domain_scale $scale -3 sphere.geo -o sphere-x$scale.msh > "$folder/log-mesh-x$scale"
done

# Launch the matrix-based and the matrix-free simulations at both Reynolds numbers
for case in sphere-0.1 sphere-150; do
  { time mpirun -np $n_proc lethe-fluid $case.prm ; } &> "$folder/log-$case"
  { time mpirun -np $n_proc lethe-fluid-matrix-free $case-matrix-free.prm ; } &> "$folder/log-$case-matrix-free"
done

# Launch the matrix-free simulations on the enlarged domains at Re = 0.1
for scale in 2 4; do
  { time mpirun -np $n_proc lethe-fluid-matrix-free sphere-0.1-matrix-free-x$scale.prm ; } &> "$folder/log-sphere-0.1-matrix-free-x$scale"
done

# Process the simulations
python3 ./postprocess_drag.py --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Add the plots to the report
add_to_report "$folder" $plots
