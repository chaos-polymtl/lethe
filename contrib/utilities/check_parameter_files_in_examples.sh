#!/usr/bin/env bash
# This is a script that is used by the continuous integration actions
# to ensure that the parameter files used in the examples are all
# valid parameter files which do not use any deprecated parameter
# or parameter value. This file requires a valid installation
# of lethe and the presence of the lethe-parameter-check
# in the environment of the user.
# This script should run from the root of the lethe directory
#

# Force exit as soon as an error is encountered
set -e

# Check a parameter file with lethe-parameter-check. The check fails if the
# file is not valid for the solver family, or if its parsing reports the use of
# a deprecated parameter or of a deprecated parameter value.
#   $1: path of the parameter file
#   $2: solver family (lethe-fluid, lethe-particles or lethe-fluid-particles)
check_parameter_file() {
  local file_path="$1"
  local solver_family="$2"
  local output
  local status=0

  echo "$file_path"
  # The files read while parsing (e.g. the composite shapes of the sharp-edge
  # solver) are given relative to the parameter file, so the check is run from
  # its folder to prevent crashes
  output=$(cd "$(dirname "$file_path")" &&
    lethe-parameter-check "$(basename "$file_path")" "$solver_family" 2>&1) || status=$?
  if [ -n "$output" ]; then
    echo "$output"
  fi

  if [ "$status" -ne 0 ]; then
    echo "Error: $file_path is not a valid $solver_family parameter file"
    exit 1
  fi
  if echo "$output" | grep -qi "deprecated"; then
    echo "Error: $file_path uses a deprecated parameter or parameter value"
    exit 1
  fi
}

# We check the syntax of all examples by family of solver.

# Folders with the lethe-fluid syntax. This includes lethe-fluid, lethe-fluid-block, lethe-fluid-nitsche and lethe-fluid-sharp
lethe_fluid=("examples/incompressible-flow" "examples/multiphysics" "examples/sharp-immersed-boundary")

for folder in ${lethe_fluid[@]}; do
  for file in $(find "$folder" -type f -name "*.prm"); do
    check_parameter_file "$file" lethe-fluid
  done
done

# Folders with the lethe-particles syntax
lethe_particles=("examples/dem")

for folder in ${lethe_particles[@]}; do
  for file in $(find "$folder" -type f -name "*.prm"); do
    check_parameter_file "$file" lethe-particles
  done
done

# Folders with the lethe-fluid-particles syntax
lethe_fluid_particles=("examples/unresolved-cfd-dem")

for folder in ${lethe_fluid_particles[@]}; do
  for file in $(find "$folder" -type f -name "*.prm"); do
    # If "particles" is in the filename this is a
    # lethe-particles file
    if [[ "$file" == *"particles"* ]];then
        check_parameter_file "$file" lethe-particles
    #else we assume it's a CFD-DEM file
    else
        check_parameter_file "$file" lethe-fluid-particles
    fi
  done
done

