#!/bin/sh
# This is a script that is used by the continuous integration actios
# to ensure that the parameter files used in the examples are all
# valid parameter files. This file requires a valid installation
# of lethe and the presence of the lethe-parameter-check
# in the environment of the user.
# This script should run from the root of the lethe directory
#

# Force exit as soon as an error is encountered
set -e

# We check the syntax of all examples by family of solver. 

# Folders with the lethe-fluid syntax. This includes lethe-fluid, lethe-fluid-block, lethe-fluid-nitsche and lethe-fluid-sharp
lethe_fluid=("examples/incompressible-flow" "examples/multiphysics" "examples/sharp-immersed-boundary")

for folder in ${lethe_fluid[@]}; do
  for file_path in $(find "$folder" -type f -name "*.prm"); do
    # Because the sharp-edge solver actually parse the composite in the parameter file
    # we need to move to the folder where we are doing the parsing itself to prevent crashes
    echo "$file_path"
    current=$(pwd)
    DIR="$(dirname "${file_path}")" 
    FILE="$(basename "${file_path}")"
    cd "$DIR"
    lethe-parameter-check $FILE lethe-fluid
    cd "$current"
  done
done

# Folders with the lethe-particles syntax
lethe_particles=("examples/dem")

for folder in ${lethe_particles[@]}; do
  for file in $(find "$folder" -type f -name "*.prm"); do
    echo $file
    lethe-parameter-check $file lethe-particles
  done
done

# Folders with the lethe-fluid-particles syntax
lethe_fluid_particles=("examples/unresolved-cfd-dem")

for folder in ${lethe_fluid_particles[@]}; do
  for file in $(find "$folder" -type f -name "*.prm"); do
    echo $file
    # If generator is in the file name this is a lethe-particles file
    if [[ "$file" == *"generator"* ]];then
        lethe-parameter-check $file lethe-particles
    #else we assume it's a CFD-DEM file
    else
        lethe-parameter-check $file lethe-fluid-particles
    fi
  done
done

