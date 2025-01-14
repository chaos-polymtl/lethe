# SPDX-FileCopyrightText: Copyright (c) 2020-2023 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the rectangular hopper
# in Lethe.

# Function to recreate the folder
recreate_folder() {
    local folder_path="$1"

    if [ -d "$folder_path" ]; then
        rm -rf "$folder_path"  # Deletes the folder and its contents
    fi
    mkdir -p "$folder_path"  # Creates the folder
}

# Store filenames of all plots in a variable (space-seperated)
plots="hopper-flow-rate.pdf"

# Store filenames of all data files in a variable (space-seperated)
data="mass-and-discharge-rate.txt solution.dat"

# Default path
default_value="./"

# Default number of cores
n_proc=1

# Parse command-line arguments
output_root=$default_value
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--path)
            if [[ -n $2 && $2 != -* ]]; then
                output_root="$2"
                shift 2
            else
                echo "Warning: Missing value for $1. Using default value: $default_value"
                shift
            fi
            ;;
        *)
    esac
    case $1 in
        -np|--nproc)
            if [[ -n $2 && $2 != -* ]]; then
                n_proc="$2"
                shift 2
            else
                echo "Warning: Missing value for $1. Using default value: $n_proc"
                shift
            fi
            ;;
        *)
    esac
done

folder="$output_root/3d-rectangular-hopper"

action="mpirun -np $n_proc lethe-particles hopper.prm" 
recreate_folder "$folder"

# Recreate the mesh
gmsh -3 hopper_structured.geo > "$folder/log-mesh"

{ time $action ; } &> "$folder/log"
python3 hopper_post_processing.py -f . --prm hopper.prm --validate

# Copy the information to the log folder
cp $plots $folder
cp $data $folder

# Append the information to the report
magick -density 300 -pointsize 12 text:mass-and-discharge-rate.txt mass-and-discharge-rate.pdf

magick -density 300  $output_root/report.pdf mass-and-discharge-rate.pdf $plots  -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf


