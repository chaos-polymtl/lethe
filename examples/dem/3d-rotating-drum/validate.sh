# SPDX-FileCopyrightText: Copyright (c) 2020-2023 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# This script automates the validation process for the 2D Taylor-Couette benchmark
# in Lethe.

# Function to recreate the folder
recreate_folder() {
    local folder_path="$1"

    if [ -d "$folder_path" ]; then
        rm -rf "$folder_path"  # Deletes the folder and its contents
    fi
    mkdir -p "$folder_path"  # Creates the folder
}

# Function to display ... animation
show_animation() {
    local message=$2
    local dots=("")
    while kill -0 $1 2>/dev/null; do
        for i in ".  " ".. " "..." "   "; do
            printf "\r%s%s " "$message" "$i"
            sleep 0.5
        done
    done
    printf "\n" # Ensure a new line after the loop
}

# Store filenames of all plots in a variable (space-seperated)
plots="lethe-rotating-drum-comparison.pdf"

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

folder="$output_root/3d-rotating-drum"
action="mpirun -np $n_proc lethe-particles load-rotating-drum.prm"
recreate_folder "$folder"

{ time $action ; } &> "$folder/log_loading" & pid=$!

show_animation $pid "---> Loading rotating drum"
wait $pid
# Loading done

action="mpirun -np $n_proc lethe-particles rotating-drum.prm"

{ time $action ; } &> "$folder/log_rotating" & pid=$!
show_animation $pid "---> Running simulation"
wait $pid

python3 post_processing_rotating_drum.py -i data_exp.csv -f ./

# Copy the information to the log folder
cp $plots $folder

# Append the information to the report
magick -density 300  $output_root/report.pdf $plots -quality 100 $output_root/temporary.pdf
cp $output_root/temporary.pdf $output_root/report.pdf



