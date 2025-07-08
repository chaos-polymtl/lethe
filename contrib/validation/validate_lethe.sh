#!/bin/bash
# SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# -----------------------------------------------------------------------------
# Script Name: validate_lethe.sh
# 
# Description:
# This script automates the validation of the Computational Fluid Dynamics (CFD) 
# and Discrete Element Method (DEM) software Lethe. It runs a series of predefined
# test cases and generates a validation report.  The script ensures that new updates 
# or configurations do not introduce unintended errors or deviations in 
# the software's behavior using validatoin benchmarks
#
# Usage:
# ./validate_lethe.sh [options]
# 
# Options:
#   -h, --help       Display this help message
#   -o, --output     Specify the output directory for validation results
#   -v, --verbose    Enable verbose mode for detailed output
#
# Requirements:
# - Lethe must be installed and accessible from the environment.
# - Required dependencies: [list dependencies, if any].
# - Reference results must be available in the specified directory.
#
# Notes:
# - The program needs to be launched from the main folder of lethe
#
# -----------------------------------------------------------------------------

# Function to display the help message
show_help() {
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  -h, --help       Display this help message"
    echo "  -i, --input      Specify the input file with the list of cases to validate (default: validation_cases.txt)"
    echo "  -o, --output     Specify the absolute path directory for validation results"
    echo "  -f, --force      Do not ask for confirmation when erasing the contents of the output folder or proceeding with the validation" 
    echo 
    echo "Description:"
    echo "This script automates the validation of the Computational Fluid Dynamics (CFD)"
    echo "and Discrete Element Method (DEM) software Lethe. It runs a series of predefined"
    echo "test cases and generates a validation report. The script ensures that new updates"
    echo "or configurations do not introduce unintended errors or deviations in the software's"
    echo "behavior using validation benchmarks."
    echo
    exit 0
}

# Function to display a splash screen
show_splash() {
    echo "#############################################"
    echo "#                                           #"
    echo "#      Lethe Validation Script Started      #"
    echo "#                                           #"
    echo "#############################################"
    echo
}

# Function to verify a folder exists, create it if necessary, or erase its contents if specified
verify_or_create_folder() {
    local folder_path="$1"

    if [ -d "$folder_path" ]; then
        echo "Folder '$folder_path' exists."
        echo -n "Do you want to erase its contents? [y/N]: "
        read -r response
        case "$response" in
            [yY][eE][sS]|[yY])
                echo "Erasing contents of '$folder_path'..."
                if rm -rf "${folder_path:?}"/*; then
                    echo "Contents of '$folder_path' erased successfully."
                else
                    echo "Error: Failed to erase contents of '$folder_path'. Exiting."
                    exit 1
                fi
                ;;
            *)
                echo "Preserving existing contents of '$folder_path'."
                ;;
        esac
    else
        echo "Folder '$folder_path' does not exist. Attempting to create it..."
        if mkdir -p "$folder_path"; then
            echo "Successfully created folder '$folder_path'."
        else
            echo "Error: Failed to create folder '$folder_path'. Exiting."
            exit 1
        fi
    fi
}

show_splash

# Default values for options
input_file="$lethe_path/contrib/validation/validation_cases.txt"
force=false # Flag to skip confirmation prompts

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            ;;
        -f|--force)
            force=true  # Set a flag to skip confirmation prompts
            shift  # Skip the flag
            ;;
        -i|--input)
            if [[ -n $2 && $2 != -* ]]; then
                input_file="$2"
                shift 2  # Skip the flag and its argument
            else
                echo "Error: Missing argument for -i|--input."
                exit 1
            fi
            ;;
        -o|--output)
            if [[ -n $2 && $2 != -* ]]; then
                output_path="$2"
                shift 2  # Skip the flag and its argument
            else
                echo "Error: Missing argument for -o|--output."
                exit 1
            fi
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if the output path was provided
if [[ -z $output_path ]]; then
    echo "Output path not specified. Use the -o or --output flag to set it."
    exit 1
fi

if [[ "$output_path" != /* ]]; then
    echo "Error: The path '$output_path' is not an absolute path." >&2
    exit 1
fi

echo "The output path has been set to $output_path"
verify_or_create_folder "$output_path"
lethe_path=$(pwd)
echo "The current path is: $lethe_path"
echo "Please ensure that the current path is the root folder of lethe"

# Check if magick is available. Otherwise reports will not be generated.
# If magick is not available, the script will exit with an error status.
# Check if the application exists
if ! command -v "magick" &> /dev/null; then
  echo "Error: magick is not installed or not in the PATH." >&2
  exit 1
fi

# Check if git is available. Otherwise it is not possible to have a git hash
if ! command -v "git" &> /dev/null; then
  error exit echo "Error: git is not installed or not in the PATH." >&2
  exit 1
fi

# Check if python3 is available. Otherwise it is not possible to post-process the simulations
if ! command -v "python3" &> /dev/null; then
  error exit echo "Error: python3 is not installed or not in the PATH." >&2
  exit 1
fi

# Output file
hash_file="$output_path/current_git_hash.txt"

# Get the current Git hash
if git rev-parse --is-inside-work-tree > /dev/null 2>&1; then
    git_hash=$(git rev-parse HEAD)
    echo "Current Git hash: $git_hash"
    echo "$git_hash" > "$hash_file"
    echo "Git hash saved to '$hash_file'."
else
    echo "Error: Not inside a Git repository."
    exit 1
fi

if ($force); then
    echo "Skipping confirmation prompts as --force was specified."
else
    echo "Press Enter to continue..."
    read -r  # Wait for user to press Enter
fi

cases=()
n_procs=()

# Loop through each line of the file and store the cases and number of procs
while IFS= read -r line; do

  # Skip lines that are empty or start with a comment
  if [[ -z "$line" || "$line" =~ ^# ]]; then
      continue
  fi

  # Split the line into two parts
  part1=$(echo "$line" | awk '{print $1}')
  part2=$(echo "$line" | awk '{print $2}')
  
  # Append the parts to the respective arrays
  cases+=("$part1")
  n_procs+=("$part2")

done < "$input_file"

#Create the first page of the PDF report
echo -e "Lethe validation \nDate: $(date) \nGit hash: $git_hash  \nComputer: $(hostname -f) \n" | \
magick -density 300 -pointsize 12 text:- $output_path/report.pdf


echo "Proceeding with validation..."
echo "-----------------------------"
date +"%Y-%m-%d: %H:%M"
for i in "${!cases[@]}"; do
  echo "---> Processing ${cases[i]} using ${n_procs[i]} cores"
  cd $lethe_path/${cases[i]}
  bash validate.sh -p $output_path -np ${n_procs[i]}
  echo "     Finished processing ${cases[i]}"
  cd $lethe_path
  date +"%Y-%m-%d: %H:%M"
done

