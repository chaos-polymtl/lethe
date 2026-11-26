#!/bin/bash
# SPDX-FileCopyrightText: Copyright (c) 2024-2026 The Lethe Authors
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
# - Required dependencies: git, and python3 with matplotlib and pypdf
#   (python3 -m pip install --user pypdf).
# - Reference results must be available in the specified directory.
#
# Notes:
# - The program needs to be launched from the main folder of lethe
# - The cases are listed in the input file with their number of cores and the
#   solver they validate. After every case, generate_validation_report.py
#   updates report.pdf in the output folder, from the files listed in
#   report_contents.txt by the validate.sh script of every case.
# - Once all the cases are done, compare_validation_runs.py compares the results
#   with the previous runs stored next to the output folder. Set the environment
#   variable LETHE_VALIDATION_SLACK_WEBHOOK to post the summary to Slack. Also
#   set LETHE_VALIDATION_SLACK_TOKEN (token of a Slack bot with the files:write
#   scope, member of the channel) and LETHE_VALIDATION_SLACK_CHANNEL (ID of the
#   channel) to share report.pdf right after the summary.
#
# -----------------------------------------------------------------------------

# Function to display the help message
show_help() {
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  -h, --help       Display this help message"
    echo "  -i, --input      Specify the input file containing the list of cases to validate (default: validation_cases.txt)"
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

lethe_path=$(pwd)

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
echo "The current path is: $lethe_path"
echo "Please ensure that the current path is the root folder of lethe"

# Check if git is available. Otherwise it is not possible to have a git hash
if ! command -v "git" &> /dev/null; then
  echo "Error: git is not installed or not in the PATH." >&2
  exit 1
fi

# Check if python3 is available. Otherwise it is not possible to post-process the simulations
if ! command -v "python3" &> /dev/null; then
  echo "Error: python3 is not installed or not in the PATH." >&2
  exit 1
fi

# Check if the Python packages used to generate the report are available
if ! python3 -c "import matplotlib, pypdf" &> /dev/null; then
  echo "Error: the Python packages matplotlib and pypdf are required to generate the report." >&2
  echo "Install them with: python3 -m pip install --user matplotlib pypdf" >&2
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
solvers=()

# Loop through each line of the file and store the cases, number of procs and solvers
line_number=0
while IFS= read -r line; do
  line_number=$((line_number + 1))

  # Skip lines that are empty or start with a comment
  if [[ -z "$line" || "$line" =~ ^# ]]; then
      continue
  fi

  # Split the line into three parts
  part1=$(echo "$line" | awk '{print $1}')
  part2=$(echo "$line" | awk '{print $2}')
  part3=$(echo "$line" | awk '{print $3}')

  if [[ -z "$part3" ]]; then
      echo "Error: line $line_number of $input_file does not give the solver of the case: '$line'." >&2
      echo "Every case is given by its path, its number of cores and the solver it validates." >&2
      exit 1
  fi

  # Append the parts to the respective arrays
  cases+=("$part1")
  n_procs+=("$part2")
  solvers+=("$part3")

done < "$input_file"

# Start the list of the content of the report, which the cases complete
report_contents="$output_path/report_contents.txt"
echo "computer $(hostname -f)" > "$report_contents"


echo "Proceeding with validation..."
echo "-----------------------------"
date +"%Y-%m-%d: %H:%M"
for i in "${!cases[@]}"; do
  echo "---> Processing ${cases[i]} using ${n_procs[i]} cores"
  echo "case ${cases[i]} ${n_procs[i]} ${solvers[i]}" >> "$report_contents"
  cd $lethe_path/${cases[i]}
  # Remove any PDF artifacts to ensure that previous post-processing outputs are not used in the generation of the report
  rm -f *.pdf
  # Proceed with the validation case
  bash validate.sh -p $output_path -np ${n_procs[i]}
  echo "     Finished processing ${cases[i]}"
  cd $lethe_path
  # Update the report, so that it is complete up to this case even if the run is interrupted
  python3 "$lethe_path/contrib/validation/generate_validation_report.py" --run "$output_path"
  date +"%Y-%m-%d: %H:%M"
done

# Compare the results with the previous runs stored next to the output folder.
# The summary is posted to Slack if LETHE_VALIDATION_SLACK_WEBHOOK is set, and
# the report is shared with it if LETHE_VALIDATION_SLACK_TOKEN and
# LETHE_VALIDATION_SLACK_CHANNEL are set too.
echo "-----------------------------"
echo "Comparing with the previous validation runs..."
python3 "$lethe_path/contrib/validation/compare_validation_runs.py" --current "$output_path" --repo "$lethe_path"

