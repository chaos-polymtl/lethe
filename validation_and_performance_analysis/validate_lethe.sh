#!/bin/bash
# SPDX-FileCopyrightText: Copyright (c) 2020-2023 The Lethe Authors
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
# - This program but be launched from the main folder of lethe
#
# -----------------------------------------------------------------------------


# Function to display the help message
show_help() {
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  -h, --help       Display this help message"
    echo "  -o, --output     Specify the output directory for validation results"
    echo "  -v, --verbose    Enable verbose mode for detailed output"
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

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            ;;
        -o|--output)
            if [[ -n $2 && $2 != -* ]]; then
                output_path="$2"
                shift 2  # Skip the flag and its argument
            else
                error_exit "Missing argument for -o|--output."
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

echo "The output path has been set to $output_path"
verify_or_create_folder "$output_path"
lethe_path=$(pwd)
echo "Please ensure that the current path is the root folder of lethe"
echo "The current path is: $lethe_path"

echo "Press Enter to continue..."
read -r  # Wait for user to press Enter
echo "Proceeding with validation..."


# Input file
input_file="validation_and_performance_analysis/validation_cases.txt"

# Loop through each line of the file
while IFS= read -r line; do
    # Skip lines that are empty or start with a comment
    if [[ -z "$line" || "$line" =~ ^# ]]; then
        continue
    fi

    # Do something with the line
    echo "Processing: $line"
    cd $lethe_path/$line
    bash validate.sh -p $output_path
    cd $lethe_path
done < "$input_file"


