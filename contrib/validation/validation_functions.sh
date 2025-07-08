# SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

# Function to recreate the folder
recreate_folder() {
    local folder_path="$1"

    if [ -d "$folder_path" ]; then
        rm -rf "$folder_path"  # Deletes the folder and its contents
    fi
    mkdir -p "$folder_path"  # Creates the folder
}

# Parse command-line arguments
# Stores both the file name and the number of processors
parse_command_line() {
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
    shift
  done
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