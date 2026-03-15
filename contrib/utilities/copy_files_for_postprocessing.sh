#!/usr/bin/env bash
# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

################################################################################
# For help the user can call the script with the flag (-h|--help).
################################################################################

USAGE="
 Bash script for copying files of interest in a specific folder.

 This script can be used for identifying and copying files of importance before
 transferring to another machine or simply to clean up irrelevant files. This
 can be especially useful for parametric studies with more than necessary
 outputs for postprocessing.

 This script identifies .vtu and .pvtu files corresponding to requested
 time-step indices and copies them. If none is requested, by default, it
 copies the last time-step files.
 This script also copies files of the output folder with the following
 extensions: .dat and .pvd
 This script copies files from the simulation folder with the following
 extensions: .out, .csv, .sh, and .prm
 This script copies files from the \$base folder with the following
 extensions: .sh, .tpl, .csv, and .txt

 * IMPORTANT *
  - When specifying \$base and \$destination use absolute paths to avoid path
    issues.
  - The parameter files associated with each simulation must be present in their
    respective simulation folders.


 How to call the script:
 - Default call; it copies the .vtu and .pvtu files of the simulation:
   ./copy_files_for_postprocessing.sh -b <path_to_base_folder> -d <path_to_destination_folder>
 - Call with specific time-step indices to be copied (0 to 20, 100 and last):
   ./copy_files_for_postprocessing.sh -b <path_to_base_folder> -d <path_to_destination_folder> -si 0-20,100,end
 
 
 Required arguments:
 -b, --base <path>
 Path to the base directory containing all simulation folders (see Tree example
 below).
 Each simulation must be contained in its own subdirectory.
 
 -d, --destination <path>
 Path to the directory where selected files will be copied.
 If the directory does not exist, it will be created.

 Optional arguments:
 -si, --step-index <list>
 Specify which VTU time-step indices to copy. This corresponds to the digits
 that appear at the end of the .pvtu filename (before the extension).
 
  Supported formats with examples:
    2309            → single step
    100,500,2309    → list of steps
    1000-1010       → range of steps
    0,500,1000-1010 → mixed list
    end             → last available step
    0-20,end        → range plus last step
 
  *** If this option is not provided, the script copies the last
  time step available in the .pvd file. ***
 
 -rmd, --delete_destination_folder
 Delete the destination directory before copying files.
 The directory will be recreated automatically.
 
 -h, --help
 Print this help message and exit.


 Assumed folder structure:
 - Each simulation case must be contained in its own folder.
 - The simulation folder must contain:
   - A parameter file with the extension '.prm'. There can only be one parameter
     file in the folder. The name of the parameter files can be different from
     one simulation folder to another.

 Tree example of \$base folder:
  base_directory
  ├── simulation_00
  │   ├── log.out
  │   ├── parameter_file.prm
  │   ├── job_script.sh
  │   ├── ...
  │   └── output_folder_00
  │       ├── simulation_00.00000.00000.vtu
  │       ├── simulation_00.00000.00001.vtu
  │       ├── ...
  │       ├── simulation_00.02309.00000.vtu
  │       ├── simulation_00.02309.00000.vtu
  │       ├── simulation_00.02309.00000.vtu
  │       ├── ...
  │       ├── simulation_00.02309.pvtu
  │       ├── simulation_00.pvd
  │       ├── ...
  │       ├── mass_conservation.dat
  │       ├── ...
  ├── simulation_01
  │   ├── log.out
  │   ├── parameter_file.prm
  │   ├── job_script.sh
  │   ├── ...
  │   └── output_folder_01
  │       ├── simulation_01.00000.00000.vtu
  │       ├── ...
  ├── ...
  ├── parameter_file.tpl
  ├── job_script.tpl
  ├── parametric_sweep.sh
  ├── ...


 How the script works:
 - The script loops over subdirectories of \$base to find simulation cases and
   copies files of relevance.
 - The name of the output directory for VTK and .dat files is taken from the
   .prm file.
 - For each simulation, the script goes through the .pvd file to identify the
   requested time-step indices or by default the last time step and copies the 
   corresponding and .pvtu files.

"
################################################################################

# Parse script arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -b|--base)
      base="$2"
      shift 2
      ;;
    -d|--destination)
      destination="$2"
      shift 2
      ;;
    -si|--step-index)
      requested_step_indices="$2"
      shift 2
      ;;
    -rmd|--delete_destination_folder)
      delete_destination_folder="true"
      shift
      ;;
    -h|--help)
      echo "$USAGE"
      exit 0
      ;;
    *)
      echo "Unknown parameter: $1"
      echo "$USAGE"
      exit 1
      ;;
  esac
done

# Assign default value if not set
delete_destination_folder="${delete_destination_folder:-"false"}"
requested_step_indices="${requested_step_indices:-""}"

# Check that required arguments are set
if [[ -z "$base" ]]; then
    echo "Error: Base directory is not set. Use -b or --base to specify it."
    echo "$USAGE"
    exit 1
fi

if [[ -z "$destination" ]]; then
    echo "Error: Destination directory is not set. Use -d or --destination to specify it."
    echo "$USAGE"
    exit 1
fi

# Check if the base directory exists
if [[ ! -d "$base" ]]; then
    echo "The base directory specified ($base) does not exist"
    exit 1
fi

# Handle deletion if requested
if [[ "$delete_destination_folder" == "true" ]]; then
    if [[ -e "$destination" ]]; then
        echo "Deleting existing destination folder: $destination"
        rm -rf "$destination"
    fi
    echo "Recreating destination folder: $destination"
    mkdir -p "$destination"
else
    # Create destination folder if it does not exist
    if [[ ! -e "$destination" ]]; then
        echo "Creating destination folder: $destination"
        mkdir -p "$destination"
    fi
fi

################################################################################
# Handles one simulation subdirectory
# Arguments: $1: Simulation name, $2: Simulation directory,
#            $3: Destination directory
################################################################################
process_simulation() {
    local sim_name="$1"
    local dir="$2"
    local dst="$3"

    #----------------------
    # PRINT VALUES TO USER
    #----------------------
    echo "************************************************************"
    echo " Simulation name:       $sim_name"
    echo " Simulation directory:  $dir"
    echo " Destination directory: $dst/$sim_name"

    # Create copy simulation directory
    mkdir -p "${dst}/${sim_name}/"

    # Ensure that globbing produces an empty list instead of a literal pattern
    shopt -s nullglob

    #-------------------------------------
    # READ OUTPUT DIRECTORY FROM PRM FILE
    #-------------------------------------
    local prm_file
    prm_file=$(find "$dir" -maxdepth 1 -name "*.prm" -print -quit)

    if [[ -n "$prm_file" ]]; then
        local output_file_name pvd_base_name
        output_file_name=$(grep -i "set output path" "$prm_file" | awk -F= '{print $2}' | xargs)
        output_file_name=$(basename "$output_file_name")
        # If "output path" not set in prm, give default value
        if [[ -z "$output_file_name" ]]; then
          output_file_name="./"
        fi
        pvd_base_name=$(grep -i "set output name" "$prm_file" | awk -F= '{print $2}' | xargs)
        pvd_base_name=$(basename "$pvd_base_name")
        # If "output name" not set in prm, give default value
        if [[ -z "$pvd_base_name" ]]; then
          pvd_base_name="out"
        fi
    else
        echo " Warning: Cannot find output path. Skipping to next simulation folder."
        return 0
    fi
    echo " Output directory:      $output_file_name"

    #---------------
    # FIND PVD FILE
    #---------------
    # Find the .pvd file in the simulation output directory
    local pvd_file=$(echo "$dir/$output_file_name/$pvd_base_name.pvd")

    if [[ ! -f "$pvd_file" ]]; then
        echo "Error: No PVD file found in $dir/$output_file_name"
        return 1
    fi
    echo " Found PVD file: $pvd_file"

    #--------------------------------
    # DETERMINE LAST TIME-STEP INDEX
    #--------------------------------
    # Extract last time step and corresponding pvtu file from the .pvd
    last_line=$(grep "<DataSet " "$pvd_file" | tail -n1)
    last_file=$(echo "$last_line" | grep -oP '(?<=file=")[^"]+')
    last_step=$(echo "$last_file" | awk -F. '{print $(NF-1)}')

    echo " Last time-step index: $last_step"

    #----------------------------
    # BUILD TIME-STEP INDEX LIST
    #----------------------------
    local step_list=()

    if [[ -n "$requested_step_indices" ]]; then
      IFS=',' read -ra step_indices <<< "$requested_step_indices"
        for step_index in "${step_indices[@]}"; do
          if [[ "$step_index" == "end" ]]; then
            step_list+=("$last_step")
          elif [[ "$step_index" =~ ^[0-9]+-[0-9]+$ ]]; then
            start=${step_index%-*}
            end=${step_index#*-}
            for ((i=start;i<=end;i++)); do
              step_list+=("$i")
            done
          elif [[ "$step_index" =~ ^[0-9]+$ ]]; then
            step_list+=("$step_index")
          else
            echo " Warning: invalid time-step index '$step_index'"
          fi
        done
      else
        # If not specified, get last-time step index
        step_list=("$last_step")
    fi

    # Remove duplicates
    step_list=($(printf "%s\n" "${step_list[@]}" | sort -n | uniq))

    #----------------
    # COPY VTK FILES
    #----------------
    # Create output directory in destination
    mkdir -p "${dst}/${sim_name}/${output_file_name}"

    # Loop through steps
    for step in "${step_list[@]}"; do
      printf -v step_padded "%05d" "$((10#$step))"
      echo ""
      echo " Copying step $step_padded: "
      pvtu_files=( "$dir/$output_file_name/"*".${step_padded}.pvtu" )

      if (( ${#pvtu_files[@]} == 0 )); then
        echo "  Step not found"
        continue
      fi

      pvtu=$(basename "${pvtu_files[0]}")

      echo "  Corresponding pvtu file: $pvtu"

      cp "$dir/$output_file_name/$pvtu" \
         "$dst/$sim_name/$output_file_name/"
      cp "$dir/$output_file_name/"*".${step_padded}."*.vtu \
         "$dst/$sim_name/$output_file_name/" 2>/dev/null || true
    done

    #-----------------------------------------------------
    # COPY OTHER OUTPUT FILES OF INTEREST (.dat and .pvd)
    #-----------------------------------------------------
    if [[ -d "$dir/${output_file_name}" ]]; then
        cp "$dir/${output_file_name}"/*.dat \
           "${dst}/${sim_name}/${output_file_name}/" 2>/dev/null || true
        cp "$dir/${output_file_name}"/*.pvd \
           "${dst}/${sim_name}/${output_file_name}/" 2>/dev/null || true
    fi

    #----------------------------------
    # COPY DIRECTORY FILES OF INTEREST
    #----------------------------------
    cp "$dir"/*.out "${dst}/${sim_name}/" 2>/dev/null || true
    cp "$dir"/*.csv "${dst}/${sim_name}/" 2>/dev/null || true
    cp "$dir"/*.sh  "${dst}/${sim_name}/" 2>/dev/null || true
    cp "$dir"/*.prm "${dst}/${sim_name}/" 2>/dev/null || true
}

#----------------------------
# LOOP OVER SIMULATION CASES
#----------------------------
for simulation in "$base"/*/; do
    simulation=${simulation%/}  # remove trailing slash
    sim_name=$(basename "$simulation")
    process_simulation "$sim_name" "$simulation" "$destination"
done

#----------------------------------------------
# COPY FILES OF INTEREST FROM PARENT DIRECTORY
#----------------------------------------------
cp "$base"/*.sh  "${destination}" 2>/dev/null || true
cp "$base"/*.tpl "${destination}" 2>/dev/null || true
cp "$base"/*.csv "${destination}" 2>/dev/null || true
cp "$base"/*.txt "${destination}" 2>/dev/null || true

echo "************************************************************"
echo "************************************************************"
echo " Done :D"
echo "************************************************************"
echo "************************************************************"