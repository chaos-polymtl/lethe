#!/bin/bash
# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

################################################################################
# Bash script for checking simulations on the cluster that have timed-out
# (TIMEOUT) or failed (FAILED) and relaunching them if wished.
#
# * IMPORTANT *
#   - Make sure that you have the proper permissions to execute the file before
#     calling it. If not, you can use "chmod" to change the files permissions.
#
# Assumed structure and prerequisites:
# - Each simulation must be contained in its own folder.
# - The simulation folder must contain:
#   - A parameter file with the suffix ".prm". There can only be one parameter
#     file in the folder. The name of the parameter files can be different from
#     one simulation folder to another.
#   - A job script for launching the simulation on the cluster.
# - All parameter files must contain already the "restart" subsection if they
#   have to be restarted from a checkpoint.
# - All job scripts should have the same name even if their content is different
#   (e.g. launch.sh or job.sh)
# - All job output log files must have the suffix ".out" and have the job ID
#   as part of their name.
# - The cluster must have the "seff" command. Alliance Canada clusters have it.
#
# Tree example 1:
#   .
#   ├── relaunch_simulations.sh
#   ├── <search_file>
#   ├── <simulation_00_folder>
#   │   ├── <job_file>
#   │   └── <parameter_file>
#   ├── <simulation_01_folder>
#   │   ├── <job_file>
#   │   └── <parameter_file>
#   ├── <simulation_02_folder>
#   │   ├── <job_file>
#   │   └── <parameter_file>
#   ...
#
#   The content of the search file would look like (1 folder = 1 line in
#   the search file):
#     simulation_00_folder
#     simulation_01_folder
#     simulation_02_folder
#     ...
#
# The script can also be called with nested folders as long the relative paths
# (originating from the call location) in the search file is well defined.
#
# Tree example 2:
#   .
#   ├── relaunch_simulations.sh
#   ├── <search_file>
#   ├── <subfolder_00>
#   │   ├── <simulation_00_folder>
#   │   │   ├── <job_file>
#   │   │   └── <parameter_file>
#   │   ├── <simulation_01_folder>
#   │   │   ├── <job_file>
#   │   │   └── <parameter_file>
#   │   ├── <simulation_02_folder>
#   │   │   ├── <job_file>
#   │   │   └── <parameter_file>
#   │   ...
#   ├── <subfolder_01>
#   │   ├── <simulation_00_folder>
#   │   │   ├── <job_file>
#   │   │   └── <parameter_file>
#   │   ├── <simulation_01_folder>
#   │   │   ├── <job_file>
#   │   │   └── <parameter_file>
#   │   ├── <simulation_02_folder>
#   │   │   ├── <job_file>
#   │   │   └── <parameter_file>
#   │   ...
#   ...
#
#   The content of the search file would look like (1 folder = 1 line in
#   the search file):
#     ./subfolder_00/simulation_00_folder
#     ./subfolder_00/simulation_01_folder
#     ./subfolder_00/simulation_02_folder
#     ...
#     ./subfolder_01/simulation_00_folder
#     ./subfolder_01/simulation_01_folder
#     ./subfolder_01/simulation_02_folder
#     ...
#
#  * Note that in the examples above, "..." is added to signify the presence
#    of more folders.
#
# How the script works:
# - A search file (-sf|--search_file) and job file name (-j|--job_file) must be
#   specified when calling the script. The search file contains the path to the
#   simulation folders to be checked.
# - The script checks the status of the latest launched job of the simulations
#   in the search file.
#   - If the status is "CANCELLED", check the previous job until a
#     non-"CANCELLED" simulation is found.
#     This way, user-stopped simulations are ignored.
# - If the job was found to be "TIMEOUT":
#   - The simulation folder name is saved in a file named
#     "timeout_simulations.txt".
#   - If the user requested that the simulation should be rerun from checkpoint
#     with (-lr|--launch_restart), then the scripts searches for the "set restart"
#     in the "subsection restart" of the parameter file and changes its value to
#     "true". Then, it launches the simulation with the specified job file.
#     The launched simulation folder name is saved into a file named
#     "relaunched_simulations.txt". The idea is to reuse this file as the
#     "search_file" when relaunching again. At the second iteration, to avoid
#     overwriting "relaunched_simulations.txt" and losing track of the previous
#     relaunch, a copy can be made by adding the argument (-cp|--copy_search_file)
#     to the script call.
# - If the job was found to be "FAILED":
#   - The simulation folder name is saved in a file named
#     "failed_simulations.txt".
#   - If the user requested that the failed simulation should be rerun with
#     (-rf|--rerun_failed), then the simulation is rerun without restart.
#   - If the user wishes to rerun the simulation from the last checkpoint, the
#     (-st|--set_restart_true) flag must be added to the call.
#   - All relaunched simulations, from a checkpoint or not, are listed in the
#     "relaunched_simulations.txt" file.
#
# For help on the different arguments and their default values, the user can
# call the script with the flag (-h|--help).
#
# How to call the script:
# - Call example for only finding the simulations that are incomplete
#   (TIMOUT or FAILED):
#     ./relaunch_simulations -sf "./simulations_list.txt" -j "launch.sh"
# - Call example for only finding the simulations that are incomplete
#   (TIMOUT or FAILED) and launching TIMEOUT ones:
#     ./relaunch_simulations -sf "./simulations_list.txt" -j "launch.sh" -lr
# - Call example for only finding the simulations that are incomplete
#   (TIMOUT or FAILED) and launching TIMEOUT and FAILED ones. The FAILED ones
#   are launched without restart:
#     ./relaunch_simulations -sf "./simulations_list.txt" -j "launch.sh" -lr -rf
# - Call example for only finding the simulations that are incomplete
#   (TIMOUT or FAILED) and launching TIMEOUT and FAILED ones. Where, the FAILED
#   ones are launched with restart enabled :
#     ./relaunch_simulations -sf "./simulations_list.txt" -j "launch.sh" -lr -rf -st
# - Call example for relaunching from the previously relaunched list. Additionally,
#   make a copy of the search file for comparing changes between launches:
#     ./relaunch_simulations -sf "./relaunched_simulations.txt" -j "launch.sh" -lr -rf -st -cp
#     or
#     ./relaunch_simulations -j "launch.sh" -lr -rf -st -cp
#   since "./relaunched_simulations.txt" is the default value for (-sf|--search_file).
#
################################################################################

# Function that prints the message when no restart subsection is found
print_no_restart_message() {
  local prm_file="$1"
  echo " No \"subsection restart\" was found in $prm_file"
  echo " The simulation cannot be restarted."
  echo " Please add the subsection and relaunch the simulation with "
  echo " checkpointing enabled."
  echo " Aborting script execution."
  echo "************************************************************"
}

# Function that prints the message when checkpoint is disabled
print_disabled_checkpoint_message() {
  local prm_file="$1"
  local simulation_folder="$2"
  echo " Checkpoint was not enabled in ${simulation_folder}/${prm_file}"
  echo " Cannot restart without checkpoint"
  echo " Aborting script execution."
  echo "************************************************************"
}

# Function that enables the restart feature in the parameter file
enable_restart() {
  local prm_file="$1"
  local simulation_folder="$2"

  # Check if there is a restart subsection in the prm file
  if grep -q "subsection.*restart" "$prm_file"
  then
    # Check if the simulation was previously checkpointed
    if grep -q "set checkpoint\s*=\s*true" "$prm_file"
    then
      # Check if there is already a "set restart" line
      if grep -q "set restart" "$prm_file"
      then
        # Set restart to "true"
        sed -i 's/\(set restart\s*=\s*\)false/\1true/' "$prm_file"
      else
        # Add "set restart" line with argument "true"
        sed -i '/set checkpoint\s*=\s*true/ {h;p;g;s/checkpoint/restart/}' "$prm_file"
      fi
    else
      # Checkpoint is not enabled
      print_disabled_checkpoint_message "$prm_file" "$simulation_folder"
      exit 1
    fi
  else
    # No restart subsection was found
    print_no_restart_message "$prm_file" "$simulation_folder"
    exit 1
  fi
}

################################################################################
USAGE="Usage: \n
             $0 [-sf <path_to_files_with_folders_to_be_searched>] \n
             [-j <slurm_job_file_filename>] [-cp] [-lr] [-rf] [-st] \n
            \n
             -sf|--search_file: Path to the file with the folders to be searched (default: ./relaunched_simulations.txt)\n
                                This file must be a \".txt\" text file containing only a list of folder names of the simulations to check.\n
             -j|--job_file: Name of the slurm job file within each simulation folder (default: job.sh)\n
             -cp|--copy_search_file: Flag indicates to copy the search file with the suffix \"_old\" (deactivated by default)\n

             For TIMEOUT simulations:
             -lr|--launch_restart: Flag indicates to launch the identified timeout simulations (deactivated by default)\n

             For FAILED simulations:
             -rf|--rerun_failed: Flag that indicates to rerun failed simulations too (deactivated by default)\n
             -st|--set_restart_true: Flag that indicates to the restart argument of FAILED simulations to \"true\" \n
                                     (deactivated by default, meaning, restart is set to \"false\" in the parameter files.\n
             "

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -sf|--search_file)
      search_file="$2"
      shift 2
      ;;
    -j|--job_file)
      job_file="$2"
      shift 2
      ;;
    -cp|--copy_search_file)
      copy_search_file="true"
      shift
      ;;
    -lr|--launch_restart)
      launch_restart="true"
      shift
      ;;
    -rf|--rerun_failed)
      rerun_failed="true"
      shift
      ;;
    -st|--set_restart_true)
      set_restart_true="true"
      shift
      ;;
    -h|--help)
      echo -e "$USAGE"
      exit 0
      ;;
    *)
      echo "Unknown parameter: $1"
      echo -e "$USAGE"
      exit 1
      ;;
  esac
done

# Assign default values if not set
search_file="${search_file:-./relaunched_simulations.txt}"
job_file="${job_file:-job.sh}"
copy_search_file="${copy_search_file:-"false"}"
launch_restart="${launch_restart:-"false"}"
rerun_failed="${rerun_failed:-"false"}"
set_restart_true="${set_restart_true:-"false"}"

# Check if the search file exists
if [ ! -f "$search_file" ]; then
  echo "File not found: $search_file"
  exit 1
fi

# Print values to user
echo "************************************************************"
echo -e " Running $0 with parameters: \n"
echo " Search file:                           $search_file"
echo " Slurm job file:                        $job_file"
echo " Copy search file:                      $copy_search_file"
echo " Launch timeout simulations:            $launch_restart"
echo " Rerun failed simulations:              $rerun_failed"
echo " Enable restart of failed simulations:  $set_restart_true"
echo "************************************************************"
echo "************************************************************"

# Copy search file if requested. Otherwise, keep it as a temporary copy.
# The temporary copy is needed since "relaunched_simulations.txt" is overwritten.
search_file_prefix=$(basename "$search_file" .txt)
search_file_copy=$(echo "${search_file_prefix}_old.txt")
cp "$search_file" "$search_file_copy"
if [ "$copy_search_file" == "true" ]
then
  echo " Search file copied as $search_file_copy"
  echo "************************************************************"
fi

################################################################################
# List all folders to be checked
echo -e " SIMULATIONS TO BE CHECKED: \n"
for simulation_folder in $(cat $search_file)
do
  echo " $simulation_folder"
done
echo ""
################################################################################
# Logging files
relaunched_simulations_file="relaunched_simulations.txt"
found_timeout_file="timeout_simulations.txt"
found_failed_file="failed_simulations.txt"

# Clear logging files
rm -f "$found_timeout_file" "$found_failed_file" "$relaunched_simulations_file"

while IFS= read -r simulation_folder
do
  simulation_folder=$(echo "$simulation_folder" | tr -d '\r'| sed 's/[[:space:]]*$//')
  echo "------------------------------------------------------------"
  echo " Checking simulation:"
  echo "  $simulation_folder"

  # Move into folder if it exists
  if [ -d "$simulation_folder" ]
  then
    cd "$simulation_folder"

    # Loop over job output files
    while IFS= read -r file
    do
      job_id=$(echo "$file" | grep -o '[0-9]\+' | tail -n 1)
      echo  " Job ID: $job_id"
      cancelled_string=$(seff $job_id|grep CANCELLED)

      # If the job was not cancelled, check if it is TIMEOUT or FAILED
      if [ -z "$cancelled_string" ]
      then
        timeout_string=$(seff $job_id|grep TIMEOUT)
        failed_string=$(seff $job_id|grep FAILED)

        # Check if the status is TIMEOUT
        if [ -n "$timeout_string" ]
        then
          echo "  Simulation job status: TIMEOUT"
          echo "$simulation_folder" >> "../$found_timeout_file"

          # Restart simulation if requested
          if [ "$launch_restart" == "true" ]
          then
            echo " Launching simulation restart"

            # Get prm file
            prm_file=$(ls | grep "\.prm$")
            echo " Found prm file: ${simulation_folder}/${prm_file}"

            # Enable the restart feature in the parameter file
            enable_restart "$prm_file" "$simulation_folder"

            # Check if the job file exists
            if [ ! -f "$job_file" ]; then
              echo "File not found: ${simulation_folder}/${job_file}"
              exit 1
            fi

            # Launch job
            sbatch -J $simulation_folder $job_file
            echo "$simulation_folder" >> "../$relaunched_simulations_file"
          fi

        # Check if the simulation FAILED
        elif [ -n "$failed_string" ]
        then
          echo "  Simulation job status: FAILED"
          echo "$simulation_folder" >> "../$found_failed_file"

          # Rerun or restart simulation if requested
          if [ "$rerun_failed" == "true" ]
          then
            echo " Re-running failed simulation"

            # Check if restart parameter should be set to "true"
            if [ "$set_restart_true" == "true" ]
            then
              # Get prm file
              prm_file=$(ls | grep "\.prm$")
            echo " Found prm file: ${simulation_folder}/${prm_file}"

              # Enable the restart feature in the parameter file
              enable_restart "$prm_file" "$simulation_folder"
            fi

            # Launch job
	          sbatch -J $simulation_folder $job_file
	          echo "$simulation_folder" >> "../$relaunched_simulations_file"
          fi
        else
          echo "  The job status is not CANCELLED, TIMEOUT or FAILED."
          echo "  No action taken."
        fi
        break
      else
        echo "  The job was cancelled. Checking previous job..."
        continue
      fi
    done  < <(ls *.out 2>/dev/null | sort -r) # Sort logging files

    # Go back to parent folder
    cd - >/dev/null || exit 1
  else
    echo " Directory not found: $simulation_folder"
    echo " Aborting script execution."
    echo "************************************************************"
    exit
  fi
done < "$search_file_copy"

# Delete temporary copy if not requested
if [ "$copy_search_file" == "false" ]
then
  rm -f "$search_file_copy"
fi

echo "************************************************************"
echo "************************************************************"
echo " Done :D"
echo "************************************************************"
echo "************************************************************"
