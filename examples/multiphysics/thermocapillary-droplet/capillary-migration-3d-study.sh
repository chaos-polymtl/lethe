# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#!/bin/bash

###############################################################################################
# Bash script for generating capillary migration 3d cases and launching them
###############################################################################################

USAGE="Usage: $0 [-t <path_to_template_parameter_file>] \n
            \t[-j <path_to_template_job_file>] \n
            \t[-p <path_to_the_csv_file_with_case_parameters>] \n
            \t[-c <\"sequence_of_case_numbers\">] \n
            \t[-rd <true/false>] \n
            \n
             -t: Path to template parameter file (default: ./capillary-migration-3d-geo.tpl) \n
             -j: Path to template job file for launching simulations on clusters (default: ./launch.tpl) \n
             -p: Path to CSV file with case parameters (default: ./capillary-migration-3d-cases.csv) \n
             -c: Case numbers (default: \"6-8\") \n
             -rd: Delete previously generated folders for the specified cases before running simulations <true/false> (default: false) \n
             -h: Show help \n
             \n
            The  <\"sequence_of_case_numbers\"> argument corresponds to the case numbers in the csv file.\n
            \t The numbers should be separated with a comma (',') delimiter. \n
            \t Ranges can also be specified using a hyphen ('-'). \n
            \t For example, \"0-5, 9, 30-32\" includes cases 0, 1, 2, 3, 4, 5, 9, 30, 31, and 32."

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -t|--template_prm)
      template_parameter_file="$2"
      shift 2
      ;;
    -j|--template_job)
      template_job_file="$2"
      shift 2
      ;;
    -p|--parameters_csv_file)
      cases_csv_file="$2"
      shift 2
      ;;
    -c|--cases)
      case_sequence="$2"
      shift 2
      ;;
    -rd|--remove_directories)
      remove_directories="$2"
      shift 2
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
template_parameter_file="${template_parameter_file:-./capillary-migration-3d.tpl}"
template_job_file="${template_job_file:-./launch.tpl}"
cases_csv_file="${cases_csv_file:-./capillary-migration-3d-cases.csv}"
case_sequence="${case_sequence:-"6-8"}"
remove_directories="${remove_directories:-"false"}"

# Check if the files exists
if [ ! -f "$template_parameter_file" ]; then
  echo "File not found: $template_parameter_file"
  exit 1
fi

if [ ! -f "$template_job_file" ]; then
  echo "File not found: $template_job_file"
  exit 1
fi

if [ ! -f "$cases_csv_file" ]; then
  echo "File not found: $cases_csv_file"
  exit 1
fi

# Print values to user
echo "Template parameter file:  $template_parameter_file"
echo "Template job file:        $template_job_file"
echo "CSV case parameter files: $cases_csv_file"
echo "Cases:                    $case_sequence"
echo "Remove directories:       $remove_directories"

# Job file for launching on cluster
job_file=$(basename "$template_job_file" .tpl)
job_file=$job_file".sh"

# Function to expand ranges specified in the <"case_sequence"> argument
expand_case_numbers() {
  local input="$1"
  local result=()
  input="${input//[\{\}]/}"  # Remove curly braces
  input="${input//[[:space:]]/}" # Remove whitespaces

  IFS=',' read -ra tokens <<< "$input"
  for token in "${tokens[@]}"; do
    if [[ "$token" =~ ^[0-9]+-[0-9]+$ ]]; then
      IFS='-' read -r start end <<< "$token"
      for ((i=start; i<=end; i++)); do
        result+=("$i")
      done
    elif [[ "$token" =~ ^[0-9]+$ ]]; then
      result+=("$token")
    else
      echo "Invalid range or number: $token"
      exit 1
    fi
  done

  echo "${result[@]}"
}

# Make an array of case numbers
IFS=',' read -ra case_numbers_array <<< "$(expand_case_numbers "$case_sequence")"

# Check if directories have to be deleted
if [ "$#" -eq 5 ]; then
  if [ "$5" = "-rd" ]; then
    echo "************************************************************"
    echo "************************************************************"
    echo " DELETING EXISTING CASE DIRECTORIES"
    echo "************************************************************"
    echo "************************************************************"

    # Loop over specified cases
    for case_number in ${case_numbers_array[@]}
    do
      # Get case name
      padded_case_number=$(printf "%02d" $case_number)
      line=$(grep "^${padded_case_number}," "$cases_csv_file")
      IFS=',' read -r -a case_parameters <<< "$line"
      case_name=${case_parameters[1]}

      # Delete directory if it exists
      if [ -d "$case_name" ]; then
          rm -rf "$case_name"
          echo " Directory deleted:"
          echo "   $case_name"
      fi
    done
  else
    echo " Incorrect arguments "
    echo -e $USAGE
    exit 1
  fi
fi


echo "************************************************************"
echo "************************************************************"
echo " GENERATE SIMULATION FILES AND LAUNCH"
echo "************************************************************"

# Loop over specified cases
for case_number in ${case_numbers_array[@]}
do
    # Get parameters of the case and store in appropriate variables
    padded_case_number=$(printf "%02d" "${case_number#0}")
    line=$(grep "^${padded_case_number}," "$cases_csv_file")
    IFS=',' read -r -a case_parameters <<< "$line"

    case_name=${case_parameters[1]}
    reinitialization_type=${case_parameters[2]}
    reinitialization_frequency=${case_parameters[3]}
    epsilon=${case_parameters[4]}
    reinitialization_distance=${case_parameters[5]}
    tanh_thickness=${case_parameters[6]}
    tanh_thickness="${tanh_thickness//[$'\t\r\n ']/}"

  echo "************************************************************"
  echo " Writing files and launching for:"
  echo "  \"$case_name\" (case: $padded_case_number)"
  echo ""

  # Check if data directory of the case already exists before creating
  if [ ! -d "$case_name" ]; then
      # If the directory doesn't exist, create it
      mkdir -p "$case_name"
      echo " Directory created:"
      echo "   $case_name"
  else
      echo " Directory already exists:"
      echo "   $case_name"
  fi
  echo ""

  echo " Parameter values:"
  echo "  REINITIALIZATION_TYPE:      $reinitialization_type"
  echo "  REINITIALIZATION_FREQUENCY: $reinitialization_frequency"
  echo "  EPSILON:                  $epsilon"
  echo "  REINITIALIZATION_DISTANCE:  $reinitialization_distance"
  echo "  TANH_THICKNESS:           $tanh_thickness"
  echo ""

  # Create parameter file for the case
  case_parameter_file=$(echo $case_name".prm")
  sed "s/CASE_NAME/$case_name/g" $template_parameter_file > "$case_name/$case_parameter_file"
  sed -i "s/REINITIALIZATION_TYPE/$reinitialization_type/g" "$case_name/$case_parameter_file"
  sed -i "s/REINITIALIZATION_FREQUENCY/$reinitialization_frequency/g" "$case_name/$case_parameter_file"
  sed -i "s/EPSILON/$epsilon/g" "$case_name/$case_parameter_file"
  sed -i "s/REINITIALIZATION_DISTANCE/$reinitialization_distance/g" "$case_name/$case_parameter_file"
  sed -i "s/TANH_THICKNESS/$tanh_thickness/g" "$case_name/$case_parameter_file"


  # Generate job file
  sed "s/CASE_PRM_FILE/$case_parameter_file/g" $template_job_file > "$case_name/$job_file"

  # Move into case directory
  cd $case_name

  # Launch simulation (feel free to use a different number of CPUs)
  # sbatch -J $case_name $job_file # On cluster
  mpirun -np 12  ~/work/lethe/build_dev/applications/lethe-fluid/lethe-fluid $case_parameter_file

  # Return to parent directory
  cd ..
done

echo "************************************************************"
echo "************************************************************"
echo " Done :D"
echo "************************************************************"
echo "************************************************************"



