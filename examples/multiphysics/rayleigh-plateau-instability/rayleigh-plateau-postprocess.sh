#!/bin/bash

################################################################################
# Bash script for postprocessing rayleigh-plateau cases
#
# Run with:
#   ./rayleigh-plateau-postprocess.sh <path_to_comparison_data_csv_file>
#   <"{sequence_of_delta_values}"> <-ne>
################################################################################

USAGE=" Usage: $0 <"path_to_comparison_data_csv_file">
        <\"{sequence_of_delta_values}\"> <-ne> \n
        Add the last argument '-ne' only if you wish to avoid extracting again breakup lengths."

# Check the number of input arguments
if [ "$#" -lt 2 ]; then
    echo " Incorrect number of arguments "
    echo -e $USAGE
    exit 1
fi

# Get We and Oh values
case_folder_name=$(ls | grep _delta | head -1)
We_value=$(echo $case_folder_name | grep -oP "We0*\K\d+")
Oh_value=$(echo $case_folder_name | grep -oP "Oh\K\d+_\d+")
Oh_value="${Oh_value//_/.}"

echo "************************************************************"
echo "************************************************************"
echo " POSTPROCESS"
echo "************************************************************"
echo "************************************************************"
echo " Weber number (We):     $We_value"
echo " Ohnesorge number (Oh): $Oh_value"
echo "************************************************************"

# Name of reference data file
reference_data_path="$1"
echo " The reference solution will be taken from: "
echo "    \""$reference_data_path"\""

# Make an array of delta_values
IFS=',' read -ra delta_values_array <<< "$(echo "$2" | tr -d '{}')"

# Check if the lengths have to be extracted
if [ "$#" -eq 3 ]; then
  if [ "$3" = "-ne" ]; then
    echo " Breakup lengths will not be extracted"
  else
    echo " Incorrect number of arguments "
    echo -e $USAGE
    exit 1
  fi
else # only 2 arguments
  # Loop over delta values to evaluate breakup lengths
  for delta_value in ${delta_values_array[@]}
  do
    padded_delta_value=$(printf "%.2f" $delta_value)
    delta_value_without_period=$(echo $padded_delta_value | sed "s/\./_/g")

    echo "************************************************************"
    echo " Postprocessing for:"
    echo " Excitation amplitude: $padded_delta_value"

    # Get current case folder
    current_case_folder=$(ls | grep _delta$delta_value_without_period)
    echo " Case directory: $current_case_folder"

    # Go to current case folder
    cd $current_case_folder

    # Get prm file of the case
    case_parameter_file=$(ls | grep .prm)

    # Launch single case postprocessing to get breakup lengths
    python3 ../rayleigh-plateau-postprocess.py . $case_parameter_file

    # Return to parent directory
    cd ..
  done
fi

# Generate comparison figure
echo "************************************************************"
echo " Generating comparison figure"
python3 rayleigh-plateau-compare.py $reference_data_path ${delta_values_array[@]}

echo "************************************************************"
echo "************************************************************"
echo " Done :D"
echo "************************************************************"
echo "************************************************************"
