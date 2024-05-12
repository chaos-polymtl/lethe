#!/bin/bash
#set -ex
################################################################################
# Bash script for generating Rayleigh-Plateau cases and launching them
#
# Run with:
#   ./rayleigh-plateau-case-generation.sh <"path_to_template_parameter_file">
#   <"{sequence_of_delta_values}">
################################################################################

USAGE="Usage: $0 <"path_to_template_parameter_file">
        <\"{sequence_of_delta_values}\">"

# Check the number of input arguments
if [ "$#" -ne 2 ]; then
    echo "Incorrect number of arguments "
    echo $USAGE
    exit 1
fi


# Parameter filenames
template_parameter_file=$1
template_parameter_file_prefix=$(basename "$template_parameter_file" .tpl)
We_value=$(echo $template_parameter_file | grep -oP "We0*\K\d+")
Oh_value=$(echo $template_parameter_file | grep -oP "Oh\K\d+_\d+")
Oh_value="${Oh_value//_/.}"


echo "************************************************************"
echo "************************************************************"
echo " GENERATE SIMULATION FILES AND LAUNCH"
echo "************************************************************"

# Make an array of delta_values
IFS=',' read -ra delta_values_array <<< "$(echo "$2" | tr -d '{}')"

# Loop over delta values to generate cases
for delta_value in ${delta_values_array[@]}
do
  padded_delta_value=$(printf "%.2f" $delta_value)
  delta_value_without_period=$(echo $padded_delta_value | sed "s/\./_/g")
  case_name=$(echo $template_parameter_file_prefix"_delta"$delta_value_without_period)
  case_parameter_file=$(echo $case_name".prm")

  echo "************************************************************"
  echo " Writing files and launching for:"
  echo "  \"$case_name\""
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
  echo " Dimensionless velocity pertubation (delta_0): $padded_delta_value"
  echo " Weber number (We):                            $We_value"
  echo " Ohnesorge number (Oh):                        $Oh_value"
  echo ""
  echo " Parameter file of the case:"
  echo "   $case_parameter_file"

  # Generate prm file for each case
  sed "s/DELTA_VALUE_OUTPUT/$delta_value_without_period/g" $template_parameter_file > "$case_name/$case_parameter_file"
  sed -i "s/DELTA_VALUE/$delta_value/g" "$case_name/$case_parameter_file"
  echo ""

  # Launch simulation (feel free to use a different number of CPUs)
  cd  $case_name
  mpirun -np 14 lethe-fluid $case_parameter_file

  # Return to parent directory
  cd ..
done

echo "************************************************************"
echo "************************************************************"
echo " Done :D"
echo "************************************************************"
echo "************************************************************"
