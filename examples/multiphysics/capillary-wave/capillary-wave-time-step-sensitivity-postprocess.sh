#!/bin/bash

##################################################################################################################################
# Bash script for computing capillary wave analytical solution and postprocessing cases
#
# Run with: ./capillary-wave-time-step-sensitivity.sh <path_to_analytical_data_csv_file> <"{sequence_of_time-step_multipliers}">
##################################################################################################################################

# Check the number of input arguments
if [ "$#" -lt 2 ]; then
    echo " Incorrect number of arguments "
    echo " Usage: $0 <path_to_analytical_data_csv_file> <"{sequence_of_time-step_multipliers}"> <-sa>"
    echo " Add the last argument '-sa' only if you wish to solve the analytical expression."
    exit 1
fi

# Capillary time-step constraint and simulation end time
calculation_output_file=calculations.output
capillary_time_step=$(grep "The capillary constraint is:" "$calculation_output_file" | awk '{printf $5}')
capillary_time_step=$(awk "BEGIN { printf \"%.10f\", $capillary_time_step }")
end_time=$(grep "End time:" "$calculation_output_file" | awk '{printf $3}')
end_time=$(awk "BEGIN { printf \"%.6f\", $end_time }")
echo "*****************************************************"
echo "*****************************************************"
echo " POSTPROCESS"
echo "*****************************************************"
echo "*****************************************************"
echo " Simulation end time: $end_time s"
echo " Capillary time-step constraint: $capillary_time_step s"
echo "*****************************************************"


# Name of analytical data file
analytical_data_path="$1"

# Generate analytical data if wanted
if [ "$#" -eq 3 ]; then
  if [ "$3" = "-sa" ]; then
    echo " The analytical solution will be saved under: "
    echo "    \""$analytical_data_path"\""
    python3 capillary-wave-prosperetti-solution.py $analytical_data_path
  else
    echo " Incorrect usage"
    echo " Usage: $0 <path_to_analytical_data_csv_file> <"{sequence_of_time-step_multipliers}"> <-sa>"
    echo " Add the last argument '-sa' only if you wish to solve the analytical expression."
    exit 1
  fi
elif [ -e $analytical_data_path ]; then
  echo " The analytical solution will be taken from: "
  echo "    \""$analytical_data_path"\""
else
  echo " The analytical solution file doesn't exist, to generate it, run with:"
  echo " $0 <path_to_analytical_data_csv_file> <"{sequence_of_time-step_multipliers}"> <-sa>"
  echo " The last argument '-sa' indicates that you wish to solve the analytical expression."
  exit 1
fi

# Make an array of time-step_multiplier
IFS=',' read -ra time_step_multiplier_array<<< "$(echo "$2" | tr -d '{}')"

# Loop over the different cases
for time_step_multiplier in ${time_step_multiplier_array[@]}
do
  echo "*****************************************************"
  echo " Postprocessing for:"
  echo " Time-step constraint multiplier: $time_step_multiplier"
  time_step=$(echo "$capillary_time_step * $time_step_multiplier" | bc)

  echo " Time-step:                       $time_step s"

  # Postprocess with Python
  time_step_multiplier_without_period=$(echo $time_step_multiplier | sed "s/\./_/g")
  python3 capillary-wave-postprocess.py . capillary-wave-TSM-$time_step_multiplier_without_period.prm $analytical_data_path
done

# Comparison figure
echo "*****************************************************"
echo " Generating comparison figure"
python3 capillary-wave-combined.py $analytical_data_path ${time_step_multiplier_array[@]}

echo "*****************************************************"
echo "*****************************************************"
echo " Done postprocessing cases!"
echo "*****************************************************"
echo "*****************************************************"
