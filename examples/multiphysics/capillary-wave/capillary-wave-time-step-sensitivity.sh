#!/bin/bash

###############################################################################################
# Bash script for generating capillary wave cases and launching them
#
# Run with: ./capillary-wave-time-step-sensitivity.sh <"{sequence_of_time_step-multipliers}">
###############################################################################################

# Check the number of input arguments
if [ "$#" -ne 1 ]; then
    echo "Incorrect number of arguments "
    echo "Usage: $0 <"{sequence_of_time_step-multipliers}">"
    exit 1
fi

echo "*****************************************************"
echo "*****************************************************"
echo " GENERATE SIMULATION FILES AND LAUNCH"

# Capillary time-step constraint and simulation end time
echo "*****************************************************"
echo "*****************************************************"
echo " Capillary wave calculations"
echo "*****************************************************"
calculation_output_file=calculations.output # Name of the file containing the calculation results
python3 capillary-wave-calculation.py "./$calculation_output_file"
capillary_time_step=$(grep "The capillary constraint is:" "$calculation_output_file" | awk '{printf $5}')
capillary_time_step=$(awk "BEGIN { printf \"%.10f\", $capillary_time_step }")
end_time=$(grep "End time:" "$calculation_output_file" | awk '{printf $3}')
end_time=$(awk "BEGIN { printf \"%.6f\", $end_time }")

echo "*****************************************************"
echo "*****************************************************"
echo " Simulation end time: $end_time s"
echo " Capillary time-step constraint: $capillary_time_step s"
echo "*****************************************************"

# Make an array of time-step_multiplier
IFS=',' read -ra time_step_multiplier_array <<< "$(echo "$1" | tr -d '{}')"

for time_step_multiplier in ${time_step_multiplier_array[@]}
do
  echo "*****************************************************"
  echo " Writing files and launching for:"
  echo " Time-step constraint multiplier:  $time_step_multiplier"
  number_of_outputs=100 # You may change the number of outputs if you wish
  time_step=$(echo "$capillary_time_step * $time_step_multiplier" | bc)
  output_frequency=$(echo "$end_time / $time_step / $number_of_outputs" | bc)

  echo " Time-step:                        $time_step s"
  echo " Output frequency:                 $output_frequency"

  # Generate prm file for every case
  time_step_multiplier_without_period=$(echo $time_step_multiplier | sed "s/\./_/g") # If you wish to use decimal multipliers, this will avoid problems when postprocessing
  sed "s/TIMESTEPMULTIPLIER/$time_step_multiplier_without_period/g" capillary-wave.prm > capillary-wave-TSM-$time_step_multiplier_without_period.prm
  sed -i "s/TIMESTEP/$time_step/g" capillary-wave-TSM-$time_step_multiplier_without_period.prm
  sed -i "s/ENDTIME/$end_time/g" capillary-wave-TSM-$time_step_multiplier_without_period.prm
  sed -i "s/OUTPUTFREQUENCY/$output_frequency/g" capillary-wave-TSM-$time_step_multiplier_without_period.prm

  # Launch simulation (feel free to use a different number of CPUs)
  mpirun -np 4 lethe-fluid capillary-wave-TSM-$time_step_multiplier_without_period.prm
done

echo "*****************************************************"
echo "*****************************************************"
echo " Done!"
echo "*****************************************************"
echo "*****************************************************"
