#!/bin/bash

# Check if number of cases provided (and positive integer)
if [ -z "$1" ] || ! [[ "$1" =~ ^[0-9]+$ ]] || [[ "$1" -eq 0 ]]; then
    echo "Run: bash $0 number_of_cases [--plot] [--clean]"
    echo "  number_of_cases : Required strictly positive integer"
    echo "  --plot          : Optional flag to plot or not the results"
    echo "  --clean         : Optional flag to delete case folders after running"
    exit 1
fi

N_CASES="$1"
CLEANUP=false
PLOT=false

# Check optional arguments
for arg in "$@"; do
    if [ $arg == "--plot" ]; then
        PLOT=true
    elif [ $arg == "--clean" ]; then
        CLEANUP=true
    fi
done

# Generate the cases
echo "Generating $N_CASES cases"
python3 ./generate_cases.py -n "$N_CASES"

# Run all simulation cases
for dir in single_particle_flow_case_*/; do

    echo "Processing $dir"
    
    # Run lethe-particles
    echo "Running lethe-particles for ${dir}"
    lethe-particles "${dir}initial-particle.prm"

    # Run lethe-fluid-vans
    echo "Running lethe-fluid-vans for ${dir}"
    lethe-fluid-vans "${dir}single-particle-flow.prm"

done

# Post-processing
echo "Running post-processing for all cases"
python3 ./postprocessing_all_cases.py -f ./

# Optional plot data
if [ "$PLOT" = true ]; then
    python3 ./plot_data.py -f ./
fi

# Optional cleanup of case folders
if [ "$CLEANUP" = true ]; then
    echo "Cleaning up case directories"
    rm -rf single_particle_flow_case_*/
fi