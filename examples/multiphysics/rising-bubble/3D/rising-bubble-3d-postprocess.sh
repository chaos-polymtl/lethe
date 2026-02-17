#!/bin/bash
mkdir -p figures
simulation_folder="."
lethe_path="/home/hepap/lethe/lethe/"

python3 postprocess-rising-bubble-3d.py -p ${simulation_folder}/rising-bubble-3d-AMR-pro-f100-eps06/output/ \
-a ${simulation_folder}/rising-bubble-3d-AMR-pde-f100-eps02/output/ -g ${simulation_folder}/rising-bubble-3d-AMR-geo-f100-eps06/output/ \
-sf ./figures -l $lethe_path -cfl False