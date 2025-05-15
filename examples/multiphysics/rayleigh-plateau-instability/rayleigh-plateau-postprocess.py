# SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#############################################################################
"""
Postprocessing code for the Rayleigh-Plateau example.
This code extracts breakup lengths of a continuous jet for a given case and
saves the data in a csv file.
"""
#############################################################################

#############################################################################
'''Importing Libraries'''
import numpy as np
import pandas as pd
import sys
sys.path.append("$LETHE_PATH/contrib/postprocessing/")
from lethe_pyvista_tools import *
#############################################################################

#############################################################################
# Run script: python3 path_to_rayleigh-plateau-postprocess.py
#             path_to_case prm_filename
#############################################################################
# Check the number of input arguments
if len(sys.argv)!= 3:
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
          "Incorrect number of arguments\n"
          "Run script with: \n"
          "\t python3 path_to_rayleigh-plateau-postprocess.py path_to_case prm_filename\n "
          "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    exit(1)

#############################################################################

''' Function definitions '''

'''
Function that determines breakup times and lengths through the jet length 
variation.

Arguments:
    - fluids: lethe_pyvista_tools object containing all postprocessing 
               information and simulation results
    - phase_limit: FLOAT phase fraction value representing the interface
    
Returns lists of breakup times and lengths.
'''
def get_breakup_times_and_lengths(fluids, phase_limit):
    # Initialize breakup time and length lists
    tb_list = []
    Lb_list = []
    # Get active times
    time_list = fluids.time_list

    # Number of breakups
    n_breakup = 0

    # Initial domain length
    x_max = 0.0916

    # Jet radius
    r_jet = 1.145e-3

    # Initialize variables for breakup identification
    sample_box_area = 2*r_jet**2
    previous_x_jet = 1e-10
    tol = r_jet*10


    for i in range(len(fluids.list_vtu)):
        # Store results in 'df'
        df = fluids.get_df(i)

        # Identify breakups
        interface = df.contour([phase_limit], scalars="filtered_phase")
        interface_blocks = interface.connectivity()
        region_ids = np.unique(interface_blocks["RegionId"])

        # Isolate continuous jet region
        jet = interface_blocks.clip_scalar(scalars="RegionId",
                                           invert=True, value=0)

        # Find breakup length
        Lb = np.max(jet.points[:, 0])

        # Compute jet length variation to identify if a breakup occurs
        x_jet_variation = (Lb-previous_x_jet)/previous_x_jet

        # Only evaluate breakup length when there is a new breakup
        if ((x_jet_variation + tol < 0) and (abs(Lb-x_max)/x_max > 0.05)):

            # Find RegionId corresponding to the first drop
            drop_index = find_drop_region_id(interface_blocks, region_ids)

            drop = interface_blocks.threshold(value=(drop_index,drop_index),
                                              scalars="RegionId",
                                              invert=False)

            # Store breakup time and length value only if the drop is not
            # satellite
            if (drop_is_not_satellite(df, drop, sample_box_area, phase_limit)):
                n_breakup += 1
                tb = time_list[i]
                Lb_list.append(Lb)
                tb_list.append(tb)
                print("------------------------------------------------------------")
                print(" breakup number:  ", n_breakup)
                print(" time at breakup: ", f"{tb:.4f} s")
                print(" Lb:              ", f"{Lb:.4f} m")

        previous_x_jet = Lb

    return tb_list, Lb_list


'''
Function that identifies the RegionID of the first drop from left.

Arguments:
    - interface_blocks: DATAFRAME object of the current time step
    - region_ids: LIST of INT representing the ids of closed loops
    
Returns an INT corresponding to the RegionID of the first breakup
'''
def find_drop_region_id(interface_blocks, region_ids):

    region_min_list = []
    for j in range(len(region_ids)):
        exec(f"drop_{j} = interface_blocks.threshold(value=({region_ids[j]},{region_ids[j]}), scalars='RegionId', invert=False)")
        exec(f"region_min_list.append(np.min(drop_{j}.points[:,0]))")

    # Remove jet value to identify the first drop from left
    jet_value = np.min(region_min_list)
    truncated_region_min_list = region_min_list.copy()
    truncated_region_min_list.remove(jet_value)
    drop_min_value = np.min(truncated_region_min_list)

    return region_min_list.index(drop_min_value)


'''
Function that identifies if the drop is a satellite drop by checking the area 
ratio of a drop inside a sample box. This sample box is chosen to be a square 
with sides of 2*r_jet m.

Arguments:
    - df: DATAFRAME object of the current time step
    - drop_df: DATAFRAME object of the drop contour
    - sample_box_area: FLOAT area of the sample box
    - phase_limit: FLOAT phase fraction value representing the interface

Returns a BOOLEAN. If TRUE, the drop is a satellite drop and the breakup
    length will not be saved.
'''
def drop_is_not_satellite(df, drop_df, sample_box_area, phase_limit):
    # Make a search box around the drop
    x, y = drop_df.points[:,0], drop_df.points[:,1]
    search_box_dimensions = [np.min(x), np.max(x), np.min(y), np.max(y), 0, 0]
    box = df.clip_box(search_box_dimensions, invert=False)

    # Evaluate drop area
    fluid1_volume = box.clip_scalar(scalars="filtered_phase", invert = False,
                                    value = phase_limit)
    drop_integration_data = fluid1_volume.integrate_data()
    current_drop_area = drop_integration_data["Area"]

    # Evaluate area ratio
    area_ratio = current_drop_area/sample_box_area
    return (area_ratio > 1)


#############################################################################
#----------------------------------
# Read, extract and compute
# quantities of interest
#----------------------------------
# Parse arguments
simulation_path = sys.argv[1]
prm_file_name = sys.argv[2]

# Name of the pvd file and extract delta value
delta_value = prm_file_name.split("delta", 1)
delta_string = delta_value[-1].split(".prm", 1)[0]
delta_value = delta_string.replace('_', '.')
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Postprocessing for excitation amplitude =",delta_value)
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
pvd_file_name = f"{prm_file_name.split('-J1', 1)[0]}.pvd"

# Create the fluids object
fluids = lethe_pyvista_tools(simulation_path, prm_file_name, pvd_file_name)

# Set phase_limit to search for height values
phase_limit = 0.5

# Get active times
time_list = fluids.time_list

# Extract breakup times and lengths
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Extracting breakup lengths (Lb)")

tb_list, Lb_list = get_breakup_times_and_lengths(fluids, phase_limit)


#############################################################################
#----------------------------------
# Write csv file with values
#----------------------------------
csv_filename = f"lethe-delta{delta_string}.csv"
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print(f"Writing data into {csv_filename}")
lethe_df = pd.DataFrame({'t': tb_list, 'Lb': Lb_list})
lethe_df.to_csv(f"../{csv_filename}", index=False)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Breakup lengths extraction complete for delta =", delta_value)
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

#############################################################################
