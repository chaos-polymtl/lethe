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
Function that calculates the length of the 1st gap from the continuous jet 
after a breakup.

Arguments:
    - fluid0_points : LIST OF LIST OS FLOATS - list of fluid 0 points along 
                      the center axis
    - tol : FLOAT - tolerance to identify the first gap's end

Returns:
    - FLOAT representing the first gap's length
    - FLOAT representing the end of the gap 
'''
def calculate_1st_gap_length(fluid0_points, tol):
    gap_start = fluid0_points[0][0]
    gap_end = fluid0_points[-1][0]

    # Identify the first gap when there are many
    for i, point in enumerate(fluid0_points[:-1]):
        dx = fluid0_points[i+1][0] - point[0]
        if dx > tol:
            gap_end = point[0]
            break
    return gap_end-gap_start, gap_end

'''
Function that identifies if the drop is a satellite drop.

Arguments:
    - df : DATAFRAME object of the current time step
    - dx : FLOAT - length of one side of the sampling box
    - gap_end : FLOAT - x-position of the end of the gap
    - delta : FLOAT - jet perturbation amplitude at the inlet

Returns a boolean. If TRUE, the drop is a satellite drop and the breakup 
    length will not be saved.
'''
def check_drop_is_satellite(df, dx, gap_end, delta):
    half_dx = 0.5 * dx
    box_area = dx**2
    if delta <= 0.5:
        search_box_dimensions = [gap_end, gap_end+dx, -half_dx, half_dx, 0, 0]
    else:
        search_box_dimensions = [gap_end-dx, gap_end, -half_dx, half_dx, 0, 0]
    box = df.clip_box(search_box_dimensions, invert=False)
    fluid1_volume = box.clip_scalar(scalars="filtered_phase", invert = False,
                                    value = 0.5)
    volume_integration_data = fluid1_volume.integrate_data()
    fluid1_volume_value = volume_integration_data["Area"]
    area_ratio = fluid1_volume_value/box_area
    # print("Area ratio: ", area_ratio)

    return (area_ratio < 0.35)

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

# Create a list that holds breakup lengths
Lb_list = []
t_list = []
dt_list = []

# Create beginning and end points for jet length evaluations
x_max = 0.0916
Lb = x_max
previous_gap_length = np.inf
P1 = [0, 0, 0]
P2 = [x_max, 0, 0]
res = 5000

# Jet radius
r_jet = 1.145e-3

# Parameters for identifying breakup
dx_sample = x_max/(res)
tol = dx_sample*10
n_breakup = 1
last_breakup_time = 0

# For satellite drop identification
drop_is_satellite = False
n_satellite = 0
satellite_lb = 1e-10
search_box_dx = r_jet


print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Extracting breakup lengths (Lb)")

# Read PVTU files
for i in range(len(fluids.list_vtu)):
    # Store results in 'df'
    df = fluids.get_df(i)

    # Extract phase values and points over the axis line
    sampled_data = df.sample_over_line(P1, P2, resolution=res)
    phase_over_line = pd.DataFrame(sampled_data["filtered_phase"])
    points_over_line = pd.DataFrame(sampled_data.points)

    # Find the furthest point of fluid 1 before a first gap
    fluid0_points = points_over_line[phase_over_line[0] <= phase_limit].values
    if len(fluid0_points) > 0:
        gap_length, gap_end = calculate_1st_gap_length(fluid0_points, tol)
        # If the first gap reduces in size considerably, it is considered that
        # a droplet is formed.
        if gap_length < previous_gap_length-tol:
            Lb = fluid0_points[0][0]
            t = time_list[i]
            dt = t - last_breakup_time

            # Check if it is a satellite drop
            if n_breakup > 1:
                drop_is_satellite = check_drop_is_satellite(df, search_box_dx, gap_end, float(delta_value))
                if drop_is_satellite:
                    n_satellite += 1
                    if n_satellite == 1:
                        satellite_lb = Lb

            # If the drop is not satellite
            if not drop_is_satellite:
                print("------------------------------------------------------------")
                print(" breakup number:  ", n_breakup)
                print(" time at breakup: ", f"{t:.4f} s")
                print(" Lb:              ", f"{Lb:.4f} m")
                Lb_list.append(Lb)
                t_list.append(t)
                dt_list.append(dt)
                last_breakup_time = t_list[-1]
                n_breakup += 1
        previous_gap_length = gap_length

print("------------------------------------------------------------")
print(f" There were {n_satellite} satellite drops.")
print(f" Their breakup length was approximately {satellite_lb:.4f} m.")
print("------------------------------------------------------------")


#############################################################################
#----------------------------------
# Write csv file with values
#----------------------------------
csv_filename = f"lethe-delta{delta_string}.csv"
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print(f"Writing data into {csv_filename}")
lethe_df = pd.DataFrame({'t': t_list, 'Lb': Lb_list, "dt":dt_list})
lethe_df.to_csv(f"../{csv_filename}", index=False)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Breakup lengths extraction complete for delta =", delta_value)
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

#############################################################################