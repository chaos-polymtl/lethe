# SPDX-FileCopyrightText: Copyright (c) 2023 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
IMPORTS
"""

from lethe_pyvista_tools import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
import multiprocessing
import pyvista as pv

"""
FUNCTION
"""

def process_file(prm_name):
    
    example = lethe_pyvista_tools('.', prm_name, 'out.pvd', n_procs=1, ignore_data = ['velocity_divergence', 'q_criterion', 'subdomain'])

    radius = 0.5
    pos_x  = 8
    pos_y  = 0
    pos_z  = 0

    number_of_times = len(example.list_vtu)
    time_step       = example.step
    time_array      = np.zeros(number_of_times)
    temperature     = np.zeros(number_of_times)

    for i_t in range(number_of_times):
        current_slice = example.get_df(i_t)
        x_positions = current_slice.points[:, 0]
        y_positions = current_slice.points[:, 1]
        z_positions = current_slice.points[:, 2]
        
        # We select only the solid cylinder geometry
        circle_sfd = ((x_positions-pos_x)**2 + (y_positions-pos_y)**2 + (z_positions-pos_z)**2)**0.5 - radius
        current_slice["circle_sdf"] = circle_sfd
        current_slice_inside = current_slice.clip_scalar(scalars="circle_sdf",invert=True, value=0.0)
        
        # Average the temperature on the sphere
        integrated_data  = current_slice_inside.integrate_data()
        surface          = integrated_data["Area"]
        temperature_mean = integrated_data["temperature"]/surface
        
        print(example.time_list[i_t])
        time_array [i_t] = example.time_list[i_t]
        temperature[i_t] = temperature_mean

    return time_array, temperature

"""
POSTPROCESSING
"""

# We construct the data frame then export to csv
prm_names = ['cylinder-xi-0_01.prm', 'cylinder-xi-1.prm', 'cylinder-xi-100.prm']
csv_names = ['cylinder-xi-0_01.csv', 'cylinder-xi-1.csv', 'cylinder-xi-100.csv']

i = 0
for prm_name in prm_names:

    df = pd.DataFrame()
    time_array, temperature = process_file(prm_name)

    df["time_array: " + prm_name[:-4]]      = time_array
    df["temperature: " + prm_name[:-4]]     = temperature

    # Export to csv
    df.to_csv('postprocessed-avg-temperature-' + csv_names[i], sep='\t', encoding='utf-8')
    i += 1

print("End")


