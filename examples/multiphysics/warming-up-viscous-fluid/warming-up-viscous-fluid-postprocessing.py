# SPDX-FileCopyrightText: Copyright (c) 2023 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
IMPORTS
"""
from lethe_pyvista_tools import *
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

"""
FUNCTION
"""
def process_temperature(params, prm_folder, prm_file_name, output_folder, pvd_file_name):
    
    temperature = lethe_pyvista_tools(prm_folder, prm_file_name, pvd_file_name)
    line = pv.Line(np.array([0, 0.5, 0]), np.array([0.5, 0.5, 0]), resolution = 100)

    time_list = temperature.time_list
    
    max_temp_list  = np.zeros(len(temperature.list_vtu))
    min_temp_list  = np.zeros(len(temperature.list_vtu))
    mean_temp_list = np.zeros(len(temperature.list_vtu))

    sim_error = np.zeros(len(temperature.list_vtu))

    for i in range(len(temperature.list_vtu)):
        sim = pv.read(output_folder + "/" + temperature.list_vtu[i])
        
        max_temp_list[i]  = max(sim["temperature"])
        min_temp_list[i]  = min(sim["temperature"])
        mean_temp_list[i] = np.mean(sim["temperature"])

        position_on_sampled_line    = line.sample(sim)["Distance"]
        temperature_on_sampled_line = line.sample(sim)["temperature"]

        # Calculate the analytical solution
        analytical_temperature = params.Tw+(((params.rho*params.nu)*params.v*params.v)/(2*params.K))*(1-(position_on_sampled_line/params.B)*(position_on_sampled_line/params.B))
        sim_error[i] = np.sqrt(np.sum((analytical_temperature - temperature_on_sampled_line)**2))
    
    return time_list, max_temp_list, min_temp_list, mean_temp_list, sim_error

"""
PARAMETERS
"""
prm_folder = "."
prm_file_name = "warming-up-viscous-fluid"

output_folder = "output"
pvd_file_name = "warming-up.pvd"

class params:
    rho = 0.9
    nu  = 0.5
    K   = 0.12
    Tw  = 80
    v   = 2
    B   = 0.5

"""
MAIN
"""
time_list, max_temp_list, min_temp_list, mean_temp_list, sim_error = process_temperature(params, prm_folder, prm_file_name, output_folder, pvd_file_name)

"""
PLOTS FOR TEMPERATURE
"""
plt.figure(figsize=(8, 6))
plt.plot(time_list, mean_temp_list, color='green', label='Mean Temperature')
plt.fill_between(time_list, min_temp_list, max_temp_list, color='green', alpha=0.3, label='Temperature envelope')

# Labeling
plt.title('Temperature envelope and mean value', fontsize=14)
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Temperature (C)', fontsize=12)
plt.ylim(0, 90)
plt.xlim(0, 7)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()

"""
PLOTS FOR ERROR
"""
plt.figure(figsize=(8, 6))
plt.title('Error on the temperature according to time', fontsize=14)
plt.plot(time_list, sim_error, color='red', label='L2 Error')
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('L2 norm of the error on temperature (-)', fontsize=12)
plt.xlim(0, 7)
plt.grid(True, linestyle='--', alpha=0.5)

plt.show()