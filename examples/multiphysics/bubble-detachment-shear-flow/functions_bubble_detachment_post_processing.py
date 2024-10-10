# SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing code for controlling the volume of the phases in the 3d bubble detachment case.
"""
# -------------------------------------------
# Modules
# -------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
from natsort import os_sorted
from alive_progress import alive_bar
from alive_progress.styles import showtime

import os

# Returns the analytical bubble volume with a given theoretical flux of air at bubble inlet.
def analytical_volume(t, flux, v_init):
    return t * flux + v_init

# Returns the computed bubble volume using the simulation data on the flux at the bubble inlet. It is different from the volume obtained by integrating the phase field in the domain.
def computed_volume(t, flux, v_init):
    volume = np.empty(len(t))
    volume[0] = v_init
    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]
        volume[i] = -(flux[i] + flux[i - 1]) * dt * 0.5 + volume[i - 1]
    return volume


# Returns the list of bubble volume for each time-step (analytical and computed and absolute) using the simulation's post-processing files.
def bubble_volume(output_path):

    # This volume is the one obtained by integrating the phase field in the domain
    time, absolute_volume = np.loadtxt(output_path + "/      phase_statistics.dat",skiprows=1,usecols=(0,6),unpack=True)
    
    volume_flux = np.loadtxt(output_path + "/flow_rate.dat",skiprows=1,usecols=3,unpack=True)
    
        
    # This volume is the one obtained by integrating the numerical inlet flux with respect to time
    time_integrated_flux_volume = computed_volume(
        time, volume_flux,absolute_volume[0])
        
    # This volume is the one obtained by integrating the theoretical inlet flux with respect to time
    analytical_volume_array = analytical_volume(time,
                                                5.0e2,
                                                (2 * np.pi * (0.5) ** 3) / 3)  


    return absolute_volume, time_integrated_flux_volume, analytical_volume_array, time

# Returns the volume of the bubble at a given time
def get_volume_at_time(volume_list,time_list,time):
    time_list = np.array(time_list)
    time_index = (np.abs(time_list - time)).argmin()

    return volume_list[time_index]

# Returns the time-step list, the detachment time and x and y coordinates of the contour of the bubble at detachment time.
def get_numerical_detachment_time(output_path):
    list_vtu = os.listdir(output_path)

    if os.path.isfile(
            output_path + '/detachment_index.npy') and os.path.isfile(
            output_path + '/time.npy'):
        print('Reading the contour area file from an external source')
        detachment_index = np.load(output_path + '/detachment_index.npy')
        pvd_file = [x for x in list_vtu if (x.endswith('.pvd'))]
        reader = pv.get_reader(output_path + "/" + pvd_file[0])
        time = np.array(reader.time_values)
        detachment_time = time[detachment_index]
        time = np.load(output_path + '/time.npy')

    else:
        print('No existing contour area file')
        pvd_file = [x for x in list_vtu if (x.endswith('.pvd'))]
        if len(pvd_file) == 0:
            raise Exception(f"Folder {output_path} is empty!")
        list_vtu = [x for x in list_vtu if ("pvtu" in x)]
        list_vtu = os_sorted(list_vtu)
        reader = pv.get_reader(output_path + "/" + pvd_file[0])
        time = np.array(reader.time_values)

        detachment_index = np.zeros(1)
        with alive_bar(len(time), refresh_secs=0) as bar:
            detachment_time = 0.
            detachment_index = int(0)
            for i, vtu_file in enumerate(list_vtu):
                # print(f'Processing file {i+1} out of {len(list_vtu)} files')
                # Sort VTU files to ensure they are in the same order as the time step
                sim = pv.read(f"{output_path}/{vtu_file}")
                # Extract pressure field
                if "phase_order" not in sim.point_data.keys():
                    print(f"Skipping {vtu_file} : corrupted")
                    continue
                sim.set_active_scalars("phase_order")
                contour_val = np.array([0.0])
                contour = sim.contour(contour_val, scalars='phase_order')
                contour_connectivity = contour.connectivity()
                regions_ids = np.unique(contour_connectivity['RegionId'])
                if len(regions_ids) > 1 and detachment_index == int(0) and detachment_time<1e-9 :
                    detachment_index = i-1
                    detachment_time = time[detachment_index]
                    print(f'Detachment time found! t_det={detachment_time}')
                    np.save(output_path + '/time.npy', time)
                    np.save(output_path + '/detachment_index.npy',
                            detachment_index)
                    break
                bar()
                                     
    x, y = get_contour_at_fixed_time(detachment_time, output_path)

    return time, detachment_time, x, y

# Returns the x and y coordinates of the bubble contour in the z=0 plane at a given time.
def get_contour_at_fixed_time(time_value, output_path):
    list_vtu = os.listdir(output_path)
    pvd_file = [x for x in list_vtu if (x.endswith('.pvd'))]
    
    if len(pvd_file) == 0:
        return [], []
        
    list_vtu = [x for x in list_vtu if ("pvtu" in x)]
    list_vtu = os_sorted(list_vtu)
    reader = pv.get_reader(output_path + "/" + pvd_file[0])

    time = np.array(reader.time_values)
    time_index = (np.abs(time - time_value)).argmin()

    vtu_file = list_vtu[time_index]
    
    sim = pv.read(f"{output_path}/{vtu_file}")
    sim.set_active_scalars("phase_order")
    slice_single = sim.slice(normal="z", origin=(0, 0, 0))
    contour_val = np.array([0.0])
    contour = slice_single.contour(contour_val, scalars="phase_order")
    x, y = contour.points[:, 0], contour.points[:, 1]
    return x, y
    
