"""
Postprocessing code for controlling the volume of the phases in the 3d bubble detachment case.
"""
import numpy
# -------------------------------------------
# Modules
# -------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv
from natsort import os_sorted
from alive_progress import alive_bar
from alive_progress.styles import showtime
import multiprocessing

import os

number_of_procs = 5


def analytical_volume(t, flux, v_init):
    return t * flux + v_init


def computed_volume(t, flux, v_init):
    volume = np.empty(len(t))
    volume[0] = v_init
    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]
        volume[i] = -(flux[i] + flux[i - 1]) * dt * 0.5 + volume[i - 1]
    return volume


# Functions to turn .dat data in numpy array
# Credits to Lucka Barbeau for the two functions below !
def is_number(s):
    try:
        float(s)
    except ValueError:
        return False
    return True


def read_my_data(results_path):
    force_file = open(results_path, 'r')
    list_of_list_of_vars_name = [[]];
    list_of_list_of_vars = [[]];
    fix_vars_name = False

    nb_set_of_vars = 0;
    for line in force_file.readlines():
        i = 0;
        for word in line.split():
            if is_number(word):
                fix_vars_name = True
                list_of_list_of_vars[nb_set_of_vars][i] = np.append(
                    list_of_list_of_vars[nb_set_of_vars][i],
                    float(word))
                i += 1
            else:
                if word != 'particle':
                    if fix_vars_name:
                        nb_set_of_vars += 1
                        list_of_list_of_vars_name.append([])
                        list_of_list_of_vars.append([])
                        fix_vars_name = False
                    list_of_list_of_vars_name[nb_set_of_vars].append(word)
                    list_of_list_of_vars[nb_set_of_vars].append(np.array([]))
    return list_of_list_of_vars_name, list_of_list_of_vars


def bubble_volume(output_path):

    list_of_list_of_vars_name, list_of_list_of_vars = read_my_data(
        output_path + "/phase_statistics.dat")
    list_of_list_of_vars_name_flux, list_of_list_of_vars_flux = read_my_data(
        output_path + "/flow_rate.dat")

    absolute_volume = list_of_list_of_vars[0][6]
    time_integrated_flux_volume = computed_volume(
        list_of_list_of_vars_flux[0][0], list_of_list_of_vars_flux[0][3],
        list_of_list_of_vars[0][6][0])
    analytical_volume_array = analytical_volume(list_of_list_of_vars[0][0],
                                                5.0e2,
                                                (2 * np.pi * (0.5) ** 3) / 3)
    time = list_of_list_of_vars[0][0]

    return absolute_volume, time_integrated_flux_volume, analytical_volume_array, time


def get_volume_at_time(volume_list,time_list,time):
    time_list = np.array(time_list)
    time_index = (np.abs(time_list - time)).argmin()

    return volume_list[time_index]

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

    time_clip = -1
    x, y = get_contour_at_detachment(detachment_index,
                                     list_vtu, output_path)

    return time, detachment_time, x, y


def get_contour_at_detachment(index_detachment, list_vtu, output_path):
    vtu_file = list_vtu[index_detachment]
    sim = pv.read(f"{output_path}/{vtu_file}")
    sim.set_active_scalars("phase_order")
    slice_single = sim.slice(normal="z", origin=(0, 0, 0))
    contour_val = np.array([0.0])
    contour = slice_single.contour(contour_val, scalars="phase_order")
    x, y = contour.points[:, 0], contour.points[:, 1]
    return x, y


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



def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
    
