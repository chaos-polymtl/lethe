"""
Definition of the functions used in static-bubble-multiple-folders.py
"""
# -------------------------------------------
# Modules
# -------------------------------------------

import numpy as np
import pyvista as pv
from natsort import os_sorted
import os
import sys


# --------------------------------------------
# Main
# --------------------------------------------

# This is the function called by static-bubble-multiple-folders.py for every
# folder encountered

def get_pressure_difference(output_path, prm):
    list_vtu = os.listdir(output_path)
    pvd_file = [x for x in list_vtu if (x.endswith('.pvd'))]
    list_vtu = [x for x in list_vtu if ("pvtu" in x)]
    list_vtu = os_sorted(list_vtu)
    if len(pvd_file) == 0:
        raise Exception(f"Folder {output_path} does not contain any .pvd file!")
    last_filename = list_vtu[-1]
    reader = pv.get_reader(output_path + "/" + pvd_file[0])
    # Sort VTU files to ensure they are in the same order as the time step
    sim = pv.read(output_path + "/" + last_filename)
    # Extract pressure field
    pressure = sim["pressure"]

    # Extract coordinates
    x = sim.points[:, 0]
    y = sim.points[:, 1]

    # Compute distance from the center of the bubble
    distance = np.sqrt(x ** 2 + y ** 2)

    # Find indices of points inside and outside the bubble
    inside_indices = np.where(distance < prm.radius - 0.05)[0]
    outside_indices = np.where(distance > prm.radius + 0.05)[0]

    # Calculate pressure inside and outside
    avg_pressure_inside = np.mean(pressure[inside_indices])
    avg_pressure_outside = np.mean(pressure[outside_indices])

    # Compute the pressure difference
    pressure_difference = avg_pressure_inside - avg_pressure_outside
    return pressure_difference
    
# Returns the pressure along the [-2.5,0,0]->[2.5,0,0] line segment and the line segment itself to allow for an easy plot.
def get_pressure_slice(output_path):
    list_vtu = os.listdir(output_path)
    pvd_file = [x for x in list_vtu if (x.endswith('.pvd'))]
    list_vtu = [x for x in list_vtu if ("pvtu" in x)]
    list_vtu = os_sorted(list_vtu)
    last_filename = list_vtu[-1]
    reader = pv.get_reader(output_path + "/" + pvd_file[0])
    # Sort VTU files to ensure they are in the same order as the time step
    sim = pv.read(f"{output_path}/{last_filename}")
    origin = [-2.5,0,0]
    end = [2.5,0,0]
    line = sim.sample_over_line(origin, end, resolution=1000)
    pressure_line = line["pressure"]
    x_line = line.points[:,0]
    
    return pressure_line, x_line

# Returns the analytical pressure difference for a given set of parameters (surface tension coefficient)
def analytical_solution(prm):
    radius_array = np.linspace(0.1, 0.5, 500)
    sigma = prm.sigma
    return radius_array, sigma / radius_array




