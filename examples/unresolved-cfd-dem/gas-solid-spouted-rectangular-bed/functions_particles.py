# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import numpy as np
import pandas as pd
from matplotlib.patches import Patch

######################################################################
# Functions for Post-Processing

def group_velocities_by_x(H, df_positions, velocity_magnitude, x_centers, dx, dy):
    """
    Select particles in a window around height H and group them by x-location.
    
    Parameters:
    - H: target y height
    - df_positions: DataFrame with columns ['x', 'y', 'z']
    - velocity_magnitude: 1D array of velocity magnitudes
    - x_centers: array of target x values
    - dx: half-width of the x-window
    - dy: half-width of the y-window
    
    Returns:
    - Dictionary mapping (H, x_center) -> list of velocity magnitudes in z for that height and x-center
    """
    grouped = {}
    mask_y = (df_positions['y'] >= H - dy) & (df_positions['y'] <= H + dy)

    # Filter positions based on y height
    filtered_df = df_positions[mask_y]
    filtered_vel = velocity_magnitude[mask_y.to_numpy()]

    for x in x_centers:
        # Select particles within dx of the current x center
        mask_x = (filtered_df['x'] >= x - dx) & (filtered_df['x'] <= x + dx)
        selected = filtered_vel[mask_x.to_numpy()]
        grouped[(H, x)] = selected.tolist()

    return grouped 

def welford(H, grouped_data, stats_dict):
    """
    Incrementally compute mean and variance using Welford's algorithm.
    
    Parameters:
    - H: height being processed
    - grouped_data: dict (H, x) -> list of velocity magnitudes
    - stats_dict: defaultdict to store running statistics
    """
    for (h, x), values in grouped_data.items():
        for value in values:
            stats = stats_dict[(h, x)]
            
            # Welford's algorithm for mean and variance
            stats['count'] += 1
            delta = value - stats['mean']
            stats['mean'] += delta / stats['count']
            delta2 = value - stats['mean']
            stats['sq_diff_accumulator'] += delta * delta2

def standard_deviation(stats):
    """
    Compute the standard deviation from Welford's statistics.
    
    Parameters:
    - stats: dict with keys 'count', 'mean', 'sq_diff_accumulator'
    
    Returns:
    - Standard deviation or 0 if count < 2
    """

    # Returns the standard deviation if there is more than one data point used
    if stats['count'] > 1:
        return np.sqrt(stats['sq_diff_accumulator'] / (stats['count'] - 1))

    else:
        return np.nan

def plot_data(keys, index, simulation_data, paper_data, axs, color_map, custom_lines):
    """
    Plot the experimental and simulation data for a given index.

    """
    for key in keys[index]:
        x = paper_data[key]["x"] / 1000 # Divided by 1000 to convert x from mm to m
        y = paper_data[key][paper_data[key].columns[1]]
        custom_lines.append(Patch(color = color_map[key], label = f"H/H₀ = 0.{key[-2:]}"))

        # Plot the experimental data
        axs[index].scatter(x, y, label = key, color = color_map[key])

        sim_data = simulation_data.get('simulation_data_' + str(key[-3:]), np.empty((0, 4)))
        if sim_data.size > 0:
            x_sim = sim_data[:, 0]
            y_sim = sim_data[:, 1]
            
            # Plot the simulation data
            axs[index].plot(x_sim, y_sim, color = color_map[key], linestyle = '-')

            # Plot the standard deviation
            std_sim = sim_data[:, 2]
            axs[index].fill_between(x_sim, y_sim - std_sim, y_sim + std_sim, color = color_map[key], alpha = 0.2)

    axs[index].set_xlabel("X (m)")
    axs[index].set_ylabel(r"$V_p$ (m/s)")
    axs[index].legend(handles = custom_lines, loc = "best")