# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import numpy as np
import pandas as pd
from matplotlib.patches import Patch

######################################################################
# Functions for Post-Processing

def get_y_slices(H, y_grid):
    """
    Get the y-slices for the given height and y-grid.
    
    Parameters:
    - H (float): target height
    - y_grid (numpy.ndarray): array of y-coordinates from the mesh
    
    Returns:
    - numpy.ndarray: indices of the y-slices
    """
    H = round(H, 5)

    # If the target height is exactly in the y-grid, return it directly
    if np.any(np.isclose(y_grid, H)):
        return [H]
        
    # If not, find the two closest values in the y-grid
    else:
        y_below = np.max(y_grid[y_grid < H])
        y_above = np.min(y_grid[y_grid > H])
        return [y_below, y_above]

def vertical_spacing_y(particle):
    """
    Calculate the vertical spacing of the mesh based on the y-coordinates.
    
    Parameters:
    - particle (lethe_pyvista_tools): instance of the lethe_pyvista_tools class
    
    Returns:
    - y_grid (array): mesh coordinates in the y-direction
    - dy (float): vertical spacing between the y-coordinates
    """
    first_df = pd.DataFrame(np.copy(particle.get_df(0).points), columns = ['x', 'y','z'])
    y_grid = np.sort(first_df['y'].unique())
    dy = round(np.min(np.diff(y_grid[y_grid >= 0])), 3) # Only consider positive y-coordinates (in the bed), exclude the channel
    
    return y_grid, dy

def standard_deviation(profile_stats):
    """
    Calculate the standard deviation from the Welford's algorithm statistics.
    
    Parameters:
    - profile_stats (dict): dictionary containing count, mean, and sq_diff_accumulator
    
    Returns:
    - float: standard deviation or NaN if count < 2
    """
    # Returns the standard deviation if there is more than one data point used
    if profile_stats['count'] > 1:
        return np.sqrt(profile_stats['sq_diff_accumulator'] / (profile_stats['count'] - 1))
    else:
        return np.nan

def welford(H, grouped_by_x, stats):
    """
    Apply Welford's algorithm to calculate mean and standard deviation for each x position at height H.
    
    Parameters:
    - H (float): height at which to calculate statistics
    - grouped_by_x (pandas.DataFrame): DataFrame containing x positions and their corresponding velocity magnitudes
    - stats (defaultdict): dictionary to store the statistics
    
    """
    for _, row in grouped_by_x.iterrows():
        x, v = row['x'], row['velocity_mag']
        key = (H, round(x, 5))
        profile_stats = stats[key]
        
        # Welford's algorithm
        profile_stats['count'] += 1
        delta = v - profile_stats['mean']
        profile_stats['mean'] += delta / profile_stats['count']
        delta2 = v - profile_stats['mean']
        profile_stats['sq_diff_accumulator'] += delta * delta2

def group_velocities_by_x(H, y_grid, df_positions, df_velocity_mag):
    """
    Group velocities by x position for a given height H.
    
    Parameters:
    - H (float): height at which to group velocities
    - y_grid (numpy.ndarray): array of y-coordinates from the mesh
    - df_positions (pandas.DataFrame): DataFrame containing particle positions
    - df_velocity_mag (numpy.ndarray): array of velocity magnitudes
    
    Returns:
    pandas.DataFrame: DataFrame with x positions and their corresponding mean velocity magnitudes
    """
    y_slices = get_y_slices(H, y_grid)
    
    mask = np.isclose(df_positions['y'].values, y_slices[0])
    
    if len(y_slices) == 2:
        mask |= np.isclose(df_positions['y'].values, y_slices[1])
    
    # Filter positions and velocity magnitudes based on the mask
    selected = df_positions[mask].copy()
    selected['velocity_mag'] = df_velocity_mag[mask]
    
    return selected.groupby('x')['velocity_mag'].mean().reset_index() # Returns the mean velocity magnitude (in depth of bed z) for each x position at height H

def plot_data(keys, index, simulation_data, paper_data, axs, color_map, custom_lines):
    """
    Plot the experimental and simulation data for a given index.
    
    Parameters:
    - keys (list): list of keys for the experimental datasets (paper data) to plot
    - index (int): index to select the specific set of keys to plot, corresponds to the index of the subgroup of dataset keys in list 'keys'
    - simulation_data (dict): dictionary containing simulation data
    - paper_data (dict): dictionary containing experimental data
    - axs (list): list of matplotlib axes to plot on
    - color_map (dict): dictionary assigning keys to colors for plotting
    - custom_lines (list): list of custom legend lines to add to the plot

    Returns:
    None: updates the axes with the plotted data

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