# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import argparse
from cycler import cycler
from collections import defaultdict
import math

def collision_frequency_by_diameter(csv_file, output_csv = 'collision_frequency.csv'):
    """
    Computes collision frequency from a CSV file and produces a csv file with the collision frequency.
    Args:
        csv_file (str): Path to the CSV file containing collision statistics.
        output_csv (str): Path to save the collision frequency data.
    Returns:
        pd.DataFrame: DataFrame containing collision frequency by particle diameter.
    """

    # Sort diameters from smallest to largest
    particle_diameters = sorted(csv_file['diameter'].unique())

    # Count collisions for each particle by diameter
    collision_counts = {diameter: defaultdict(int) for diameter in particle_diameters}
    collision_frequency = {diameter: defaultdict(int) for diameter in particle_diameters}

    # Count collisions for each particle diameter
    for diameter in particle_diameters:
        diameter_data = csv_file[csv_file['diameter'] == diameter]
        # Count collisions per particle
        for particle_id in diameter_data['particle_id']:
            collision_counts[diameter][particle_id] += 1
        # Count how many particles have each collision count
        for count in collision_counts[diameter].values():
            collision_frequency[diameter][count] += 1
    
    # Create a DataFrame to store the results
    rows = []
    for diameter in particle_diameters:
        for count, num_particles in collision_frequency[diameter].items():
            rows.append({'diameter': diameter, 'collision_count': count, 'num_particles': num_particles})
    df_collision_frequency = pd.DataFrame(rows)

    # Save the DataFrame to a CSV file
    df_collision_frequency.to_csv(output_csv, index=False)

    return df_collision_frequency

def plot_collision_frequency(df_collision_frequency, output_png = 'collision_frequency.png'):
    """
    Plots collision frequency by particle diameter.
    Args:
        df_collision_frequency (pd.DataFrame): DataFrame containing collision frequency data.
        output_png (str): Path to save the plot.
    """
    
    particle_diameters = sorted(df_collision_frequency['diameter'].unique())

    bins_labels = [1, 2, 3, 4, '>4']
    collision_frequency_binned = {diameter: {label:0 for label in bins_labels} for diameter in particle_diameters}
    for _, row in df_collision_frequency.iterrows():
        diameter = row['diameter']
        count = row['collision_count']
        num_particles = row['num_particles']
        if count > 4:
            collision_frequency_binned[diameter]['>4'] += num_particles
        else:
            collision_frequency_binned[diameter][count] += num_particles
    
    # Plot
    x_positions = np.arange(len(bins_labels))
    bar_width = 0.2

    fig, ax = plt.subplots(figsize=(12, 10))
    max_count = 0
    
    for i, diameter in enumerate(particle_diameters):
        y = [collision_frequency_binned[diameter][label] for label in bins_labels]
        max_count = max(max_count, max(y))
        label = diameter / 1e-6  # Convert to micrometers for legend
        ax.bar(x_positions + i * bar_width, y, width=bar_width, label=f'Diameter: {label:.0f} μm')
    
    top_tick = 10 ** math.ceil(math.log10(max_count))
    ax.set_ylim(0.5, top_tick)
    ax.set_yscale('log')
    
    # Set x-ticks
    ax.set_xticks(x_positions + bar_width * (len(particle_diameters)-1)/2)
    ax.set_xticklabels(bins_labels)

    # Labels & title
    ax.set_xlabel('Number of Collisions')
    ax.set_ylabel('Number of Particles')
    ax.set_title('Collision Frequency by Particle Diameter')
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    plt.show()

def mean_speed_difference(csv_file):
    """
    Analyzes mean delta speed from a CSV file.
    Args:
        csv_file (str): Path to the CSV file containing collision statistics.
    Returns:
        mean_delta_velocity (dict): Dictionary with particle diameters as keys and mean speed differences as values.
    """

    particle_diameter = csv_file['diameter'].unique()

    mean_delta_velocity = defaultdict(float)


    # Analyze mean delta velocity for each particle diameter
    for diameter in particle_diameter:
        diameter_data = csv_file[csv_file['diameter'] == diameter]
        start_velocities = np.array([diameter_data['start_particle_velocity_x'], diameter_data['start_particle_velocity_y'], diameter_data['start_particle_velocity_z']]).T
        end_velocities = np.array([diameter_data['end_particle_velocity_x'], diameter_data['end_particle_velocity_y'], diameter_data['end_particle_velocity_z']]).T
        
        # Calculate delta speed
        delta_velocity = np.linalg.norm(end_velocities, axis=1) - np.linalg.norm(start_velocities, axis=1)
        mean_delta_velocity[diameter] = np.mean(delta_velocity)

    return mean_delta_velocity