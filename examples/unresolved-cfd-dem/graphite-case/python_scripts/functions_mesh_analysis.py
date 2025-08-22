# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

######################################################################
# Import Libraries

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from cycler import cycler
from collections import defaultdict

######################################################################

def get_positions_and_velocities(flow_mesh):
    """
    Extracts positions and velocities from the flow mesh, avoiding duplicates.
    
    Parameters:
    - flow_mesh: Lethe PyVista object containing the mesh data.
    
    Returns:
    - x_positions: List of unique y positions.
    - avg_velocity: List of average velocities corresponding to y positions.
    """
    x_positions = []
    avg_velocity = []
    seen_set = set()

    # Extract the first DataFrame and average velocity from the flow mesh
    first_df = pd.DataFrame(np.copy(flow_mesh.get_df(0).points), columns=['x', 'y', 'z'])
    avg_velocity_mesh = np.copy(flow_mesh.get_df(0)['average_velocity'])

    for i in range(len(first_df)):
        if abs(first_df['x'][i] - 0) < 1e-10 and abs(first_df['z'][i] - 0.01) < 1e-8:
            if first_df['y'][i] not in seen_set:
                seen_set.add(first_df['y'][i])
                x_positions.append(first_df['y'][i])
                avg_velocity.append(np.linalg.norm(avg_velocity_mesh[i]))

    # Sort by x position
    x_positions, avg_velocity = zip(*sorted(zip(x_positions, avg_velocity)))

    return list(x_positions), list(avg_velocity)
