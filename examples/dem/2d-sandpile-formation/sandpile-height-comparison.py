# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


rollingmethod = ["constant", "viscous", "epsd"]
color = ['blue', 'orange', 'red']


for i, x in enumerate(rollingmethod):

    # Paths to files
    sim_file = 'data/height_' + x + '.csv'
    paper_file = 'data/reference/extraction_model_' + x + '.csv'
    
    # Check if files exist
    if not os.path.exists(sim_file):
        print("Before launching this code, make sure you have \n"    
              "- launched the simulation for rolling resistance models constant, viscous and epsd \n"    
              "- launched sandpile-postprocessing.py after each simulation and using the right command\n" 
              "See Lethe-Documentation for more details\n")
        raise FileNotFoundError(f"Error: The file '{sim_file}' does not exist.")

    # Read the CSV files into Pandas DataFrames
    df_sim = pd.read_csv(sim_file)
    df_paper = pd.read_csv(paper_file)

    # Convert the DataFrames to NumPy arrays
    time_sim = df_sim['time'].to_numpy()
    time_paper = df_paper['x'].to_numpy()
    height_sim = df_sim['height'].to_numpy()
    height_paper = df_paper['Curve1'].to_numpy()

    # Plot the heights of the simulation and of the paper
    plt.plot(time_sim, height_sim, label="Lethe-DEM " + x, color=color[i])
    plt.plot(time_paper, height_paper, '--', label="Ai2010 " + x, color=color[i])


plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Height of the pile (m)')
plt.yticks(np.arange(0.1, 0.32, 0.02))
plt.xticks(np.arange(0, 60, 10))
plt.grid()
plt.title("Evolution of the height of the pile with different rolling resistance models")
plt.savefig("data/height-comparison.png")
plt.show()
