# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


######################################################################
# Import Libraries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# Set plot parameters
plt.rcParams['figure.figsize'] = (10,8)
plt.rcParams['lines.markersize'] = '11'
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.fancybox'] = False
plt.rcParams['font.size'] = 20
plt.rcParams['font.family']='DejaVu Serif'
plt.rcParams['legend.handlelength']=2
plt.rcParams['lines.linewidth'] = 3

######################################################################


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
        raise FileNotFoundError('Error: The file ' + sim_file + ' does not exist.')

    # Read the CSV files into Pandas DataFrames
    df_sim = pd.read_csv(sim_file)
    df_paper = pd.read_csv(paper_file)

    # Plot the heights of the simulation and of the paper
    plt.plot(df_sim['time'], df_sim['height'], label="Lethe-DEM " + x, color=color[i])
    plt.plot(df_paper['x'], df_paper['Curve1'], '--', label="Ai2010 " + x, color=color[i])


plt.subplots_adjust(top=0.8, left=0.2, right=0.9)
plt.legend(fontsize=15 , loc="lower left")
plt.xlabel('Time (s)', labelpad=7)
plt.ylabel('Height of the pile (m)', labelpad=20)
plt.yticks(np.arange(0.1, 0.32, 0.02))
plt.xticks(np.arange(0, 60, 10))
plt.grid()
plt.title("Evolution of the height of the pile \n with different rolling resistance models" , pad=25)
plt.savefig("data/figure-height-comparison.png")
plt.show()
