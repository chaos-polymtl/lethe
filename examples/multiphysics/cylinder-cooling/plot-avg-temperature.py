"""
IMPORTS
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from matplotlib.lines import Line2D
import numpy as np
from cycler import cycler

"""
FILES
"""

csv_files_list = ['postprocessed-avg-temperature-cylinder-xi-0_01.csv','postprocessed-avg-temperature-cylinder-xi-1.csv', 'postprocessed-avg-temperature-cylinder-xi-100.csv']

"""
PLOT
"""

# Plot style
colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
plt.rcParams['axes.prop_cycle'] = cycler(color = colors)
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize'] = (10,9)
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['lines.markersize'] = '11'
plt.rcParams['markers.fillstyle'] = "none"
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.columnspacing'] = 2
plt.rcParams['legend.handlelength'] = 3
plt.rcParams['legend.handletextpad'] = 0.2
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fancybox'] = False
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['font.size'] = '25'
plt.rcParams['font.family']='DejaVu Serif'
plt.rcParams['font.serif']='cm'
plt.rcParams['savefig.bbox']='tight'
plt.rcParams['legend.handlelength']=1

# Separate the figure into subplots
fig,ax = plt.subplots(1,1)

# Loop through the CSV files and plot the avg temperature
for i, csv_file in enumerate(csv_files_list):
    
    # Load the CSV file
    data = pd.read_csv(csv_file, delimiter='\t')

    # Check if there are enough columns
    if data.shape[0] < 2:
        raise ValueError("The CSV file must contain at least two columns.")

    # Plot
    time = data.iloc[:, 1]
    avg_temperature = data.iloc[:, 2]
    ax.plot(time, avg_temperature)

# Reference data
x_ref_01_stars = np.array([0.000010079075541254116, 0.00009947627933938772, 0.0003507841123301386, 
              0.0010026289340245849, 0.003498643420573148, 0.010105572765880899, 
              0.10079067465433475])
y_ref_01_stars = np.array([0.989573438491937, 0.8938388510862751, 0.7327014176110445, 
              0.5507108972779542, 0.3127961943930949, 0.13554504649695417, 
              0.018009492322302888])

x_ref_1_stars = np.array([0.00010158776376630683, 0.0010239107450676843, 0.009999991987538346, 
              0.03489469834597206, 0.10079067465433475, 0.3480317513545161, 
              1.0052647793432752, 10.02628130671995])

y_ref_1_stars = np.array([0.9924170601461935, 0.9668246370096506, 0.8834123257364785, 
              0.7895734529552435, 0.6502369591303025, 0.3819905329884835, 
              0.10521332323999436, 0])

x_ref_100 = np.array([0.0010347204540048808, 0.009999991987538346, 
              0.09973779591478213, 1.0052647793432752, 3.471194760127587, 
              10.02628130671995, 34.986434205731406, 98.95522222452317])

y_ref_100 = np.array([1.0009478618319967, 0.9981042582568734, 
              0.9933649219781902, 0.9658767661380873, 0.900473929108085, 
              0.7573459879551572, 0.40853084507572274, 0.0881516837101216])

# Plot the reference data
ax.plot(x_ref_01_stars, y_ref_01_stars, 'o', color=colors[0], label=r'Juncu, $\beta$ = 0.01')
ax.plot(x_ref_1_stars, y_ref_1_stars, 'o', color=colors[1], label=r'Juncu, $\beta$ = 1')
ax.plot(x_ref_100, y_ref_100, 'o', color=colors[2], label=r'Juncu, $\beta$ = 100')

# Add labels and legend
legend_elements = [Line2D([0], [0], marker='o', color='k', markersize=0, lw=3, markeredgecolor=colors[0], markeredgewidth=2), Line2D([0], [0], marker='o', color=colors[0], markersize=12, lw=0, markeredgecolor='k', markeredgewidth=2), Line2D([0], [0], marker='o', color=colors[0], markersize=12, lw=3, markeredgecolor=colors[0], markeredgewidth=2),
                          Line2D([0], [0], marker='o', color=colors[1], markersize=12, lw=3, markeredgecolor=colors[1], markeredgewidth=2),
                          Line2D([0], [0], marker='o', color=colors[2], markersize=12, lw=3, markeredgecolor=colors[2], markeredgewidth=2)]
labels = ['Lethe', 'Juncu', r'$\beta$ = 0.01', r'$\beta$ = 1', r'$\beta$ = 100']
fig.legend(legend_elements, labels, loc='upper center', ncol = 6,  facecolor = 'white', framealpha = 0.75,  edgecolor = 'black', fancybox = False, shadow = False, fontsize=20, bbox_to_anchor=(0.5, 0.98))

# Set the axis labels and limits
ax.set_xlabel('Dimensionless time [-]')
ax.set_ylabel('Dimensionless solid temperature [-]')
ax.set_xscale("log")
ax.set_xlim(0.000001, 1000)
ax.set_ylim(-0.01, 1)

plt.show()