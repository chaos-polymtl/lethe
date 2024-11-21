"""
IMPORTS
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

"""
MAIN
"""

# Read data
found_positions = pd.read_csv("found_positions.csv", header=1).values 
real_positions = np.loadtxt("real-positions.txt", skiprows=1)

# Extract X and Y coordinates
found_x, found_y = found_positions[:, 0], found_positions[:, 1]
real_x, real_y = real_positions[:, 0], real_positions[:, 1]

"""
PLOTS
"""

plt.figure(figsize=(6, 6))
plt.scatter(real_x, real_y, c='red', label='Experimental position', s=10, alpha=0.8)
plt.scatter(found_x, found_y, c='black', label='Reconstructed position', s=10, alpha=0.8)

plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.legend(loc="upper center", bbox_to_anchor=(0.5, 1.1), ncol=2)
plt.axis('equal') # ensures the grid is square

plt.show()
