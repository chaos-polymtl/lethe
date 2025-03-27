# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
This code takes the output from the simulations at the different velocity and pressure shape function degrees and mesh resolutions and plots the corresponding error
This script should be run after running organize_output.py which rearranges the output of the simplex mesh simulations to match the hierarchy and formatting of that obtained for the quad mesh
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import font_manager
from matplotlib import rcParams
import numpy as np
from pathlib import Path
import math
import re
import os

# Define the repository
base_dir = os.getcwd()
repo = base_dir + '/' + 'mms_quad'

# Define the base filename and possible suffixes
path=Path(repo)
pattern = r"degu_(\d+)_degp_(\d+)"
num_columns = 7
col_widths = [6, 12, 5, 11, 5, 8, 7]  # Adjust based on your file structure

m_uL2E = np.empty((0, 0))
m_pL2E = np.empty((0, 0))
m_O_u = np.empty((0, 0))
m_O_p = np.empty((0, 0))
v_O_u_mean = np.empty((0))
v_O_p_mean = np.empty((0))
v_degu = np.empty((0))
v_degp = np.empty((0))

for folders_MMS in path.glob("mms*"):
  if folders_MMS.is_dir():
    print(folders_MMS.name)
    match = re.search(pattern, folders_MMS.name)
    degu= match.group(1)  # The degree of the velocity shape functions
    degp = match.group(2)  # The degree of the pressure shape functions
    file = "L2Error.dat"
    print(folders_MMS, file, f"degu:{degu}", f"degp:{degp}")

    file = f"{folders_MMS}/L2Error.dat"
    print ("Filename-> %s" %file)

    df = pd.read_fwf(file, widths=col_widths, skiprows=0)
    n = df.iloc[:, 0].to_numpy()
    uL2E = df.iloc[:, 1].to_numpy()
    O_u = df.iloc[:, 2].to_numpy()
    O_u_mean = np.mean(np.asarray(O_u[1:-1], dtype=np.float64))
    pL2E = df.iloc[:, 3].to_numpy()
    O_p = df.iloc[:, 4].to_numpy()
    O_p_mean = np.mean(np.asarray(O_p[1:-1], dtype=np.float64))
    # print(n, uL2E, O_u, pL2E, O_p)
    # The cell volume is Vtotal/n where Vtot=4. h is defined here as the diameter of a circle with the same area as the cell
    # h = (4*4/(np.pi*n))**0.5
    # The cell volume is taken here to be equal to the sqaure side
    h = 2/(n**0.5)

    if m_uL2E.size == 0:  # First iteration, initialize the matrix
        m_uL2E = uL2E.reshape(-1, 1)  # Convert vector to column vector
        m_pL2E = pL2E.reshape(-1, 1) 
        m_O_u = O_u.reshape(-1, 1)  
        m_O_p = O_p.reshape(-1, 1)  
        v_O_u_mean = O_u_mean
        v_O_p_mean = O_p_mean
        v_degu = degu
        v_degp = degp
    else:
        m_uL2E= np.column_stack(( m_uL2E, uL2E))
        m_pL2E= np.column_stack(( m_pL2E, pL2E))
        m_O_u= np.column_stack(( m_O_u, O_u))
        m_O_p= np.column_stack(( m_O_p, O_p))
        v_O_u_mean = np.append(v_O_u_mean, O_u_mean)
        v_O_p_mean = np.append(v_O_p_mean, O_p_mean)
        v_degu = np.append(v_degu, degu)
        v_degp = np.append(v_degp, degp)

# Orders to display on the plot
v_O_u_mean = np.round(v_O_u_mean)
v_O_p_mean = np.round(v_O_p_mean)
# Modify font of the graphic
font = {'weight' : 'normal',
        'size'   : 18}
plt.rc('font', **font)
plt.rcParams['legend.numpoints'] = 1
params = {'backend': 'ps',
             'axes.labelsize': 24,
             'font.size': 28,
             'legend.fontsize': 17,
             'xtick.labelsize': 15,
             'ytick.labelsize': 15,
             'text.usetex': False,
             }
plt.rcParams.update(params)

# Plot velocity L2 error
plt.figure(num=1)
# Plot the different lines
lines = plt.plot(h, m_uL2E, label='m_uL2E Line')

# Shaded triangular area to show the slope
highlight_x = np.log10(np.array([5e-3, 6.5e-3]))  # X-limits of the triangle

# Define a function to apply np.polyfit to each column
def fit_log(log_h, log_uL2E):
    return np.polyfit(log_h, log_uL2E, 1)

# We use the logarithm of the x and y coordinates to avoid having to fit polynomials of varying degrees
log_h = np.log10(h)  # log(h)
log_m_uL2E = np.log10(m_uL2E)  # log(uL2E)

# Apply polyfit to each column of log_uL2E
coeffs = np.apply_along_axis(fit_log, 0, log_h, log_m_uL2E)
a = coeffs[0, :]  # Slope of each column
b = coeffs[1, :]  # Intercept of each column
a = a.reshape(1, -1) # This is a row vector
b = b.reshape(1, -1) # This is a row vector
highlight_x = highlight_x.reshape(-1, 1) # This is a column vector
# We do the previous reshaping to be able to do the following broadcasting
# Each column in highlight_y corresponds to a column in highlight_x, and therefore to a line in the plot
highlight_y = highlight_x*a + b  # Corresponding y-values

# Add third vertex of each triangle
highlight_x =  np.tile(highlight_x, (1, highlight_y.shape[1])) 

# Linearize values again
highlight_y = np.power(10,highlight_y)
highlight_x = np.power(10,highlight_x)

row_to_add = highlight_x[1,:].reshape(1,-1)
highlight_x = np.vstack((highlight_x, row_to_add))
row_to_add = highlight_y[0,:].reshape(1,-1)
highlight_y = np.vstack((highlight_y, row_to_add))

legend_labels = []

for i in range(highlight_x.shape[1]):  # Iterate over the columns (over the plot lines)
    # Fill legend labels for each line
    legend_labels.append(f"Q{v_degu[i]} - Q{v_degp[i]}")

    if int(v_degp[i]) != 3:
      x_coords = highlight_x[:, i]  # Get the x-coordinates for the triangle
      y_coords = highlight_y[:, i]  # Get the y-coordinates for the triangle
    
      # Plot the shaded triangle
      plt.fill(x_coords, y_coords, 'gray', alpha=0.3)  # 'b' is blue, alpha controls transparency
      # Write the order of convergence next to the filled triangle
      alpha = 1.075
      x_center = highlight_x[1,i]*alpha  # You can calculate the center of the triangle to place the text
      y_center = 5*(highlight_y[1,i] + highlight_y[2,i])/12
      # y_center = highlight_y[2,i] 
      plt.text(x_center, y_center, f"{v_O_u_mean[i]}", fontsize=6, ha='center', va='center')
   

# Define 6 distinct colors for the lines
colors = ["#d73027","#fdae61", "#fee08b", "#d9ef8b", "#1591EA", "#1a9850"]
markers = ['o', 's', '^', 'D', 'x', '*'] # Define markers for each line

# We do this because the files are not read in increasing order of degu and then degp
# Combine degu and degp into a list of tuples (degu, degp)
my_list = list(zip(v_degu, v_degp))

# Sort the list first by degu and then by degp (for each degu)
sorted_list = sorted(my_list, key=lambda x: (x[0], x[1]))

# Separate the sorted list into degu and degp again
sorted_degu, sorted_degps = zip(*sorted_list)

# Get the sorted indices based on the sorting of degu and degp
sorted_indices = [my_list.index(t) for t in sorted_list]

# Now, reorder another_list based on the same sorted indices
legend_labels_sorted = [legend_labels[i] for i in sorted_indices]

handles, legend_labels_1 = plt.gca().get_legend_handles_labels()

quad_handles = [handles[i] for i in sorted_indices]
quad_labels = legend_labels_sorted

lines = [lines[i] for i in sorted_indices]

for i, (line, color, marker, label) in enumerate(zip(lines, colors, markers, legend_labels_sorted)):
    line.set_color(color)
    line.set_marker(marker)
    line.set_label(label)

# Apply the function along each column of uL2E
plt.xscale('log') 
plt.yscale('log') 

# This has to be after setting the log-scale
tick_values = [2*10**(-3), 10**(-2), 7*10**(-2)]  # Tick positions
tick_labels = [r"$2\times 10^{-3}$", r"$10^{-2}$", r"$7\times 10^{-2}$"]  # LaTeX labels
plt.xticks(tick_values, tick_labels)  # Set x-ticks and labels
plt.xlabel('h')       # label for x-axis
plt.ylabel('$\Vert e\Vert_{2}$')  # label for y-axis

####################################################################################################################################################
# Plot pressure L2 error
plt.figure(num=2)
lines = plt.plot(h, m_pL2E)
# Shaded triangular area to show the slope
highlight_x = np.log10(np.array([5e-3, 6.5e-3]))  # X-coordinates of the triangle
# We use the logarithm of the x and y coordinates to avoid having to fit polynomials of varying degrees
log_h = np.log10(h)  # log(h)
log_m_pL2E = np.log10(m_pL2E)  # log(uL2E)

# Apply polyfit to each column of log_uL2E
coeffs = np.apply_along_axis(fit_log, 0, log_h, log_m_pL2E)
a = coeffs[0, :]  # Slope of each column
b = coeffs[1, :]  # Intercept of each column
a = a.reshape(1, -1) # This is a row vector
b = b.reshape(1, -1) # This is a row vector
highlight_x = highlight_x.reshape(-1, 1) # This is a column vector
# Each column in highlight_y corresponds to a column in highlight_x, and therefore to a line in the plot
highlight_y = highlight_x*a+ b  # Corresponding y-values

# Add third vertex of each triangle
highlight_x =  np.tile(highlight_x, (1, highlight_y.shape[1])) 

# Linearize values again
highlight_y = np.power(10,highlight_y)
highlight_x = np.power(10,highlight_x)

row_to_add = highlight_x[1,:].reshape(1,-1)
highlight_x = np.vstack((highlight_x, row_to_add))
row_to_add = highlight_y[0,:].reshape(1,-1)
highlight_y = np.vstack((highlight_y, row_to_add))

for i in range(highlight_x.shape[1]):  # Iterate over the columns
    x_coords = highlight_x[:, i]  # Get the x-coordinates for the triangle
    y_coords = highlight_y[:, i]  # Get the y-coordinates for the triangle
    
    # Plot the shaded triangle
    plt.fill(x_coords, y_coords, 'gray', alpha=0.3)  # 'b' is blue, alpha controls transparency
    # Write the order of convergence next to the filled triangle
    alpha = 1.075
    x_center = highlight_x[1,i]*alpha  # You can calculate the center of the triangle to place the text
    y_center = 5*(highlight_y[1,i] + highlight_y[2,i])/12
    # y_center = highlight_y[2,i] 
    plt.text(x_center, y_center, f"{v_O_p_mean[i]}", fontsize=6, ha='center', va='center')

plt.xscale('log') 
plt.yscale('log') 
tick_values = [2*10**(-3), 10**(-2), 7*10**(-2)]  # Tick positions
tick_labels = [r"$2\times 10^{-3}$", r"$10^{-2}$", r"$7\times 10^{-2}$"]  # LaTeX labels
plt.xticks(tick_values, tick_labels)  # Set x-ticks and labels
plt.xlabel('h')       # label for x-axis
plt.ylabel('$\Vert e\Vert_{2}$')  # label for y-axis

lines = [lines[i] for i in sorted_indices]
for line, color, marker, label in zip(lines, colors, markers, legend_labels_sorted):
    line.set_color(color)
    line.set_marker(marker)
    line.set_label(label)

#########################################################################################################################################################################################################
#                                                                                                        SIMPLEX mesh                                                                                   #
#########################################################################################################################################################################################################
# Define the repository
repo = base_dir + '/' + 'mms_simplex'
# Define the base filename 
path=Path(repo)

lines=[]
m_uL2E = np.empty((0, 0))
uL2E = np.empty((0))
m_pL2E = np.empty((0, 0))
pL2E = np.empty((0))
m_O_u = np.empty((0, 0))
O_u = np.empty((0))
m_O_p = np.empty((0, 0))
O_p = np.empty((0))
v_O_u_mean = np.empty((0))
O_u_mean = 0
v_O_p_mean = np.empty((0))
O_p_mean = 0
v_degu = np.empty((0))
degu = 0
v_degp = np.empty((0))
degp = 0


for folders_MMS in path.glob("mms*"):
  if folders_MMS.is_dir():
    print(folders_MMS.name)
    match = re.search(pattern, folders_MMS.name)
    degu= match.group(1)  # The degree of the velocity shape functions
    degp = match.group(2)  # The degree of the pressure shape functions
    file = "L2Error.dat"
    print(folders_MMS, file, f"degu:{degu}", f"degp:{degp}")

    file = f"{folders_MMS}/L2Error.dat"
    print ("Filename-> %s" %file)

    #mat = np.genfromtxt(file, skip_header=1, dtype=float, delimiter='none', invalid_raise= = False)
    df = pd.read_fwf(file, widths=col_widths, skiprows=0)
    n = df.iloc[:, 0].to_numpy()
    uL2E = df.iloc[:, 1].to_numpy()
    O_u = df.iloc[:, 2].to_numpy()
    O_u_mean = np.mean(np.asarray(O_u[1:-1], dtype=np.float64))
    pL2E = df.iloc[:, 3].to_numpy()
    O_p = df.iloc[:, 4].to_numpy()
    O_p_mean = np.mean(np.asarray(O_p[1:-1], dtype=np.float64))
    # print(n, uL2E, O_u, pL2E, O_p)
    # The cell volume is Vtotal/n where Vtot=4. h is defined here as the diameter of a circle with the same area as the cell
    # h = (4/(np.pi*n))**0.5
    # The side of one triangle is equal to half of the side of the square it belongs to: https://www.dealii.org/current/doxygen/deal.II/namespaceGridGenerator.html#ac7515d2b17c025dddc0e37286fb8d216
    h = 2/((n/8)**0.5)*0.5

    if m_uL2E.size == 0:  # First iteration, initialize the matrix
        m_uL2E = uL2E.reshape(-1, 1)  # Convert vector to column vector
        m_pL2E = pL2E.reshape(-1, 1) 
        m_O_u = O_u.reshape(-1, 1)  
        m_O_p = O_p.reshape(-1, 1)  
        v_O_u_mean = O_u_mean
        v_O_p_mean = O_p_mean
        v_degu = degu
        v_degp = degp
    else:
        m_uL2E= np.column_stack(( m_uL2E, uL2E))
        m_pL2E= np.column_stack(( m_pL2E, pL2E))
        m_O_u= np.column_stack(( m_O_u, O_u))
        m_O_p= np.column_stack(( m_O_p, O_p))
        v_O_u_mean = np.append(v_O_u_mean, O_u_mean)
        v_O_p_mean = np.append(v_O_p_mean, O_p_mean)
        v_degu = np.append(v_degu, degu)
        v_degp = np.append(v_degp, degp)

# Orders to display on the plot
v_O_u_mean = np.round(v_O_u_mean)
v_O_p_mean = np.round(v_O_p_mean)


# Plot velocity L2 error
plt.figure(num=1)

lines = plt.plot(h, m_uL2E)

# Shaded triangular area to show the slope
highlight_x = np.log10(np.array([3.5e-2, 4.5e-2]))  # X-coordinates of the triangle

# Define a function to apply np.polyfit to each column
def fit_log(log_h, log_uL2E):
    return np.polyfit(log_h, log_uL2E, 1)

# We use the logarithm of the x and y coordinates to avoid having to fit polynomials of varying degrees
log_h = np.log10(h)  # log(h)
log_m_uL2E = np.log10(m_uL2E)  # log(uL2E)

# Apply polyfit to each column of log_uL2E
coeffs = np.apply_along_axis(fit_log, 0, log_h, log_m_uL2E)
a = coeffs[0, :]  # Slope of each column
b = coeffs[1, :]  # Intercept of each column
a = a.reshape(1, -1) # This is a row vector
b = b.reshape(1, -1) # This is a row vector
highlight_x = highlight_x.reshape(-1, 1) # This is a column vector
# We do the previous reshaping to be able to do the following broadcasting
# Each column in highlight_y corresponds to a column in highlight_x, and therefore to a line in the plot
highlight_y = highlight_x*a+ b  # Corresponding y-values

# Add third vertex of each triangle
highlight_x =  np.tile(highlight_x, (1, highlight_y.shape[1])) 

# Linearize values again
highlight_y = np.power(10,highlight_y)
highlight_x = np.power(10,highlight_x)

row_to_add = highlight_x[1,:].reshape(1,-1)
highlight_x = np.vstack((highlight_x, row_to_add))
row_to_add = highlight_y[0,:].reshape(1,-1)
highlight_y = np.vstack((highlight_y, row_to_add))

legend_labels = []

for i in range(highlight_x.shape[1]):  # Iterate over the columns
    # Fill legend labels for each line
    legend_labels.append(f"P{v_degu[i]} - P{v_degp[i]}")
    if (int(v_degu[i]) != 3):
      x_coords = highlight_x[:, i]  # Get the x-coordinates for the triangle
      y_coords = highlight_y[:, i]  # Get the y-coordinates for the triangle
    
      # Plot the shaded triangle
      plt.fill(x_coords, y_coords, 'gray', alpha=0.3)  # 'b' is blue, alpha controls transparency
      # Write the order of convergence next to the filled triangle
      alpha = 1.075
      x_center = highlight_x[1,i]*alpha  # You can calculate the center of the triangle to place the text
      y_center = 5*(highlight_y[1,i] + highlight_y[2,i])/12
      plt.text(x_center, y_center, f"{v_O_u_mean[i]}", fontsize=6, ha='center', va='center')



# We do this because the files are not read in increasing order of degu and then degp
# Combine degu and degp into a list of tuples (degu, degp)
my_list = list(zip(v_degu, v_degp))

# Sort the list first by degu and then by degp (for each degu)
sorted_list = sorted(my_list, key=lambda x: (x[0], x[1]))

# Separate the sorted list into degu and degp again
sorted_degu, sorted_degps = zip(*sorted_list)

# Get the sorted indices based on the sorting of degu and degp
sorted_indices = [my_list.index(t) for t in sorted_list]


legend_labels_sorted = [legend_labels[i] for i in sorted_indices]
simplex_labels = legend_labels_sorted

lines = [lines[i] for i in sorted_indices]
simplex_handles = []

for line, color, marker, label in zip(lines, colors, markers, legend_labels_sorted):
    line.set_color(color)
    line.set_marker(marker)
    line.set_label(label)
    line.set_linestyle('--') 
    simplex_handles.append(line)

# Apply the function along each column of uL2E
plt.xscale('log') 
plt.yscale('log') 

# This has to be after setting the log-scale
tick_values = [2*10**(-3), 10**(-2), 7*10**(-2)]  # Tick positions
tick_labels = [r"$2\times 10^{-3}$", r"$10^{-2}$", r"$7\times 10^{-2}$"]  # LaTeX labels
plt.xticks(tick_values, tick_labels)  # Set x-ticks and labels
plt.xlabel('h')       # label for x-axis
plt.ylabel('$\Vert e\Vert_{2}$')  # label for y-axis

legend_handles = quad_handles + simplex_handles
legend_labels = quad_labels + simplex_labels

# Define the font properties explicitly
font_properties = font_manager.FontProperties(size=9.5)
plt.legend(legend_handles, legend_labels, prop=font_properties , loc='lower right', ncol=2)
plt.savefig("order_of_convergence_velocity.png", dpi=300, bbox_inches='tight')
# plt.show()
plt.clf() 

################################################################################################################################################################
# Plot pressure L2 error
lines=[]
plt.figure(num=2)
lines = plt.plot(h, m_pL2E)
# Shaded triangular area to show the slope
highlight_x = np.log10(np.array([3.5e-2, 4.5e-2]))  # X-coordinates of the triangle
# We use the logarithm of the x and y coordinates to avoid having to fit polynomials of varying degrees
log_h = np.log10(h)  # log(h)
log_m_pL2E = np.log10(m_pL2E)  # log(uL2E)

# Apply polyfit to each column of log_uL2E
coeffs = np.apply_along_axis(fit_log, 0, log_h, log_m_pL2E)
a = coeffs[0, :]  # Slope of each column
b = coeffs[1, :]  # Intercept of each column
a = a.reshape(1, -1) # This is a row vector
b = b.reshape(1, -1) # This is a row vector
highlight_x = highlight_x.reshape(-1, 1) # This is a column vector
# Each column in highlight_y corresponds to a column in highlight_x, and therefore to a line in the plot
highlight_y = highlight_x*a+ b  # Corresponding y-values

# Add third vertex of each triangle
highlight_x =  np.tile(highlight_x, (1, highlight_y.shape[1])) 

# Linearize values again
highlight_y = np.power(10,highlight_y)
highlight_x = np.power(10,highlight_x)

row_to_add = highlight_x[1,:].reshape(1,-1)
highlight_x = np.vstack((highlight_x, row_to_add))
row_to_add = highlight_y[0,:].reshape(1,-1)
highlight_y = np.vstack((highlight_y, row_to_add))

for i in range(highlight_x.shape[1]):  # Iterate over the columns
    if (int(v_degu[i]) != 3):
      x_coords = highlight_x[:, i]  # Get the x-coordinates for the triangle
      y_coords = highlight_y[:, i]  # Get the y-coordinates for the triangle
    
      # Plot the shaded triangle
      plt.fill(x_coords, y_coords, 'gray', alpha=0.3)  # 'b' is blue, alpha controls transparency
      # Write the order of convergence next to the filled triangle
      alpha = 1.075
      x_center = highlight_x[1,i]*alpha  # You can calculate the center of the triangle to place the text
      y_center = 5*(highlight_y[1,i] + highlight_y[2,i])/12
      # y_center = highlight_y[2,i] 
      plt.text(x_center, y_center, f"{v_O_p_mean[i]}", fontsize=6, ha='center', va='center')

plt.xscale('log') 
plt.yscale('log') 
tick_values = [2*10**(-3), 10**(-2), 7*10**(-2)]  # Tick positions
tick_labels = [r"$2\times 10^{-3}$", r"$10^{-2}$", r"$7\times 10^{-2}$"]  # LaTeX labels
plt.xticks(tick_values, tick_labels)  # Set x-ticks and labels
plt.xlabel('h')       # label for x-axis
plt.ylabel('$\Vert e\Vert_{2}$')  # label for y-axis

lines = [lines[i] for i in sorted_indices]

for line, color, marker, label in zip(lines, colors, markers, legend_labels_sorted):
    line.set_color(color)
    line.set_marker(marker)
    line.set_label(label)
    line.set_linestyle('--') 

# Define the font properties explicitly
font_properties = font_manager.FontProperties(size=9.5) 
plt.legend(legend_handles, legend_labels, prop=font_properties , loc='lower right', ncol=2)
# Reusing the same legend handles and labels   
plt.savefig("order_of_convergence_pressure.png", dpi=300, bbox_inches='tight')
plt.clf() 

