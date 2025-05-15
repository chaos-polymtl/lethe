# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
This code takes the output from the simulations at the different velocity and pressure shape function degrees and mesh resolutions and plots the corresponding error
This script should be run after running organize_output.py which rearranges the output of the simplex mesh simulations to match the hierarchy and formatting of that obtained for the quad mesh
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
from matplotlib import font_manager
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
from cycler import cycler
import numpy as np
from pathlib import Path
import math
import re
import os

# Modify font and color of the plots
colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
markers = ['o', 's', '^', 'D', 'x', '*'] # Define markers for each line

font = {'weight' : 'normal',
        'size'   : 18}
plt.rc('font', **font)
plt.rcParams['legend.numpoints'] = 1
params = {#'backend': 'ps',
             'axes.labelsize': 24,
             'font.size': 28,
             'legend.fontsize': 17,
             'xtick.labelsize': 15,
             'ytick.labelsize': 15,
             'text.usetex': False,
             }

# Set the font for the order of convergence indication
order_font=8

# Additional plot parameters to manually control plot formatting
# Tick values and labels for the x axis (mesh size)
x_tick_values = [10**(-2), 1*10**(-1)]  # Tick positions
x_tick_labels = [r"$10^{-2}$", r"$1\times 10^{-1}$"]  # LaTeX labels

#Tick values and label for the pressure error
p_tick_values = [10**(-8), 10**(-6), 10**(-4), 10**(-2), 10**(0)]  # Tick positions
p_tick_labels = [r"$10^{-8}$", r"$10^{-6}$", r"$10^{-4}$", r"$10^{-2}$", r"$10^{0}$"]  # LaTeX labels
p_range = [10**(-8), 10**(0)]  # Range for the y-axis ressure error

#Tick values and label for the velocity error
u_tick_values = [10**(-10), 10**(-8), 10**(-6), 10**(-4), 10**(-2), 10**(0)]  # Tick positions
u_tick_labels = [r"$10^{-10}$",r"$10^{-8}$", r"$10^{-6}$", r"$10^{-4}$", r"$10^{-2}$", r"$10^{0}$"]  # LaTeX labels
u_range = [10**(-10), 10**(0)]  # Range for the y-axis for velocity error


# Number of decimals to use when rounding the order of convergence
n_round=1

def load_a_file(file):
    """
    Load a file and return its contents. The content is the L2 error of the velocity and pressure as well as the mean order of convergence 
    The contents are the
        
        Inputs:
        - file : A string representing the name of the file to be loaded. This string, in general, points to an L2Error.dat file
    """
    print ("Filename-> %s" %file)

    n,uL2E,pL2E = np.loadtxt(file, skiprows=1, unpack=True,usecols=(0,1,3), dtype=float)
    O_u,O_p = np.loadtxt(file, skiprows=2, unpack=True,usecols=(2,4), dtype=float)

    O_u_mean = np.mean(O_u)
    O_p_mean = np.mean(O_p)
    # Round the order of convergence to n_round decimals
    O_u_mean = round(O_u_mean, n_round)
    O_p_mean = round(O_p_mean, n_round)

    return (n, uL2E, pL2E, O_u, O_p, O_u_mean, O_p_mean)

def process_folder(folder):
    """
    Process a folder to gather all of the data from the different simulations at various orders of accuracy.
    returns the number of cells (n)

        Inputs:
        - folder : A string representing the name of the folder containing all simulation results.
          This string should either point to mms_quad or mms_simplex
           
        Outputs:
        - n : The number of cells in the mesh
        - m_uL2E : The L2 error of the velocity
        - m_pL2E : The L2 error of the pressure
        - m_O_u : The order of convergence of the velocity
        - m_O_p : The order of convergence of the pressure
        - v_O_u_mean : The mean order of convergence of the velocity
        - v_O_p_mean : The mean order of convergence of the pressure
        - v_degu : The degree of the velocity shape functions
        - v_degp : The degree of the pressure shape functions
    """

    base_dir = os.getcwd()
    repo = base_dir + '/' + folder + '/'

    # Define the base filename and possible suffixes
    path=Path(repo)
    pattern = r"degu_(\d+)_degp_(\d+)"

    # For the first file we need to allocate storage array to which we will append
    # when reading subsequent files
    first_file=True
    for folders_MMS in path.glob("mms*"):
        if folders_MMS.is_dir():
            match = re.search(pattern, folders_MMS.name)
            degu = match.group(1)  # The degree of the velocity shape functions
            degp = match.group(2)  # The degree of the pressure shape functions
            file = "L2Error.dat"
            print(folders_MMS.name, file, f"degu:{degu}", f"degp:{degp}")
            print(os.getcwd())

            file = f"{repo+folders_MMS.name}/L2Error.dat"

            n,uL2E,pL2E,O_u,O_p,O_u_mean,O_p_mean = load_a_file(file)

            if first_file:  # First iteration, initialize the matrix
                first_file = False
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

    return n, m_uL2E, m_pL2E, m_O_u, m_O_p, v_O_u_mean, v_O_p_mean, v_degu, v_degp

# Define a function to apply np.polyfit to each column
def fit_log(log_h, log_err):
    return np.polyfit(log_h, log_err, 1)

### Load the data for both the quad (_q) and simplex (_s) meshes

n_q, m_uL2E_q, m_pL2E_q, m_O_u_q, m_O_p_q, v_O_u_mean_q, v_O_p_mean_q, v_degu_q, v_degp_q = process_folder("mms_quad")
# The cell volume is taken here to be equal to the square side
h_q = 2/(n_q**0.5)

n_s, m_uL2E_s, m_pL2E_s, m_O_u_s, m_O_p_s, v_O_u_mean_s, v_O_p_mean_s, v_degu_s, v_degp_s = process_folder("mms_simplex")
# The side of one triangle is equal to half of the side of the square it belongs to: https://www.dealii.org/current/doxygen/deal.II/namespaceGridGenerator.html#ac7515d2b17c025dddc0e37286fb8d216
h_s = 2/((n_s/8)**0.5)*0.5


### Plot for the L2 error on the velocity

### Quad mesh
# We first begin my plotting the quad mesh results
# Plot the different lines with dummy labels that we will modify a posteriori
lines = plt.plot(h_q, m_uL2E_q, label='Q Line - Velocity')

# Shaded triangular area to show the slope
highlight_x = np.log10(np.array([1e-2, 1.3e-2]))  # X-limits of the triangle

# We use the logarithm of the x and y coordinates to obtain the regression through the error
log_h = np.log10(h_q)  # log(h)
log_m_uL2E = np.log10(m_uL2E_q)  # log(uL2E)

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
    legend_labels.append(f"Q{v_degu_q[i]} - Q{v_degp_q[i]}")

    if int(v_degp_q[i]) != 3:
      x_coords = highlight_x[:, i]  # Get the x-coordinates for the triangle
      y_coords = highlight_y[:, i]  # Get the y-coordinates for the triangle
    
      # Plot the shaded triangle
      plt.fill(x_coords, y_coords, 'gray', alpha=0.3)  # 'b' is blue, alpha controls transparency
      # Write the order of convergence next to the filled triangle
      alpha = 1.075
      x_center = highlight_x[1,i]*alpha  # You can calculate the center of the triangle to place the text
      y_center = 5*(highlight_y[1,i] + highlight_y[2,i])/12
      # y_center = highlight_y[2,i] 
      plt.text(x_center, y_center, f"{v_O_u_mean_q[i]}", fontsize=order_font, ha='center', va='center')
   
# We do this because the files are not read in increasing order of degu and then degp
# Combine degu and degp into a list of tuples (degu, degp)
my_list_q = list(zip(v_degu_q, v_degp_q))

# Sort the list first by degu and then by degp (for each degu)
sorted_list_q = sorted(my_list_q, key=lambda x: (x[0], x[1]))

# Separate the sorted list into degu and degp again
sorted_degu_q, sorted_degps_q = zip(*sorted_list_q)

# Get the sorted indices based on the sorting of degu and degp
sorted_indices_q = [my_list_q.index(t) for t in sorted_list_q]

# Now, reorder another_list based on the same sorted indices
legend_labels_sorted_q = [legend_labels[i] for i in sorted_indices_q]

handles, legend_labels_1 = plt.gca().get_legend_handles_labels()

quad_handles = [handles[i] for i in sorted_indices_q]
quad_labels = legend_labels_sorted_q

lines = [lines[i] for i in sorted_indices_q]
for i, (line, color, marker, label) in enumerate(zip(lines, colors, markers, legend_labels_sorted_q)):
    line.set_color(color)
    line.set_marker(marker)
    line.set_label(label)

### Simplex mesh
lines = plt.plot(h_s, m_uL2E_s, label='S Line - Velocity')

# Shaded triangular area to show the slope
highlight_x = np.log10(np.array([3.5e-2, 4.5e-2]))  # X-coordinates of the triangle

# We use the logarithm of the x and y coordinates to avoid having to fit polynomials of varying degrees
log_h = np.log10(h_s)  # log(h)
log_m_uL2E_s = np.log10(m_uL2E_s)  # log(uL2E)

# Apply polyfit to each column of log_uL2E
coeffs = np.apply_along_axis(fit_log, 0, log_h, log_m_uL2E_s)
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
    legend_labels.append(f"P{v_degu_s[i]} - P{v_degp_s[i]}")
    if (int(v_degu_s[i]) != 3):
      x_coords = highlight_x[:, i]  # Get the x-coordinates for the triangle
      y_coords = highlight_y[:, i]  # Get the y-coordinates for the triangle
    
      # Plot the shaded triangle
      plt.fill(x_coords, y_coords, 'gray', alpha=0.3)  # 'b' is blue, alpha controls transparency
      # Write the order of convergence next to the filled triangle
      alpha = 1.075
      x_center = highlight_x[1,i]*alpha  # You can calculate the center of the triangle to place the text
      y_center = 5*(highlight_y[1,i] + highlight_y[2,i])/12
      plt.text(x_center, y_center, f"{v_O_u_mean_s[i]}", fontsize=order_font, ha='center', va='center')

# We do this because the files are not read in increasing order of degu and then degp
# Combine degu and degp into a list of tuples (degu, degp)
my_list_s = list(zip(v_degu_s, v_degp_s))

# Sort the list first by degu and then by degp (for each degu)
sorted_list_s = sorted(my_list_s, key=lambda x: (x[0], x[1]))

# Separate the sorted list into degu and degp again
sorted_degu_s, sorted_degp_s = zip(*sorted_list_s)

# Get the sorted indices based on the sorting of degu and degp
sorted_indices_s = [my_list_s.index(t) for t in sorted_list_s]

legend_labels_sorted_s = [legend_labels[i] for i in sorted_indices_s]
simplex_labels = legend_labels_sorted_s

lines = [lines[i] for i in sorted_indices_s]
simplex_handles = []

for line, color, marker, label in zip(lines, colors, markers, legend_labels_sorted_s):
    line.set_color(color)
    line.set_marker(marker)
    line.set_label(label)
    line.set_linestyle('--') 
    simplex_handles.append(line)

legend_handles = quad_handles + simplex_handles
legend_labels = quad_labels + simplex_labels

# Define the font properties explicitly
font_properties = font_manager.FontProperties(size=9.5)

# Apply the function along each column of uL2E
plt.xscale('log') 
plt.yscale('log') 

# This has to be after setting the log-scale
plt.xticks(x_tick_values, x_tick_labels)  # Set x-ticks and labels
plt.yticks(u_tick_values, u_tick_labels)  # Set y-ticks and labels

plt.xlabel('h')       # label for x-axis
plt.ylabel("$\\Vert e_{\\mathbf{u}}\\Vert_{2}$")  # label for y-axis

plt.ylim(u_range)

plt.legend(legend_handles, legend_labels, prop=font_properties , loc='lower right', ncol=2)
plt.savefig("order_of_convergence_velocity.png", dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()

### Plot for the L2 error on the pressure

### Quad mesh
lines = plt.plot(h_q, m_pL2E_q, label='Q Line - Pressure')

# Shaded triangular area to show the slope
highlight_x = np.log10(np.array([1e-2, 1.3e-2]))  # X-limits of the triangle

# We use the logarithm of the x and y coordinates to avoid having to fit polynomials of varying degrees
log_h = np.log10(h_q)  # log(h)
log_m_pL2E = np.log10(m_pL2E_q)  # log(uL2E)

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
    plt.text(x_center, y_center, f"{v_O_p_mean_q[i]}", fontsize=order_font, ha='center', va='center')

lines = [lines[i] for i in sorted_indices_q]
for line, color, marker, label in zip(lines, colors, markers, legend_labels_sorted_q):
    line.set_color(color)
    line.set_marker(marker)
    line.set_label(label)


### Simplex mesh
lines = plt.plot(h_s, m_pL2E_s, label='S Line - Pressure')
# Shaded triangular area to show the slope
highlight_x = np.log10(np.array([3.5e-2, 4.5e-2]))  # X-coordinates of the triangle
# We use the logarithm of the x and y coordinates to avoid having to fit polynomials of varying degrees
log_h = np.log10(h_q)  # log(h)
log_m_pL2E_s = np.log10(m_pL2E_s)  # log(uL2E)

# Apply polyfit to each column of log_uL2E
coeffs = np.apply_along_axis(fit_log, 0, log_h, log_m_pL2E_s)
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
    if (int(v_degu_s[i]) != 3):
      x_coords = highlight_x[:, i]  # Get the x-coordinates for the triangle
      y_coords = highlight_y[:, i]  # Get the y-coordinates for the triangle
    
      # Plot the shaded triangle
      plt.fill(x_coords, y_coords, 'gray', alpha=0.3)  # 'b' is blue, alpha controls transparency
      # Write the order of convergence next to the filled triangle
      alpha = 1.075
      x_center = highlight_x[1,i]*alpha  # You can calculate the center of the triangle to place the text
      y_center = 5*(highlight_y[1,i] + highlight_y[2,i])/12
      # y_center = highlight_y[2,i] 
      plt.text(x_center, y_center, f"{v_O_p_mean_s[i]}", fontsize=order_font, ha='center', va='center')


lines = [lines[i] for i in sorted_indices_s]

for line, color, marker, label in zip(lines, colors, markers, legend_labels_sorted_s):
    line.set_color(color)
    line.set_marker(marker)
    line.set_label(label)
    line.set_linestyle('--') 

plt.xscale('log') 
plt.yscale('log') 

plt.xticks(x_tick_values, x_tick_labels)  # Set x-ticks and labels
plt.yticks(p_tick_values, p_tick_labels)  # Set y-ticks and labels

plt.xlabel('h')       # label for x-axis
plt.ylabel('$\\Vert e_p\\Vert_{2}$')  # label for y-axis

plt.ylim(p_range)

# Define the font properties explicitly
font_properties = font_manager.FontProperties(size=9.5) 
plt.legend(legend_handles, legend_labels, prop=font_properties , loc='lower right', ncol=2)
plt.savefig("order_of_convergence_pressure.png", dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()

