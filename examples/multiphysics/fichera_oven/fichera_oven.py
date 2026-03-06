# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#############################################################################
"""
Postprocessing code for the fichera-oven example.

Reads the DPG error norm from each adaptive refinement iteration stored in
a PVD/VTU output and produces a convergence plot (error vs. iteration and
error vs. number of cells).

"""
#############################################################################

'''Importing Libraries'''
import numpy as np
import sys
import matplotlib.pyplot as plt
import pyvista as pv
import argparse

#############################################################################

'''Plot formatting'''

from cycler import cycler

colors = ['#008c66', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']

plt.rcParams['axes.prop_cycle'] = cycler(color=colors)
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['lines.markersize'] = '11'
plt.rcParams['markers.fillstyle'] = 'none'
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.columnspacing'] = 2
plt.rcParams['legend.handlelength'] = 2.8
plt.rcParams['legend.handletextpad'] = 0.2
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.fontsize'] = '18'
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['font.size'] = '25'
plt.rcParams['font.family'] = 'DejaVu Serif'
plt.rcParams['font.serif'] = 'cm'
plt.rcParams['savefig.bbox'] = 'tight'

plt.rcParams.update({
    'text.usetex': True,
    'text.latex.preamble': r'\usepackage{amsfonts}'
})

#############################################################################

parser = argparse.ArgumentParser(
    description='Arguments for the validation of the fichera-oven benchmark')
parser.add_argument(
    "--validate", action="store_true",
    help="Launches the script in validation mode. This will log the content "
         "of the graph and prevent the display of figures",
    default=False)
parser.add_argument(
    "-f", "--folder", type=str,
    help="Path to the output folder for the fichera-oven results. This is "
         "the folder that contains the results of the simulation "
         "(.vtu, .pvtu, .dat and .pvd files)",
    required=False)

args, leftovers = parser.parse_known_args()

#############################################################################
# Read the PVD file and extract error data at each refinement iteration
#############################################################################

# Resolve output path
output_path = args.folder if args.folder else "./output"
pvd_file = f"{output_path}/out.pvd"

reader = pv.get_reader(pvd_file)
time_values = reader.time_values

iterations = []
n_cells_list = []
n_dofs_list = []
l2_error_list = []
max_error_list = []

for t in time_values:
    reader.set_active_time_value(t)
    data = reader.read()

    # Handle MultiBlock datasets (parallel output)
    if hasattr(data, 'combine'):
        data = data.combine()

    iteration = int(t)
    n_cells = data.n_cells
    n_dofs = data.n_points

    # Extract the DPG error norm
    dpg_error = data['dpg_error_norm']

    # Skip the initial iteration if no error has been computed yet
    if dpg_error.max() == 0.0:
        continue

    # Compute L2-like error norm: sqrt( sum_K ( e_K^2 * |K| ) )
    sized = data.compute_cell_sizes(length=False, area=False, volume=True)
    cell_volumes = sized.cell_data['Volume']

    cell_data = data.point_data_to_cell_data()
    cell_err = cell_data['dpg_error_norm']

    l2_error = np.sqrt(np.sum(cell_err**2 * cell_volumes))
    max_error = dpg_error.max()

    iterations.append(iteration)
    n_cells_list.append(n_cells)
    n_dofs_list.append(n_dofs)
    l2_error_list.append(l2_error)
    max_error_list.append(max_error)

iterations = np.array(iterations)
n_cells_arr = np.array(n_cells_list)
n_dofs_arr = np.array(n_dofs_list)
l2_error_arr = np.array(l2_error_list)
max_error_arr = np.array(max_error_list)

#############################################################################
# Convergence plot — DPG error norm vs. number of cells
#############################################################################

fig, ax = plt.subplots()

ax.semilogy(n_cells_arr, l2_error_arr,
            's-', color=colors[0], markerfacecolor='none',
            label=r'$\| e \|_{L^2(\Omega)}$')
ax.semilogy(n_cells_arr, max_error_arr,
            'o-', color=colors[1], markerfacecolor='none',
            label=r'$\| e \|_{L^\infty(\Omega)}$')

ax.set_xlabel(r'Number of cells')
ax.set_ylabel(r'DPG error norm')
ax.legend()
ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

plt.tight_layout()
plt.savefig('fichera_oven_convergence.pdf')

if not args.validate:
    plt.show()

#############################################################################
# Save data for validation
#############################################################################

if args.validate:
    solution_data = np.column_stack([
        iterations, n_cells_arr, n_dofs_arr, l2_error_arr, max_error_arr
    ])
    #The number of dofs is for the interior solution and is only for one field component. This is not the number of dofs in the linear system, which only involves the trace dofs which are not outputed.
    header = "iteration n_cells n_dofs_per_component l2_error max_error"
    np.savetxt("solution-fichera-oven.dat", solution_data, header=header)


