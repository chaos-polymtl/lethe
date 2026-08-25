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
n_face_list = []
n_dofs_interior_list = []
n_dofs_skeleton_list = []
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

    # Extract the DPG error norm
    dpg_error = data['dpg_error_norm']

    # Skip the initial iteration if no error has been computed yet
    if dpg_error.max() == 0.0:
        continue

    # N_total_face = N_face_boundary + N_face_interior -> 6*N_cells = N_face_boundary + 2*N_face_interior -> N_face_interior = (6*N_cells - N_face_boundary)/2
    n_face_boundary = data.clean().extract_surface(algorithm="dataset_surface").n_cells
    n_face_interior = ( 6 * n_cells - n_face_boundary) /2 

    #We also need the edge information to compute the number of dofs in the skeleton space. 
    n_edge_boundary = data.clean().extract_surface(algorithm='dataset_surface').extract_all_edges().n_lines
    n_edge_all = data.clean().extract_all_edges().n_lines

    #Assuming degree 4 nedelec face (degree 3 DG for interior), the number of dofs for the skeleton space is 24 for the interior of a face and 4 per edges. There are 4 fields (E_r, E_i, H_r, H_i), so the number of dofs for the skeleton space is 4*((n_face_boundary + n_face_interior)* 24 + n_edge_all * 4)
    n_dofs_skeleton_list.append(4*((n_face_boundary + n_face_interior)* 24 + n_edge_all * 4))

    #Assuming degree 3 DG for interior, so 64 dofs per cell, the number of dofs for the interior solution is 64 dofs * number of field * number of components. We have E_r, E_i, H_r, H_i, so 4 field. For the interior solution there are 3 components in 3D, so 4*3*64 = 768 dofs.
    n_dofs_interior_list.append(n_cells * 768)

    # Compute L2-like error norm: sqrt( sum_K ( e_K^2 * |K| ) ). The output of DPG error is already multiplied by the cell volume, so we can just sum the squares of the error and take the square root.
    sized = data.compute_cell_sizes(length=False, area=False, volume=True)

    cell_data = data.point_data_to_cell_data()
    cell_err = cell_data['dpg_error_norm']

    l2_error = np.sqrt(np.sum(cell_err**2))
    max_error = dpg_error.max()

    iterations.append(iteration)
    n_cells_list.append(n_cells)
    n_face_list.append(n_face_interior + n_face_boundary)
    l2_error_list.append(l2_error)
    max_error_list.append(max_error)

iterations = np.array(iterations)
n_cells_arr = np.array(n_cells_list)
n_face_arr = np.array(n_face_list)
n_dofs_interior_arr = np.array(n_dofs_interior_list)
n_dofs_face_arr = np.array(n_dofs_skeleton_list)
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
plt.savefig('fichera_oven_convergence_cell.pdf')

if not args.validate:
    plt.show()

#############################################################################
# Convergence plot — DPG error norm vs. number of dofs
#############################################################################

# Data from the reference https://doi.org/10.1016/j.camwa.2021.01.017,.
n_dofs_ref= np.array([2076, 9070, 37956, 89400, 184874, 324918, 469785, 616082, 1071094, 2338390, 3344512, 4482201, 4941714, 5942106, 9169030])
residuals_ref = np.array([0.2912234144941983, 0.2596028668687208, 0.20471425561217332, 0.1769783464666799, 0.14063276723208087, 0.10755088210362346, 0.08288371676687098, 0.07056440607367706, 0.05479850606955148, 0.04095555466309438, 0.033815680884222044, 0.029010979756341898, 0.027495896657863637, 0.02451043853017416, 0.020083001420030503])


# Plotting the convergence of the residuals
fig, ax = plt.subplots(figsize=(5.5, 4.5))

ax.loglog(n_dofs_ref, residuals_ref, linestyle='-', color="#0e26b1", linewidth=4, label=r'\texttt{hp3d}')
ax.scatter(n_dofs_ref, residuals_ref, facecolor='k', color='none', zorder=5, s=75, linewidths=2)
ax.loglog(n_dofs_face_arr, l2_error_arr, linestyle='-', color='#1b9e77', linewidth=4, label=r'\texttt{lethe}')
ax.scatter(n_dofs_face_arr, l2_error_arr, facecolor='k', color='none', zorder=5, s=75, linewidths=2)


ax.grid(True, which="major", ls="-", lw=0.4, alpha=0.35)
ax.grid(True, which="minor", ls=":", lw=0.3, alpha=0.25)
ax.minorticks_on()
ax.legend(loc='upper right', frameon=True, fontsize=14)
ax.set_ylim(1e-2, 1)
ax.set_xlabel(r'Number of degrees of freedom')
ax.set_ylabel(r'$\|u_h-u\|_E$')
fig.savefig("fichera_oven_convergence_dofs.pdf",  bbox_inches="tight")

if not args.validate:
    plt.show()

#############################################################################
# Save data for validation
#############################################################################

if args.validate:
    solution_data = np.column_stack([
        iterations, n_cells_arr, n_face_arr, n_dofs_interior_arr, n_dofs_face_arr, l2_error_arr, max_error_arr
    ])
    #The number of dofs is for the interior solution and is only for one field component. This is not the number of dofs in the linear system, which only involves the trace dofs which are not outputed.
    header = "iteration n_cells  n_face n_dofs_interior n_dofs_face l2_error max_error"
    np.savetxt("solution-fichera-oven.dat", solution_data, header=header)


