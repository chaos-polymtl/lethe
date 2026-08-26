# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Postprocessing script for the fichera-oven example.

This script walks through every adaptive refinement iteration stored in a
PVD/VTU time-series output, extracts the DPG (Discontinuous Petrov-Galerkin)
error indicator computed by the solver, and builds two convergence plots:

  1. DPG error norm vs. number of cells.
  2. DPG error norm vs. number of degrees of freedom (compared against a
     reference hp3d solution from the literature).

If run with ``--validate``, the script skips the interactive plot windows
and instead dumps the raw convergence data to a ``.dat`` file so it can be
compared against a reference solution in an automated test.

Usage
-----
    python fichera_oven_postprocess.py -f /path/to/output
    python fichera_oven_postprocess.py --validate
"""

#############################################################################
# Imports
#############################################################################

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from cycler import cycler

#############################################################################
# Physical / discretization constants
#
# These assumptions describe the finite element spaces used by the solver
# and are needed to convert cell/face/edge counts into degree-of-freedom
# (DoF) counts. They must be kept in sync with the actual simulation
# parameters (polynomial degree, number of fields) used to produce the
# output being postprocessed.
#############################################################################

# Number of solution fields being solved for: E_r, E_i, H_r, H_i
# (real and imaginary parts of the electric and magnetic fields).
N_FIELDS = 4

# Interior (volumetric) DG space: degree-3 DG -> 64 shape functions per
# cell, each field has 3 vector components in 3D.
N_DOFS_PER_CELL_PER_FIELD = 27
N_COMPONENTS_3D = 3
N_DOFS_INTERIOR_PER_CELL = N_DOFS_PER_CELL_PER_FIELD * N_COMPONENTS_3D  # = 81

# Skeleton (trace) space: degree-4 Nedelec face element -> 24 DoFs
# associated with the interior of a face, plus 4 DoFs per edge.
N_DOFS_PER_FACE = 12
N_DOFS_PER_EDGE = 3

#############################################################################
# Reference data
#
# Digitized DoF-vs-error convergence data for the hp3d code, taken from the
# reference paper https://doi.org/10.1016/j.camwa.2021.01.017. Used only to
# provide a comparison curve on the second convergence plot.
#############################################################################

N_DOFS_REF = np.array([
    2028, 8942, 37802, 88782, 183840, 321810, 469604, 612680,
    1072494, 2348662, 3332700, 4471528, 4931841, 5916114, 9003416
])
RESIDUALS_REF = np.array([
    0.2929607385716628, 0.26155250662763546, 0.2064311160765089,
    0.17630014791054022, 0.14191694747369588, 0.10874331535504288,
    0.08332414729113355, 0.07046374783289007, 0.05561372856114809,
    0.04096563585421059, 0.03430298439717367, 0.029296073857166302,
    0.027342035013820435, 0.02477446488053953, 0.020340010731841072
])


#############################################################################
# Helper functions
#############################################################################

def parse_args():
    """Parse command-line arguments for the fichera-oven postprocessing."""
    parser = argparse.ArgumentParser(
        description="Arguments for the validation of the fichera-oven benchmark")
    parser.add_argument(
        "--validate", action="store_true", default=False,
        help="Launch the script in validation mode. This logs the content "
             "of the graph to a .dat file and suppresses the display of "
             "figures (used for automated regression testing).")
    parser.add_argument(
        "-f", "--folder", type=str, required=False,
        help="Path to the output folder for the fichera-oven results. This "
             "is the folder that contains the results of the simulation "
             "(.vtu, .pvtu, .dat and .pvd files). Defaults to './output'.")
    args, _leftovers = parser.parse_known_args()
    return args


def setup_plot_style():
    """Configure global matplotlib rcParams for consistent, publication-style figures."""
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
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'cm'
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams.update({
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsfonts}',
    })

    return colors


def count_faces_and_edges(mesh):
    """Derive interior/boundary face and edge counts for a hexahedral mesh.

    PyVista does not directly expose the number of *interior* faces or
    edges, so we derive them from the boundary surface and simple counting
    relations for a mesh of hexahedra (6 faces per cell). This assumes that the
    mesh is a conforming hexahedral mesh and that there is no subdivision of 
    the cell to improve the rendering of high-order elements in Paraview.

    Parameters
    ----------
    mesh : pyvista.DataSet
        The (already combined) mesh for the current refinement iteration.

    Returns
    -------
    n_face_boundary : int
        Number of faces lying on the domain boundary.
    n_face_interior : int
        Number of faces shared between two cells (interior faces).
    n_edge_boundary : int
        Number of edges lying on the domain boundary surface.
    n_edge_interior : int
        Number of edges shared between four cells (interior edges).
    """
    # Extracting the outer surface once and reusing it avoids recomputing
    # the same boundary extraction twice (once for faces, once for edges).
    n_cells = mesh.n_cells
    cleaned = mesh.clean()
    boundary_surface = cleaned.extract_surface(algorithm="dataset_surface")

    n_face_boundary = boundary_surface.n_cells

    # Each hexahedral cell has 6 faces. Every interior face is shared by
    # exactly two cells, every boundary face by exactly one:
    #   6 * n_cells = n_face_boundary + 2 * n_face_interior
    # Solving for the number of interior faces:
    n_face_interior = (6 * n_cells - n_face_boundary) / 2
    n_edge_boundary = boundary_surface.extract_all_edges().n_lines
    n_edge_interior = cleaned.extract_all_edges().n_lines - n_edge_boundary

    print(f"n_cells={n_cells}, n_face_boundary={n_face_boundary}, "
          f"n_face_interior={n_face_interior}, n_edge_boundary={n_edge_boundary}, "
          f"n_edge_interior={n_edge_interior}")

    return n_face_boundary, n_face_interior, n_edge_boundary, n_edge_interior





def compute_dof_counts(n_cells, n_face_boundary, n_faces_interior, n_edge_boundary, n_edge_interior):
    """Compute interior and skeleton (trace) DoF counts for one iteration.

    Assumes a degree-4 Nedelec face element for the skeleton space (24
    DoFs per face interior + 4 DoFs per edge) and a degree-3 DG space for
    the interior (64 DoFs per cell, 3 vector components), each replicated
    across the 4 solved fields (E_r, E_i, H_r, H_i).

    Returns
    -------
    n_dofs_interior : int
        Total interior-solution DoFs (all fields, all components).
    n_dofs_skeleton : int
        Total skeleton/trace-space DoFs (all fields).
    """
    n_faces = n_face_boundary + n_faces_interior
    n_edge_all = n_edge_boundary + n_edge_interior
    n_dofs_skeleton = N_FIELDS * (
        n_faces * N_DOFS_PER_FACE + n_edge_all * N_DOFS_PER_EDGE
    )
    #Only the boundary dofs of E_r are constrained, so we need to subtract them from the total dofs to get the number of free dofs
    n_dofs_skeleton_constrained = n_face_boundary * N_DOFS_PER_FACE + n_edge_boundary * N_DOFS_PER_EDGE
    n_dofs_skeleton -= n_dofs_skeleton_constrained 

    n_dofs_interior = n_cells * N_FIELDS * N_DOFS_INTERIOR_PER_CELL

    print(f"n_dofs_interior={n_dofs_interior}, n_dofs_skeleton={n_dofs_skeleton}")

    return n_dofs_interior, n_dofs_skeleton


def compute_error_norms(mesh, dpg_error_point_data):
    """Compute the global L2-like and max DPG error norms for one iteration.

    The point-data array ``dpg_error_norm`` already stores, per cell (after
    conversion), a quantity that has been pre-multiplied by the cell
    volume by the solver. Consequently the global L2-like norm can be
    obtained simply as the square root of the sum of squared per-cell
    contributions, i.e. sqrt( sum_K e_K^2 ), without an extra volume
    weighting step.

    Returns
    -------
    l2_error : float
        sqrt( sum_K e_K^2 ), a global L2-like error norm.
    max_error : float
        The maximum pointwise error indicator over the mesh.
    """
    cell_data = mesh.point_data_to_cell_data()
    cell_err = cell_data['dpg_error_norm']

    l2_error = np.sqrt(np.sum(cell_err ** 2))
    max_error = dpg_error_point_data.max()

    return l2_error, max_error


def read_convergence_data(pvd_file):
    """Read every timestep of a PVD file and extract convergence metrics.

    Each timestep in the PVD corresponds to one adaptive refinement
    iteration. The very first iteration is skipped, since nothing is 
    computed at that iteration as it is the initial condition.

    Returns
    -------
    dict of str -> np.ndarray
        Arrays (all the same length, one entry per retained iteration):
        'iterations', 'n_cells', 'n_faces', 'n_dofs_interior',
        'n_dofs_skeleton', 'l2_error', 'max_error'.
    """
    reader = pv.get_reader(pvd_file)

    iterations = []
    n_cells_list = []
    n_face_list = []
    n_dofs_interior_list = []
    n_dofs_skeleton_list = []
    l2_error_list = []
    max_error_list = []

    for t in reader.time_values:
        reader.set_active_time_value(t)
        data = reader.read()

        # Parallel runs produce a MultiBlock dataset (one block per rank);
        # combine them into a single mesh before doing any geometric or
        # field-based computation.
        if hasattr(data, 'combine'):
            data = data.combine()

        n_cells = data.n_cells
        dpg_error = data['dpg_error_norm']

        # Skip the initial iteration: since no solution is computed at that 
        # iteration, the DPG error indicator is zero so we can use it to filter 
        # out the first iteration. This avoids a spurious zero in the 
        # convergence plots.
        if dpg_error.max() == 0.0:
            continue

        n_face_boundary, n_face_interior, n_edge_boundary, n_edge_interior = \
            count_faces_and_edges(data)

        n_dofs_interior, n_dofs_skeleton = compute_dof_counts(
            n_cells, n_face_boundary, n_face_interior, n_edge_boundary, n_edge_interior
        )

        l2_error, max_error = compute_error_norms(data, dpg_error)

        iterations.append(int(t))
        n_cells_list.append(n_cells)
        n_face_list.append(n_face_interior + n_face_boundary)
        n_dofs_interior_list.append(n_dofs_interior)
        n_dofs_skeleton_list.append(n_dofs_skeleton)
        l2_error_list.append(l2_error)
        max_error_list.append(max_error)

    return {
        'iterations': np.array(iterations),
        'n_cells': np.array(n_cells_list),
        'n_faces': np.array(n_face_list),
        'n_dofs_interior': np.array(n_dofs_interior_list),
        'n_dofs_skeleton': np.array(n_dofs_skeleton_list),
        'l2_error': np.array(l2_error_list),
        'max_error': np.array(max_error_list),
    }


def plot_convergence_vs_cells(results, colors, show):
    """Plot DPG L2 and max error norms as a function of the number of cells."""
    fig, ax = plt.subplots()

    ax.loglog(results['n_cells'], results['l2_error'],
                's-', color=colors[0], markerfacecolor='none',
                label=r'$\| e \|_{L^2(\Omega)}$')
    ax.loglog(results['n_cells'], results['max_error'],
                'o-', color=colors[1], markerfacecolor='none',
                label=r'$\| e \|_{L^\infty(\Omega)}$')

    ax.set_xlabel(r'Number of cells')
    ax.set_ylabel(r'DPG error norm')
    ax.legend()
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

    fig.tight_layout()
    fig.savefig('fichera_oven_convergence_cell.pdf')

    if show:
        fig.show()

    fig.clear()


def plot_convergence_vs_dofs(results, show):
    """Plot the DPG L2 error norm vs. DoF count, compared against hp3d reference data."""
    fig, ax = plt.subplots()

    ax.loglog(N_DOFS_REF, RESIDUALS_REF,
              linestyle='-', color="#0e26b1", linewidth=4, label=r'\texttt{hp3d}',zorder=3)
    ax.scatter(N_DOFS_REF, RESIDUALS_REF,
               facecolor='k', color='none', zorder=4, s=75, linewidths=2)

    ax.loglog(results['n_dofs_skeleton'], results['l2_error'],
              linestyle='-', color='#1b9e77', linewidth=4, label=r'\texttt{lethe}',zorder=5)
    ax.scatter(results['n_dofs_skeleton'], results['l2_error'],
               facecolor='k', color='none', zorder=6, s=75, linewidths=2)

    ax.grid(True, which="major", ls="-", lw=0.4, alpha=0.35)
    ax.grid(True, which="minor", ls=":", lw=0.3, alpha=0.25)
    ax.minorticks_on()
    ax.legend(loc='upper right', frameon=True, fontsize=14)
    ax.set_ylim(1e-2, 1)
    ax.set_xlabel(r'Number of degrees of freedom')
    ax.set_ylabel(r'$\|u_h-u\|_E$')

    fig.savefig("fichera_oven_convergence_dofs.pdf", bbox_inches="tight")

    if show:
        fig.show()

    fig.clear()


def save_validation_data(results, filename="solution-fichera-oven.dat"):
    """Dump the raw convergence data to a text file for automated validation.
    """
    solution_data = np.column_stack([
        results['iterations'],
        results['n_cells'],
        results['n_faces'],
        results['n_dofs_interior'],
        results['n_dofs_skeleton'],
        results['l2_error'],
        results['max_error'],
    ])
    header = "iteration n_cells n_faces n_dofs_interior n_dofs_skeleton l2_error max_error"
    np.savetxt(filename, solution_data, header=header)


#############################################################################
# Main
#############################################################################

def main():
    args = parse_args()
    colors = setup_plot_style()

    output_path = args.folder if args.folder else "./output"
    pvd_file = f"{output_path}/out.pvd"

    results = read_convergence_data(pvd_file)

    show_figures = not args.validate
    plot_convergence_vs_cells(results, colors, show=show_figures)
    plot_convergence_vs_dofs(results, show=show_figures)

    if args.validate:
        save_validation_data(results)


if __name__ == "__main__":
    main()