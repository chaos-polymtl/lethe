// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_read_mesh_h
#define lethe_read_mesh_h

#include <dem/dem_solver_parameters.h>

/**
 * @brief Define or read the mesh based on the information provided by the
 * user. Gmsh file can also be read in this function
 *
 * @param mesh_parameters Input DEM parameters in the parameter handler file
 * @param pcout Printing in parallel
 * @param restart If in restart situation
 * @param triangulation Triangulation
 * @param bc_params Boundary conditions parameters for DEM
 */
template <int dim, int spacedim = dim>
void
read_mesh(const Parameters::Mesh   &mesh_parameters,
          const bool                restart,
          const ConditionalOStream &pcout,
          parallel::DistributedTriangulationBase<dim, spacedim> &triangulation,
          const Parameters::Lagrangian::BCDEM                   &bc_params);

/**
 * @brief Allow periodic faces mapping in triangulation of mesh with boundary
 * ids parameters.
 *
 * @param triangulation Triangulation
 * @param bc_params Boundary conditions parameters for DEM
 */
template <int dim, int spacedim>
void
match_periodic_boundaries(Triangulation<dim, spacedim>        &triangulation,
                          const Parameters::Lagrangian::BCDEM &bc_params);

#endif
