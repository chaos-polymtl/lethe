/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

#ifndef lethe_read_mesh_h
#define lethe_read_mesh_h

#include <dem/dem_solver_parameters.h>

#include <deal.II/distributed/tria.h>

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
