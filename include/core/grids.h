/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 -  by the Lethe authors
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

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2020 -
 */

#ifndef lethe_grids_h
#define lethe_grids_h

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <core/boundary_conditions.h>
#include <core/manifolds.h>
#include <core/parameters.h>



using namespace dealii;


/**
 * @brief Attaches a grid to a triangulation using mesh parameters
 *
 * @param triangulation The triangulation to which a grid is attached
 *
 * @param mesh_parameters The mesh parameters used to decide what type of mesh or primitive is  used
 *
 * @param boundary_conditions The information about the boundary conditions id. This is used to set-up the periodicity of the domain
 */
template <int dim, int spacedim = dim>
void
attach_grid_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  const Parameters::Mesh &                                     mesh_parameters,
  const BoundaryConditions::BoundaryConditions<dim> &boundary_conditions);

/**
 * @brief Completely set-up a mesh and its manifold
 *
 * @param triangulation The triangulation to which a grid is attached
 *
 * @param mesh_parameters The mesh parameters used to decide what type of mesh or primitive is  used
 *
 * @param manifolds_parameters The information about the type of manifolds attached to the boundary conditions
 *
 * @param boundary_conditions The information about the boundary conditions id. This is used to set-up the periodicity of the domain
 */
template <int dim, int spacedim = dim>
void
read_mesh_and_manifolds(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  const Parameters::Mesh &                                     mesh_parameters,
  const Parameters::Manifolds &                      manifolds_parameters,
  const BoundaryConditions::BoundaryConditions<dim> &boundary_conditions);


#endif
