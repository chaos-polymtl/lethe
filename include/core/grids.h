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

template <int dim>
void
attach_grid_to_triangulation(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  const Parameters::Mesh &                                     mesh_parameters,
  const BoundaryConditions::BoundaryConditions<dim> &boundary_conditions);

template <int dim>
void
read_mesh_and_manifolds(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  const Parameters::Mesh &                                     mesh_parameters,
  const Parameters::Manifolds &                      manifolds_parameters,
  const BoundaryConditions::BoundaryConditions<dim> &boundary_conditions);


#endif
