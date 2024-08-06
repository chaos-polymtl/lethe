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

#ifndef lethe_boundary_cells_info_struct_h
#define lethe_boundary_cells_info_struct_h

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

using namespace dealii;

/**
 * @brief Handle the information related to the boundary cells
 * which will be used for particle-wall contact detection.
 */
template <int dim>
struct boundary_cells_info_struct
{
  // The boundary cell
  typename Triangulation<dim>::active_cell_iterator cell;

  // ID of boundary
  types::boundary_id boundary_id;

  // ID of boundary face in the boundary cell
  unsigned int global_face_id;

  // Normal vector of the boundary face
  Tensor<1, dim> normal_vector;

  // A point on the boundary face
  Point<dim> point_on_face;
};

/**
 * @brief Handle the information related to the periodic boundary cells
 * which will be used for particle-wall contact detection.
 */
template <int dim>
struct periodic_boundaries_cells_info_struct
{
  // The boundary cell
  typename Triangulation<dim>::active_cell_iterator cell;

  // ID of boundary
  types::boundary_id boundary_id;

  // Normal vector of the boundary face
  Tensor<1, dim> normal_vector;

  // A point on the boundary face
  Point<dim> point_on_face;

  // The periodic boundary cell
  typename Triangulation<dim>::active_cell_iterator periodic_cell;

  // Normal vector of the periodic boundary face
  Tensor<1, dim> periodic_normal_vector;

  // A point on the peridic boundary face
  Point<dim> point_on_periodic_face;
};

#endif
