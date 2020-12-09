/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */

#ifndef boundary_cells_info_struct_h
#define boundary_cells_info_struct_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/grid/tria.h>


/**
 * This struct handles the information related to the boundary cells which will
 * be used for particle-wall contact detection
 */

using namespace dealii;

template <int dim>
struct boundary_cells_info_struct
{
  // The boundary cell
  typename Triangulation<dim>::active_cell_iterator cell;

  // ID of boundary
  unsigned int boundary_id;

  // ID of boundary face in the boundary cell
  unsigned int boundary_face_id;

  // Normal vector of the boundary face
  Tensor<1, dim> normal_vector;

  // A point on the boundary face
  Point<dim> point_on_face;
};

#endif /* boundary_cells_info_struct_h */
