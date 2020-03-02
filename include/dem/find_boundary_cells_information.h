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

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <iostream>
#include <vector>

#include "dem/boundary_cells_info_struct.h"

using namespace dealii;

#ifndef FINDBOUNDARYCELLSINFORMATION_H_
#  define FINDBOUNDARYCELLSINFORMATION_H_

/**
 * Finds all the boundary cells and faces in the triangulation, for each cell
 * the boundary face is specified and the normal vector as well as a point on
 * the boundary face are obtained
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class FindBoundaryCellsInformation
{
public:
  FindBoundaryCellsInformation<dim>();

  /**
   * Loops over all the cells to find boundary cells, find the boundary faces of
   * boundary cells and for the boundary faces the normal vector and a point
   * locating on the face are obtained
   *
   * @param triangulation Triangulation to access the information of the cells
   * @return A vector of structs (boundary_cells_info_struct), this struct
   * contains 1, an integer which shows the number of boundary face 2, the
   * corresponding boundary cell 3, the face number of the cell 4, normal vector
   * of the face 5, a point on the face which will be used to obtain the
   * distance between center of particles and the face
   */

  std::vector<boundary_cells_info_struct<dim>>
  find_boundary_cells_information(
    const parallel::distributed::Triangulation<dim> &triangulation);
};

#endif /* FINDBOUNDARYCELLSINFORMATION_H_ */
