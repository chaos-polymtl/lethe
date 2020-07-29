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

#include <deal.II/grid/grid_tools.h>

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
   * @param boundary_cells_with_faces A vector which contains all the boundary
   * cells which has atleast one boundary face
   * @param triangulation Triangulation to access the information of the cells
   * @return A map of structs (boundary_cells_info_struct), this struct
   * contains 1, an integer which shows the number of boundary face 2, the
   * corresponding boundary cell 3, the face number of the cell 4, normal vector
   * of the face 5, a point on the face which will be used to obtain the
   * distance between center of particles and the face
   */

  std::map<int, boundary_cells_info_struct<dim>>
  find_boundary_cells_information(
    std::vector<typename Triangulation<dim>::active_cell_iterator>
      &                                              boundary_cells_with_faces,
    const parallel::distributed::Triangulation<dim> &triangulation);

  /**
   * Loops over all the cells to find cells which should be searched for
   * particle-line contact, boundary lines and a point locating on each line are
   * obtained
   *
   * @param boundary_cells_with_faces A vector which contains all the boundary
   * cells which has atleast one boundary face
   * @param triangulation Triangulation to access the information of the cells
   * @param boundary_cells_with_lines A vector of tuples which contains the
   * cells with boundary lines and locations of the beginning and ending
   * vertices of the boundary line
   * @param boundary_cells_with_points A vector of pairs which contains the
   * cells with boundary points, and the location of the point
   */

  void
  find_particle_point_and_line_contact_cells(
    const std::vector<typename Triangulation<dim>::active_cell_iterator>
      &                                              boundary_cells_with_faces,
    const parallel::distributed::Triangulation<dim> &triangulation,
    std::vector<std::tuple<typename Triangulation<dim>::active_cell_iterator,
                           Point<dim>,
                           Point<dim>>> &            boundary_cells_with_lines,
    std::vector<std::pair<typename Triangulation<dim>::active_cell_iterator,
                          Point<dim>>> &boundary_cells_with_points);

private:
  // This vector stores the location of vertices on boundaries for each cell.
  // The size of this vector can be 1 or 2, since cells with points have one
  // boundary vertex and cells with lines have two boundary vertices
  std::vector<Point<dim>> boundary_points;

  // This integer is used to count the number of vertices on boundaries for each
  // cell
  unsigned int number_of_boundary_vertices = 0;

  // This vector stores both the cells with boundary lines and cells with
  // boundary points
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    boundary_cells_with_lines_or_points;

  // This map stores the vertex index and position of boundary vertices
  std::map<int, Point<dim>> boundary_vertices;
};

#endif /* FINDBOUNDARYCELLSINFORMATION_H_ */
