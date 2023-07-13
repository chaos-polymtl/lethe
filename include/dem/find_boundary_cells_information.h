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

#include <dem/boundary_cells_info_struct.h>
#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>

#include <iostream>
#include <unordered_set>
#include <vector>

using namespace dealii;

#ifndef find_boundary_cells_information_h
#  define find_boundary_cells_information_h

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
class BoundaryCellsInformation
{
public:
  BoundaryCellsInformation<dim>();

  /**
   * @brief The build function builds all the boundary cell information
   * structure that finds which cell has a boundary face, boundary line or a
   * boundary point. Additionally, it identifies which cells are located in the
   * vicinity (controlled by the max cell diameter of floating  walls
   * @param triangulation Triangulation to access the information of the cells
   * @param floating_wall_properties Properties of floating walls specified in
   * the parameter handler file
   * @param outlet_boundaries A vector which contains the outlet boundary IDs
   * @param check_diamond_cells If true, the diamond shaped cells are found
   * and added to the particle-wall contact search cells
   * @param expand_particle_wall_contact_search If true, expands the
   * particle-wall contact search by adding the contact search between
   * particles in the boundary cells with the boundary neighbors of these
   * boundary cells
   *
   * @param pcout
   */
  void
  build(
    const parallel::distributed::Triangulation<dim> & triangulation,
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
    const std::vector<unsigned int> &                 outlet_boundaries,
    const bool &                                      check_diamond_cells,
    const bool &              expand_particle_wall_contact_search,
    const ConditionalOStream &pcout);

  void
  build(const parallel::distributed::Triangulation<dim> &triangulation,
        const std::vector<unsigned int> &                outlet_boundaries,
        const bool &                                     check_diamond_cells,
        const ConditionalOStream &                       pcout);

  std::map<int, boundary_cells_info_struct<dim>> &
  get_boundary_cells_information()
  {
    return boundary_cells_information;
  }

  std::unordered_map<
    std::string,
    std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>> &
  get_boundary_cells_with_points()
  {
    return boundary_cells_with_points;
  }

  std::unordered_map<
    std::string,
    std::tuple<typename Triangulation<dim>::active_cell_iterator,
               Point<dim>,
               Point<dim>>> &
  get_boundary_cells_with_lines()
  {
    return boundary_cells_with_lines;
  }

  std::unordered_map<
    types::particle_index,
    std::set<typename Triangulation<dim>::active_cell_iterator>> &
  get_boundary_cells_with_floating_walls()
  {
    return boundary_cells_for_floating_walls;
  }

  /**
   * Carries out updating the boundary information (point on boundary face and
   * its normal vector) after moving the grid.
   *
   * @param updated_boundary_points_and_normal_vectors Normal vector and a point
   * on the boundary faces which are updated after motion of the grid
   */
  void
  update_boundary_info_after_grid_motion(
    typename DEM::dem_data_structures<dim>::boundary_points_and_normal_vectors
      &updated_boundary_points_and_normal_vectors);

private:
  /**
   * Loops over all the cells to find boundary cells, find the boundary faces of
   * boundary cells and for the boundary faces the normal vector and a point
   * locating on the face are obtained
   *
   * @param triangulation Triangulation to access the information of the cells
   * @param outlet_boundaries A vector which contains the outlet boundary IDs
   */
  void
  find_boundary_cells_information(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const std::vector<unsigned int> &                outlet_boundaries);

  /**
   * Loops over all the cells to find boundary cells, find the global boundary
   * cells with faces. Note that the boundary_cells_with_faces container only
   * stores the local boundary cells, while global_boundary_cells_with_faces
   * stores all the boundary cells with faces.
   *
   * @param triangulation Triangulation to access the information of the cells
   * @param outlet_boundaries A vector which contains the outlet boundary IDs
   */
  void
  find_global_boundary_cells_with_faces(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const std::vector<unsigned int> &                outlet_boundaries);

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
    const parallel::distributed::Triangulation<dim> &triangulation);

  /**
   * Adds the cells with boundary lines to the boundary cells information
   * (boundary_cells_information), First , it loops through the
   * boundary_cells_with_lines vector, then it loops through the neighbors of
   * boundary cells with lines, if a neighbor is a member of
   * boundary_cells_with_faces, then it adds the neighbor cell to the
   * boundary_cells_information with the normal vector and point of the neighbor
   * cell.
   *
   * @param triangulation Triangulation to access the information of the cells
   * @param outlet_boundaries A vector which contains the outlet boundary IDs
   * @param check_diamond_cells If true, the diamond shaped cells are found and added to the particle-wall contact search cells
   * @param pcout
   */
  void
  add_cells_with_boundary_lines_to_boundary_cells(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const std::vector<unsigned int> &                outlet_boundaries,
    const bool &                                     check_diamond_cells,
    const ConditionalOStream &                       pcout);

  /**
   * Loops over all the cells to find cells which should be searched for
   * particle-floating wall contact, normal vector of each cell is also obtained
   *
   * @param triangulation Triangulation to access the information of the cells
   * @param floating_wall_properties Properties of floating walls specified in
   * the parameter handler file
   * @param boundary_cells_for_floating_walls An unordered_set which contains
   * all the boundary cells of floating walls
   * @param maximum_cell_diameter Maximum cell length in the triangulation
   */
  void
  find_boundary_cells_for_floating_walls(
    const parallel::distributed::Triangulation<dim> & triangulation,
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
    const double                                      maximum_cell_diameter);

  /**
   * Carries out adding new elements with the boundary cells (cells with
   * boundary faces) and information (point and normal vector) of their neighbor
   * boundary cells to the boundary_cells_information. By adding these
   * elements, the broad and fine searches will search for potential contacts
   * of the particles in each boundary cell with the boundary neighbors of these
   * cells. This function is used to solve the problem of late contact detection
   * in triangulations with curved boundaries.
   *
   * @param boundary_cells_information A container that contains the information
   * of all the boundary cells with boundary faces
   * @param outlet_boundaries A vector which contains the outlet boundary IDs
   * @param global_boundary_cells_information A vector that contains the geometrical
   * information of all (global) boundary cells
   */
  void
  add_boundary_neighbors_of_boundary_cells(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const std::vector<unsigned int> &                outlet_boundaries,
    std::map<int, boundary_cells_info_struct<dim>> & boundary_cells_information,
    const std::map<int, boundary_cells_info_struct<dim>>
      &global_boundary_cells_information);

  // Structure that contains the necessary information for boundaries
  std::map<int, boundary_cells_info_struct<dim>> boundary_cells_information;

  // A vector that contains the geometrical information of all (global) boundary
  // cells. This vector is used in
  // add_cells_with_boundary_lines_to_boundary_cells function
  std::map<int, boundary_cells_info_struct<dim>>
    global_boundary_cells_information;

  // Structure that contains the boundary cells which have a line
  std::unordered_map<
    std::string,
    std::tuple<typename Triangulation<dim>::active_cell_iterator,
               Point<dim>,
               Point<dim>>>
    boundary_cells_with_lines;

  // A vector of all the local cells with boundary lines. This vector is used in
  // add_cells_with_boundary_lines_to_boundary_cells function.
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    local_cells_with_boundary_lines;

  // Structure that contains the boundary cells which have a point
  std::unordered_map<
    std::string,
    std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>>
    boundary_cells_with_points;

  // Structure that contains the boundary cells with floating walls
  std::unordered_map<
    unsigned int,
    std::set<typename Triangulation<dim>::active_cell_iterator>>
    boundary_cells_for_floating_walls;

  // Structure that contains the local boundary cells with faces
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    boundary_cells_with_faces;

  // Structure that contains the global boundary cells with faces
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    global_boundary_cells_with_faces;

  // A map that contains updated points on boundaries and normal vectors of the
  // boundary faces. This container is used when the grid is moving. We use this
  // vector to update the boundary points and normal vectors in particle-wall
  // contact list.
  std::map<unsigned int, std::pair<Tensor<1, dim>, Point<dim>>>
    updated_boundary_points_and_normal_vectors;

  // Boolean variables used for printing a warning only once
  bool display_diamond_cells_warning        = true;
  bool display_pw_contact_expansion_warning = true;
};

#endif /* find_boundary_cells_information_h */
