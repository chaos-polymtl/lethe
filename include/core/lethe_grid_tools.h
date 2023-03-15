/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 -  by the Lethe authors
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
 * -------------------------------------------------------------------*/

#ifndef lethe_lethegridtools_h
#define lethe_lethegridtools_h

#include <core/serial_solid.h>

#include <deal.II/base/table_handler.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_q1.h>

#include <unordered_set>

using namespace dealii;
namespace LetheGridTools
{
  /**
   * @brief
   * Map the vertex index to the cells that includes that vertex.
   *
   * @param dof_handler DofHandler of the triangulation on which to create the map
   *
   * @param vertices_cell_map Object to which the data will be written
   */

  template <int dim>
  void
  vertices_cell_mapping(
    const DoFHandler<dim> &dof_handler,
    std::map<unsigned int,
             std::set<typename DoFHandler<dim>::active_cell_iterator>>
      &vertices_cell_map);

  /**
   * @brief
   * Return the cell around a point by starting from the least refined level and
   * iterating over children.
   *
   * @param point The point around which we want to find the cell.
   *
   */
  template <int dim>
  typename DoFHandler<dim>::active_cell_iterator
  find_cell_around_point_with_tree(const DoFHandler<dim> &dof_handler,
                                   const Point<dim> &     point);

  /**
   * @brief
   * Return the cell around a point based on an initial guess of a nearby cell
   *
   * @param cell The initial cell. We suspect that the point is in one of the neighbours of this cell.
   *
   * @param point The point around which we want to find the cell.
   *
   * @param vertices_cell_map see vertices_cell_mapping function description.
   */
  template <int dim>
  typename DoFHandler<dim>::active_cell_iterator
  find_cell_around_point_with_neighbors(
    const DoFHandler<dim> &dof_handler,
    std::map<unsigned int,
             std::set<typename DoFHandler<dim>::active_cell_iterator>>
      &                                                   vertices_cell_map,
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const Point<dim> &                                    point);

  /**
   * @brief
   * Return a vector of cells around a cell including vertex neighbors
   *
   * @param cell The initial cell. We want to know all the cells that share a vertex with this cell.
   *
   * @param vertices_cell_map see vertices_cell_mapping function description.
   */
  template <int dim>
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
  find_cells_around_cell(
    std::map<unsigned int,
             std::set<typename DoFHandler<dim>::active_cell_iterator>>
      &                                                   vertices_cell_map,
    const typename DoFHandler<dim>::active_cell_iterator &cell);

  /**
   * @brief
   * Return a vector of cells around a flat cell in spacedim. A flat cell has a
   * dim of spacedim-1.
   *
   * @param cell The cell that describes the flat. We want to know all the cells that are cut by the flat.
   *
   * @param vertices_cell_map see vertices_cell_mapping function description.
   */

  template <int dim>
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
  find_cells_around_flat_cell(
    const DoFHandler<dim> &                                        dof_handler,
    const typename DoFHandler<dim - 1, dim>::active_cell_iterator &cell,
    std::map<unsigned int,
             std::set<typename DoFHandler<dim>::active_cell_iterator>>
      &vertices_cell_map);

  /**
   * @brief
   * Return a vector of cells around an edge cell in spacedim. An edge cell as a
   * dim of 1.
   *
   * @param point_1 The first point that describes the edge
   *
   * @param point_2 The second point that describes the edge
   *
   * @param vertices_cell_map see vertices_cell_mapping function description.
   */
  template <int dim>
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
  find_cells_around_edge(
    const DoFHandler<dim> &dof_handler,
    std::map<unsigned int,
             std::set<typename DoFHandler<dim>::active_cell_iterator>>
      &                                                   vertices_cell_map,
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    Point<dim> &                                          point_1,
    Point<dim> &                                          point_2);

  /**
   * @brief
   * Return a vector of cells from the dof_handler that are inside a cell from
   * another mesh.
   *
   * @param cell Describes the volume for which we want to obtain the enclosed cells.
   */
  template <int dim>
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
  find_cells_in_cells(
    const DoFHandler<dim> &                               dof_handler,
    const typename DoFHandler<dim>::active_cell_iterator &cell);

  /**
   * @brief
   * Return a bool that indicates if a cell is cut by a flat cell
   *
   * @param cell The cell that we want to check.
   *
   * @param cell_flat The cell that describes the flat.
   */
  template <int dim>
  bool
  cell_cut_by_flat(
    const typename DoFHandler<dim>::active_cell_iterator &         cell,
    const typename DoFHandler<dim - 1, dim>::active_cell_iterator &cell_flat);

  /**
   * @brief
   * Return a bool that say if a cell is cut by a edge cell
   *
   * @param cell The cell that we want to check.
   *
   * @param cell_edge The cell that describes the edge.
   */
  template <int dim>
  bool
  cell_pierced_by_edge(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const TriaIterator<CellAccessor<1, dim>> &            cell_edge);

  /**
   * @brief
   * Return a bool that say if a cell is cut by a edge cell
   *
   * @param cell The cell that we want to check.
   *
   * @param point_1 The first point that describes the edge
   *
   * @param point_2 The second point that describes the edge
   */
  template <int dim>
  bool
  cell_pierced_by_edge(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    Point<dim>                                            point_1,
    Point<dim>                                            point_2);


  /**
   * @brief
   * A copy of the project_to_d_linear from dealii that also outputs the normal
   * of the flat
   *
   * @param object The flat or the edge.
   *
   * @param trial_point The point we project on the flat or edge
   *
   */
  template <int spacedim, int structdim>
  std::pair<std::pair<Point<spacedim>, bool>, Tensor<1, spacedim>>
  project_to_d_linear_object(
    const typename DoFHandler<structdim, spacedim>::active_cell_iterator
      &                    object,
    const Point<spacedim> &trial_point);

  /**
   * @brief
   * Function returns all the boundary cells with at least one vertex in a
   * sphere.
   *
   * @param dof_handler the dof handler containing all the elements.
   *
   * @param center The center of the sphere.
   *
   * @param radius The radius of the sphere.
   *
   */
  template <int dim>
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
  find_boundary_cells_in_sphere(const DoFHandler<dim> &dof_handler,
                                const Point<dim> &     center,
                                const double           radius);


  /**
   * @brief
   * Function returns all the cells cut by a list of object defined by a mesh
   *
   * @param dof_handler the dof handler containing all the elements.
   * @param vertices_cell_map An objects that maps the vertices to cells.
   * @param list_of_objects Solid objects
   *
   *
   */
  template <int spacedim, int structdim>
  std::map<
    typename DoFHandler<spacedim>::active_cell_iterator,
    std::map<unsigned int,
             typename DoFHandler<structdim, spacedim>::active_cell_iterator>>
  find_cells_cut_by_object(
    const DoFHandler<spacedim> &dof_handler,
    std::map<unsigned int,
             std::set<typename DoFHandler<spacedim>::active_cell_iterator>>
      &                                            vertices_cell_map,
    std::vector<SerialSolid<structdim, spacedim>> &list_of_objects);

  /**
   * @brief Calculates the distance between particles and a triangle (defined using
   * three vertices)
   *
   * @return A tuple in which 0. a vector of bools to determine if the particle is
   * close to the triangle plane, 1. a vector of projected location of particles
   * on the triangle, 2. a vector of normal vectors of the triangles
   *
   * @param triangle A vector of points that defines a triangle
   * @param particles A particle_iterator_range that refers to all the particles
   * located in the background (base) cell
   * @param n_particles_in_base_cell Number of particles in the base cell
   *
   */
  template <int dim>
  std::
    tuple<std::vector<bool>, std::vector<Point<3>>, std::vector<Tensor<1, 3>>>
    find_particle_triangle_projection(
      const std::vector<Point<dim>> &                      triangle,
      const std::vector<Particles::ParticleIterator<dim>> &particles,
      const unsigned int &n_particles_in_base_cell);


  /**
   * @brief Calculates the distance between points and a triangle (defined using
   * three vertices)
   *
   * @return A distance
   *
   * @param triangle_vertices A vector of points that defines a triangle
   * @param point A point for which we want to find the distance to the triangle
   *
   */
  template <int dim>
  double
  find_point_triangle_distance(const std::vector<Point<dim>> &triangle_vertices,
                               const Point<dim> &             point);

  /**
   * @brief
   * A functor that provides a unique and uniformly distributed hash for a
   * a cell. Used to store cells in hash sets.
   */
  template <int dim>
  struct hash_cell
  {
    std::size_t
    operator()(const typename DoFHandler<dim>::active_cell_iterator &cell) const
      noexcept
    {
      return std::hash<int>()(cell->global_active_cell_index());
    }
  };

  /**
   * @brief
   * A functor that provides a way to check if two cells are the same cell.
   * Used to store cells in hash sets.
   */
  template <int dim>
  struct equal_cell
  {
    std::size_t
    operator()(const typename DoFHandler<dim>::active_cell_iterator &cell1,
               const typename DoFHandler<dim>::active_cell_iterator &cell2)
      const noexcept
    {
      return cell1->global_active_cell_index() ==
             cell2->global_active_cell_index();
    }
  };

} // namespace LetheGridTools


#endif // lethe_lethegridtools_h
