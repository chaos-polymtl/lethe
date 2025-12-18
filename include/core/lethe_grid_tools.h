// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_lethe_grid_tools_h
#define lethe_lethe_grid_tools_h

#include <core/serial_solid.h>

#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q_cache.h>

using namespace dealii;
namespace LetheGridTools
{
  /**
   * @brief
   * Map the vertex index to the cells that includes that vertex.
   *
   * @param dof_handler DofHandler of the triangulation on which to create the map
   *
   * @param vertices_cell_map The map container of vertex ids with its set of neighbor cells
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
   * Map the vertex index to the cells on periodic boundaries that includes the
   * periodic coinciding vertex.
   *
   * @param dof_handler DofHandler of the triangulation on which to create the map
   *
   * @param vertices_cell_map The map of vertex ids at periodic boundary to the
   * set of periodic neighbor cells
   */

  template <int dim>
  void
  vertices_cell_mapping_with_periodic_boundaries(
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
                                   const Point<dim>      &point);

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
                                                         &vertices_cell_map,
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const Point<dim>                                     &point);

  /**
   * @brief
   * Return a vector of cells around a cell. The vector of cells includes all
   * the cells that share a vertex with the initial cell, including the initial
   * cell.
   *
   * @param cell The initial cell. We want to know all the cells that share a vertex with this cell.
   * @param vertices_cell_map see vertices_cell_mapping function description.
   */
  template <int dim>
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
  find_cells_around_cell(
    std::map<unsigned int,
             std::set<typename DoFHandler<dim>::active_cell_iterator>>
                                                         &vertices_cell_map,
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
    const DoFHandler<dim>                                         &dof_handler,
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
                                                         &vertices_cell_map,
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    Point<dim>                                           &point_1,
    Point<dim>                                           &point_2);

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
    const DoFHandler<dim>                                &dof_handler,
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
    const typename DoFHandler<dim>::active_cell_iterator          &cell,
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
    const TriaIterator<CellAccessor<1, dim>>             &cell_edge);

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
                          &object,
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
                                const Point<dim>      &center,
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
                                                  &vertices_cell_map,
    std::vector<SerialSolid<structdim, spacedim>> &list_of_objects);


  /**
   *  @enum ParticleTriangleContactIndicator
   *  @brief Indicates the type of contact between a particle and a triangle.
   */
  enum class ParticleTriangleContactIndicator : std::uint8_t
  {
    face_contact,
    edge_contact,
    vertex_contact,
    no_contact
  };

  /**
   * @brief Calculates the distance between particles and a triangle (defined using
   * three vertices). The full calculation is taken from  Geometric Tools for
   * Computer Graphics, Eberly 2003 Chapter 10.3.2 - Point to triangle.
   * The entire reference is available at:
   * https://www.sciencedirect.com/science/article/pii/B9781558605947500138
   * A full implementation of the  above reference is also available here:
   * https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
   * Variables name are taken straight from this reference to ensure a better
   * readability.
   *
   * @tparam dim An integer that denotes the number of spatial dimensions.
   * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
   *
   * @param triangle A vector of points that defines a triangle
   * @param particle A particle_iterator that refers to the particle
   * located in the background (base) cell
   *
   * @return A tuple in which: 0. a bool to determine if the particle is
   * less than a radius away from the triangle plane, 1. a projected
   * location of the particle on the triangle, 2. a normal contact
   * vector between the triangle and the particle, 3. a contact
   * type between the particle and the triangle.
   */
  template <int dim, typename PropertiesIndex>
  std::tuple<bool, Point<3>, Tensor<1, 3>, ParticleTriangleContactIndicator>
  find_particle_triangle_projection(
    const std::vector<Point<dim>>          &triangle,
    const Particles::ParticleIterator<dim> &particle);

  /**
   * @brief Calculates the distance between points and a triangle (defined using
   * three vertices). The full calculation is taken from  Geometric Tools for
   * Computer Graphics, Eberly 2003 Chapter 10.3.2 - Point to triangle.
   * The entire reference is available at:
   * https://www.sciencedirect.com/science/article/pii/B9781558605947500138
   * A full implementation of the  above reference is also available here:
   * https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
   * Variables name are taken straight from this reference to ensure a better
   * readability.
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
                               const Point<dim>              &point);


  /**
   * @brief Calculates the intersection points between a line and a sphere. The
   * full calculation is taken from  Geometric Tools for Computer Graphics,
   * Eberly 2003 Chapter 11.3.2 - LINEAR COMPONENTS AND A SPHERE. The entire
   * reference is available at:
   * https://www.sciencedirect.com/science/article/pii/B978155860594750014X
   *
   * @param[in] line_start The starting point of the line
   * @param[in] line_direction The direction of the line
   * @param[in] sphere_center The center of the sphere
   * @param[in] sphere_diameter The diameter of the sphere
   * @return A vector of points that are the intersection points between the line
   * and the sphere. The vector can be empty if there is no intersection, or
   * contain one or two points if there is an intersection.
   */
  template <int dim>
  std::vector<Point<dim>>
  find_line_sphere_intersection(const Point<dim>     &line_start,
                                const Tensor<1, dim> &line_direction,
                                const Point<dim>     &sphere_center,
                                const double         &sphere_diameter);


  /**
   * @brief Calculate the minimum distance between a point and a line (defined using
   * its orign and direction). The full calculation is taken from
   * Geometric Tools for Computer Graphics, Eberly 2003 Chapter 10.2 - Point
   * to linear component. The entire reference is available at:
   * https://www.sciencedirect.com/science/article/pii/B9781558605947500138
   * A full implementation of the  above reference is also available here:
   * https://www.geometrictools.com/Documentation/DistancePointLine.pdf
   * Variables name are taken straight from this reference to ensure a better
   * readability.
   *
   * @return A distance
   *
   * @param[in] line_origin Origin of the line (any point on the line)
   * @param[in] line_direction Vector tangent to the line
   * @param[in] point A point for which we want to find the minimum distance to
   * the line
   *
   * return distance
   *
   */
  template <int dim>
  inline double
  find_point_line_distance(const Point<dim>     &line_origin,
                           const Tensor<1, dim> &line_direction,
                           const Point<dim>     &point)
  {
    Tensor<1, dim> diff = point - line_origin;

    double t_0 = line_direction * diff;
    t_0 /= line_direction.norm();

    const Point<dim> closest_point_on_line = line_origin + t_0 * line_direction;

    return point.distance(closest_point_on_line);
  }

  /**
   * @brief
   * A functor that provides a unique and uniformly distributed hash for a
   * a cell. Used to store cells in hash sets.
   */
  template <int dim>
  struct hash_cell
  {
    std::size_t
    operator()(
      const typename DoFHandler<dim>::active_cell_iterator &cell) const noexcept
    {
      return std::hash<int>()(cell->global_active_cell_index());
    }
  };

  /**
   * @brief Rotate the mapping according to the rotation angle in a rotor-stator configuration
   *
   * @param[in] dof_handler DofHandler of the triangulation
   * @param[in] mapping Current mapping
   * @param[in] radius Radius or the rotor domain which is going to be rotated
   * @param[in] rotation_angle Rotation angle of the rotor domain
   * @param[in] center_of_rotation Center of rotation of the rotor domain.
   * Default is the origin
   * @param[in] rotation_axis Rotation axis of the rotor domain for 3D case
   */
  template <int dim>
  void
  rotate_mapping(const DoFHandler<dim> &dof_handler,
                 MappingQCache<dim>    &mapping_cache,
                 const Mapping<dim>    &mapping,
                 const double          &radius,
                 const double          &rotation_angle,
                 const Point<dim>      &center_of_rotation = Point<dim>(),
                 const Tensor<1, dim>  &rotation_axis      = Tensor<1, dim>());

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

  /**
   * @brief Map the vertex index to the cells which include that vertex.
   *
   * @param tria Triangulation
   * @param vertices_cell_map The map container of vertex ids with its set of
   * neighbor cells.
   */
  template <int dim, int spacedim>
  void
  vertices_cell_mapping(
    const Triangulation<dim, spacedim> &tria,
    std::map<
      unsigned int,
      std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
      &vertices_cell_map);


  /**
   * @brief Return a vector of cells around a cell. This vector includes
   * all the cells that share a vertex with the initial cell, including the
   * initial cell.
   *
   * @param cell The initial cell. We want to know all the cells that share a
   * vertex with this cell.
   * @param vertices_cell_map Map of vertex ids to the set of neighbor cells.
   */
  template <int dim, int spacedim = dim>
  std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>
  find_cells_around_cell(
    std::map<
      unsigned int,
      std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
      &vertices_cell_map,
    const typename Triangulation<dim, spacedim>::active_cell_iterator &cell);

} // namespace LetheGridTools
#endif // lethe_lethegridtools_h
