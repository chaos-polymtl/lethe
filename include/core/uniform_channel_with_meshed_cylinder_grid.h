// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_uniform_channel_with_meshed_cylinder_grid_h
#define lethe_uniform_channel_with_meshed_cylinder_grid_h

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <sstream>

using namespace dealii;

/**
 * @brief Class that creates a custom uniform channel with meshed cylinder geometry. The size of the cylinder, of the transition region and of the channel can be specified in the constructor. the number of subdivisions in the x and y direction can be adjusted by the padding parameters. The height of the channel and the number of subdivisions in the z direction can also be specified for the extruded mesh.
 * The geometry manifold is constructed with SphericalManifold for the inner
 * circle and TransfiniteManifold for the rest of the domain.
 */

template <int dim, int spacedim>
class UniformChannelWithMeshedCylinderGrid
{
public:
  UniformChannelWithMeshedCylinderGrid(const std::string &grid_arguments);

  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  std::string  grid_arguments;
  Point<dim>   bottom_left;
  Point<dim>   top_right;
  Point<dim>   center;
  double       inner_radius;
  double       outer_radius;
  unsigned int pad_bottom;
  unsigned int pad_top;
  unsigned int pad_left;
  unsigned int pad_right;
  double       height;
  unsigned int n_slices;
  bool         colorize;
};


/**
 * @brief Constructor for the UniformChannelWithMeshedCylinderGrid.
 *
 * @param grid_arguments. A string with 11 parameters
 * @param bottom_left : coordinates of the bottom left corner of the channel
 * @param top_right : coordinates of the top right corner of the channel
 * @param center : coordinates of the center of the cylinder
 * @param inner_radius : radius of the cylinder
 * @param outer_radius : radius of the transition region between the cylinder and the channel
 * @param pad_bottom : number of subdivisions in the bottom padding region
 * @param pad_top : number of subdivisions in the top padding region
 * @param pad_left : number of subdivisions in the left padding region
 * @param pad_right : number of subdivisions in the right padding region
 * @param height : height of the channel used for the extruded mesh
 * @param n_slices : number of subdivisions in the z direction for the extruded mesh
 * @param colorize : whether to set different boundary ids for the different sides of the channel (true/false)
 */

template <int dim, int spacedim>
UniformChannelWithMeshedCylinderGrid<dim, spacedim>::
  UniformChannelWithMeshedCylinderGrid(const std::string &grid_arguments)
{
  if constexpr (dim == 1 || spacedim == 1)
    {
      AssertThrow(
        false,
        std::runtime_error(
          "The uniform channel with meshed cylinder mesh is only supported in 2d and 3d space with 2d and 3d elements."));
    }
  else if constexpr (dim == 2 && spacedim == 3)
    {
      AssertThrow(
        false,
        std::runtime_error(
          "The uniform channel with meshed cylinder mesh is only supported in 3d space with 3d elements."));
    }

  this->grid_arguments = grid_arguments;
  const std::vector<std::string> arguments =
    Utilities::split_string_list(grid_arguments, ':');

  if (arguments.size() < 5)
    {
      AssertThrow(
        false,
        std::runtime_error(
          "Mandatory uniform channel with meshed cylinder parameters are (bottom left point : top right point : center of the cylinder : inner radius : outer radius). The points should be given as x,y and the radii should be given as a single number. The optional parameters are (pad bottom : pad top : pad left : pad right : height : n_slices : colorize). The padding parameters should be given as a single number, the height should be given as a single number, the n_slices should be given as a single integer and the colorize parameter should be given as true/false."));
    }

  // Parse bottom_left point
  std::stringstream   bottom_left_stream(arguments[0]);
  std::vector<double> bottom_left_coords;
  std::string         coord_str;
  while (getline(bottom_left_stream, coord_str, ','))
    {
      bottom_left_coords.push_back(std::stod(coord_str));
    }

  if (bottom_left_coords.size() != dim)
    {
      AssertThrow(false,
                  std::runtime_error("The bottom left point should have " +
                                     std::to_string(dim) + " coordinates."));
    }

  this->bottom_left = Point<dim>(bottom_left_coords);

  // Parse top_right point
  std::stringstream   top_right_stream(arguments[1]);
  std::vector<double> top_right_coords;
  while (getline(top_right_stream, coord_str, ','))
    {
      top_right_coords.push_back(std::stod(coord_str));
    }

  if (top_right_coords.size() != dim)
    {
      AssertThrow(false,
                  std::runtime_error("The top right point should have " +
                                     std::to_string(dim) + " coordinates."));
    }

  this->top_right = Point<dim>(top_right_coords);

  // Parse center point
  std::stringstream   center_stream(arguments[2]);
  std::vector<double> center_coords;
  while (getline(center_stream, coord_str, ','))
    {
      center_coords.push_back(std::stod(coord_str));
    }

  if (center_coords.size() != dim)
    {
      AssertThrow(false,
                  std::runtime_error("The center point should have " +
                                     std::to_string(dim) + " coordinates."));
    }

  this->center = Point<dim>(center_coords);

  // Parse inner_radius
  this->inner_radius = std::stod(arguments[3]);

  // Parse outer_radius
  this->outer_radius = std::stod(arguments[4]);

  // Check if the cylinder fits in the rectangle
  AssertThrow((center[0] - outer_radius >= bottom_left[0]) &&
                (center[0] + outer_radius <= top_right[0]) &&
                (center[1] - outer_radius >= bottom_left[1]) &&
                (center[1] + outer_radius <= top_right[1]),
              ExcMessage(
                "The cylinder outer radius does not fit in the channel."));

  // Parse optional parameters if provided. If not provided, they will be set to
  // default values.
  this->pad_bottom = (arguments.size() > 5) ? std::stoi(arguments[5]) : 0;
  this->pad_top    = (arguments.size() > 6) ? std::stoi(arguments[6]) : 0;
  this->pad_left   = (arguments.size() > 7) ? std::stoi(arguments[7]) : 0;
  this->pad_right  = (arguments.size() > 8) ? std::stoi(arguments[8]) : 0;

  // Check if the padding is required and consistent with the geometry.
  AssertThrow(!((pad_bottom == 0) &&
                (center[1] - outer_radius == bottom_left[1])),
              ExcMessage("Padding at the bottom is required (the outer radius "
                         "doesn't fill the whole channel bottom part)."));
  AssertThrow(!((pad_top == 0) && (center[1] + outer_radius == top_right[1])),
              ExcMessage("Padding at the top is required (the outer radius "
                         "doesn't fill the whole channel top part)."));
  AssertThrow(!((pad_left == 0) &&
                (center[0] - outer_radius == bottom_left[0])),
              ExcMessage("Padding at the left is required (the outer radius "
                         "doesn't fill the whole channel left part)."));
  AssertThrow(!((pad_right == 0) && (center[0] + outer_radius == top_right[0])),
              ExcMessage("Padding at the right is required (the outer "
                         "radius doesn't fill the whole channel right part)."));

  this->height   = (arguments.size() > 9) ? std::stod(arguments[9]) : 1.0;
  this->n_slices = (arguments.size() > 10) ? std::stoi(arguments[10]) : 0;
  this->colorize =
    (arguments.size() > 11 && arguments[11] == "true") ? true : false;
}

/**
 * @brief make_grid. The make_grid function generates a balanced circle at position center and with radius inner_radius. The transition is obtained with an hyper_shell mesh between inner_radius and outer_radius which we redefine the corners to create a square. The rest of the domain is meshed with a hyper_rectangle to fit the channel size. The geometry manifold is constructed with SphericalManifold for the inner circle and TransfiniteManifold for the rest of the domain.
 *
 * @param triangulation. The triangulation object on which the grid is generated
 */
template <>
void
UniformChannelWithMeshedCylinderGrid<2, 2>::make_grid(
  Triangulation<2, 2> &triangulation)
{
  const types::manifold_id tfi_manifold_id   = 0;
  const types::manifold_id polar_manifold_id = 1;

  // Create the inner cylinder
  Triangulation<2> cylinder_tria;
  GridGenerator::hyper_ball_balanced(cylinder_tria, center, inner_radius);

  // Create the outer box that will contain the cylinder. We use a hyper shell
  // because it enables to have a center point offset of (0,0). We will need to
  // adjust the vertices that lie on the diagonals to form a square.
  Triangulation<2> box_tria;
  GridGenerator::hyper_shell(box_tria, center, inner_radius, outer_radius, 8);
  for (const auto &cell : box_tria.active_cell_iterators())
    {
      for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
        {
          Point<2>    &vertex = cell->vertex(v);
          const double dist   = vertex.distance(center);
          if (/* is the vertex on the outer ring? */
              (std::abs(dist - outer_radius) < 1e-12 * outer_radius)
              /* and */
              &&
              /* is the vertex on one of the two diagonals? */
              (std::abs(std::abs(vertex[0] - center[0]) -
                        std::abs(vertex[1] - center[1]))) <
                1e-12 * outer_radius)
            {
              vertex = center + std::sqrt(2.) * (vertex - center);
            }
        }
    }

  // Create the the padding with 8 rectangles. The points used to generate the
  // rectangles follow the convention that the lower left is always number 1 and
  // the top right is number 2.
  Triangulation<2> pad_bottom_tria;
  if (pad_bottom > 0)
    {
      Point<2> pad_bottom_point_1(center[0] - outer_radius, bottom_left[1]);
      Point<2> pad_bottom_point_2(center[0] + outer_radius,
                                  center[1] - outer_radius);
      GridGenerator::subdivided_hyper_rectangle(pad_bottom_tria,
                                                {2, pad_bottom},
                                                pad_bottom_point_1,
                                                pad_bottom_point_2);
    }

  Triangulation<2> pad_top_tria;
  if (pad_top > 0)
    {
      Point<2> pad_top_point_1(center[0] - outer_radius,
                               center[1] + outer_radius);
      Point<2> pad_top_point_2(center[0] + outer_radius, top_right[1]);
      GridGenerator::subdivided_hyper_rectangle(pad_top_tria,
                                                {2, pad_top},
                                                pad_top_point_1,
                                                pad_top_point_2);
    }

  Triangulation<2> pad_left_tria;
  if (pad_left > 0)
    {
      Point<2> pad_left_point_1(bottom_left[0], center[1] - outer_radius);
      Point<2> pad_left_point_2(center[0] - outer_radius,
                                center[1] + outer_radius);
      GridGenerator::subdivided_hyper_rectangle(pad_left_tria,
                                                {pad_left, 2},
                                                pad_left_point_1,
                                                pad_left_point_2);
    }

  Triangulation<2> pad_right_tria;
  if (pad_right > 0)
    {
      Point<2> pad_right_point_1(center[0] + outer_radius,
                                 center[1] - outer_radius);
      Point<2> pad_right_point_2(top_right[0], center[1] + outer_radius);
      GridGenerator::subdivided_hyper_rectangle(pad_right_tria,
                                                {pad_right, 2},
                                                pad_right_point_1,
                                                pad_right_point_2);
    }

  Triangulation<2> pad_bottom_left_corner_tria;
  if (pad_bottom > 0 && pad_left > 0)
    {
      Point<2> pad_bottom_left_corner_point_1(bottom_left[0], bottom_left[1]);
      Point<2> pad_bottom_left_corner_point_2(center[0] - outer_radius,
                                              center[1] - outer_radius);
      GridGenerator::subdivided_hyper_rectangle(pad_bottom_left_corner_tria,
                                                {pad_left, pad_bottom},
                                                pad_bottom_left_corner_point_1,
                                                pad_bottom_left_corner_point_2);
    }

  Triangulation<2> pad_bottom_right_corner_tria;
  if (pad_bottom > 0 && pad_right > 0)
    {
      Point<2> pad_bottom_right_corner_point_1(center[0] + outer_radius,
                                               bottom_left[1]);
      Point<2> pad_bottom_right_corner_point_2(top_right[0],
                                               center[1] - outer_radius);
      ;
      GridGenerator::subdivided_hyper_rectangle(
        pad_bottom_right_corner_tria,
        {pad_right, pad_bottom},
        pad_bottom_right_corner_point_1,
        pad_bottom_right_corner_point_2);
    }

  Triangulation<2> pad_top_left_corner_tria;
  if (pad_top > 0 && pad_left > 0)
    {
      Point<2> pad_top_left_corner_point_1(bottom_left[0],
                                           center[1] + outer_radius);
      Point<2> pad_top_left_corner_point_2(center[0] - outer_radius,
                                           top_right[1]);
      GridGenerator::subdivided_hyper_rectangle(pad_top_left_corner_tria,
                                                {pad_left, pad_top},
                                                pad_top_left_corner_point_1,
                                                pad_top_left_corner_point_2);
    }

  Triangulation<2> pad_top_right_corner_tria;
  if (pad_top > 0 && pad_right > 0)
    {
      Point<2> pad_top_right_corner_point_1(center[0] + outer_radius,
                                            center[1] + outer_radius);
      Point<2> pad_top_right_corner_point_2(top_right[0], top_right[1]);
      GridGenerator::subdivided_hyper_rectangle(pad_top_right_corner_tria,
                                                {pad_right, pad_top},
                                                pad_top_right_corner_point_1,
                                                pad_top_right_corner_point_2);
    }

  GridGenerator::merge_triangulations({&cylinder_tria,
                                       &box_tria,
                                       &pad_bottom_tria,
                                       &pad_top_tria,
                                       &pad_left_tria,
                                       &pad_right_tria,
                                       &pad_bottom_left_corner_tria,
                                       &pad_bottom_right_corner_tria,
                                       &pad_top_left_corner_tria,
                                       &pad_top_right_corner_tria},
                                      triangulation);

  triangulation.reset_all_manifolds();
  triangulation.set_all_manifold_ids(0);

  // Set the polar manifold for the inner cylinder
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      const Point<2> cell_center = cell->center();
      const double   cell_dist   = cell_center.distance(center);

      // Cells that are part of the inner disk
      if (cell_dist < inner_radius)
        {
          for (const auto &face : cell->face_iterators())
            {
              // Only set manifold on faces that are at the inner_radius
              // boundary
              bool all_vertices_on_circle = true;
              for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_face;
                   ++v)
                {
                  const double vertex_dist = face->vertex(v).distance(center);
                  if (std::abs(vertex_dist - inner_radius) >
                      1e-10 * inner_radius)
                    {
                      all_vertices_on_circle = false;
                      break;
                    }
                }
              if (all_vertices_on_circle)
                {
                  face->set_all_manifold_ids(polar_manifold_id);
                }
            }
        }
    }

  PolarManifold<2, 2> polar_manifold(center);
  triangulation.set_manifold(polar_manifold_id, polar_manifold);

  // Set the TFI manifold for the outer box
  TransfiniteInterpolationManifold<2> tfi_manifold;
  tfi_manifold.initialize(triangulation);
  triangulation.set_manifold(tfi_manifold_id, tfi_manifold);

  // Set the boundary ids. we follow the same convention as for the
  // subdivided_hyper_rectangle
  if (colorize)
    {
      for (const auto &face : triangulation.active_face_iterators())
        {
          if (!face->at_boundary())
            continue;
          const Point<2> face_center = face->center();

          // Left boundary (-x direction)
          if (std::abs(face_center[0] - bottom_left[0]) <
              1e-10 * (top_right[0] - bottom_left[0]))
            {
              face->set_boundary_id(0);
              continue;
            }
          // Right boundary (+ x direction)
          if (std::abs(face_center[0] - top_right[0]) <
              1e-10 * (top_right[0] - bottom_left[0]))
            {
              face->set_boundary_id(1);
              continue;
            }
          // Bottom boundary (- y direction)
          if (std::abs(face_center[1] - bottom_left[1]) <
              1e-10 * (top_right[1] - bottom_left[1]))
            {
              face->set_boundary_id(2);
              continue;
            }
          // Top boundary (+ y direction)
          if (std::abs(face_center[1] - top_right[1]) <
              1e-10 * (top_right[1] - bottom_left[1]))
            {
              face->set_boundary_id(3);
              continue;
            }
        }
    }
}

template <>
void
UniformChannelWithMeshedCylinderGrid<3, 3>::make_grid(
  Triangulation<3, 3> &triangulation)
{
  // Create the 2D domain
  Triangulation<2> tria_2d;
  Point<2>         center_2d(center[0], center[1]);
  Point<2>         bottom_left_2d(bottom_left[0], bottom_left[1]);
  Point<2>         top_right_2d(top_right[0], top_right[1]);
  channel_with_filled_cylinder<2>(tria_2d,
                                  bottom_left_2d,
                                  top_right_2d,
                                  center_2d,
                                  inner_radius,
                                  outer_radius,
                                  pad_bottom,
                                  pad_top,
                                  pad_left,
                                  pad_right,
                                  L,
                                  n_slices,
                                  colorize);

  // Extrude in the 3rd dimension
  GridGenerator::extrude_triangulation(
    tria_2d, n_slices, L, triangulation, true);

  // Set the 3D manifolds
  const types::manifold_id tfi_manifold_id         = 0;
  const types::manifold_id cylindrical_manifold_id = 1;
  const Tensor<1, 3>       direction{{0.0, 0.0, 1.0}};

  triangulation.set_manifold(cylindrical_manifold_id, FlatManifold<3>());
  triangulation.set_manifold(tfi_manifold_id, FlatManifold<3>());
  const CylindricalManifold<3> cylindrical_manifold(direction, center);
  triangulation.set_manifold(cylindrical_manifold_id, cylindrical_manifold);

  TransfiniteInterpolationManifold<3> inner_manifold;
  inner_manifold.initialize(triangulation);
  triangulation.set_manifold(tfi_manifold_id, inner_manifold);

  // Set the boundary ids is not necessary since they are extruded from 2D and
  // give automatically the id 4 for the new bottom and 5 for the new top.
}

#endif
