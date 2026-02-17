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
 * @brief Generates a rectangular channel mesh with a meshed cylinder obstacle.
 *
 * The mesh consists of three concentric regions around the cylinder:
 * -# An inner balanced disk of radius @p inner_radius (hyper_ball_balanced).
 * -# A transition shell between @p inner_radius and @p outer_radius, whose
 *    diagonal vertices are adjusted to form a square.
 * -# Rectangular padding regions that fill the gap between the transition
 *    square and the channel boundary defined by @p bottom_left / @p top_right.
 *
 * The number of subdivisions in each padding direction (bottom, top, left,
 * right) can be controlled independently. In 3D, the 2D cross-section is
 * extruded along the z-axis with a configurable height and number of slices.
 *
 * Manifold objects attached for proper mesh refinement:
 * - PolarManifold (2D) or CylindricalManifold (3D) on the inner cylinder
 *   boundary faces (manifold_id = 1).
 * - TransfiniteInterpolationManifold on the remaining domain (manifold_id = 0).
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
 * @brief Constructor that parses geometry parameters from a colon-separated
 * string.
 *
 * @param[in] grid_arguments A colon-separated string with the following fields
 * (coordinates are comma-separated):
 *
 * | #  | Field        | Format     | Required | Description                                    |
 * |----|--------------|------------|----------|------------------------------------------------|
 * |  0 | bottom_left  | x,y[,z]   | yes      | Bottom-left corner of the channel              |
 * |  1 | top_right    | x,y[,z]   | yes      | Top-right corner of the channel                |
 * |  2 | center       | x,y[,z]   | yes      | Center of the cylinder                         |
 * |  3 | inner_radius | double     | yes      | Radius of the cylinder                         |
 * |  4 | outer_radius | double     | yes      | Radius of the transition region                |
 * |  5 | pad_bottom   | int        | no       | Subdivisions below the transition (default: 0) |
 * |  6 | pad_top      | int        | no       | Subdivisions above the transition (default: 0) |
 * |  7 | pad_left     | int        | no       | Subdivisions left of the transition (default: 0)|
 * |  8 | pad_right    | int        | no       | Subdivisions right of the transition (default: 0)|
 * |  9 | height       | double     | no       | Extrusion height in z, 3D only (default: 1.0)  |
 * | 10 | n_slices     | int        | no       | Number of z-layers, 3D only (default: 0)       |
 * | 11 | colorize     | true/false | no       | Assign distinct boundary IDs (default: false)  |
 *
 * Example: @code "0,0 : 10,2 : 5,1 : 0.1 : 0.3 : 2 : 2 : 5 : 5 : 2. : 2 : True"
 * @endcode
 *
 * When @p colorize is true, boundary IDs follow the
 * subdivided_hyper_rectangle convention: 0 (-x), 1 (+x), 2 (-y), 3 (+y).
 * In 3D, the extruded bottom (z=0) gets id 4 and top (z=height) gets id 5.
 */

template <int dim, int spacedim>
UniformChannelWithMeshedCylinderGrid<dim, spacedim>::
  UniformChannelWithMeshedCylinderGrid(const std::string &grid_arguments)
{
  if constexpr (dim == 1 || spacedim == 1)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The uniform channel with meshed cylinder mesh is only supported in 2d and 3d space with 2d and 3d elements."));
    }
  else if constexpr (dim == 2 && spacedim == 3)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The uniform channel with meshed cylinder mesh is only supported in 3d space with 3d elements."));
    }

  this->grid_arguments = grid_arguments;
  const std::vector<std::string> arguments =
    Utilities::split_string_list(grid_arguments, ':');

  if (arguments.size() < 5)
    {
      AssertThrow(
        false,
        ExcMessage(
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
                  ExcMessage("The bottom left point should have " +
                             std::to_string(dim) + " coordinates."));
    }

  this->bottom_left =
    (dim == 2) ? Point<dim>(bottom_left_coords[0], bottom_left_coords[1]) :
                 Point<dim>(bottom_left_coords[0],
                            bottom_left_coords[1],
                            bottom_left_coords[2]);

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
                  ExcMessage("The top right point should have " +
                             std::to_string(dim) + " coordinates."));
    }

  this->top_right =
    (dim == 2) ?
      Point<dim>(top_right_coords[0], top_right_coords[1]) :
      Point<dim>(top_right_coords[0], top_right_coords[1], top_right_coords[2]);

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
                  ExcMessage("The center point should have " +
                             std::to_string(dim) + " coordinates."));
    }

  this->center =
    (dim == 2) ?
      Point<dim>(center_coords[0], center_coords[1]) :
      Point<dim>(center_coords[0], center_coords[1], center_coords[2]);

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

  // Padding is required in each direction where the transition square does not
  // reach the channel boundary (i.e., there is a gap to fill with rectangles).
  AssertThrow(!((pad_bottom == 0) && (std::abs((center[1] - outer_radius) -
                                               bottom_left[1]) > 1e-12)),
              ExcMessage("Bottom padding is required because the transition "
                         "region does not reach the bottom channel boundary."));
  AssertThrow(!((pad_top == 0) &&
                (std::abs((center[1] + outer_radius) - top_right[1]) > 1e-12)),
              ExcMessage("Top padding is required because the transition "
                         "region does not reach the top channel boundary."));
  AssertThrow(!((pad_left == 0) && (std::abs((center[0] - outer_radius) -
                                             bottom_left[0]) > 1e-12)),
              ExcMessage("Left padding is required because the transition "
                         "region does not reach the left channel boundary."));
  AssertThrow(!((pad_right == 0) &&
                (std::abs((center[0] + outer_radius) - top_right[0]) > 1e-12)),
              ExcMessage("Right padding is required because the transition "
                         "region does not reach the right channel boundary."));

  this->height   = (arguments.size() > 9) ? std::stod(arguments[9]) : 1.0;
  this->n_slices = (arguments.size() > 10) ? std::stoi(arguments[10]) : 0;
  this->colorize =
    (arguments.size() > 11 && arguments[11] == "true") ? true : false;
}


/**
 * @brief Generate the 2D channel mesh geometry with a meshed cylinder.
 *
 * This inline helper creates the composite 2D triangulation by merging:
 * - A balanced disk of radius @p inner_radius (the cylinder cross-section).
 * - A hyper_shell transition region with squared diagonal vertices between
 *   @p inner_radius and @p outer_radius.
 * - Up to 8 rectangular padding sub-meshes filling the gap between the
 *   transition square and the channel boundary.
 *
 * Manifold integer IDs are assigned (0 = TFI, 1 = polar) but no manifold
 * objects are attached — the caller must do that after this function returns.
 *
 * If @p colorize is true, boundary IDs are set following the
 * subdivided_hyper_rectangle convention: 0 (-x), 1 (+x), 2 (-y), 3 (+y).
 *
 * @param[out] triangulation  The triangulation to fill.
 * @param[in]  bottom_left    Bottom-left corner of the channel.
 * @param[in]  top_right      Top-right corner of the channel.
 * @param[in]  center         Center of the cylinder.
 * @param[in]  inner_radius   Radius of the inner cylinder.
 * @param[in]  outer_radius   Radius of the transition region.
 * @param[in]  pad_bottom     Number of subdivisions in the bottom padding.
 * @param[in]  pad_top        Number of subdivisions in the top padding.
 * @param[in]  pad_left       Number of subdivisions in the left padding.
 * @param[in]  pad_right      Number of subdivisions in the right padding.
 * @param[in]  colorize       Whether to assign distinct boundary IDs.
 */
inline void
generate_2d_channel_mesh(Triangulation<2>  &triangulation,
                         const Point<2>    &bottom_left,
                         const Point<2>    &top_right,
                         const Point<2>    &center,
                         const double       inner_radius,
                         const double       outer_radius,
                         const unsigned int pad_bottom,
                         const unsigned int pad_top,
                         const unsigned int pad_left,
                         const unsigned int pad_right,
                         const bool         colorize)
{
  const types::manifold_id tfi_manifold_id   = 0;
  const types::manifold_id polar_manifold_id = 1;

  // Create the inner cylinder (balanced disk)
  Triangulation<2> cylinder_tria;
  GridGenerator::hyper_ball_balanced(cylinder_tria, center, inner_radius);

  // Create the transition region between the cylinder and the rectangular
  // channel. A hyper_shell with 8 segments produces vertices at 45-degree
  // intervals. The diagonal vertices are moved outward by a factor of sqrt(2)
  // so that the outer boundary forms a square of half-side outer_radius.
  Triangulation<2> box_tria;
  GridGenerator::hyper_shell(box_tria, center, inner_radius, outer_radius, 8);
  for (const auto &cell : box_tria.active_cell_iterators())
    {
      for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
        {
          Point<2>    &vertex = cell->vertex(v);
          const double dist   = vertex.distance(center);
          // Check if the vertex is on the outer ring AND on a diagonal
          if ((std::abs(dist - outer_radius) < 1e-12 * outer_radius) &&
              (std::abs(std::abs(vertex[0] - center[0]) -
                        std::abs(vertex[1] - center[1]))) <
                1e-12 * outer_radius)
            {
              vertex = center + std::sqrt(2.) * (vertex - center);
            }
        }
    }

  // Create the padding with up to 8 rectangles that fill the gap between the
  // transition square and the channel boundary. The four sides and four
  // corners are meshed independently.
  Triangulation<2> pad_bottom_tria;
  if (pad_bottom > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_bottom_tria,
        {2, pad_bottom},
        Point<2>(center[0] - outer_radius, bottom_left[1]),
        Point<2>(center[0] + outer_radius, center[1] - outer_radius));
    }

  Triangulation<2> pad_top_tria;
  if (pad_top > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_top_tria,
        {2, pad_top},
        Point<2>(center[0] - outer_radius, center[1] + outer_radius),
        Point<2>(center[0] + outer_radius, top_right[1]));
    }

  Triangulation<2> pad_left_tria;
  if (pad_left > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_left_tria,
        {pad_left, 2},
        Point<2>(bottom_left[0], center[1] - outer_radius),
        Point<2>(center[0] - outer_radius, center[1] + outer_radius));
    }

  Triangulation<2> pad_right_tria;
  if (pad_right > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_right_tria,
        {pad_right, 2},
        Point<2>(center[0] + outer_radius, center[1] - outer_radius),
        Point<2>(top_right[0], center[1] + outer_radius));
    }

  Triangulation<2> pad_bottom_left_corner_tria;
  if (pad_bottom > 0 && pad_left > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_bottom_left_corner_tria,
        {pad_left, pad_bottom},
        bottom_left,
        Point<2>(center[0] - outer_radius, center[1] - outer_radius));
    }

  Triangulation<2> pad_bottom_right_corner_tria;
  if (pad_bottom > 0 && pad_right > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_bottom_right_corner_tria,
        {pad_right, pad_bottom},
        Point<2>(center[0] + outer_radius, bottom_left[1]),
        Point<2>(top_right[0], center[1] - outer_radius));
    }

  Triangulation<2> pad_top_left_corner_tria;
  if (pad_top > 0 && pad_left > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_top_left_corner_tria,
        {pad_left, pad_top},
        Point<2>(bottom_left[0], center[1] + outer_radius),
        Point<2>(center[0] - outer_radius, top_right[1]));
    }

  Triangulation<2> pad_top_right_corner_tria;
  if (pad_top > 0 && pad_right > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_top_right_corner_tria,
        {pad_right, pad_top},
        Point<2>(center[0] + outer_radius, center[1] + outer_radius),
        top_right);
    }

  // Merge all sub-triangulations into the final channel mesh
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

  // Assign manifold IDs:
  //  - id 0 (TFI) on all cells and faces (default)
  //  - id 1 (polar/cylindrical) on inner-cylinder boundary faces
  triangulation.reset_all_manifolds();
  triangulation.set_all_manifold_ids(tfi_manifold_id);

  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->center().distance(center) < inner_radius)
        {
          for (const auto &face : cell->face_iterators())
            {
              bool all_vertices_on_circle = true;
              for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_face;
                   ++v)
                {
                  if (std::abs(face->vertex(v).distance(center) -
                               inner_radius) > 1e-10 * inner_radius)
                    {
                      all_vertices_on_circle = false;
                      break;
                    }
                }
              if (all_vertices_on_circle)
                face->set_all_manifold_ids(polar_manifold_id);
            }
        }
    }

  // Assign boundary IDs following the subdivided_hyper_rectangle convention:
  //   0: left (-x),  1: right (+x),  2: bottom (-y),  3: top (+y)
  if (colorize)
    {
      const double tol_x = 1e-10 * (top_right[0] - bottom_left[0]);
      const double tol_y = 1e-10 * (top_right[1] - bottom_left[1]);

      for (const auto &face : triangulation.active_face_iterators())
        {
          if (!face->at_boundary())
            continue;

          const Point<2> face_center = face->center();

          if (std::abs(face_center[0] - bottom_left[0]) < tol_x)
            face->set_boundary_id(0);
          else if (std::abs(face_center[0] - top_right[0]) < tol_x)
            face->set_boundary_id(1);
          else if (std::abs(face_center[1] - bottom_left[1]) < tol_y)
            face->set_boundary_id(2);
          else if (std::abs(face_center[1] - top_right[1]) < tol_y)
            face->set_boundary_id(3);
        }
    }
}


/**
 * @brief Generate the 2D channel mesh with a meshed cylinder and attach
 * manifold objects for proper mesh refinement.
 *
 * Delegates geometry construction to generate_2d_channel_mesh(), then attaches:
 * - PolarManifold on the inner-cylinder boundary faces (manifold_id = 1).
 * - TransfiniteInterpolationManifold on the remaining domain (manifold_id = 0).
 *
 * @param[out] triangulation The triangulation to fill with the channel mesh.
 */
template <>
void
UniformChannelWithMeshedCylinderGrid<2, 2>::make_grid(
  Triangulation<2, 2> &triangulation)
{
  generate_2d_channel_mesh(triangulation,
                           bottom_left,
                           top_right,
                           center,
                           inner_radius,
                           outer_radius,
                           pad_bottom,
                           pad_top,
                           pad_left,
                           pad_right,
                           colorize);

  // Attach manifold objects for proper refinement behavior
  PolarManifold<2, 2> polar_manifold(center);
  triangulation.set_manifold(1, polar_manifold);

  TransfiniteInterpolationManifold<2> tfi_manifold;
  tfi_manifold.initialize(triangulation);
  triangulation.set_manifold(0, tfi_manifold);
}

/**
 * @brief Generate the 3D channel mesh by extruding the 2D cross-section and
 * attaching manifold objects for proper mesh refinement.
 *
 * The 2D cross-section is generated by generate_2d_channel_mesh() and then
 * extruded along the z-axis. Manifold IDs are inherited from the 2D mesh
 * during extrusion (copy_manifold_ids = true). 3D manifold objects are then
 * attached:
 * - CylindricalManifold on the cylinder surface (manifold_id = 1).
 * - TransfiniteInterpolationManifold on the remaining domain (manifold_id = 0).
 *
 * Lateral boundary IDs (0-3) are inherited from the 2D mesh when colorize is
 * enabled. The extruded bottom (z = 0) gets boundary_id = 4 and the top
 * (z = height) gets boundary_id = 5.
 *
 * @param[out] triangulation The triangulation to fill with the channel mesh.
 */
template <>
void
UniformChannelWithMeshedCylinderGrid<3, 3>::make_grid(
  Triangulation<3, 3> &triangulation)
{
  // Generate the 2D cross-section (geometry + manifold IDs + boundary IDs)
  Triangulation<2> tria_2d;
  generate_2d_channel_mesh(tria_2d,
                           Point<2>(bottom_left[0], bottom_left[1]),
                           Point<2>(top_right[0], top_right[1]),
                           Point<2>(center[0], center[1]),
                           inner_radius,
                           outer_radius,
                           pad_bottom,
                           pad_top,
                           pad_left,
                           pad_right,
                           colorize);

  // Extrude the 2D cross-section along the z-axis. Manifold IDs from the
  // 2D mesh are copied to the lateral faces of the 3D mesh.
  GridGenerator::extrude_triangulation(
    tria_2d, n_slices, height, triangulation, true);

  // Attach 3D manifold objects. The manifold IDs (0 for TFI, 1 for
  // cylindrical) were inherited from the 2D mesh during extrusion.
  const Tensor<1, 3>           direction{{0.0, 0.0, 1.0}};
  const CylindricalManifold<3> cylindrical_manifold(direction, center);
  triangulation.set_manifold(1, cylindrical_manifold);

  TransfiniteInterpolationManifold<3> tfi_manifold;
  tfi_manifold.initialize(triangulation);
  triangulation.set_manifold(0, tfi_manifold);
}

#endif
