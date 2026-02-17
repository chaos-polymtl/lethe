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
 * The geometry manifold is constructed with SphericalManifold for the inner circle and TransfiniteManifold for the rest of the domain.
 */

template <int dim, int spacedim>
class UniformChannelWithMeshedCylinderGrid
{
public:
  UniformChannelWithMeshedCylinderGrid(const std::string &grid_arguments);

  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  std::string grid_arguments;
  Point<dim>  bottom_left;
  Point<dim>  top_right;
  Point<dim>  center;
  double      inner_radius;
  double      outer_radius;
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
 */

template <int dim, int spacedim>
UniformChannelWithMeshedCylinderGrid<dim, spacedim>::UniformChannelWithMeshedCylinderGrid(const std::string &grid_arguments)
{

  if constexpr (dim == 1 || spacedim == 1)
    {
      throw std::runtime_error(
        "The uniform channel with meshed cylinder mesh is only supported in 2d and 3d space with 2d and 3d elements.");
    }
  else if constexpr (dim == 2 && spacedim == 3)
    {
      throw std::runtime_error(
        "The uniform channel with meshed cylinder mesh is only supported in 3d space with 3d elements.");
    }

  this->grid_arguments = grid_arguments;
  const std::vector<std::string> arguments = Utilities::split_string_list(grid_arguments, ':');
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

}

#endif
