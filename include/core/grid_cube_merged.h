// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_cube_merged_grid_h
#define lethe_cube_merged_grid_h

#include <deal.II/base/utilities.h>

#include <deal.II/grid/manifold_lib.h>

#include <sstream>

using namespace dealii;

/**
 * @brief Class that creates a custom geometry that combines an extruded plate
 * with a hole and an inner cylinder to form a cube.
 *
 * @tparam dim The dimension of the mesh, must be 3.
 * @tparam spacedim The dimension of the space where the mesh is defined, must be 3.
 */

template <int dim, int spacedim>
class GridCubeMerged
{
public:
  /**
   * @brief Constructor for the GridCubeMerged
   *
   * @param[in] grid_arguments A string with 4 parameters, `R_cylinder :
   * L_plate : H : n_subdivisions`, separated by colons. `R_cylinder` is the
   * radius of the cylinder mesh and of the hole in the extruded square plate;
   * `L_plate` is the half of the edge length of the square; `H` is the height
   * in the z direction; `n_subdivisions` is the number of subdivisions along
   * the z direction.
   */
  GridCubeMerged(const std::string &grid_arguments);

  /**
   * @brief Generate the merged cube mesh.
   *
   * @param[out] triangulation The triangulation to fill with the mesh.
   */
  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  /// Arguments used to generate the grid
  std::string grid_arguments;
  /// Radius of the cylindrical mesh
  double cylinder_radius;
  /// Half edge length of the square in the extruded plate with hole mesh
  double plate_half_length;
  /// Height in the extruded direction (z axis)
  double height;
  /// Number of subdivisions in the extruded direction (z axis)
  unsigned int n_subdivisions;
};

#endif
