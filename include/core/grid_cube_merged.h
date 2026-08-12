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
   * @param[in] grid_arguments A string with 4 parameters, 'plate file name :
   * cylinder file name : cylindrical BID plate : cylindrical BID cylinder',
   * separated by colons. The 'plate file name' and 'cylinder file name' are
   * the .msh file names of the two initial geometries, and 'cylindrical BID
   * plate', 'cylindrical BID cylinder' are the boundary IDs of the cylindrical
   * boundaries in such initial geometries (needed to properly apply the
   * manifolds)
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
  /// File name of the extruded plate with hole mesh
  std::string plate_file_name;
  /// ID corresponding to the cylindrical boundary of the extruded plate with
  /// hole mesh
  unsigned int cylindrical_bid_plate;
  /// File name of the cylinder mesh
  std::string cylinder_file_name;
  /// ID corresponding to the cylindrical boundary of the cylinder mesh
  unsigned int cylindrical_bid_cylinder;
};

#endif
