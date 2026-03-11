// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_grid_birmingham_fluidized_bed_h
#define lethe_grid_birmingham_fluidized_bed_h

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <sstream>

using namespace dealii;

/**
 * @brief Generates a 3D mesh for the Birmingham fluidized bed geometry.
 *
 * The geometry consists of three sections joined along the x-axis:
 * - A bottom cylinder of small radius (0.077 m),
 * - A truncated cone expanding from the small to a larger radius (0.127 m),
 * - A top cylinder of the larger radius.
 *
 * Boundary IDs are assigned as follows:
 * - 0: lateral (curved) surfaces,
 * - 1: inlet end cap (x = 0),
 * - 2: outlet end cap (x = total length).
 *
 * A CylindricalManifold is attached to the lateral surfaces so that mesh
 * refinement preserves the circular cross-sections and the cone taper.
 *
 * @tparam dim The dimension of the mesh, must be 3.
 * @tparam spacedim The dimension of the space, must be 3.
 */
template <int dim, int spacedim>
class BirminghamFluidizedBedGrid
{
public:
  /**
   * @brief Constructor for BirminghamFluidizedBedGrid.
   *
   * @param[in] grid_arguments A colon-separated string of optional
   * arguments. Currently unused because the geometry is fully defined by
   * fixed physical dimensions, but accepted for interface consistency with
   * other Lethe grid classes.
   */
  BirminghamFluidizedBedGrid(const std::string &grid_arguments);

  /**
   * @brief Generate the Birmingham fluidized bed mesh.
   *
   * Three sub-triangulations (bottom cylinder, truncated cone, top cylinder)
   * are created and merged. Manifold IDs and boundary IDs are then reassigned
   * so that the curved surfaces receive a CylindricalManifold along the
   * x-axis and the flat end caps are identified as inlet (boundary 1) and
   * outlet (boundary 2).
   *
   * @param[out] triangulation The triangulation to fill with the mesh.
   */
  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  std::string grid_arguments;
};


#endif
