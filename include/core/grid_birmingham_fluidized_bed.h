// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_grid_birmingham_fluidized_bed_h
#define lethe_grid_birmingham_fluidized_bed_h

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <array>
#include <sstream>

using namespace dealii;

/**
 * @brief Generates a 3D mesh for the Birmingham fluidized bed geometry.
 *
 * The geometry is taken from:
 * P. Fede, O. Simonin, and A. Ingram, "3D numerical simulation of a
 * lab-scale pressurized dense fluidized bed focussing on the effect of
 * the particle-particle restitution coefficient and particle-wall boundary
 * conditions," Chemical Engineering Science, vol. 142, pp. 215-235, 2016.
 * DOI: 10.1016/j.ces.2015.11.016
 *
 * The geometry consists of three or four sections joined along the x-axis:
 * - A bottom cylinder of small radius (0.077 m),
 * - A truncated cone expanding from the small to a larger radius (0.127 m),
 * - A top cylinder of the larger radius,
 * - (optional) A rectangular chimney pipe extending from an off-center cell
 *   of the top end cap, serving as the outlet.
 *
 * The chimney is controlled by the @p enable_chimney constructor argument.
 *
 * Boundary IDs are assigned as follows:
 * - 0: wall surfaces (curved cylinder/cone laterals; when the chimney is
 *   enabled this also includes the chimney walls and the annular top wall
 *   around the chimney opening),
 * - 1: inlet end cap (x = 0),
 * - 2: outlet end cap (x = total length), or chimney top face when the
 *   chimney is enabled.
 *
 * A CylindricalManifold is attached to the curved lateral surfaces so that
 * mesh refinement preserves the circular cross-sections and the cone taper.
 * When the chimney is enabled, its walls and the annular top wall remain
 * flat.
 *
 * @tparam dim The dimension of the mesh, must be 3.
 * @tparam spacedim The dimension of the space, must be 3.
 */
template <int dim, int spacedim>
class BirminghamFluidizedBedGrid
{
public:
  /**
   * @brief Constructor that parses the chimney flag from a colon-separated
   * string.
   *
   * @param[in] grid_arguments A string with the following format:
   * @code "enable_chimney" @endcode
   *
   * | #  | Field          | Format     | Required | Description                          |
   * |----|----------------|------------|----------|--------------------------------------|
   * |  0 | enable_chimney | true/false | no       | Add chimney outlet (default: true)   |
   *
   * Example: @code "true" @endcode
   */
  BirminghamFluidizedBedGrid(const std::string &grid_arguments);

  /**
   * @brief Generate the Birmingham fluidized bed mesh.
   *
   * The cylindrical geometry (bottom cylinder, truncated cone, top cylinder)
   * is first merged and refined with a CylindricalManifold to obtain a
   * finer mesh on the top end cap. The refined mesh is then flattened and
   * an off-center face on the top end cap is selected for the chimney
   * attachment. A hexahedral chimney pipe is built from that face's exact
   * vertices and merged with the flattened mesh. The number of pre-chimney
   * refinements is controlled by the hardcoded parameter
   * @p n_chimney_refinements.
   *
   * @param[out] triangulation The triangulation to fill with the mesh.
   */
  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  std::string grid_arguments;
  bool        enable_chimney;
};


#endif
