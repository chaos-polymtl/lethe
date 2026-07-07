// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_grid_impinging_jet_mixer_h
#define lethe_grid_impinging_jet_mixer_h

#include <deal.II/base/utilities.h>

#include <deal.II/grid/manifold_lib.h>

#include <string>

using namespace dealii;

/**
 * @brief Generates a 3D mesh for an impinging-jet mixer geometry.
 *
 * The geometry (axis of the mixing chamber = z, vertical) consists of:
 * - A large vertical mixing chamber (cylinder of radius R_chamber) whose top
 *   is closed by a hemispherical dome;
 * - A conical reduction (R_chamber -> r_outlet) followed by a straight outlet
 *   pipe (radius r_outlet) below the chamber;
 * - Two opposing, smaller horizontal inlet pipes (radius r_inlet) entering the
 *   side wall of the chamber from +x and -x, whose jets impinge on the chamber
 *   axis. The mesh is datumed so that the common axis of the two inlet pipes
 *   lies on the z = 0 plane (the impingement plane).
 *
 * How the inlets are connected to the chamber
 * --------------------------------------------
 * A round pipe entering the curved wall of a larger cylinder is a saddle
 * intersection that cannot be produced conformingly by merging library
 * primitives (a naive merge leaves the pipes floating against the wall). The
 * inlet pipes are therefore grown out of the chamber wall itself:
 *   1. The vessel (chamber + dome + reduction + outlet) is built as one
 *      conforming triangulation and refined so the wall is finely resolved;
 *   2. at each inlet location, the small patch of chamber-wall faces that falls
 *      inside the inlet footprint is selected;
 *   3. that patch is extruded outward in +/- x to build the inlet pipe. Because
 *      the new pipe cells reuse the wall-patch vertices, the patch faces become
 *      interior faces shared between chamber and pipe -- the inlet opens into
 *      the chamber and the mesh stays watertight and conforming.
 *
 * Boundary IDs are assigned as follows:
 * - 0: inlet 1 (+x),
 * - 1: inlet 2 (-x),
 * - 2: outlet,
 * - 3: walls (everything else).
 *
 * Manifold IDs are assigned as follows:
 * - 0: cylindrical about the z-axis (chamber, reduction, outlet laterals),
 * - 1: spherical (hemispherical dome),
 * - 2: cylindrical about the x-axis (both inlet pipes).
 *
 * @tparam dim The dimension of the mesh, must be 3.
 * @tparam spacedim The dimension of the space, must be 3.
 */
template <int dim, int spacedim>
class GridImpingingJetMixer
{
public:
  /**
   * @brief Constructor that parses the geometry dimensions from a
   * colon-separated string.
   *
   * The eight geometric dimensions (all in metres) are, in order:
   *   R_chamber : r_inlet : r_outlet : H_chamber : L_cone : L_outlet :
   *   L_inlet : z_inlet
   * An empty string keeps the default dimensions. The mesh discretisation is
   * not exposed and stays hardcoded.
   *
   * @param[in] grid_arguments A colon-separated string with the geometry
   * dimensions, or an empty string to use the defaults.
   */
  GridImpingingJetMixer(const std::string &grid_arguments);

  /**
   * @brief Generate the impinging-jet mixer mesh.
   *
   * @param[out] triangulation The triangulation to fill with the mesh.
   */
  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  /// Arguments used to generate the grid
  std::string grid_arguments;

  /// Geometry dimensions (SI units, metres). Defaults below are overridden by
  /// the constructor when a non-empty argument string is provided.
  double R_chamber = 0.05;  ///< mixing-chamber radius
  double r_inlet   = 0.02;  ///< inlet-pipe radius
  double r_outlet  = 0.025; ///< outlet-pipe radius
  double H_chamber = 0.16;  ///< chamber height (axial extent)
  double L_cone    = 0.06;  ///< axial length of the conical reduction
  double L_outlet  = 0.05;  ///< length of the straight outlet pipe
  double L_inlet   = 0.08;  ///< length of each inlet pipe
  double z_inlet   = 0.08;  ///< inlet-axis height above the chamber floor
};

#endif
