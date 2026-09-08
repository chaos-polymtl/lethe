// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_grid_rushton_mixer_h
#define lethe_grid_rushton_mixer_h

#include <deal.II/base/numbers.h>

#include <deal.II/grid/tria.h>

#include <string>

using namespace dealii;

/**
 * @brief Geometry of a baffled stirred tank agitated by a Rushton turbine.
 *
 * This is the single container for the dimensions of the vessel and of the
 * impeller, together with the discretisation controls. The defaults reproduce
 * the standard Rushton configuration normalised by the tank diameter
 * @f$T = 1@f$: liquid height @f$H = T@f$, impeller diameter @f$D = T/3@f$, disc
 * diameter @f$0.75 D@f$, blades spanning radially from @f$D/4@f$ to @f$D/2@f$
 * and of height @f$D/5@f$, four baffles of width @f$T/10@f$ and an off-bottom
 * clearance of @f$T/3@f$.
 *
 * The tank axis is the z axis, its floor lies at @f$z = 0@f$ and the free
 * surface at @f$z = H@f$. All dimensions are diameters or lengths; the member
 * functions below expose the radii and the axial levels derived from them, and
 * these derived values are the ones the mesh is built from.
 */
struct RushtonGeometry
{
  double tank_diameter       = 1.0;        ///< vessel inner diameter, T
  double liquid_height       = 1.0;        ///< liquid height, H
  double impeller_diameter   = 1.0 / 3.0;  ///< impeller diameter, D
  double blade_root_diameter = 1.0 / 6.0;  ///< inner edge of the blades, D/2
  double disc_diameter       = 0.25;       ///< Rushton disc diameter
  double disc_thickness      = 0.02;       ///< axial thickness of the disc
  double blade_height        = 1.0 / 15.0; ///< axial height of a blade
  double shaft_diameter      = 0.1;        ///< agitation shaft diameter
  double impeller_clearance  = 1.0 / 3.0;  ///< disc mid-plane above floor
  double baffle_width        = 0.1;        ///< radial width of a baffle
  double interface_diameter  = 0.5;        ///< mortar interface diameter

  unsigned int n_blades       = 6;  ///< number of blades on the disc
  unsigned int n_baffles      = 4;  ///< number of baffles on the tank wall
  unsigned int n_theta        = 48; ///< angular cells over the full circle
  unsigned int n_blade_cells  = 1;  ///< angular cells occupied by one blade
  unsigned int n_baffle_cells = 1;  ///< angular cells occupied by one baffle

  /// Target cell size used to subdivide the radial and axial segments. A value
  /// of zero (the default) means that it is derived from the circumferential
  /// cell size at the mortar interface, which yields a near-isotropic mesh.
  double target_cell_size = 0.0;

  /// @name Derived radii.
  /// @{
  double
  tank_radius() const
  {
    return 0.5 * tank_diameter;
  }
  double
  shaft_radius() const
  {
    return 0.5 * shaft_diameter;
  }
  double
  disc_radius() const
  {
    return 0.5 * disc_diameter;
  }
  double
  blade_root_radius() const
  {
    return 0.5 * blade_root_diameter;
  }
  double
  blade_tip_radius() const
  {
    return 0.5 * impeller_diameter;
  }
  double
  interface_radius() const
  {
    return 0.5 * interface_diameter;
  }
  double
  baffle_root_radius() const
  {
    return tank_radius() - baffle_width;
  }
  /// @}

  /// @name Derived axial levels.
  /// @{
  double
  blade_bottom() const
  {
    return impeller_clearance - 0.5 * blade_height;
  }
  double
  blade_top() const
  {
    return impeller_clearance + 0.5 * blade_height;
  }
  double
  disc_bottom() const
  {
    return impeller_clearance - 0.5 * disc_thickness;
  }
  double
  disc_top() const
  {
    return impeller_clearance + 0.5 * disc_thickness;
  }
  /// @}

  /**
   * @brief Cell size used to subdivide the radial and axial segments.
   *
   * @return The user-provided target cell size when it is strictly positive,
   * and otherwise the circumferential size of an interface cell,
   * @f$\pi d_\mathrm{interface} / n_\theta@f$.
   */
  double
  cell_size() const
  {
    return (target_cell_size > 0.0) ?
             target_cell_size :
             numbers::PI * interface_diameter / static_cast<double>(n_theta);
  }
};

/**
 * @brief Generates the rotor or the stator half of a baffled stirred tank
 * agitated by a Rushton turbine, for use with the mortar element method.
 *
 * The two halves are separated by a cylindrical mortar interface of diameter
 * @p interface_diameter, coaxial with the tank:
 * - `rushton_rotor` meshes the fluid inside that cylinder, i.e. the annular
 *   region between the shaft and the interface, from which the disc and the
 *   blades are carved out;
 * - `rushton_stator` meshes the fluid outside it, i.e. the annular region
 *   between the interface and the tank wall, from which the baffles are carved
 *   out.
 *
 * Both are generated from the *same* grid-argument string and differ only by
 * the grid type, so that they necessarily share the number of circumferential
 * cells and the axial levels at the interface. This is what the mortar
 * implementation requires: the two sides must carry the same number of
 * interface faces, those faces must subdivide the full circle uniformly, and
 * the axial levels must coincide (the circumferential count is read from the
 * rotor while the axial levels are read from the stator).
 *
 * How the geometry is built
 * -------------------------
 * Each half is one structured grid in the cylindrical index space
 * @f$(i_r, i_\theta, i_z)@f$. The radial and axial node lists are built by
 * concatenating segments so that a grid line falls exactly on every geometric
 * feature (shaft, disc rim, blade tip, interface, baffle root and tank wall;
 * floor, blade bottom, disc bottom, disc top, blade top and free surface), and
 * each solid part is then a range of cells simply left out of the grid:
 * - the shaft removes every cell inside @p shaft_diameter, over the full
 *   height, so the domain is an annulus and no cell touches the axis;
 * - the disc removes the cells between the shaft and the disc rim over the
 *   thickness of the disc;
 * - each blade removes an angular wedge between the blade root and the blade
 *   tip over the height of the blade. The blade root normally sits inside the
 *   disc rim, as on a real Rushton turbine, so that immediately above and below
 *   the disc the blade reaches further inwards than the disc itself;
 * - each baffle removes an angular wedge between the baffle root and the tank
 *   wall over the full height.
 *
 * Blades and baffles therefore occupy whole cells of the polar grid: they are
 * bounded by surfaces of constant @f$\theta@f$ and are wedges rather than
 * parallel-sided plates, so they taper slightly from root to tip. In exchange
 * every cell is a polar-extruded hexahedron, for which a single
 * CylindricalManifold about the tank axis is exact -- radial edges keep a
 * constant angle, circumferential edges a constant radius and axial edges both
 * -- so the geometry is preserved under refinement.
 *
 * Boundary IDs, before the shift applied when the two halves are merged, are:
 * - rotor: 0 shaft, 1 disc, 2 blades, 3 floor, 4 free surface, 5 mortar
 *   interface;
 * - stator: 0 mortar interface, 1 tank wall, 2 baffles, 3 floor, 4 free
 *   surface.
 *
 * Since the stator carries five distinct boundary IDs, the merge performed by
 * read_mesh_and_manifolds_for_stator_and_rotor() shifts the rotor IDs by five,
 * so that in the merged triangulation used by the solver they read: 0 stator
 * interface, 1 tank wall, 2 baffles, 3 stator floor, 4 stator free surface,
 * 5 shaft, 6 disc, 7 blades, 8 rotor floor, 9 rotor free surface and 10 rotor
 * interface.
 *
 * @tparam dim The dimension of the mesh, must be 3.
 * @tparam spacedim The dimension of the space, must be 3.
 */
template <int dim, int spacedim>
class GridRushtonMixer
{
public:
  /**
   * @brief Constructor that selects the half to generate and parses the
   * geometry and discretisation from a colon-separated string.
   *
   * @param[in] grid_type Either `rushton_rotor` or `rushton_stator`.
   *
   * @param[in] grid_arguments A colon-separated string with the following
   * fields. An empty string keeps every default; otherwise all seventeen
   * fields
   * must be given.
   *
   * | #  | Field              | Default   | Description                                     |
   * |----|--------------------|-----------|-------------------------------------------------|
   * |  0 | tank_diameter      | 1         | Vessel inner diameter, T                        |
   * |  1 | liquid_height      | 1         | Liquid height, H                                |
   * |  2 | impeller_diameter  | 1/3       | Impeller diameter D, blade tip at D/2           |
   * |  3 | blade_root_diameter| 1/6       | Inner edge of the blades, D/2                   |
   * |  4 | disc_diameter      | 0.25      | Diameter of the Rushton disc                    |
   * |  5 | disc_thickness     | 0.02      | Axial thickness of the disc                     |
   * |  6 | blade_height       | 1/15      | Axial height of a blade                         |
   * |  7 | shaft_diameter     | 0.1       | Diameter of the agitation shaft                 |
   * |  8 | impeller_clearance | 1/3       | Height of the disc mid-plane above the floor    |
   * |  9 | baffle_width       | 0.1       | Radial width of a baffle                        |
   * | 10 | interface_diameter | 0.5       | Diameter of the mortar interface                |
   * | 11 | n_blades           | 6         | Number of blades                                |
   * | 12 | n_baffles          | 4         | Number of baffles                               |
   * | 13 | n_theta            | 48        | Angular cells over the full circle              |
   * | 14 | n_blade_cells      | 1         | Angular cells occupied by one blade             |
   * | 15 | n_baffle_cells     | 1         | Angular cells occupied by one baffle            |
   * | 16 | target_cell_size   | 0         | Radial/axial cell size, 0 to derive it          |
   *
   * Example: @code "1 : 1 : 0.3333333333333333 : 0.16666666666666666 : 0.25 :
   * 0.02 : 0.06666666666666667 : 0.1 : 0.3333333333333333 : 0.1 : 0.5 : 6 : 4 :
   * 48 : 1 : 1 : 0"
   * @endcode
   */
  GridRushtonMixer(const std::string &grid_type,
                   const std::string &grid_arguments);

  /**
   * @brief Generate the requested half of the Rushton mixer and attach a
   * CylindricalManifold about the tank axis to the whole triangulation.
   *
   * @param[out] triangulation The triangulation to fill with the mesh.
   */
  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  /// Whether the rotor (inside the mortar interface) or the stator (outside it)
  /// is generated.
  bool generate_rotor;

  /// All the dimensions and discretisation controls of the mixer, with their
  /// defaults. Populated from the grid-argument string by the constructor when
  /// a non-empty one is provided.
  RushtonGeometry geometry;
};

#endif
