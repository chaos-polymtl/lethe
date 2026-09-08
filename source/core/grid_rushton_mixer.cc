// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grid_rushton_mixer.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/cell_data.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

namespace
{
  // --- Boundary ids of the rotor half, before the shift applied at merge time.
  constexpr types::boundary_id rotor_shaft_boundary_id        = 0;
  constexpr types::boundary_id rotor_disc_boundary_id         = 1;
  constexpr types::boundary_id rotor_blade_boundary_id        = 2;
  constexpr types::boundary_id rotor_floor_boundary_id        = 3;
  constexpr types::boundary_id rotor_free_surface_boundary_id = 4;
  constexpr types::boundary_id rotor_interface_boundary_id    = 5;

  // --- Boundary ids of the stator half. -----------------------------------
  constexpr types::boundary_id stator_interface_boundary_id    = 0;
  constexpr types::boundary_id stator_wall_boundary_id         = 1;
  constexpr types::boundary_id stator_baffle_boundary_id       = 2;
  constexpr types::boundary_id stator_floor_boundary_id        = 3;
  constexpr types::boundary_id stator_free_surface_boundary_id = 4;

  // The whole triangulation carries a single cylindrical manifold about the
  // tank axis. Using one manifold id per half is also what the rotor/stator
  // merge in read_mesh_and_manifolds_for_stator_and_rotor() handles cleanly.
  constexpr types::manifold_id axis_manifold_id = 0;

  // Relative tolerance used both to merge coincident grid lines and to decide,
  // when assigning boundary ids, which geometric feature a face belongs to. It
  // is applied to the smallest spacing of the grid, which is many orders of
  // magnitude above the round-off of the vertex coordinates.
  constexpr double relative_tolerance = 1e-6;

  /**
   * @brief Sort a list of grid-line coordinates and merge the coincident ones.
   *
   * Two levels are considered coincident when they differ by less than the
   * relative tolerance times the overall extent of the list. This is what makes
   * degenerate configurations well defined, e.g. a disc as thick as the blades,
   * for which the disc and blade levels collapse onto each other.
   *
   * @param[in] levels The grid-line coordinates, in any order.
   *
   * @return The sorted coordinates, without duplicates.
   */
  std::vector<double>
  merge_levels(std::vector<double> levels)
  {
    std::ranges::sort(levels);

    const double tolerance =
      relative_tolerance * (levels.back() - levels.front());
    const auto duplicates =
      std::ranges::unique(levels, [tolerance](const double a, const double b) {
        return std::abs(a - b) <= tolerance;
      });
    levels.erase(duplicates.begin(), duplicates.end());

    return levels;
  }

  /**
   * @brief Subdivide each segment between consecutive grid lines uniformly.
   *
   * The number of cells in a segment is the one that brings the cell size
   * closest to @p cell_size, with at least one cell per segment so that no
   * geometric feature is lost.
   *
   * @param[in] levels The sorted, duplicate-free segment boundaries.
   * @param[in] cell_size The target cell size.
   *
   * @return The complete list of node coordinates, segment boundaries included.
   */
  std::vector<double>
  subdivide_levels(const std::vector<double> &levels, const double cell_size)
  {
    std::vector<double> nodes;
    nodes.push_back(levels.front());

    for (unsigned int i = 0; i + 1 < levels.size(); ++i)
      {
        const double       length = levels[i + 1] - levels[i];
        const unsigned int n_cells =
          std::max(1u,
                   static_cast<unsigned int>(std::round(length / cell_size)));

        for (unsigned int k = 1; k < n_cells; ++k)
          nodes.push_back(levels[i] + length * static_cast<double>(k) /
                                        static_cast<double>(n_cells));

        // Close the segment on the level itself rather than on the accumulated
        // sum, so that the grid lines sit exactly on the geometric features.
        // In particular the rotor and the stator then share the interface
        // radius bit for bit.
        nodes.push_back(levels[i + 1]);
      }

    return nodes;
  }

  /**
   * @brief Assign the boundary ids of one half of the mixer.
   *
   * Every boundary face of the structured polar grid is of one of three kinds,
   * which are told apart by the coordinate that does not vary over the face: a
   * face of constant height (the floor, the free surface, the flat faces of the
   * disc and of the blades), a face of constant radius (the shaft, the disc
   * rim, the blade tips, the mortar interface, the baffle roots and the tank
   * wall) or a face of constant angle (the sides of the blades and of the
   * baffles). Within each kind the feature is then identified by the value of
   * that coordinate.
   *
   * @param[in,out] triangulation The triangulation whose boundary faces are
   * labelled.
   * @param[in] geometry The dimensions of the mixer.
   * @param[in] build_rotor Whether the rotor or the stator half is labelled.
   * @param[in] tolerance The tolerance used to compare coordinates.
   */
  void
  set_boundary_ids(Triangulation<3, 3>   &triangulation,
                   const RushtonGeometry &geometry,
                   const bool             build_rotor,
                   const double           tolerance)
  {
    const auto is_at = [tolerance](const double a, const double b) {
      return std::abs(a - b) < tolerance;
    };

    for (const auto &cell : triangulation.active_cell_iterators())
      for (const auto &face : cell->face_iterators())
        {
          if (!face->at_boundary())
            continue;

          // Radial and axial extent of the face.
          double radius_min = std::numeric_limits<double>::max();
          double radius_max = std::numeric_limits<double>::lowest();
          double height_min = std::numeric_limits<double>::max();
          double height_max = std::numeric_limits<double>::lowest();
          for (const auto vertex_no : face->vertex_indices())
            {
              const Point<3> &vertex = face->vertex(vertex_no);
              const double    radius = std::hypot(vertex[0], vertex[1]);

              radius_min = std::min(radius_min, radius);
              radius_max = std::max(radius_max, radius);
              height_min = std::min(height_min, vertex[2]);
              height_max = std::max(height_max, vertex[2]);
            }
          const double radius = 0.5 * (radius_min + radius_max);
          const double height = 0.5 * (height_min + height_max);

          types::boundary_id boundary_id = numbers::invalid_boundary_id;

          if (height_max - height_min < tolerance)
            {
              // Face of constant height.
              if (is_at(height, 0.))
                boundary_id = build_rotor ? rotor_floor_boundary_id :
                                            stator_floor_boundary_id;
              else if (is_at(height, geometry.liquid_height))
                boundary_id = build_rotor ? rotor_free_surface_boundary_id :
                                            stator_free_surface_boundary_id;
              else
                {
                  // The only remaining horizontal faces are the top and
                  // bottom faces of the disc and of the blades. Since the
                  // blades reach inside the disc rim, the two are told apart by
                  // height as well as by radius: a face lying inside the rim at
                  // one of the two disc levels belongs to the disc, and every
                  // other one to a blade. Testing the radius as well keeps a
                  // disc as thick as its blades, for which the disc and blade
                  // levels coincide, splitting correctly at the rim.
                  AssertThrow(
                    build_rotor,
                    ExcMessage(
                      "Unexpected horizontal boundary face in the stator of the Rushton mixer mesh."));
                  const bool at_disc_level =
                    is_at(height, geometry.disc_bottom()) ||
                    is_at(height, geometry.disc_top());
                  boundary_id =
                    (at_disc_level && radius < geometry.disc_radius()) ?
                      rotor_disc_boundary_id :
                      rotor_blade_boundary_id;
                }
            }
          else if (radius_max - radius_min < tolerance)
            {
              // Face of constant radius.
              if (is_at(radius, geometry.interface_radius()))
                boundary_id = build_rotor ? rotor_interface_boundary_id :
                                            stator_interface_boundary_id;
              else if (build_rotor && is_at(radius, geometry.shaft_radius()))
                boundary_id = rotor_shaft_boundary_id;
              else if (build_rotor && is_at(radius, geometry.disc_radius()))
                // The disc rim, over the thickness of the disc. Above and below
                // it this radius is normally interior, since the blades reach
                // further in; it is only exposed when the blade root is placed
                // on the rim or outside it.
                boundary_id = (height > geometry.disc_bottom() &&
                               height < geometry.disc_top()) ?
                                rotor_disc_boundary_id :
                                rotor_blade_boundary_id;
              else if (build_rotor &&
                       is_at(radius, geometry.blade_root_radius()))
                boundary_id = rotor_blade_boundary_id;
              else if (build_rotor &&
                       is_at(radius, geometry.blade_tip_radius()))
                boundary_id = rotor_blade_boundary_id;
              else if (!build_rotor && is_at(radius, geometry.tank_radius()))
                boundary_id = stator_wall_boundary_id;
              else if (!build_rotor &&
                       is_at(radius, geometry.baffle_root_radius()))
                boundary_id = stator_baffle_boundary_id;
              else
                AssertThrow(
                  false,
                  ExcMessage(
                    "A cylindrical boundary face of the Rushton mixer mesh could not be matched to any geometric feature."));
            }
          else
            // Face of constant angle: the side of a blade or of a baffle.
            boundary_id =
              build_rotor ? rotor_blade_boundary_id : stator_baffle_boundary_id;

          face->set_boundary_id(boundary_id);
        }
  }

  /**
   * @brief Build one half of the Rushton mixer as a structured polar grid.
   *
   * The radial and axial node lists are built by concatenating the segments
   * delimited by the geometric features, and the cells covering a solid part
   * are simply left out of the grid.
   *
   * @param[out] triangulation The triangulation to fill.
   * @param[in] geometry The dimensions of the mixer.
   * @param[in] build_rotor Whether the rotor (inside the mortar interface) or
   * the stator (outside it) is built.
   */
  void
  build_part(Triangulation<3, 3>   &triangulation,
             const RushtonGeometry &geometry,
             const bool             build_rotor)
  {
    // Radial grid lines: the features of the half being built.
    const std::vector<double> radial_levels =
      build_rotor ? merge_levels({geometry.shaft_radius(),
                                  geometry.blade_root_radius(),
                                  geometry.disc_radius(),
                                  geometry.blade_tip_radius(),
                                  geometry.interface_radius()}) :
                    merge_levels({geometry.interface_radius(),
                                  geometry.baffle_root_radius(),
                                  geometry.tank_radius()});

    // Axial grid lines. Both halves deliberately use the same list, even though
    // only the rotor carries the disc and the blades: the mortar implementation
    // reads the number of circumferential cells from the rotor side of the
    // interface and the axial levels from the stator side, so the two must
    // discretize the height identically.
    const std::vector<double> axial_levels =
      merge_levels({0.,
                    geometry.blade_bottom(),
                    geometry.disc_bottom(),
                    geometry.disc_top(),
                    geometry.blade_top(),
                    geometry.liquid_height});

    const std::vector<double> radii =
      subdivide_levels(radial_levels, geometry.cell_size());
    const std::vector<double> heights =
      subdivide_levels(axial_levels, geometry.cell_size());

    const unsigned int n_radial      = radii.size() - 1;
    const unsigned int n_axial       = heights.size() - 1;
    const unsigned int n_angular     = geometry.n_theta;
    const unsigned int n_axial_nodes = heights.size();

    // Vertices of the (r, theta, z) lattice. The angular direction wraps
    // around, so there are as many angular nodes as angular cells.
    std::vector<Point<3>> vertices;
    vertices.reserve(radii.size() * n_angular * n_axial_nodes);
    for (const double radius : radii)
      for (unsigned int i_theta = 0; i_theta < n_angular; ++i_theta)
        {
          const double angle = 2. * numbers::PI * static_cast<double>(i_theta) /
                               static_cast<double>(n_angular);
          const double x = radius * std::cos(angle);
          const double y = radius * std::sin(angle);

          for (const double height : heights)
            vertices.emplace_back(x, y, height);
        }

    const auto vertex_index = [n_angular,
                               n_axial_nodes](const unsigned int i_radius,
                                              const unsigned int i_theta,
                                              const unsigned int i_height) {
      return (i_radius * n_angular + i_theta % n_angular) * n_axial_nodes +
             i_height;
    };

    // A blade (respectively a baffle) occupies the first n_blade_cells
    // (n_baffle_cells) angular cells of each of the n_blades (n_baffles)
    // identical angular sectors.
    const unsigned int blade_period  = n_angular / geometry.n_blades;
    const unsigned int baffle_period = n_angular / geometry.n_baffles;

    std::vector<CellData<3>> cells;
    for (unsigned int i_radius = 0; i_radius < n_radial; ++i_radius)
      {
        const double radius = 0.5 * (radii[i_radius] + radii[i_radius + 1]);

        for (unsigned int i_theta = 0; i_theta < n_angular; ++i_theta)
          {
            const bool in_blade =
              (i_theta % blade_period) < geometry.n_blade_cells;
            const bool in_baffle =
              (i_theta % baffle_period) < geometry.n_baffle_cells;

            for (unsigned int i_height = 0; i_height < n_axial; ++i_height)
              {
                const double height =
                  0.5 * (heights[i_height] + heights[i_height + 1]);

                // Cells covered by a solid part are left out of the grid.
                bool is_solid = false;
                if (build_rotor)
                  {
                    const bool in_disc_band = height > geometry.disc_bottom() &&
                                              height < geometry.disc_top();
                    const bool in_blade_band =
                      height > geometry.blade_bottom() &&
                      height < geometry.blade_top();

                    is_solid =
                      (radius < geometry.disc_radius() && in_disc_band) ||
                      (radius > geometry.blade_root_radius() &&
                       radius < geometry.blade_tip_radius() && in_blade_band &&
                       in_blade);
                  }
                else
                  is_solid =
                    radius > geometry.baffle_root_radius() && in_baffle;

                if (is_solid)
                  continue;

                CellData<3> cell;
                for (unsigned int k = 0; k < 2; ++k)
                  for (unsigned int j = 0; j < 2; ++j)
                    for (unsigned int i = 0; i < 2; ++i)
                      cell.vertices[i + 2 * j + 4 * k] =
                        vertex_index(i_radius + i, i_theta + j, i_height + k);
                cell.material_id = 0;

                cells.push_back(cell);
              }
          }
      }

    SubCellData subcelldata;
    GridTools::delete_unused_vertices(vertices, cells, subcelldata);
    triangulation.create_triangulation(vertices, cells, subcelldata);

    // Tolerance for the identification of the boundary faces, based on the
    // smallest spacing of the grid.
    double minimum_spacing = std::numeric_limits<double>::max();
    for (unsigned int i = 0; i < n_radial; ++i)
      minimum_spacing = std::min(minimum_spacing, radii[i + 1] - radii[i]);
    for (unsigned int i = 0; i < n_axial; ++i)
      minimum_spacing = std::min(minimum_spacing, heights[i + 1] - heights[i]);

    set_boundary_ids(triangulation,
                     geometry,
                     build_rotor,
                     relative_tolerance * minimum_spacing);
  }
} // namespace



template <int dim, int spacedim>
GridRushtonMixer<dim, spacedim>::GridRushtonMixer(
  const std::string &grid_type,
  const std::string &grid_arguments)
{
  if constexpr (!(dim == 3 && spacedim == 3))
    {
      AssertThrow(
        false,
        ExcMessage(
          "The Rushton mixer mesh is only supported in 3D space with 3D elements."));
      return;
    }

  if (grid_type == "rushton_rotor")
    this->generate_rotor = true;
  else if (grid_type == "rushton_stator")
    this->generate_rotor = false;
  else
    AssertThrow(
      false,
      ExcMessage(
        "Unsupported Rushton mixer grid, the choices are: rushton_rotor|rushton_stator"));

  // Separate the arguments of the string
  std::vector<std::string> arguments;
  std::stringstream        s_stream(grid_arguments);
  while (s_stream.good())
    {
      std::string substring;
      getline(s_stream, substring, ':');
      arguments.push_back(substring);
    }

  // An empty argument string keeps the default dimensions; otherwise every
  // field must be provided.
  if (!grid_arguments.empty())
    {
      AssertThrow(
        arguments.size() == 17,
        ExcMessage(
          "The Rushton mixer parameters are (tank diameter : liquid height : impeller diameter : blade root diameter : disc diameter : disc thickness : blade height : shaft diameter : impeller clearance : baffle width : interface diameter : number of blades : number of baffles : number of angular cells : number of angular cells per blade : number of angular cells per baffle : target cell size). Provide all seventeen of them, or an empty string to use the default dimensions."));

      geometry.tank_diameter       = Utilities::string_to_double(arguments[0]);
      geometry.liquid_height       = Utilities::string_to_double(arguments[1]);
      geometry.impeller_diameter   = Utilities::string_to_double(arguments[2]);
      geometry.blade_root_diameter = Utilities::string_to_double(arguments[3]);
      geometry.disc_diameter       = Utilities::string_to_double(arguments[4]);
      geometry.disc_thickness      = Utilities::string_to_double(arguments[5]);
      geometry.blade_height        = Utilities::string_to_double(arguments[6]);
      geometry.shaft_diameter      = Utilities::string_to_double(arguments[7]);
      geometry.impeller_clearance  = Utilities::string_to_double(arguments[8]);
      geometry.baffle_width        = Utilities::string_to_double(arguments[9]);
      geometry.interface_diameter  = Utilities::string_to_double(arguments[10]);
      geometry.target_cell_size    = Utilities::string_to_double(arguments[16]);

      // The counts are read as signed integers and checked before they are
      // stored, since a negative value would otherwise wrap around silently.
      const std::vector<int> counts = {Utilities::string_to_int(arguments[11]),
                                       Utilities::string_to_int(arguments[12]),
                                       Utilities::string_to_int(arguments[13]),
                                       Utilities::string_to_int(arguments[14]),
                                       Utilities::string_to_int(arguments[15])};
      AssertThrow(
        std::ranges::all_of(counts, [](const int count) { return count > 0; }),
        ExcMessage(
          "The numbers of blades, of baffles, of angular cells, of angular cells per blade and of angular cells per baffle of the Rushton mixer must all be strictly positive."));

      geometry.n_blades       = static_cast<unsigned int>(counts[0]);
      geometry.n_baffles      = static_cast<unsigned int>(counts[1]);
      geometry.n_theta        = static_cast<unsigned int>(counts[2]);
      geometry.n_blade_cells  = static_cast<unsigned int>(counts[3]);
      geometry.n_baffle_cells = static_cast<unsigned int>(counts[4]);
    }

  // Every dimension must be strictly positive before the derived radii and
  // heights below can be compared meaningfully.
  AssertThrow(
    geometry.tank_diameter > 0. && geometry.liquid_height > 0. &&
      geometry.impeller_diameter > 0. && geometry.blade_root_diameter > 0. &&
      geometry.disc_diameter > 0. && geometry.disc_thickness > 0. &&
      geometry.blade_height > 0. && geometry.shaft_diameter > 0. &&
      geometry.impeller_clearance > 0. && geometry.baffle_width > 0. &&
      geometry.interface_diameter > 0.,
    ExcMessage(
      "All the dimensions of the Rushton mixer must be strictly positive."));

  // The blades must be a non-degenerate radial band lying outside the shaft.
  // The blade root may sit inside the disc rim, which is the usual Rushton
  // configuration, on it, or outside it.
  AssertThrow(
    geometry.shaft_radius() < geometry.blade_root_radius() &&
      geometry.blade_root_radius() < geometry.blade_tip_radius(),
    ExcMessage(
      "The blades of the Rushton mixer must extend from a root strictly outside the shaft to a tip strictly beyond that root. Check the shaft, blade root and impeller diameters."));

  // The radii of the geometric features must be strictly increasing, so that
  // every part of the mixer occupies a non-degenerate region of the mesh.
  AssertThrow(
    geometry.shaft_radius() < geometry.disc_radius() &&
      geometry.disc_radius() < geometry.blade_tip_radius() &&
      geometry.blade_tip_radius() < geometry.interface_radius() &&
      geometry.interface_radius() < geometry.baffle_root_radius() &&
      geometry.baffle_root_radius() < geometry.tank_radius(),
    ExcMessage(
      "The radii of the Rushton mixer must satisfy shaft < disc < blade tip < mortar interface < baffle root < tank wall. Check the shaft, disc, impeller, interface and tank diameters and the baffle width."));

  // The impeller must be fully immersed, and the disc may not be thicker than
  // the blades it carries.
  AssertThrow(
    geometry.blade_bottom() > 0. &&
      geometry.blade_top() < geometry.liquid_height,
    ExcMessage(
      "The blades of the Rushton mixer must lie strictly between the floor and the free surface. Check the impeller clearance, the blade height and the liquid height."));
  AssertThrow(
    geometry.disc_thickness <= geometry.blade_height,
    ExcMessage(
      "The disc of the Rushton mixer may not be thicker than its blades."));

  // The blades and the baffles have to be laid out on the angular grid, which
  // requires the angular cells to be shared equally between them, and requires
  // some fluid to be left between two consecutive blades or baffles.
  AssertThrow(
    geometry.n_theta % geometry.n_blades == 0 &&
      geometry.n_theta % geometry.n_baffles == 0,
    ExcMessage(
      "The number of angular cells of the Rushton mixer must be a multiple of both the number of blades and the number of baffles."));
  AssertThrow(
    geometry.n_blade_cells > 0 &&
      geometry.n_blade_cells * geometry.n_blades < geometry.n_theta &&
      geometry.n_baffle_cells > 0 &&
      geometry.n_baffle_cells * geometry.n_baffles < geometry.n_theta,
    ExcMessage(
      "The blades and the baffles of the Rushton mixer must each span at least one angular cell, and must leave some fluid between two consecutive ones."));

  AssertThrow(
    geometry.target_cell_size >= 0.,
    ExcMessage(
      "The target cell size of the Rushton mixer cannot be negative."));
}



template <>
void
GridRushtonMixer<3, 3>::make_grid(Triangulation<3, 3> &triangulation)
{
  build_part(triangulation, this->geometry, this->generate_rotor);

  // Every cell of both halves is a polar-extruded hexahedron, for which a
  // cylindrical manifold about the tank axis is exact. Since the shaft is
  // meshed out over the full height, no cell touches the axis, where such a
  // manifold would be singular.
  triangulation.set_all_manifold_ids(axis_manifold_id);
  triangulation.set_manifold(axis_manifold_id, CylindricalManifold<3, 3>(2));
}



// Fallback make_grid definition for unsupported template parameters. This
// provides a linker-visible symbol and a clear runtime error when the
// class is instantiated for dim/spacedim combinations that are not
// specialized above.
template <int dim, int spacedim>
void
GridRushtonMixer<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> & /*triangulation*/)
{
  AssertThrow(
    false,
    ExcMessage(
      "GridRushtonMixer is only implemented for dim = 3 and spacedim = 3."));
}

// Explicit template instantiations
template class GridRushtonMixer<2, 2>;
template class GridRushtonMixer<2, 3>;
template class GridRushtonMixer<3, 3>;
