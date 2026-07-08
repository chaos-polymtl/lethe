// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grid_impinging_jet_mixer.h>

#include <deal.II/base/point.h>

#include <deal.II/grid/cell_data.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace
{
  // --- Geometry of the mixer (SI units, metres). --------------------------
  // The geometric dimensions are parsed from the grid arguments (see the
  // constructor) and gathered into this struct, which also exposes the axial
  // coordinates derived from them.
  //
  // All coordinates here are expressed in the CONSTRUCTION datum, in which the
  // chamber floor sits at z = 0 and the chamber spans z in [0, H_chamber].  The
  // vessel is built and all of its boundary ids are assigned in this frame.  As
  // a final step, make_grid() rigidly shifts the whole mesh downward by the
  // locked inlet-axis height so that, in the delivered mesh, the common axis of
  // the two inlet pipes lies on the z = 0 (impingement) plane.
  struct MixerGeometry
  {
    double R_chamber; // mixing-chamber radius
    double r_inlet;   // inlet-pipe radius
    double r_outlet;  // outlet-pipe radius
    double H_chamber; // chamber height (axial extent)
    double L_cone;    // axial length of the conical reduction
    double L_outlet;  // length of the straight outlet pipe
    double L_inlet;   // length of each inlet pipe
    double z_inlet;   // inlet-axis height above the chamber floor

    // Derived axial coordinates along z (construction datum).
    double z_chamber_top() const { return H_chamber; }
    double z_cone_bottom() const { return -L_cone; }
    double z_outlet_bottom() const { return -L_cone - L_outlet; }
    double z_dome_center() const { return H_chamber; }
  };

  // --- Discretisation parameters. -----------------------------------------
  // The axial subdivisions below are chosen so that, after global refinement,
  // every region carries axial cells close to the ~0.01 m size of the chamber
  // wall faces (which is itself set by the circumferential resolution, so that
  // the wall faces are square).  Keeping all regions near that single target
  // size makes the volumetric mesh about as uniform as this topology allows,
  // and avoids the previous strong over-refinement of the (short) cone and
  // outlet sections.
  //   chamber : 0.16 / (2 * 8) = 0.0100
  //   cone    : 0.06 / (1 * 8) = 0.0075
  //   outlet  : 0.05 / (1 * 8) = 0.0063
  //   pipe    : 0.08 / 6       = 0.0133  (matches the ~r_inlet/3 cross-section)
  constexpr unsigned int chamber_axial_subdiv = 2;
  constexpr unsigned int cone_axial_subdiv    = 1;
  constexpr unsigned int outlet_axial_subdiv  = 1;
  // global_refinements sets the circumferential/radial (and thus the wall)
  // resolution.  It is effectively fixed at 3: it is what keeps the chamber
  // wall a 32-faceted circle, so the 3x3 inlet-port patch stays smaller than
  // the inlet diameter and maps onto the inlet circle with little chamber
  // distortion.  Lowering it makes the wall coarse and blocky at the ports.
  constexpr unsigned int global_refinements = 3; // wall resolution for ports
  // Cells along each inlet pipe.  Six layers give ~0.0133 m cells, matching the
  // in-plane pipe cell size (r_inlet mapped through the 3x3 patch) so the pipe
  // cells stay roughly cubic instead of being stretched along the axis.
  constexpr unsigned int n_pipe_layers = 6; // cells along the pipe axis

  // --- Manifold ids. ------------------------------------------------------
  constexpr types::manifold_id cylinder_z_manifold_id = 0;
  constexpr types::manifold_id dome_manifold_id       = 1;
  constexpr types::manifold_id inlet_manifold_id      = 2;

  // --- Boundary ids. ------------------------------------------------------
  constexpr types::boundary_id inlet1_boundary_id = 0; // +x
  constexpr types::boundary_id inlet2_boundary_id = 1; // -x
  constexpr types::boundary_id outlet_boundary_id = 2;
  constexpr types::boundary_id wall_boundary_id   = 3;

  // --- Material ids (internal: used only to tag the inlet-pipe cells vs the
  //     vessel cells so the correct manifold can be assigned to each). ------
  constexpr types::material_id vessel_material_id = 0;
  constexpr types::material_id pipe_material_id   = 1;
  // The central pipe column straddles the pipe axis, where the pipe's
  // cylinder-about-x manifold is singular; it is tagged separately so it can be
  // left flat while the surrounding ring cells take the cylindrical manifold.
  constexpr types::material_id pipe_core_material_id = 2;

  // Rotation mapping the native cylinder axis (+x) onto the mixer axis (+z):
  // a rotation of -pi/2 about +y sends (1,0,0) -> (0,0,1).
  void
  align_x_axis_to_z(Triangulation<3> &tria)
  {
    GridTools::rotate(Tensor<1, 3>({0., 1., 0.}), -numbers::PI_2, tria);
  }

  // Relabel every boundary face (and its lines) with manifold id `from` to id
  // `to` (used to give the dome its own spherical-manifold id before merging).
  void
  relabel_manifold(Triangulation<3>        &tria,
                   const types::manifold_id from,
                   const types::manifold_id to)
  {
    for (const auto &cell : tria.active_cell_iterators())
      for (const unsigned int f : GeometryInfo<3>::face_indices())
        if (cell->face(f)->at_boundary() &&
            cell->face(f)->manifold_id() == from)
          cell->face(f)->set_all_manifold_ids(to);
  }

  // Build the vessel (chamber + dome + reduction + outlet) as one conforming,
  // manifold-equipped triangulation.
  void
  build_vessel(Triangulation<3> &vessel, const MixerGeometry &g)
  {
    // Mixing chamber: vertical cylinder of radius R_chamber,
    // z in [0, H_chamber].
    Triangulation<3> chamber;
    GridGenerator::subdivided_cylinder(chamber,
                                       chamber_axial_subdiv,
                                       g.R_chamber,
                                       g.H_chamber / 2.0);
    align_x_axis_to_z(chamber);
    GridTools::shift(Tensor<1, 3>({0., 0., g.H_chamber / 2.0}), chamber);

    // Hemispherical dome closing the top; it shares the chamber's cross
    // section.
    Triangulation<3> dome;
    GridGenerator::half_hyper_ball(dome, Point<3>(0, 0, 0), g.R_chamber);
    align_x_axis_to_z(dome);
    GridTools::shift(Tensor<1, 3>({0., 0., g.z_dome_center()}), dome);
    relabel_manifold(dome, cylinder_z_manifold_id, dome_manifold_id);

    // Conical reduction R_chamber -> r_outlet, built as a cylinder whose radius
    // is tapered linearly along the axis so it keeps the chamber's cross
    // section and merges conformingly with both the chamber and the outlet
    // pipe.
    Triangulation<3> cone;
    GridGenerator::subdivided_cylinder(cone,
                                       cone_axial_subdiv,
                                       g.R_chamber,
                                       g.L_cone / 2.0);
    GridTools::transform(
      [&g](const Point<3> &p) {
        const double s = (p[0] + g.L_cone / 2.0) / g.L_cone; // 0 at -x, 1 at +x
        const double factor =
          g.r_outlet / g.R_chamber + s * (1.0 - g.r_outlet / g.R_chamber);
        return Point<3>(p[0], factor * p[1], factor * p[2]);
      },
      cone);
    align_x_axis_to_z(cone);
    GridTools::shift(Tensor<1, 3>({0., 0., g.z_cone_bottom() + g.L_cone / 2.0}),
                     cone);

    // Straight outlet pipe of radius r_outlet.
    Triangulation<3> outlet;
    GridGenerator::subdivided_cylinder(outlet,
                                       outlet_axial_subdiv,
                                       g.r_outlet,
                                       g.L_outlet / 2.0);
    align_x_axis_to_z(outlet);
    GridTools::shift(Tensor<1, 3>(
                       {0., 0., g.z_outlet_bottom() + g.L_outlet / 2.0}),
                     outlet);

    const double tol = 1e-6 * g.r_outlet;
    GridGenerator::merge_triangulations({&chamber, &dome, &cone, &outlet},
                                        vessel,
                                        tol,
                                        /*copy_manifold_ids=*/true,
                                        /*copy_boundary_ids=*/false);

    vessel.set_manifold(cylinder_z_manifold_id,
                        CylindricalManifold<3>(Tensor<1, 3>({0., 0., 1.}),
                                               Point<3>(0., 0., 0.)));
    vessel.set_manifold(dome_manifold_id,
                        SphericalManifold<3>(Point<3>(0., 0.,
                                                      g.z_dome_center())));
  }
} // namespace


template <int dim, int spacedim>
GridImpingingJetMixer<dim, spacedim>::GridImpingingJetMixer(
  const std::string &grid_arguments)
{
  if constexpr (!(dim == 3 && spacedim == 3))
    {
      AssertThrow(
        false,
        ExcMessage(
          "The impinging-jet mixer mesh is only supported in 3D space with 3D elements."));
      return;
    }

  this->grid_arguments = grid_arguments;

  // The geometry dimensions default to the values set in the header.  An empty
  // argument string keeps those defaults; otherwise all eight dimensions must
  // be supplied, colon-separated, in this order:
  //   R_chamber:r_inlet:r_outlet:H_chamber:L_cone:L_outlet:L_inlet:z_inlet
  // (The discretisation is intentionally not exposed and stays hardcoded.)
  if (!grid_arguments.empty())
    {
      AssertThrow(
        grid_arguments.find(';') == std::string::npos,
        ExcMessage("The impinging-jet mixer grid arguments must be separated "
                   "by colons (:), not semicolons (;)."));

      const std::vector<std::string> tokens =
        Utilities::split_string_list(grid_arguments, ':');
      AssertThrow(
        tokens.size() == 8,
        ExcMessage(
          "The impinging-jet mixer expects exactly 8 colon-separated geometry "
          "arguments, in the order "
          "R_chamber:r_inlet:r_outlet:H_chamber:L_cone:L_outlet:L_inlet:"
          "z_inlet."));

      this->R_chamber = Utilities::string_to_double(tokens[0]);
      this->r_inlet   = Utilities::string_to_double(tokens[1]);
      this->r_outlet  = Utilities::string_to_double(tokens[2]);
      this->H_chamber = Utilities::string_to_double(tokens[3]);
      this->L_cone    = Utilities::string_to_double(tokens[4]);
      this->L_outlet  = Utilities::string_to_double(tokens[5]);
      this->L_inlet   = Utilities::string_to_double(tokens[6]);
      this->z_inlet   = Utilities::string_to_double(tokens[7]);
    }

  // Sanity checks: all dimensions positive, the inlet/outlet pipes fit inside
  // the chamber cross-section, and the inlet axis lies strictly inside the
  // cylindrical part of the chamber wall (where the port patch is grown).
  AssertThrow(this->R_chamber > 0.0 && this->r_inlet > 0.0 &&
                this->r_outlet > 0.0 && this->H_chamber > 0.0 &&
                this->L_cone > 0.0 && this->L_outlet > 0.0 &&
                this->L_inlet > 0.0,
              ExcMessage("All impinging-jet mixer dimensions must be strictly "
                         "positive."));
  AssertThrow(this->r_inlet < this->R_chamber &&
                this->r_outlet < this->R_chamber,
              ExcMessage("The inlet and outlet radii must be smaller than the "
                         "chamber radius."));
  AssertThrow(this->z_inlet > 0.0 && this->z_inlet < this->H_chamber,
              ExcMessage("The inlet axis height (z_inlet) must lie strictly "
                         "between the chamber floor (0) and its top "
                         "(H_chamber)."));
}


template <>
void
GridImpingingJetMixer<3, 3>::make_grid(Triangulation<3, 3> &triangulation)
{
  // Gather the parsed dimensions into the local geometry description used
  // throughout this routine.  The derived axial coordinates that the rest of
  // the code refers to by name are pulled out as locals here so the body below
  // reads exactly as it did when these were file-scope constants.
  const MixerGeometry g{R_chamber, r_inlet,  r_outlet, H_chamber,
                        L_cone,    L_outlet, L_inlet,  z_inlet};
  const double        z_chamber_top   = g.z_chamber_top();
  const double        z_outlet_bottom = g.z_outlet_bottom();
  const double        z_dome_center   = g.z_dome_center();

  // ---------------------------------------------------------------------
  // 1. Vessel, refined so the chamber wall carries fine faces, then
  //    flattened into a plain coarse mesh whose curved geometry is baked in.
  // ---------------------------------------------------------------------
  Triangulation<3> vessel;
  build_vessel(vessel, g);
  // Rotate the chamber about its axis by half a refined wall-face width so that
  // a wall face ends up centred on the +/- x meridian.  The inlet-port patch is
  // then symmetric about +/- x and each tube stays aligned with the x-axis.
  // (The z-cylinder and dome manifolds are invariant under this rotation.)
  GridTools::rotate(Tensor<1, 3>({0., 0., 1.}),
                    numbers::PI_4 / std::pow(2.0, global_refinements),
                    vessel);
  vessel.refine_global(global_refinements);

  Triangulation<3> flat_vessel;
  GridGenerator::flatten_triangulation(vessel, flat_vessel);

  // ---------------------------------------------------------------------
  // 2. Copy the vessel cells into a (vertices, cells) description we can add
  //    the inlet-pipe cells to.
  // ---------------------------------------------------------------------
  std::vector<Point<3>>    vertices = flat_vessel.get_vertices();
  std::vector<CellData<3>> cells;
  for (const auto &cell : flat_vessel.active_cell_iterators())
    {
      CellData<3> cell_data;
      for (const unsigned int v : GeometryInfo<3>::vertex_indices())
        cell_data.vertices[v] = cell->vertex_index(v);
      cell_data.material_id = vessel_material_id;
      cells.push_back(cell_data);
    }

  // ---------------------------------------------------------------------
  // 3. For each inlet (+x, -x) collect the chamber-wall faces inside the port
  //    footprint and extrude them outward to build the pipe.  New pipe cells
  //    reuse the wall-patch vertices, so the patch faces become interior and
  //    the inlet genuinely opens into the chamber.
  // ---------------------------------------------------------------------
  // Height of the inlet axis actually used: set from the selected port patch
  // (the nearest wall-face-centre row to z_inlet) and reused for the manifold.
  // The first pipe fixes this height; the second pipe then targets the very
  // same z (rather than z_inlet) so both inlets sit on the same wall-face row
  // and face one another exactly.  This matters because z_inlet can fall
  // exactly between two candidate rows, where an independent search on each
  // side could otherwise tie-break to different rows.
  double inlet_axis_z  = z_inlet;
  bool   axis_z_locked = false;

  const auto grow_inlet = [&](const double direction /* +1 : +x, -1 : -x */) {
    // A chamber lateral-wall face on the requested (+x or -x) side.
    const auto is_wall_face = [&](const Point<3> &c) {
      const double r_xy = std::sqrt(c[0] * c[0] + c[1] * c[1]);
      return std::abs(r_xy - R_chamber) <= 0.15 * R_chamber &&
             c[2] > 1e-4 * R_chamber &&
             c[2] < z_chamber_top - 1e-4 * R_chamber && direction * c[0] > 0.0;
    };

    // 3a. Find the wall face F0 nearest the ideal port centre, then take the
    //     3x3 block of wall faces around it (F0 plus every wall face sharing a
    //     vertex with it).  A 3x3 face block is a 4x4 vertex grid whose centre
    //     is a FACE, so no vertex maps onto the pipe axis: the tube is then a
    //     clean butterfly (square core + surrounding ring) with no gaps and no
    //     axis singularity.
    // Target the height fixed by the first pipe once it is locked, so the
    // second pipe lands on the same wall-face row and both inlets are coaxial.
    const double   target_z = axis_z_locked ? inlet_axis_z : z_inlet;
    const Point<3> ideal(direction * R_chamber, 0.0, target_z);
    double         best = std::numeric_limits<double>::max();
    std::array<unsigned int, 4> f0{};
    Point<3>                    f0_center;
    bool                        found = false;
    for (const auto &cell : flat_vessel.active_cell_iterators())
      for (const unsigned int f : GeometryInfo<3>::face_indices())
        {
          const auto face = cell->face(f);
          if (!face->at_boundary() || !is_wall_face(face->center()))
            continue;
          const double d = (face->center() - ideal).norm();
          if (d < best)
            {
              best      = d;
              f0_center = face->center();
              for (unsigned int i = 0; i < 4; ++i)
                f0[i] = face->vertex_index(i);
              found = true;
            }
        }
    if (!found)
      return std::size_t(0);
    inlet_axis_z  = f0_center[2];
    axis_z_locked = true;

    const std::set<unsigned int>             f0_vertices(f0.begin(), f0.end());
    std::vector<std::array<unsigned int, 4>> port_faces;
    for (const auto &cell : flat_vessel.active_cell_iterators())
      for (const unsigned int f : GeometryInfo<3>::face_indices())
        {
          const auto face = cell->face(f);
          if (!face->at_boundary() || !is_wall_face(face->center()))
            continue;
          bool shares = false;
          for (unsigned int i = 0; i < 4; ++i)
            if (f0_vertices.count(face->vertex_index(i)))
              shares = true;
          if (!shares)
            continue;
          std::array<unsigned int, 4> vids;
          for (unsigned int i = 0; i < 4; ++i)
            vids[i] = face->vertex_index(i);
          port_faces.push_back(vids);
        }

    // 3b. Half-extents of the patch (about its centre y = 0, z = inlet_axis_z),
    //     used to map the patch onto the inlet disk of radius r_inlet.
    std::set<unsigned int> patch_vertices;
    for (const auto &face : port_faces)
      for (const unsigned int vid : face)
        patch_vertices.insert(vid);

    double a_max = 0.0; // max |y|
    double b_max = 0.0; // max |z - inlet_axis_z|
    for (const unsigned int vid : patch_vertices)
      {
        a_max = std::max(a_max, std::abs(vertices[vid][1]));
        b_max = std::max(b_max, std::abs(vertices[vid][2] - inlet_axis_z));
      }

    // 3c. One column of vertices outward per patch vertex.  Layer 0 stays on
    //     the chamber wall (so the tube is connected and the mouth is a
    //     circular hole); the square-to-disk map sends the patch boundary onto
    //     the inlet circle while the four inner vertices keep a square core ->
    //     a butterfly cross section that refines cleanly (no vertex on the
    //     axis).
    const double dL = L_inlet / n_pipe_layers;
    std::map<unsigned int, std::vector<unsigned int>> column;
    for (const unsigned int vid : patch_vertices)
      {
        const double u = (a_max > 0.0) ? vertices[vid][1] / a_max : 0.0;
        const double v =
          (b_max > 0.0) ? (vertices[vid][2] - inlet_axis_z) / b_max : 0.0;
        const double disk_y = u * std::sqrt(std::max(0.0, 1.0 - 0.5 * v * v));
        const double disk_z = v * std::sqrt(std::max(0.0, 1.0 - 0.5 * u * u));
        const double ty     = r_inlet * disk_y;
        const double tz     = inlet_axis_z + r_inlet * disk_z;

        std::vector<unsigned int> col(n_pipe_layers + 1);
        const double              x_wall =
          std::sqrt(std::max(0.0, R_chamber * R_chamber - ty * ty));
        vertices[vid] = Point<3>(direction * x_wall, ty, tz);
        col[0]        = vid;
        for (unsigned int j = 1; j <= n_pipe_layers; ++j)
          {
            vertices.emplace_back(direction * (R_chamber + j * dL), ty, tz);
            col[j] = static_cast<unsigned int>(vertices.size() - 1);
          }
        column[vid] = std::move(col);
      }

    // 3d. One hexahedron per patch face per pipe layer.  The column grown from
    //     the patch's centre face f0 straddles the pipe axis, so it is tagged
    //     with pipe_core_material_id and later kept flat; the surrounding ring
    //     columns take the pipe's cylinder-about-x manifold (see the manifold
    //     assignment below).
    for (const auto &face : port_faces)
      {
        const bool is_core_column =
          std::set<unsigned int>(face.begin(), face.end()) == f0_vertices;
        for (unsigned int j = 1; j <= n_pipe_layers; ++j)
          {
            CellData<3> hex;
            for (unsigned int i = 0; i < 4; ++i)
              {
                hex.vertices[i]     = column[face[i]][j - 1];
                hex.vertices[i + 4] = column[face[i]][j];
              }
            hex.material_id =
              is_core_column ? pipe_core_material_id : pipe_material_id;
            cells.push_back(hex);
          }
      }

    return port_faces.size();
  };

  grow_inlet(+1.0);
  grow_inlet(-1.0);

  // ---------------------------------------------------------------------
  // 4. Assemble the final triangulation.  Fix cell orientation robustly.
  // ---------------------------------------------------------------------
  GridTools::invert_cells_with_negative_measure(vertices, cells);
  GridTools::consistently_order_cells(cells);

  triangulation.create_triangulation(vertices, cells, SubCellData());

  // ---------------------------------------------------------------------
  // 5. Reassign boundary ids and manifold ids per mesh region (material id).
  //    Edges shared between the pipe and the vessel -- the port rim, i.e. the
  //    junction between each tube and the chamber wall -- are written by
  //    whichever pass runs LAST (set_all_manifold_ids touches a face and all
  //    its lines).  We deliberately run the pipe pass last so the rim edges
  //    take the pipe's cylinder-about-x manifold.  The rim vertices already sit
  //    at exactly r_inlet from the pipe axis, so this keeps the circular mouth
  //    round under refinement; leaving the rim on the chamber's z-cylinder
  //    (the other pass order) instead lets the small opening drift off-circle.
  // ---------------------------------------------------------------------
  const auto is_axial_opening = [](const auto &face, const auto &cell) {
    const Tensor<1, 3> outward = face->center() - cell->center();
    return std::abs(outward[0]) > 0.6 * outward.norm();
  };

  // pass 1: vessel cells -> dome / outlet cap / cylindrical walls.
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->material_id() == vessel_material_id)
      for (const unsigned int f : GeometryInfo<3>::face_indices())
        {
          const auto face = cell->face(f);
          if (!face->at_boundary())
            continue;

          const Point<3> c = face->center();
          if (c[2] < z_outlet_bottom + 1e-4 * R_chamber)
            {
              // Bottom of the outlet pipe.
              face->set_boundary_id(outlet_boundary_id);
              face->set_all_manifold_ids(numbers::flat_manifold_id);
            }
          else if (c[2] > z_chamber_top + 1e-4 * R_chamber)
            {
              // Hemispherical dome.
              face->set_boundary_id(wall_boundary_id);
              face->set_all_manifold_ids(dome_manifold_id);
            }
          else
            {
              // Chamber / reduction / outlet lateral walls.
              face->set_boundary_id(wall_boundary_id);
              face->set_all_manifold_ids(cylinder_z_manifold_id);
            }
        }

  // pass 1.5: give the inlet-pipe cells a volume manifold so the pipes refine
  //         cylindrically (keeping the cell shape) instead of only their outer
  //         wall faces following the circle while the interior is interpolated
  //         linearly -- the latter is what slowly degrades the pipe cells under
  //         repeated refinement.  set_all_manifold_ids() tags each cell together
  //         with its faces and lines.  The ring columns take the pipe's
  //         cylinder-about-x manifold; the central core column is then reset to
  //         flat, because that manifold is singular on the pipe axis the column
  //         straddles (and a flat core keeps its square cross-section crisp).
  //         Running the ring pass before the core pass makes the shared
  //         ring/core faces flat, so the square core stays straight-sided.
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->material_id() == pipe_material_id)
      cell->set_all_manifold_ids(inlet_manifold_id);
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->material_id() == pipe_core_material_id)
      cell->set_all_manifold_ids(numbers::flat_manifold_id);

  // pass 2 (last): inlet-pipe cells -> cylindrical about x; the axial end is
  //         the inlet opening (flat), the rest is the circular pipe wall.
  //         Running last, this also claims the shared port-rim edges for the
  //         pipe's x-cylinder, so each junction stays a clean circle when
  //         refined.
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->material_id() == pipe_material_id ||
        cell->material_id() == pipe_core_material_id)
      for (const unsigned int f : GeometryInfo<3>::face_indices())
        {
          const auto face = cell->face(f);
          if (!face->at_boundary())
            continue;

          // The pipe *wall* (circumference, at radius r_inlet) takes the
          // x-cylinder manifold; the end-cap disks stay flat.  The cap must not
          // be cylindrical: the central butterfly column's cap is centred on
          // the pipe axis, where the cylindrical manifold is singular.
          if (is_axial_opening(face, cell))
            {
              face->set_boundary_id(face->center()[0] > 0.0 ?
                                      inlet1_boundary_id :
                                      inlet2_boundary_id);
              face->set_all_manifold_ids(numbers::flat_manifold_id);
            }
          else
            {
              face->set_boundary_id(wall_boundary_id);
              face->set_all_manifold_ids(inlet_manifold_id);
            }
        }

  // The material ids were only an internal device to steer the manifold
  // assignment above; reset them so the whole domain is a single material.
  for (const auto &cell : triangulation.active_cell_iterators())
    cell->set_material_id(0);

  // Re-datum the mesh vertically so that the common axis of the two inlet pipes
  // lies on the z = 0 plane (the plane where the opposing jets impinge).  The
  // whole geometry is translated rigidly by -inlet_axis_z, so the mesh that
  // already works is untouched -- only the origin moves.  Because the shift is
  // applied to the vertices after all boundary ids have been assigned (those
  // use face centres in the original datum), the assignments are unaffected;
  // the manifolds below are then defined in this shifted frame, so the inlet
  // cylinder axis passes through z = 0 and the dome centre is translated by the
  // same amount.  The z-cylinder axis passes through the origin either way.
  const double z_shift = -inlet_axis_z;
  GridTools::shift(Tensor<1, 3>({0., 0., z_shift}), triangulation);

  triangulation.set_manifold(cylinder_z_manifold_id,
                             CylindricalManifold<3>(Tensor<1, 3>({0., 0., 1.}),
                                                    Point<3>(0., 0., 0.)));
  triangulation.set_manifold(dome_manifold_id,
                             SphericalManifold<3>(
                               Point<3>(0., 0., z_dome_center + z_shift)));
  triangulation.set_manifold(inlet_manifold_id,
                             CylindricalManifold<3>(Tensor<1, 3>({1., 0., 0.}),
                                                    Point<3>(0., 0., 0.)));
}

// Fallback make_grid definition for unsupported template parameters. This
// provides a linker-visible symbol and a clear runtime error when the
// class is instantiated for dim/spacedim combinations that are not
// specialized above.
template <int dim, int spacedim>
void
GridImpingingJetMixer<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> & /*triangulation*/)
{
  AssertThrow(
    false,
    ExcMessage(
      "GridImpingingJetMixer is only implemented for dim = 3 and spacedim = 3."));
}

// Explicit template instantiations
template class GridImpingingJetMixer<2, 2>;
template class GridImpingingJetMixer<2, 3>;
template class GridImpingingJetMixer<3, 3>;
