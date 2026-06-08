// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/config.h>

#include <core/mortar_coupling_manager.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.templates.h>

#include <deal.II/fe/fe_nothing.h>

#include <algorithm>

DeclException1(
  InterfaceRadiusTolerance,
  double,
  << "The computed radius of the rotor mesh has a variation of " << arg1
  << " along the interface boundary, which exceeds the prescribed tolerance.");

/*-------------- MortarManagerBase -------------------------------*/

template <int dim>
bool
MortarManagerBase<dim>::is_mesh_aligned() const
{
  AssertThrow(dim != 1, ExcInternalError());

  // The meshes are aligned when, on every axial stage, merging the rotor and
  // stator breakpoints does not introduce any additional segments, i.e. every
  // rotor breakpoint coincides with a stator breakpoint.
  for (unsigned int s = 0; s < n_stages(); ++s)
    if (merged_breakpoints[s].size() != rotor_breakpoints[s].size() ||
        merged_breakpoints[s].size() != stator_breakpoints[s].size())
      return false;

  return true;
}

template <int dim>
unsigned int
MortarManagerBase<dim>::get_n_total_mortars() const
{
  if constexpr (dim == 1)
    return 1;

  // Total number of mortar segments summed over all axial stages. The per-stage
  // prefix sum already accounts for every stage.
  return segment_offset.back();
}

template <int dim>
unsigned int
MortarManagerBase<dim>::get_n_mortars(const Point<dim> &face_center,
                                      const bool        is_inner) const
{
  if constexpr (dim == 1)
    return 1;

  const auto [id_out_plane, start_index, arc_breaks] =
    get_face_arrangement(face_center, is_inner);
  (void)id_out_plane;
  (void)start_index;

  return static_cast<unsigned int>(arc_breaks.size() - 1);
}

template <int dim>
std::vector<unsigned int>
MortarManagerBase<dim>::get_mortar_indices(const Point<dim> &face_center,
                                           const bool        is_inner) const
{
  if constexpr (dim == 1)
    return std::vector<unsigned int>{0};

  const auto [id_out_plane, start_index, arc_breaks] =
    get_face_arrangement(face_center, is_inner);

  const unsigned int n_seg = static_cast<unsigned int>(arc_breaks.size() - 1);
  const unsigned int n_seg_stage =
    static_cast<unsigned int>(merged_breakpoints[id_out_plane].size());

  std::vector<unsigned int> indices;
  indices.reserve(n_seg);

  for (unsigned int i = 0; i < n_seg; ++i)
    {
      const unsigned int local = (start_index + i) % n_seg_stage;

      AssertIndexRange(local, n_seg_stage);

      indices.emplace_back(segment_offset[id_out_plane] + local);
    }

  return indices;
}

template <int dim>
unsigned int
MortarManagerBase<dim>::get_n_total_points() const
{
  if constexpr (dim == 1)
    return 1;

  return get_n_total_mortars() * n_quadrature_points;
}

template <int dim>
unsigned int
MortarManagerBase<dim>::get_n_points(const Point<dim> &face_center,
                                     const bool        is_inner) const
{
  if constexpr (dim == 1)
    return 1;

  return get_n_mortars(face_center, is_inner) * n_quadrature_points;
}

template <int dim>
std::vector<Point<dim>>
MortarManagerBase<dim>::get_points(const Point<dim> &face_center,
                                   const bool        is_inner) const
{
  if constexpr (dim == 1)
    return std::vector<Point<dim>>{face_center};

  const auto [id_out_plane, start_index, arc_breaks] =
    get_face_arrangement(face_center, is_inner);

  // Height of the cell in the direction of the rotation axis
  double delta_1 = 1.0;
  if constexpr (dim == 3)
    delta_1 = stage_heights[id_out_plane + 1] - stage_heights[id_out_plane];

  const unsigned int n_seg = static_cast<unsigned int>(arc_breaks.size() - 1);

  std::vector<Point<dim>> points;
  points.reserve(n_seg * n_quadrature_points);

  // Loop over the mortar segments covered by this face, emitting the
  // quadrature points of each segment in increasing-angle order.
  for (unsigned int s = 0; s < n_seg; ++s)
    {
      const double s_lo = arc_breaks[s];
      const double s_hi = arc_breaks[s + 1];

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        {
          const auto x = from_1D(s_lo + quadrature.point(q)[0] * (s_hi - s_lo));

          if constexpr (dim == 3)
            points.emplace_back(
              x[0],
              x[1],
              stage_heights[id_out_plane] +
                quadrature.point(q)[1] *
                  delta_1); // TODO Generalize for x and y directions
          else
            points.emplace_back(x);
        }
    }

  return points;
}

template <int dim>
std::vector<Point<std::max(1, dim - 1)>>
MortarManagerBase<dim>::get_points_ref(const Point<dim> &face_center,
                                       const bool        is_inner) const
{
  if (dim == 1)
    return std::vector<Point<std::max(1, dim - 1)>>{
      Point<std::max(1, dim - 1)>()};

  const auto [id_out_plane, start_index, arc_breaks] =
    get_face_arrangement(face_center, is_inner);
  (void)id_out_plane;
  (void)start_index;

  const unsigned int n_seg = static_cast<unsigned int>(arc_breaks.size() - 1);

  // The reference coordinate is expressed within the face's own arc, which
  // spans [arc_breaks.front(), arc_breaks.back()].
  const double face_lo    = arc_breaks.front();
  const double face_width = arc_breaks.back() - face_lo;

  std::vector<Point<std::max(1, dim - 1)>> points;
  points.reserve(n_seg * n_quadrature_points);

  for (unsigned int s = 0; s < n_seg; ++s)
    {
      const double s_lo = arc_breaks[s];
      const double s_hi = arc_breaks[s + 1];

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        {
          const double angle = s_lo + quadrature.point(q)[0] * (s_hi - s_lo);
          const double x      = (angle - face_lo) / face_width;

          if (dim == 2)
            points.emplace_back(x);
          else
            points.emplace_back(x, quadrature.point(q)[1]);
        }
    }

  return points;
}

template <int dim>
std::vector<double>
MortarManagerBase<dim>::get_weights(const Point<dim> &face_center,
                                    const bool        is_inner) const
{
  if (dim == 1)
    return std::vector<double>{1.0};

  const auto [id_out_plane, start_index, arc_breaks] =
    get_face_arrangement(face_center, is_inner);
  (void)start_index;

  double delta_1 = 1.0;
  if (dim == 3)
    delta_1 = stage_heights[id_out_plane + 1] - stage_heights[id_out_plane];

  const unsigned int n_seg = static_cast<unsigned int>(arc_breaks.size() - 1);

  std::vector<double> weights;
  weights.reserve(n_seg * n_quadrature_points);

  for (unsigned int s = 0; s < n_seg; ++s)
    {
      const double seg_width = arc_breaks[s + 1] - arc_breaks[s];

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        weights.emplace_back(radius[0] * quadrature.weight(q) * seg_width *
                             delta_1);
    }

  return weights;
}

template <int dim>
std::vector<Tensor<1, dim, double>>
MortarManagerBase<dim>::get_normals(const Point<dim> &face_center,
                                    const bool        is_inner) const
{
  // Coordinates of cell quadrature points
  const auto points = get_points(face_center, is_inner);

  std::vector<Tensor<1, dim, double>> result;

  result.reserve(points.size());

  for (const auto &point : points)
    result.emplace_back(get_normal(point));

  return result;
}

template <int dim>
void
MortarManagerBase<dim>::build_arrangement(
  const std::vector<std::vector<double>> &rotor_bp_per_stage,
  const std::vector<std::vector<double>> &stator_bp_per_stage)
{
  if constexpr (dim == 1)
    return;

  AssertDimension(rotor_bp_per_stage.size(), stator_bp_per_stage.size());

  const double two_pi = 2 * numbers::PI;
  const double tol    = breakpoint_tolerance;

  // Normalize angles into [0, 2*pi), sort, and remove duplicates (including the
  // periodic duplicate of a breakpoint near 0 with one near 2*pi).
  const auto normalize_sort_unique = [&](std::vector<double> &v) {
    for (auto &a : v)
      {
        a = std::fmod(a, two_pi);
        if (a < 0.0)
          a += two_pi;
      }
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(),
                        v.end(),
                        [&](double a, double b) {
                          return std::abs(a - b) < tol;
                        }),
            v.end());
    if (v.size() > 1 && std::abs((v.front() + two_pi) - v.back()) < tol)
      v.pop_back();
  };

  const unsigned int n_st = static_cast<unsigned int>(rotor_bp_per_stage.size());

  rotor_breakpoints.assign(n_st, {});
  stator_breakpoints.assign(n_st, {});
  merged_breakpoints.assign(n_st, {});
  segment_offset.assign(n_st + 1, 0);

  // Build one independent in-plane arrangement per axial stage, then accumulate
  // the per-stage segment counts into a prefix sum so that the global index of
  // local segment k in stage s is segment_offset[s] + k.
  for (unsigned int s = 0; s < n_st; ++s)
    {
      std::vector<double> rotor_bp  = rotor_bp_per_stage[s];
      std::vector<double> stator_bp = stator_bp_per_stage[s];

      normalize_sort_unique(rotor_bp);
      normalize_sort_unique(stator_bp);

      rotor_breakpoints[s]  = rotor_bp;
      stator_breakpoints[s] = stator_bp;

      // The merged, sorted, deduplicated union of both breakpoint sets defines
      // the mortar segments of this stage. Segment k spans
      // [merged[k], merged[(k + 1) % size]).
      std::vector<double> merged;
      merged.reserve(rotor_bp.size() + stator_bp.size());
      merged.insert(merged.end(), rotor_bp.begin(), rotor_bp.end());
      merged.insert(merged.end(), stator_bp.begin(), stator_bp.end());
      normalize_sort_unique(merged);

      merged_breakpoints[s] = merged;

      segment_offset[s + 1] =
        segment_offset[s] + static_cast<unsigned int>(merged.size());
    }
}

template <int dim>
void
MortarManagerBase<dim>::build_arrangement(std::vector<double> rotor_bp,
                                          std::vector<double> stator_bp)
{
  // Wrap the flat (single-stage) breakpoints used by 2D problems and unit-test
  // subclasses into the per-stage representation.
  build_arrangement(std::vector<std::vector<double>>{std::move(rotor_bp)},
                    std::vector<std::vector<double>>{std::move(stator_bp)});
}

template <int dim>
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
MortarManagerBase<dim>::compute_breakpoints_from_mesh(
  const Triangulation<dim>      &triangulation,
  const Mapping<dim>            &mapping,
  const Parameters::Mortar<dim> &mortar_parameters) const
{
  // One bucket of breakpoint angles per axial stage (a single bucket in 2D).
  const unsigned int n_st = (dim == 3) ? stage_heights.size() - 1 : 1;

  std::vector<std::vector<double>> rotor_local(n_st), stator_local(n_st);

  // Collect the in-plane angle of every interface vertex, using the manager's
  // own to_1D() transform (so this works for both circular and linear
  // interfaces, and accounts for the rotor rotation already baked into the
  // rotated mapping). Each interface face lies entirely within one axial stage
  // band, so it contributes its vertex angles to that stage's bucket only.
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto face_no : cell->face_indices())
        {
          const auto face = cell->face(face_no);

          if (!face->at_boundary())
            continue;

          const bool is_rotor =
            face->boundary_id() == mortar_parameters.rotor_boundary_id;
          const bool is_stator =
            face->boundary_id() == mortar_parameters.stator_boundary_id;

          if (!is_rotor && !is_stator)
            continue;

          const auto vertices = mapping.get_vertices(cell, face_no);

          // Determine the axial stage of this face from its center height. The
          // rotation axis is aligned with z (component 2) in 3D.
          unsigned int stage = 0;
          if constexpr (dim == 3)
            {
              double z_sum = 0.0;
              for (unsigned int v = 0; v < face->n_vertices(); ++v)
                z_sum += vertices[v][2];
              const double face_z = z_sum / face->n_vertices();

              const auto it =
                std::upper_bound(stage_heights.begin(),
                                 stage_heights.end(),
                                 face_z);
              stage = static_cast<unsigned int>(
                std::min<long>(std::distance(stage_heights.begin(), it) - 1,
                               static_cast<long>(n_st) - 1));
            }

          for (unsigned int v = 0; v < face->n_vertices(); ++v)
            (is_rotor ? rotor_local : stator_local)[stage].push_back(
              this->to_1D(vertices[v]));
        }

  // Gather the breakpoints of each stage from all processes so that every rank
  // holds the full global arrangement (mirrors compute_stage_heights).
  const auto gather = [&](const std::vector<double> &local) {
    std::vector<double> global;
    const auto          all =
      Utilities::MPI::all_gather(triangulation.get_mpi_communicator(), local);
    for (const auto &per_process : all)
      global.insert(global.end(), per_process.begin(), per_process.end());
    return global;
  };

  std::vector<std::vector<double>> rotor_global(n_st), stator_global(n_st);
  for (unsigned int s = 0; s < n_st; ++s)
    {
      rotor_global[s]  = gather(rotor_local[s]);
      stator_global[s] = gather(stator_local[s]);
    }

  return {rotor_global, stator_global};
}

template <int dim>
std::tuple<unsigned int, unsigned int, std::vector<double>>
MortarManagerBase<dim>::get_face_arrangement(const Point<dim> &face_center,
                                             const bool        is_inner) const
{
  const double two_pi = 2 * numbers::PI;
  const double tol    = breakpoint_tolerance;

  // Axial stage index (always 0 in 2D)
  unsigned int id_out_plane = 0;
  if constexpr (dim == 3)
    {
      auto upper_height_iterator =
        std::upper_bound(stage_heights.begin(),
                         stage_heights.end(),
                         face_center[2]); // TODO Generalize for x and y
      id_out_plane =
        std::distance(stage_heights.begin(), upper_height_iterator) - 1;
    }

  // Angle of the face center
  double c = std::fmod(to_1D(face_center), two_pi);
  if (c < 0.0)
    c += two_pi;

  // This stage's own and merged arrangements.
  const auto &own =
    (is_inner ? rotor_breakpoints : stator_breakpoints)[id_out_plane];
  const auto        &merged = merged_breakpoints[id_out_plane];
  const unsigned int S      = static_cast<unsigned int>(merged.size());

  // Locate, on the face's own side, the arc [lo, hi) that contains the center,
  // treating the breakpoint list cyclically (the last arc wraps past 2*pi).
  double     lo, hi;
  const auto it = std::upper_bound(own.begin(), own.end(), c);
  if (it == own.begin() || it == own.end())
    {
      lo = own.back();
      hi = own.front() + two_pi;
    }
  else
    {
      lo = *(it - 1);
      hi = *it;
    }

  // Find the merged-breakpoint index closest to lo (lo is itself a breakpoint,
  // hence present in the merged arrangement within tolerance).
  const auto find_merged = [&](const double angle) -> unsigned int {
    const auto lb   = std::lower_bound(merged.begin(), merged.end(), angle);
    const long base = lb - merged.begin();

    unsigned int best  = 0;
    double       bestd = std::numeric_limits<double>::max();
    for (long off = -1; off <= 1; ++off)
      {
        const unsigned int k =
          static_cast<unsigned int>(((base + off) % S + S) % S);
        double d = std::abs(merged[k] - angle);
        d        = std::min(d, two_pi - d);
        if (d < bestd)
          {
            bestd = d;
            best  = k;
          }
      }
    return best;
  };

  const unsigned int p = find_merged(lo);

  // Walk the merged breakpoints from p, accumulating physical (unwrapped)
  // angles until reaching hi. The face arc is a union of whole merged segments,
  // so this terminates cleanly on a breakpoint coinciding with hi.
  std::vector<double> arc_breaks;
  arc_breaks.push_back(merged[p]);

  unsigned int k       = p;
  double       current = merged[p];
  while (current < hi - tol)
    {
      const unsigned int knext = (k + 1) % S;
      double             next  = merged[knext];
      while (next <= current + tol)
        next += two_pi;
      arc_breaks.push_back(next);
      current = next;
      k       = knext;
    }

  return {id_out_plane, p, arc_breaks};
}


/*-------------- Auxiliary Functions -------------------------------*/
template <int dim>
std::vector<unsigned int>
compute_number_interface_cells(const Triangulation<dim>      &triangulation,
                               const Parameters::Mortar<dim> &mortar_parameters)
{
  // Number of subdivisions per process
  unsigned int n_subdivisions_local = 0;
  // Number of subdivisions in the radial direction per process
  unsigned int n_subdivisions_plane_local = 0;
  // Tolerance for rotor radius computation
  const double radius_tolerance = mortar_parameters.radius_tolerance;

  // Coordinate of the reference cell for computation of the number of
  // subdivisions in the radial direction. Used in 3D case
  double coord_ref_local = std::numeric_limits<double>::max();
  // Rotation axis direction. Used in 3D case
  unsigned int direction = 0;

  // Verify if rotation axis is a unit vector in x, y, or z
  if constexpr (dim == 3)
    {
      // First check if the vector has a unit norm
      bool is_unit_axis =
        mortar_parameters.rotation_axis.norm() == 1 ? true : false;

      // Now check if the vector represents the x, y, or z directions
      // specifically (we assume those are the only options for now)
      for (int i = 0; i < dim; i++)
        if (mortar_parameters.rotation_axis[i] != 0. &&
            mortar_parameters.rotation_axis[i] != 1.)
          is_unit_axis = false;

      AssertThrow(
        is_unit_axis,
        ExcMessage(
          "The rotation axis must be a unit vector in x, y, or z direction."));

      // Check if the rotation axis is aligned with the z direction
      // For now, it is the only direction supported; this throw can be removed
      // when the x and y directions are also supported
      AssertThrow(
        mortar_parameters.rotation_axis[2] == 1.,
        ExcMessage(
          "Currently, only rotation axes aligned with the z direction are supported."));

      // Find the direction of the rotation axis
      for (int d = 0; d < dim; d++)
        if (mortar_parameters.rotation_axis[d] != 0.0)
          direction = d;

      // To compute the number of subdivisions in the radial direction, we will
      // use a reference cell. Its coordinate in the direction of the rotation
      // axis is stored; then, all the cells at the same plane will be counted
      // as an "in-plane" (radial direction) subdivision
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              for (const auto face_no : cell->face_indices())
                {
                  const auto face = cell->face(face_no);

                  if (face->at_boundary() &&
                      face->boundary_id() ==
                        mortar_parameters.rotor_boundary_id)
                    coord_ref_local =
                      std::min(coord_ref_local, cell->center()[direction]);
                }
            }
        }
    }

  // Coordinate of reference cell at all processes
  const double coord_ref =
    Utilities::MPI::min(coord_ref_local, triangulation.get_mpi_communicator());

  // Check number of faces and vertices at the rotor-stator interface
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (const auto face_no : cell->face_indices())
            {
              const auto face = cell->face(face_no);

              if (face->at_boundary())
                {
                  if (face->boundary_id() ==
                      mortar_parameters.rotor_boundary_id)
                    {
                      n_subdivisions_local++;

                      // Store the number of subdivisions in the radial
                      // direction
                      if constexpr (dim == 3)
                        {
                          // Check if the current cell is contained in the same
                          // plane as the reference cell
                          if (std::abs(cell->center()[direction] - coord_ref) <
                              radius_tolerance)
                            n_subdivisions_plane_local++;
                        }
                      else
                        n_subdivisions_plane_local++;
                    }
                }
            }
        }
    }

  // Total number of faces
  const unsigned int n_subdivisions =
    Utilities::MPI::sum(n_subdivisions_local,
                        triangulation.get_mpi_communicator());

  // Total number of faces at the radial direction
  const unsigned int n_subdivisions_plane =
    Utilities::MPI::sum(n_subdivisions_plane_local,
                        triangulation.get_mpi_communicator());

  return {n_subdivisions_plane, n_subdivisions / n_subdivisions_plane};
}

template <int dim>
std::tuple<std::vector<double>, double>
compute_interface_dimensions_circular(
  const Triangulation<dim>      &triangulation,
  const Mapping<dim>            &mapping,
  const Parameters::Mortar<dim> &mortar_parameters)
{
  // Tolerance for rotor radius computation
  const double radius_tolerance = mortar_parameters.radius_tolerance;
  // Min and max values for rotor radius computation
  double radius_min = std::numeric_limits<double>::max();
  double radius_max = 0;
  // Minimum rotation angle in initial mesh configuration
  double pre_rotation_min = std::numeric_limits<double>::max();

  // Rotation axis direction. Used in 3D case
  unsigned int direction = 0;
  // Min and max vertex coordinates for length computation in the axial
  // direction. Used in 3D case
  double vertex_min = std::numeric_limits<double>::max();
  double vertex_max = std::numeric_limits<double>::lowest();

  // Verify if rotation axis is a unit vector in x, y, or z
  if constexpr (dim == 3)
    {
      Assert(mortar_parameters.rotation_axis.norm() > 0,
             ExcMessage("The rotation axis must be non-zero."));

      // First check if the vector has a unit norm
      bool is_unit_axis =
        mortar_parameters.rotation_axis.norm() == 1 ? true : false;

      // Now check if the vector represents the x, y, or z directions
      // specifically (we assume those are the only options for now)
      for (int i = 0; i < dim; i++)
        if (mortar_parameters.rotation_axis[i] != 0. &&
            mortar_parameters.rotation_axis[i] != 1.)
          is_unit_axis = false;

      AssertThrow(
        is_unit_axis,
        ExcMessage(
          "The rotation axis must be a unit vector in x, y, or z direction."));

      // Find the direction of the rotation axis
      for (int d = 0; d < dim; d++)
        if (mortar_parameters.rotation_axis[d] != 0.0)
          direction = d;
    }


  // Check number of faces and vertices at the rotor-stator interface
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (const auto face_no : cell->face_indices())
            {
              const auto face = cell->face(face_no);

              if (face->at_boundary())
                {
                  if (face->boundary_id() ==
                      mortar_parameters.rotor_boundary_id)
                    {
                      const auto vertices = mapping.get_vertices(cell, face_no);

                      for (unsigned int vertex_index = 0;
                           vertex_index < face->n_vertices();
                           vertex_index++)
                        {
                          const auto v = vertices[vertex_index];
                          double     radius_current =
                            v.distance(mortar_parameters.center_of_rotation);

                          // In 3D, the interface radius is computed as the
                          // minimum distance between the current vertex and the
                          // rotation axis, assuming that the center of rotation
                          // point is contained within the rotation axis line
                          if constexpr (dim == 3)
                            {
                              vertex_min = std::min(vertex_min, v[direction]);
                              vertex_max = std::max(vertex_max, v[direction]);

                              radius_current =
                                LetheGridTools::find_point_line_distance(
                                  mortar_parameters.center_of_rotation,
                                  mortar_parameters.rotation_axis,
                                  v);
                            }

                          radius_min = std::min(radius_min, radius_current);
                          radius_max = std::max(radius_max, radius_current);
                        }
                    }
                  // Obtain the minimum initial rotation angle based on the
                  // stator interface
                  else if (face->boundary_id() ==
                           mortar_parameters.stator_boundary_id)
                    {
                      const auto vertices = mapping.get_vertices(cell, face_no);

                      for (unsigned int vertex_index = 0;
                           vertex_index < face->n_vertices();
                           vertex_index++)
                        {
                          const auto v = vertices[vertex_index];

                          pre_rotation_min = std::min(
                            pre_rotation_min,
                            point_to_angle(mortar_parameters.center_of_rotation,
                                           v));
                        }
                    }
                }
            }
        }
    }

  // Min and max values over all processes
  radius_min =
    Utilities::MPI::min(radius_min, triangulation.get_mpi_communicator());
  radius_max =
    Utilities::MPI::max(radius_max, triangulation.get_mpi_communicator());

  pre_rotation_min =
    Utilities::MPI::min(pre_rotation_min, triangulation.get_mpi_communicator());

  vertex_min =
    Utilities::MPI::min(vertex_min, triangulation.get_mpi_communicator());
  vertex_max =
    Utilities::MPI::max(vertex_max, triangulation.get_mpi_communicator());

  // Radius variation
  const auto radius_diff = std::abs(radius_max - radius_min);

  // Check if variation is withing the prescribed tolerance
  AssertThrow(radius_diff < radius_tolerance,
              InterfaceRadiusTolerance(radius_diff));

  // Final radius value
  const double radius = radius_min;
  // Length along the axial direction
  const double length_rot_axis = dim == 3 ? vertex_max - vertex_min : 0;

  return {{radius, length_rot_axis}, pre_rotation_min};
}

template <int dim>
std::pair<double, double>
compute_interface_dimensions_linear(
  const Triangulation<dim>      &triangulation,
  const Mapping<dim>            &mapping,
  const Parameters::Mortar<dim> &mortar_parameters)
{
  // y coordinates of the interface limits
  // We assume that the mortar interface is always parallel to the y axis
  double coord_min = std::numeric_limits<double>::max();
  double coord_max = 0;

  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (const auto face_no : cell->face_indices())
            {
              const auto face = cell->face(face_no);

              if (face->at_boundary())
                {
                  if (face->boundary_id() ==
                      mortar_parameters.rotor_boundary_id)
                    {
                      const auto vertices = mapping.get_vertices(cell, face_no);

                      for (unsigned int vertex_index = 0;
                           vertex_index < face->n_vertices();
                           vertex_index++)
                        {
                          const auto v = vertices[vertex_index];
                          coord_min    = std::min(coord_min, v[1]);
                          coord_max    = std::max(coord_max, v[1]);
                        }
                    }
                }
            }
        }
    }

  // Min and max values over all processes
  coord_min =
    Utilities::MPI::min(coord_min, triangulation.get_mpi_communicator());
  coord_max =
    Utilities::MPI::max(coord_max, triangulation.get_mpi_communicator());

  return {coord_min, coord_max};
}

template <int dim>
Quadrature<dim>
construct_quadrature(const Quadrature<dim>         &quadrature,
                     const Parameters::Mortar<dim> &mortar_parameters)
{
  const int oversampling_factor =
    static_cast<int>(mortar_parameters.oversampling_factor);

  for (unsigned int i = 1; i <= 10; ++i)
    if (quadrature == QGauss<dim>(i))
      return QGauss<dim>(i * oversampling_factor);

  AssertThrow(false, ExcNotImplemented());

  return quadrature;
}

template <int dim>
std::vector<double>
compute_stage_heights(const Triangulation<dim>      &triangulation,
                      const Parameters::Mortar<dim> &mortar_parameters)
{
  if constexpr (dim == 3)
    {
      // Direction of the rotation axis
      int            direction = 0;
      Tensor<1, dim> axis      = mortar_parameters.rotation_axis;
      for (int d = 1; d < dim; ++d)
        if (std::abs(axis[d]) > std::abs(axis[direction]))
          direction = d;

      AssertThrow(
        std::abs(axis[direction]) > 0.99,
        ExcMessage(
          "Rotation axis is not aligned with a coordinate direction."));

      // Container storing all vertex coordinates in the rotation axis direction
      // for the local cells at the mortar boundary
      std::vector<double> stage_heights_local;
      // Smallest mortar cell face diameter, used to set the tolerance for
      // height comparison
      double minimum_face_diameter_local = std::numeric_limits<double>::max();

      // Loop over the cells to store the vertex coordinates in the rotation
      // axis direction
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              for (const auto face_no : cell->face_indices())
                {
                  const auto face = cell->face(face_no);

                  // Check if the face is at the boundary and belongs to the
                  // stator boundary
                  // The choice of the stator boundary is arbitrary
                  if (face->at_boundary() &&
                      face->boundary_id() ==
                        mortar_parameters.stator_boundary_id)
                    {
                      minimum_face_diameter_local =
                        std::min(minimum_face_diameter_local, face->diameter());

                      for (const auto vertex_no : face->vertex_indices())
                        stage_heights_local.push_back(
                          face->vertex(vertex_no)[direction]);
                    }
                }
            }
        }
      // Set the minimum face size over all processes
      double minimum_face_diameter =
        Utilities::MPI::min(minimum_face_diameter_local,
                            triangulation.get_mpi_communicator());

      // Store the vertex coordinates in the rotation axis direction from all
      // processes
      const auto stage_heights_all =
        Utilities::MPI::all_gather(triangulation.get_mpi_communicator(),
                                   stage_heights_local);
      std::vector<double> stage_heights;
      for (const auto &heights_per_process : stage_heights_all)
        stage_heights.insert(stage_heights.end(),
                             heights_per_process.begin(),
                             heights_per_process.end());

      // Set the tolerance for height comparison based on the minimum face size
      const double height_tolerance = minimum_face_diameter * 1e-8;

      // Remove duplicate heights within the specified tolerance
      std::ranges::sort(stage_heights);
      auto result =
        std::ranges::unique(stage_heights,
                            [height_tolerance](const double a, const double b) {
                              return std::abs(a - b) <= height_tolerance;
                            });
      stage_heights.erase(result.begin(), result.end());

      return stage_heights;
    }
  else
    return {0.0};
}

template <int dim>
void
mortar_workload_imbalance(const Triangulation<dim>      &triangulation,
                          const Parameters::Mortar<dim> &mortar_parameters,
                          const ConditionalOStream      &pcout)
{
  unsigned int n_mortar_cells = 0;

  // Identify number of cells in each process (local workload)
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto face_no : cell->face_indices())
        {
          const auto face = cell->face(face_no);

          if (face->at_boundary() &&
              (face->boundary_id() == mortar_parameters.rotor_boundary_id ||
               face->boundary_id() == mortar_parameters.stator_boundary_id))
            ++n_mortar_cells;
        }

  // Compute minimum, maximum, and summation of cells over all processes
  const auto [n_mortar_cells_sum,
              n_mortar_cells_min,
              n_mortar_cells_max,
              _,
              __,
              ___] =
    Utilities::MPI::min_max_avg(n_mortar_cells,
                                triangulation.get_mpi_communicator());

  const unsigned int n_proc =
    Utilities::MPI::n_mpi_processes(triangulation.get_mpi_communicator());

  // Ideal work: same number of cells per process
  const double ideal_work = n_mortar_cells_sum / static_cast<double>(n_proc);

  // Workload imbalance (the closest to 1.0, the better)
  const double workload_imbalance = n_mortar_cells_max / ideal_work;

  pcout << "Workload imbalance: " << workload_imbalance << std::endl;
  pcout << std::defaultfloat;
  pcout << "Number of cells per process - Min.: " << n_mortar_cells_min
        << std::endl;
  pcout << "                              Max.: " << n_mortar_cells_max
        << std::endl;
  pcout << std::scientific;
}

/*-------------- MortarManagerCircle -------------------------------*/
template <int dim>
Point<dim>
MortarManagerCircle<dim>::from_1D(const double angle_rad) const
{
  return radius_to_point<dim>(this->radius[0], angle_rad + pre_rotation_angle);
}

template <int dim>
double
MortarManagerCircle<dim>::to_1D(const Point<dim> &point) const
{
  return std::fmod(point_to_angle(point, this->center_of_rotation) -
                     pre_rotation_angle + 2 * numbers::PI,
                   2 * numbers::PI);
}

template <int dim>
Tensor<1, dim, double>
MortarManagerCircle<dim>::get_normal(const Point<dim> &point) const
{
  return point / point.norm();
}


/*-------------- MortarManagerLinear -------------------------------*/
template <int dim>
Point<dim>
MortarManagerLinear<dim>::from_1D(const double angle_rad) const
{
  return Point<dim>(0.5,
                    angle_rad / (2.0 * numbers::PI) * (coord_max - coord_min) +
                      coord_min);
}

template <int dim>
double
MortarManagerLinear<dim>::to_1D(const Point<dim> &point) const
{
  return (2.0 * numbers::PI) * (point[1] - coord_min) / (coord_max - coord_min);
}

template <int dim>
Tensor<1, dim, double>
MortarManagerLinear<dim>::get_normal(const Point<dim> &) const
{
  return Point<dim>(1.0, 0.0);
}

/*-------------- CouplingOperator -------------------------------*/
template <int dim, typename Number>
CouplingOperator<dim, Number>::CouplingOperator(
  const Mapping<dim>                                        &mapping,
  const DoFHandler<dim>                                     &dof_handler,
  const AffineConstraints<Number>                           &constraints,
  const std::shared_ptr<CouplingEvaluationBase<dim, Number>> evaluator,
  const std::shared_ptr<MortarManagerBase<dim>>              mortar_manager,
  const unsigned int                                         bid_m,
  const unsigned int                                         bid_p,
  const double                                               sip_factor)
  : mapping(mapping)
  , dof_handler(dof_handler)
  , bid_m(bid_m)
  , bid_p(bid_p)
  , evaluator(evaluator)
  , mortar_manager(mortar_manager)
{
  this->q_data_size          = evaluator->data_size();
  this->relevant_dof_indices = evaluator->get_relevant_dof_indices();
  this->n_dofs_per_cell      = this->relevant_dof_indices.size();

  data.penalty_factor =
    compute_penalty_factor(dof_handler.get_fe().degree, sip_factor);

  // Number of cells
  const unsigned int n_sub_cells = mortar_manager->get_n_total_mortars();

  std::vector<types::global_dof_index> is_local_cell;
  std::vector<types::global_dof_index> is_ghost_cell;

#ifdef DEBUG
  std::vector<double> vec_local_cells(n_sub_cells * 2, 0.0);
  std::vector<double> vec_ghost_cells(n_sub_cells * 2, 0.0);
#endif

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto face_no : cell->face_indices())
        if ((cell->face(face_no)->boundary_id() == bid_m) ||
            (cell->face(face_no)->boundary_id() == bid_p))
          {
            const auto face = cell->face(face_no);

            const auto center = get_face_center(cell, face);

            // Indices of mortars on face of cell.
            const auto indices = mortar_manager->get_mortar_indices(
              center, cell->face(face_no)->boundary_id() == bid_m);

            // Number of mortars covered by this face (1, 2, or more). Stored in
            // face-iteration order so the assembly loops can advance their
            // per-face pointers without relying on a global constant.
            this->n_mortars_per_face.emplace_back(
              static_cast<unsigned int>(indices.size()));

            const auto local_dofs = this->get_dof_indices(cell);

            // Loop over the mortar indices at the rotor-stator
            // interface. The logic of local (rotor) and ghost (stator) is the
            // same as in the previous loop.
            for (unsigned int ii = 0; ii < indices.size(); ++ii)
              {
                const unsigned int i        = indices[ii];
                unsigned int       id_local = 0, id_ghost = 0;

                if (face->boundary_id() == bid_m)
                  {
                    id_local = i;
                    id_ghost = i + n_sub_cells;
                  }
                else if (face->boundary_id() == bid_p)
                  {
                    id_local = i + n_sub_cells;
                    id_ghost = i;
                  }

#ifdef DEBUG
                vec_local_cells[id_local] += 1.0;
                vec_ghost_cells[id_ghost] += 1.0;
#endif
                is_local_cell.emplace_back(id_local);
                is_ghost_cell.emplace_back(id_ghost);

                for (const auto l_dof : local_dofs)
                  dof_indices.emplace_back(l_dof);
              }

            // Weights of quadrature points
            const auto weights =
              mortar_manager->get_weights(get_face_center(cell, face),
                                          cell->face(face_no)->boundary_id() ==
                                            bid_m);
            data.all_weights.insert(data.all_weights.end(),
                                    weights.begin(),
                                    weights.end());

            // Normals of quadrature points
            if constexpr (dim == 3)
              {
                const auto points = mortar_manager->get_points(
                  get_face_center(cell, face),
                  cell->face(face_no)->boundary_id() == bid_m);
                std::vector<Point<dim, Number>> points_ref(points.size());
                mapping.transform_points_real_to_unit_cell(cell,
                                                           points,
                                                           points_ref);
                all_points_ref.insert(all_points_ref.end(),
                                      points_ref.begin(),
                                      points_ref.end());

                std::vector<Point<dim - 1>> quad;

                for (const auto p : points_ref)
                  {
                    Point<dim - 1> temp;
                    for (int i = 0, j = 0; i < dim; ++i)
                      if ((face_no / 2) != static_cast<unsigned int>(i))
                        temp[j++] = p[i];

                    if ((dim == 3) && ((face_no / 2) == 1))
                      std::swap(temp[0], temp[1]);

                    quad.emplace_back(temp);
                  }

                FEFaceValues<dim> fe_face_values(mapping,
                                                 cell->get_fe(),
                                                 quad,
                                                 update_normal_vectors);

                fe_face_values.reinit(cell, face_no);

                data.all_normals.insert(
                  data.all_normals.end(),
                  fe_face_values.get_normal_vectors().begin(),
                  fe_face_values.get_normal_vectors().end());
              }
            else
              {
                auto normals = mortar_manager->get_normals(
                  get_face_center(cell, face),
                  cell->face(face_no)->boundary_id() == bid_m);
                if (face->boundary_id() == bid_p)
                  for (auto &normal : normals)
                    normal *= -1.0;
                data.all_normals.insert(data.all_normals.end(),
                                        normals.begin(),
                                        normals.end());

                if constexpr (dim == 1)
                  {
                    if (face_no == 0)
                      all_points_ref.emplace_back(0.0);
                    else if (face_no == 1)
                      all_points_ref.emplace_back(1.0);
                    else
                      AssertThrow(false, ExcNotImplemented());
                  }
                else if constexpr (dim == 2)
                  {
                    auto points = mortar_manager->get_points_ref(
                      get_face_center(cell, face),
                      cell->face(face_no)->boundary_id() == bid_m);

                    const bool flip =
                      (face->vertex(0)[0] * face->vertex(1)[1] -
                       face->vertex(0)[1] * face->vertex(1)[0]) < 0.0;

                    if (flip)
                      for (auto &p : points)
                        p[0] = 1.0 - p[0];

                    if (face_no / 2 == 0)
                      {
                        for (auto &p : points)
                          all_points_ref.emplace_back(face_no % 2, p[0]);
                      }
                    else if (face_no / 2 == 1)
                      {
                        for (auto &p : points)
                          all_points_ref.emplace_back(p[0], face_no % 2);
                      }
                    else
                      {
                        AssertThrow(false, ExcNotImplemented());
                      }
                  }
                else if constexpr (dim == 3)
                  {
                    AssertThrow(false, ExcNotImplemented()); // TODO
                  }
                else
                  AssertThrow(false, ExcNotImplemented());
              }

            // Penalty parameter
            const Number penalty_parameter = compute_penalty_parameter(cell);

            // Store penalty parameter for all quadrature points of this face
            const unsigned int n_face_points = mortar_manager->get_n_points(
              center, face->boundary_id() == bid_m);
            for (unsigned int i = 0; i < n_face_points; ++i)
              data.all_penalty_parameter.emplace_back(penalty_parameter);
          }

#ifdef DEBUG
  Utilities::MPI::sum(vec_local_cells,
                      dof_handler.get_mpi_communicator(),
                      vec_local_cells);
  Utilities::MPI::sum(vec_ghost_cells,
                      dof_handler.get_mpi_communicator(),
                      vec_ghost_cells);

  std::set<unsigned int> vec_local_cells_0;
  std::set<unsigned int> vec_local_cells_2;

  if (Utilities::MPI::this_mpi_process(dof_handler.get_mpi_communicator()) == 0)
    {
      for (unsigned int i = 0; i < vec_local_cells.size(); ++i)
        {
          if (vec_local_cells[i] == 0)
            vec_local_cells_0.insert(i);
          if (vec_local_cells[i] > 1)
            vec_local_cells_2.insert(i);
        }

      std::set<unsigned int> vec_ghost_cells_0;
      std::set<unsigned int> vec_ghost_cells_2;

      for (unsigned int i = 0; i < vec_ghost_cells.size(); ++i)
        {
          if (vec_ghost_cells[i] == 0)
            vec_ghost_cells_0.insert(i);
          if (vec_ghost_cells[i] > 1)
            vec_ghost_cells_2.insert(i);
        }

      if (!(vec_local_cells_0.empty() && vec_local_cells_2.empty() &&
            vec_ghost_cells_0.empty() && vec_ghost_cells_2.empty()))
        {
          std::cout << "CouplingOperator mortar matching failed:" << std::endl;

          if (!vec_local_cells_0.empty())
            {
              std::cout << " - some local cells are not owned: ";
              for (const auto &i : vec_local_cells_0)
                std::cout << i << " ";
              std::cout << std::endl;
            }
          if (!vec_local_cells_2.empty())
            {
              std::cout << " - some local cells are owned multiple times: ";
              for (const auto &i : vec_local_cells_2)
                std::cout << i << " ";
              std::cout << std::endl;
            }
          if (!vec_ghost_cells_0.empty())
            {
              std::cout << " - some ghost cells are not owned: ";
              for (const auto &i : vec_ghost_cells_0)
                std::cout << i << " ";
              std::cout << std::endl;
            }
          if (!vec_ghost_cells_2.empty())
            {
              std::cout << " - some ghost cells are owned multiple times: ";
              for (const auto &i : vec_ghost_cells_2)
                std::cout << i << " ";
              std::cout << std::endl;
            }

          MPI_Barrier(dof_handler.get_mpi_communicator());

          AssertThrow(false, ExcInternalError());
        }
    }
#endif

  // Setup communication
  partitioner.reinit(is_local_cell,
                     is_ghost_cell,
                     dof_handler.get_mpi_communicator());

  // Finalized penalty parameters
  const unsigned n_q_points = mortar_manager->get_n_quadrature_points();
  std::vector<Number> all_penalty_parameter_ghost(
    data.all_penalty_parameter.size());
  partitioner.export_to_ghosted_array<Number, 0>(data.all_penalty_parameter,
                                                 all_penalty_parameter_ghost,
                                                 n_q_points);
  for (unsigned int i = 0; i < data.all_penalty_parameter.size(); ++i)
    data.all_penalty_parameter[i] =
      std::max(data.all_penalty_parameter[i], all_penalty_parameter_ghost[i]);

  // Finalize DoF indices and update constraints
  dof_indices_ghost.resize(dof_indices.size());
  partitioner.export_to_ghosted_array<types::global_dof_index, 0>(
    dof_indices, dof_indices_ghost, n_dofs_per_cell);

  {
    auto locally_owned_dofs = constraints.get_locally_owned_indices();
    auto constraints_to_make_consistent = constraints.get_local_lines();


    for (unsigned int i = 0; i < dof_indices.size(); ++i)
      {
        constraints_to_make_consistent.add_index(dof_indices[i]);
        constraints_to_make_consistent.add_index(dof_indices_ghost[i]);
      }

    constraints_extended.reinit(locally_owned_dofs,
                                constraints_to_make_consistent);
    constraints_extended.merge(
      constraints,
      AffineConstraints<Number>::MergeConflictBehavior::no_conflicts_allowed,
      true);

    constraints_extended.make_consistent_in_parallel(
      locally_owned_dofs,
      constraints_to_make_consistent,
      dof_handler.get_mpi_communicator());

    partitioner_extended = std::make_shared<const Utilities::MPI::Partitioner>(
      locally_owned_dofs,
      constraints_extended.get_local_lines(),
      dof_handler.get_mpi_communicator());
  }
}

template <int dim, typename Number>
const AffineConstraints<Number> &
CouplingOperator<dim, Number>::get_affine_constraints() const
{
  return constraints_extended;
}

template <int dim, typename Number>
Number
CouplingOperator<dim, Number>::compute_penalty_factor(const unsigned int degree,
                                                      const Number       factor)
{
  return factor * (degree + 1.0) * (degree + 1.0);
}

template <int dim, typename Number>
Number
CouplingOperator<dim, Number>::compute_penalty_parameter(
  const typename Triangulation<dim>::cell_iterator &cell) const
{
  const unsigned int degree = dof_handler.get_fe().degree;

  FE_Nothing<dim> fe_nothing;

  dealii::QGauss<dim>   quadrature(degree + 1);
  dealii::FEValues<dim> fe_values(mapping,
                                  fe_nothing,
                                  quadrature,
                                  dealii::update_JxW_values);

  dealii::QGauss<dim - 1>   face_quadrature(degree + 1);
  dealii::FEFaceValues<dim> fe_face_values(mapping,
                                           fe_nothing,
                                           face_quadrature,
                                           dealii::update_JxW_values);

  fe_values.reinit(cell);

  Number volume = 0;
  for (unsigned int q = 0; q < quadrature.size(); ++q)
    volume += fe_values.JxW(q);

  Number surface_area = 0;
  for (const auto f : cell->face_indices())
    {
      fe_face_values.reinit(cell, f);

      const Number factor =
        (cell->at_boundary(f) && !cell->has_periodic_neighbor(f) &&
         (cell->face(f)->boundary_id() != bid_m &&
          cell->face(f)->boundary_id() != bid_p)) ?
          1. :
          0.5;

      for (unsigned int q = 0; q < face_quadrature.size(); ++q)
        surface_area += fe_face_values.JxW(q) * factor;
    }

  return surface_area / volume;
}

template <int dim, typename Number>
Point<dim>
CouplingOperator<dim, Number>::get_face_center(
  const typename Triangulation<dim>::cell_iterator &cell,
  const typename Triangulation<dim>::face_iterator &face) const
{
  return mapping.transform_unit_to_real_cell(
    cell, MappingQ1<dim>().transform_real_to_unit_cell(cell, face->center()));
}

template <int dim, typename Number>
std::vector<types::global_dof_index>
CouplingOperator<dim, Number>::get_dof_indices(
  const typename DoFHandler<dim>::active_cell_iterator &cell) const
{
  std::vector<types::global_dof_index> local_dofs_all(
    dof_handler.get_fe().n_dofs_per_cell());
  cell->get_dof_indices(local_dofs_all);

  std::vector<types::global_dof_index> local_dofs(n_dofs_per_cell);

  for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
    local_dofs[i] = local_dofs_all[relevant_dof_indices[i]];

  return local_dofs;
}

template <int dim, typename Number>
template <typename VectorType>
void
CouplingOperator<dim, Number>::vmult_add(VectorType       &dst,
                                         const VectorType &src) const
{
  VectorType dst_internal;

  if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
    {
      dst_internal.reinit(this->partitioner_extended->locally_owned_range(),
                          this->partitioner_extended->get_mpi_communicator());
    }
  else
    dst_internal.reinit(this->partitioner_extended);

  VectorType src_internal;
  if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
    src_internal.reinit(this->partitioner_extended->locally_owned_range(),
                        this->partitioner_extended->ghost_indices(),
                        this->partitioner_extended->get_mpi_communicator());
  else
    src_internal.reinit(this->partitioner_extended);
  src_internal = src;
  src_internal.update_ghost_values();

  // Number of quadrature points per mortar segment (constant; only the number
  // of mortars per face varies).
  const unsigned int n_q_per_mortar = mortar_manager->get_n_quadrature_points();

  // 1) Evaluate
  unsigned int ptr_q      = 0;
  unsigned int face_index = 0;

  Vector<Number> buffer;

  std::vector<Number> all_values_local(data.all_normals.size() * q_data_size);
  std::vector<Number> all_values_ghost(data.all_normals.size() * q_data_size);

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
          {
            // Quadrature points at the current face (all mortars together).
            const unsigned int n_q_points =
              n_mortars_per_face[face_index] * n_q_per_mortar;

            evaluator->local_reinit(
              cell,
              ArrayView<const Point<dim, Number>>(all_points_ref.data() + ptr_q,
                                                  n_q_points));

            buffer.reinit(n_dofs_per_cell);

            const auto local_dofs = this->get_dof_indices(cell);

            for (unsigned int i = 0; i < local_dofs.size(); ++i)
              buffer[i] = src_internal[local_dofs[i]];

            evaluator->local_evaluate(data,
                                      buffer,
                                      ptr_q,
                                      1,
                                      all_values_local.data() +
                                        ptr_q * q_data_size);

            ptr_q += n_q_points;
            ++face_index;
          }

  // 2) Communicate
  partitioner.export_to_ghosted_array<Number, 0>(
    ArrayView<const Number>(reinterpret_cast<Number *>(all_values_local.data()),
                            all_values_local.size()),
    ArrayView<Number>(reinterpret_cast<Number *>(all_values_ghost.data()),
                      all_values_ghost.size()),
    n_q_per_mortar * q_data_size);

  // 3) Integrate
  ptr_q      = 0;
  face_index = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
          {
            // Quadrature points at the current face (all mortars together).
            const unsigned int total_n_q_points =
              n_mortars_per_face[face_index] * n_q_per_mortar;

            evaluator->local_reinit(
              cell,
              ArrayView<const Point<dim, Number>>(all_points_ref.data() + ptr_q,
                                                  total_n_q_points));

            buffer.reinit(n_dofs_per_cell);
            evaluator->local_integrate(data,
                                       buffer,
                                       ptr_q,
                                       1,
                                       all_values_local.data() +
                                         ptr_q * q_data_size,
                                       all_values_ghost.data() +
                                         ptr_q * q_data_size);

            const auto local_dofs = this->get_dof_indices(cell);
            constraints_extended.distribute_local_to_global(buffer,
                                                            local_dofs,
                                                            dst_internal);

            ptr_q += total_n_q_points;
            ++face_index;
          }

  dst_internal.compress(VectorOperation::add);
  dst.add(1.0, dst_internal);
}

template <int dim, typename Number>
template <typename VectorType>
void
CouplingOperator<dim, Number>::add_diagonal_entries(VectorType &diagonal) const
{
  VectorType diagonal_internal;

  if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
    {
      diagonal_internal.reinit(
        this->partitioner_extended->locally_owned_range(),
        this->partitioner_extended->get_mpi_communicator());
    }
  else
    diagonal_internal.reinit(this->partitioner_extended);

  const unsigned int n_q_per_mortar = mortar_manager->get_n_quadrature_points();

  unsigned int ptr_q      = 0;
  unsigned int face_index = 0;

  Vector<Number>      buffer, diagonal_local;
  std::vector<Number> all_values_local, all_values_ghost;

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
          {
            // Quadrature points at the current face (all mortars together).
            const unsigned int n_q_points =
              n_mortars_per_face[face_index] * n_q_per_mortar;

            evaluator->local_reinit(
              cell,
              ArrayView<const Point<dim, Number>>(all_points_ref.data() + ptr_q,
                                                  n_q_points));

            buffer.reinit(n_dofs_per_cell);
            diagonal_local.reinit(n_dofs_per_cell);
            all_values_local.resize(n_q_points * q_data_size);
            all_values_ghost.resize(n_q_points * q_data_size);

            for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
              {
                // Create i-th basis vector
                for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                  buffer[j] = static_cast<Number>(i == j);

                // Interpolate i-th basis vector to the quadrature points
                evaluator->local_evaluate(
                  data, buffer, ptr_q, 1, all_values_local.data());

                buffer.reinit(n_dofs_per_cell);

                // integrate the coupling terms of the mortar using the
                // interpolated information from the cell
                evaluator->local_integrate(data,
                                           buffer,
                                           ptr_q,
                                           1,
                                           all_values_local.data(),
                                           all_values_ghost.data());

                diagonal_local[i] = buffer[i];
              }

            const auto local_dofs = this->get_dof_indices(cell);
            constraints_extended.distribute_local_to_global(diagonal_local,
                                                            local_dofs,
                                                            diagonal_internal);

            ptr_q += n_q_points;
            ++face_index;
          }

  diagonal_internal.compress(VectorOperation::add);
  diagonal.add(1.0, diagonal_internal);
}

template <int dim, typename Number>
void
CouplingOperator<dim, Number>::add_sparsity_pattern_entries(
  SparsityPatternBase &dsp) const
{
  for (unsigned int i = 0; i < dof_indices.size(); i += n_dofs_per_cell)
    {
      std::vector<types::global_dof_index> a(dof_indices.begin() + i,
                                             dof_indices.begin() + i +
                                               n_dofs_per_cell);
      std::vector<types::global_dof_index> b(dof_indices_ghost.begin() + i,
                                             dof_indices_ghost.begin() + i +
                                               n_dofs_per_cell);

      constraints_extended.add_entries_local_to_global(a, b, dsp);
      constraints_extended.add_entries_local_to_global(b, a, dsp);
    }
}

template <int dim, typename Number>
void
CouplingOperator<dim, Number>::add_system_matrix_entries(
  TrilinosWrappers::SparseMatrix &system_matrix) const
{
  std::vector<Number> all_values_local(data.all_normals.size() *
                                       n_dofs_per_cell * q_data_size);
  std::vector<Number> all_values_ghost(data.all_normals.size() *
                                       n_dofs_per_cell * q_data_size);

  const unsigned int n_q_per_mortar = mortar_manager->get_n_quadrature_points();

  unsigned int ptr_q      = 0;
  unsigned int face_index = 0;

  Vector<Number> buffer;

  // 1) Evaluate
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
          {
            // Quadrature points at the current face (all mortars together).
            const unsigned int n_q_points =
              n_mortars_per_face[face_index] * n_q_per_mortar;

            // Initialize coupling evaluator with the current cell and
            // the relevant quadrature points
            evaluator->local_reinit(
              cell,
              ArrayView<const Point<dim, Number>>(all_points_ref.data() + ptr_q,
                                                  n_q_points));

            // Initialize buffer to store information at dof level
            buffer.reinit(n_dofs_per_cell);

            for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
              {
                // Create i-th basis vector
                for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                  buffer[j] = static_cast<Number>(i == j);

                // Interpolate i-th basis vector to the quadrature points
                evaluator->local_evaluate(data,
                                          buffer,
                                          ptr_q,
                                          n_dofs_per_cell,
                                          all_values_local.data() +
                                            (ptr_q * n_dofs_per_cell + i) *
                                              q_data_size);
              }

            ptr_q += n_q_points;
            ++face_index;
          }

  // 2) Communicate: export data from local to ghost side
  partitioner.export_to_ghosted_array<Number, 0>(
    ArrayView<const Number>(reinterpret_cast<Number *>(all_values_local.data()),
                            all_values_local.size()),
    ArrayView<Number>(reinterpret_cast<Number *>(all_values_ghost.data()),
                      all_values_ghost.size()),
    n_dofs_per_cell * n_q_per_mortar * q_data_size);


  ptr_q                 = 0;
  unsigned int ptr_dofs = 0;
  face_index            = 0;

  // 3) Integrate
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
          {
            // Number of mortars covered by the current face (1 for aligned
            // meshes, 2 for the legacy non-aligned case, 3 or more for
            // non-uniform interface meshes).
            const unsigned int n_mortars = n_mortars_per_face[face_index];

            // Loop over mortars
            for (unsigned int m = 0; m < n_mortars; ++m)
              {
                evaluator->local_reinit(cell,
                                        ArrayView<const Point<dim, Number>>(
                                          all_points_ref.data() + ptr_q,
                                          n_q_per_mortar));

                // Loop over local and ghost cells attached the mortar
                for (unsigned int b = 0; b < 2; ++b)
                  {
                    FullMatrix<Number> cell_matrix(n_dofs_per_cell,
                                                   n_dofs_per_cell);

                    // Loop over cell dofs and integrate the coupling terms of
                    // the mortar using the interpolated information from the
                    // cell
                    for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                      {
                        buffer.reinit(n_dofs_per_cell);
                        if (b == 0) // local cell
                          evaluator->local_integrate(
                            data,
                            buffer,
                            ptr_q,
                            n_dofs_per_cell,
                            all_values_local.data() +
                              (ptr_q * n_dofs_per_cell + i) * q_data_size,
                            nullptr);
                        else // ghost cell
                          evaluator->local_integrate(
                            data,
                            buffer,
                            ptr_q,
                            n_dofs_per_cell,
                            nullptr,
                            all_values_ghost.data() +
                              (ptr_q * n_dofs_per_cell + i) * q_data_size);

                        // Copy data from buffer to cell matrix
                        for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                          cell_matrix[j][i] = buffer[j];
                      }

                    // Vector of local dof indices from the cell in the negative
                    // ('mortar') side
                    std::vector<types::global_dof_index> local_dof_indices(
                      dof_indices.begin() + ptr_dofs,
                      dof_indices.begin() + ptr_dofs + n_dofs_per_cell);

                    if (b == 0) // local cell -> local-local block
                      {
                        constraints_extended.distribute_local_to_global(
                          cell_matrix, local_dof_indices, system_matrix);
                      }
                    else // ghost cell -> local-ghost block
                      {
                        std::vector<types::global_dof_index>
                          local_dof_indices_ghost(dof_indices_ghost.begin() +
                                                    ptr_dofs,
                                                  dof_indices_ghost.begin() +
                                                    ptr_dofs + n_dofs_per_cell);

                        constraints_extended.distribute_local_to_global(
                          cell_matrix,
                          local_dof_indices,
                          local_dof_indices_ghost,
                          system_matrix);
                      }
                  }

                ptr_dofs += n_dofs_per_cell;

                ptr_q += n_q_per_mortar;
              }

            ++face_index;
          }

  AssertDimension(ptr_q, data.all_normals.size());
  AssertDimension(ptr_dofs, dof_indices.size());
}


/*-------------- CouplingEvaluationSIPG -------------------------------*/
template <int dim, int n_components, typename Number>
CouplingEvaluationSIPG<dim, n_components, Number>::CouplingEvaluationSIPG(
  const Mapping<dim>    &mapping,
  const DoFHandler<dim> &dof_handler,
  const unsigned int     first_selected_component)
  : fe_sub(dof_handler.get_fe().base_element(
             dof_handler.get_fe()
               .component_to_base_index(first_selected_component)
               .first),
           n_components)
  , phi_m(mapping, fe_sub, update_values | update_gradients)
{
  for (unsigned int i = 0; i < dof_handler.get_fe().n_dofs_per_cell(); ++i)
    if ((first_selected_component <=
         dof_handler.get_fe().system_to_component_index(i).first) &&
        (dof_handler.get_fe().system_to_component_index(i).first <
         first_selected_component + n_components))
      relevant_dof_indices.push_back(i);
}

template <int dim, int n_components, typename Number>
unsigned int
CouplingEvaluationSIPG<dim, n_components, Number>::data_size() const
{
  return n_components * 2;
}

template <int dim, int n_components, typename Number>
const std::vector<unsigned int> &
CouplingEvaluationSIPG<dim, n_components, Number>::get_relevant_dof_indices()
  const
{
  return relevant_dof_indices;
}

template <int dim, int n_components, typename Number>
void
CouplingEvaluationSIPG<dim, n_components, Number>::local_reinit(
  const typename Triangulation<dim>::cell_iterator &cell,
  const ArrayView<const Point<dim, Number>>        &points) const
{
  this->phi_m.reinit(cell, points);
}

template <int dim, int n_components, typename Number>
void
CouplingEvaluationSIPG<dim, n_components, Number>::local_evaluate(
  const CouplingEvaluationData<dim, Number> &data,
  const Vector<Number>                      &buffer,
  const unsigned int                         ptr_q,
  const unsigned int                         q_stride,
  Number                                    *all_values_m) const
{
  this->phi_m.evaluate(buffer,
                       EvaluationFlags::values | EvaluationFlags::gradients);

  for (const auto q : this->phi_m.quadrature_point_indices())
    {
      // Quadrature point index ('global' index within the rotor-stator
      // interface)
      const unsigned int q_index = ptr_q + q;

      // Normal, value, and gradient referring to the quadrature point
      const auto normal     = data.all_normals[q_index];
      const auto value_m    = this->phi_m.get_value(q);
      const auto gradient_m = contract(this->phi_m.get_gradient(q), normal);

      // Initialize buffer for 'negative' side of the interface (i.e. rotor),
      // where information is evaluated
      BufferRW<Number> buffer_m(all_values_m, q * 2 * n_components * q_stride);
      // Store values and gradients at the created buffer
      buffer_m.write(value_m);
      buffer_m.write(gradient_m);
    }
}

template <int dim, int n_components, typename Number>
void
CouplingEvaluationSIPG<dim, n_components, Number>::local_integrate(
  const CouplingEvaluationData<dim, Number> &data,
  Vector<Number>                            &buffer,
  const unsigned int                         ptr_q,
  const unsigned int                         q_stride,
  Number                                    *all_values_m,
  Number                                    *all_values_p) const
{
  for (const auto q : this->phi_m.quadrature_point_indices())
    {
      const unsigned int q_index = ptr_q + q;
      // Initialize buffer for both 'mortar' and 'non-mortar' sides of the
      // interface
      BufferRW<Number> buffer_m(all_values_m, q * 2 * n_components * q_stride);
      BufferRW<Number> buffer_p(all_values_p, q * 2 * n_components * q_stride);
      // Read shape functions values and gradients stored in the buffer
      const auto value_m           = buffer_m.template read<value_type>();
      const auto value_p           = buffer_p.template read<value_type>();
      const auto normal_gradient_m = buffer_m.template read<value_type>();
      const auto normal_gradient_p = buffer_p.template read<value_type>();

      const auto JxW               = data.all_weights[q_index];
      const auto penalty_parameter = data.all_penalty_parameter[q_index];
      const auto normal            = data.all_normals[q_index];

      // The expression for the jump on the mortar interface is
      // jump(u) = u_m * normal_m + u_p * normal_p. Since we are accessing only
      // the value of normal_m, we use a minus sign here because normal_p = -
      // normal_m
      const auto value_jump = outer((value_m - value_p), normal);

      // The expression for the average on the mortar interface is
      // avg(∇u).n = (∇u_m.normal_m + ∇u_p.normal_p) * 0.5. For the same reason
      // above, we include the negative sign here
      const auto gradient_normal_avg =
        (normal_gradient_m - normal_gradient_p) * 0.5;

      // SIPG penalty parameter
      const double sigma = penalty_parameter * data.penalty_factor;

      // - (n avg(∇v), jump(u))
      this->phi_m.submit_gradient(-value_jump * 0.5 * JxW, q);

      // + (jump(v), σ jump(u) - avg(∇u) n)
      this->phi_m.submit_value(
        (contract(value_jump, normal) * sigma - gradient_normal_avg) * JxW, q);
    }
  // Multiply previous terms by respective test functions values/gradients
  this->phi_m.test_and_sum(buffer,
                           EvaluationFlags::values |
                             EvaluationFlags::gradients);
}


/*----------- NavierStokesCouplingEvaluation -------------------------*/

template <int dim, typename Number>
NavierStokesCouplingEvaluation<dim, Number>::NavierStokesCouplingEvaluation(
  const Mapping<dim>    &mapping,
  const DoFHandler<dim> &dof_handler,
  const double           kinematic_viscosity)
  : fe_sub_u(dof_handler.get_fe().base_element(
               dof_handler.get_fe().component_to_base_index(0).first),
             dim)
  , fe_sub_p(dof_handler.get_fe().base_element(
               dof_handler.get_fe().component_to_base_index(dim).first),
             1)
  , phi_u_m(mapping, fe_sub_u, update_values | update_gradients)
  , phi_p_m(mapping, fe_sub_p, update_values)
  , kinematic_viscosity(kinematic_viscosity)
{
  for (unsigned int i = 0; i < dof_handler.get_fe().n_dofs_per_cell(); ++i)
    if (dof_handler.get_fe().system_to_component_index(i).first < dim)
      relevant_dof_indices.push_back(i);

  for (unsigned int i = 0; i < dof_handler.get_fe().n_dofs_per_cell(); ++i)
    if (dof_handler.get_fe().system_to_component_index(i).first == dim)
      relevant_dof_indices.push_back(i);

  AssertDimension(dof_handler.get_fe().n_dofs_per_cell(),
                  relevant_dof_indices.size());
}

template <int dim, typename Number>
unsigned int
NavierStokesCouplingEvaluation<dim, Number>::data_size() const
{
  return 4 * dim;
}

template <int dim, typename Number>
const std::vector<unsigned int> &
NavierStokesCouplingEvaluation<dim, Number>::get_relevant_dof_indices() const
{
  return relevant_dof_indices;
}

template <int dim, typename Number>
void
NavierStokesCouplingEvaluation<dim, Number>::local_reinit(
  const typename Triangulation<dim>::cell_iterator &cell,
  const ArrayView<const Point<dim, Number>>        &points) const
{
  this->phi_u_m.reinit(cell, points);
  this->phi_p_m.reinit(cell, points);
}

template <int dim, typename Number>
void
NavierStokesCouplingEvaluation<dim, Number>::local_evaluate(
  const CouplingEvaluationData<dim, Number> &data,
  const Vector<Number>                      &buffer,
  const unsigned int                         ptr_q,
  const unsigned int                         q_stride,
  Number                                    *all_values_m) const
{
  AssertDimension(buffer.size(),
                  fe_sub_u.n_dofs_per_cell() + fe_sub_p.n_dofs_per_cell());

  ArrayView<const Number> buffer_u(buffer.data() + 0,
                                   fe_sub_u.n_dofs_per_cell());
  ArrayView<const Number> buffer_p(buffer.data() + fe_sub_u.n_dofs_per_cell(),
                                   fe_sub_p.n_dofs_per_cell());

  this->phi_u_m.evaluate(buffer_u,
                         EvaluationFlags::values | EvaluationFlags::gradients);
  this->phi_p_m.evaluate(buffer_p, EvaluationFlags::values);

  for (const auto q : this->phi_u_m.quadrature_point_indices())
    {
      // Quadrature point index ('global' index within the rotor-stator
      // interface)
      const unsigned int q_index = ptr_q + q;

      // Normal, value, and gradient referring to the quadrature point
      const auto normal = data.all_normals[q_index];

      const auto u_value = this->phi_u_m.get_value(q);
      const auto u_grad_normal =
        contract(this->phi_u_m.get_gradient(q), normal);
      const auto p_value_normal = this->phi_p_m.get_value(q) * normal;

      // Initialize buffer for local side of the interface (i.e. rotor),
      // where information is evaluated
      BufferRW<Number> buffer_m(all_values_m, q * 4 * dim * q_stride);

      // Store values and gradients at the created buffer
      buffer_m.write(u_value);
      buffer_m.write(u_grad_normal);
      buffer_m.write(p_value_normal);
    }
}

template <int dim, typename Number>
void
NavierStokesCouplingEvaluation<dim, Number>::local_integrate(
  const CouplingEvaluationData<dim, Number> &data,
  Vector<Number>                            &buffer,
  const unsigned int                         ptr_q,
  const unsigned int                         q_stride,
  Number                                    *all_values_m,
  Number                                    *all_values_p) const
{
  for (const auto q : this->phi_u_m.quadrature_point_indices())
    {
      const unsigned int q_index = ptr_q + q;

      // Initialize buffer for both local and ghost sides
      BufferRW<Number> buffer_m(all_values_m, q * 4 * dim * q_stride);
      BufferRW<Number> buffer_p(all_values_p, q * 4 * dim * q_stride);

      // Read shape functions values and gradients stored in the buffer
      const auto u_value_m        = buffer_m.template read<u_value_type>();
      const auto u_value_p        = buffer_p.template read<u_value_type>();
      const auto u_grad_normal_m  = buffer_m.template read<u_value_type>();
      const auto u_grad_normal_p  = buffer_p.template read<u_value_type>();
      const auto p_value_normal_m = buffer_m.template read<u_value_type>();
      const auto p_value_normal_p = buffer_p.template read<u_value_type>();

      const auto JxW               = data.all_weights[q_index];
      const auto penalty_parameter = data.all_penalty_parameter[q_index];
      const auto normal            = data.all_normals[q_index];

      // jump(u) = u_m - u_p
      const auto u_value_jump = u_value_m - u_value_p;

      // The expression for the average on the mortar interface is
      // avg(∇u).n = (∇u_m.normal_m + ∇u_p.normal_p) * 0.5. Since we are
      // accessing only the value of normal_m, we use a minus sign here because
      // normal_p = - normal_m
      const auto u_grad_avg = (u_grad_normal_m - u_grad_normal_p) * 0.5;

      // {{p}} = (p_m + p_p)/2
      const auto p_value_avg = (p_value_normal_m - p_value_normal_p) * 0.5;

      typename FEPointIntegratorU::value_type u_grad_result  = {};
      typename FEPointIntegratorU::value_type u_value_result = {};
      typename FEPointIntegratorP::value_type p_value_result = {};

      // SIPG penalty parameter
      const double sigma = penalty_parameter * data.penalty_factor;

      /* Contributions from viscous term */
      // - (n avg(∇v), jump(u))
      u_grad_result -= u_value_jump;
      // - (jump(v), ν avg(∇δu) n)
      u_value_result -= u_grad_avg * this->kinematic_viscosity;
      // + (jump(v), ν σ jump(δu))
      u_value_result += sigma * this->kinematic_viscosity * u_value_jump;

      /* Contribution from pressure gradient term */
      // + (jump(v), avg(δp) n)
      u_value_result += p_value_avg;

      /* Contribution from velocity divergence term */
      // - (avg(q), jump(u) n)
      p_value_result -= 0.5 * contract(u_value_jump, normal);

      // - (n avg(∇v), ν/2 jump(δu))
      phi_u_m.submit_gradient(outer(u_grad_result, normal) *
                                this->kinematic_viscosity * 0.5 * JxW,
                              q);
      phi_u_m.submit_value(u_value_result * JxW, q);
      phi_p_m.submit_value(p_value_result * JxW, q);
    }

  AssertDimension(buffer.size(),
                  fe_sub_u.n_dofs_per_cell() + fe_sub_p.n_dofs_per_cell());

  ArrayView<Number> buffer_u(buffer.data() + 0, fe_sub_u.n_dofs_per_cell());
  ArrayView<Number> buffer_p(buffer.data() + fe_sub_u.n_dofs_per_cell(),
                             fe_sub_p.n_dofs_per_cell());

  this->phi_u_m.test_and_sum(buffer_u,
                             EvaluationFlags::values |
                               EvaluationFlags::gradients);
  this->phi_p_m.test_and_sum(buffer_p, EvaluationFlags::values);
}


/*-------------- Explicit Instantiations -------------------------------*/
template class MortarManagerBase<1>;
template class MortarManagerBase<2>;
template class MortarManagerBase<3>;

template class MortarManagerCircle<1>;
template class MortarManagerCircle<2>;
template class MortarManagerCircle<3>;

template class MortarManagerLinear<1>;
template class MortarManagerLinear<2>;
template class MortarManagerLinear<3>;

template class CouplingOperator<1, double>;
template class CouplingOperator<2, double>;
template class CouplingOperator<3, double>;

template void
CouplingOperator<1, double>::vmult_add(
  LinearAlgebra::distributed::Vector<double> &,
  const LinearAlgebra::distributed::Vector<double> &) const;
template void
CouplingOperator<2, double>::vmult_add(
  LinearAlgebra::distributed::Vector<double> &,
  const LinearAlgebra::distributed::Vector<double> &) const;
template void
CouplingOperator<3, double>::vmult_add(
  LinearAlgebra::distributed::Vector<double> &,
  const LinearAlgebra::distributed::Vector<double> &) const;

template void
CouplingOperator<1, double>::vmult_add(
  TrilinosWrappers::MPI::Vector &,
  const TrilinosWrappers::MPI::Vector &) const;
template void
CouplingOperator<2, double>::vmult_add(
  TrilinosWrappers::MPI::Vector &,
  const TrilinosWrappers::MPI::Vector &) const;
template void
CouplingOperator<3, double>::vmult_add(
  TrilinosWrappers::MPI::Vector &,
  const TrilinosWrappers::MPI::Vector &) const;

template void
CouplingOperator<2, double>::vmult_add(
  LinearAlgebra::distributed::Vector<float> &,
  const LinearAlgebra::distributed::Vector<float> &) const;

template void
CouplingOperator<3, double>::vmult_add(
  LinearAlgebra::distributed::Vector<float> &,
  const LinearAlgebra::distributed::Vector<float> &) const;

template void
CouplingOperator<1, double>::add_diagonal_entries(
  LinearAlgebra::distributed::Vector<double> &) const;
template void
CouplingOperator<2, double>::add_diagonal_entries(
  LinearAlgebra::distributed::Vector<double> &) const;
template void
CouplingOperator<3, double>::add_diagonal_entries(
  LinearAlgebra::distributed::Vector<double> &) const;

template void
CouplingOperator<2, double>::add_diagonal_entries(
  LinearAlgebra::distributed::Vector<float> &) const;
template void
CouplingOperator<3, double>::add_diagonal_entries(
  LinearAlgebra::distributed::Vector<float> &) const;

template class CouplingEvaluationSIPG<1, 1, double>;
template class CouplingEvaluationSIPG<1, 2, double>;
template class CouplingEvaluationSIPG<2, 1, double>;
template class CouplingEvaluationSIPG<2, 2, double>;
template class CouplingEvaluationSIPG<2, 3, double>;
template class CouplingEvaluationSIPG<3, 1, double>;
template class CouplingEvaluationSIPG<3, 3, double>;
template class CouplingEvaluationSIPG<3, 4, double>;

template class NavierStokesCouplingEvaluation<2, double>;
template class NavierStokesCouplingEvaluation<3, double>;

template std::vector<unsigned int>
compute_number_interface_cells<2>(
  const Triangulation<2>      &triangulation,
  const Parameters::Mortar<2> &mortar_parameters);

template std::vector<unsigned int>
compute_number_interface_cells<3>(
  const Triangulation<3>      &triangulation,
  const Parameters::Mortar<3> &mortar_parameters);

template std::tuple<std::vector<double>, double>
compute_interface_dimensions_circular<2>(
  const Triangulation<2>      &triangulation,
  const Mapping<2>            &mapping,
  const Parameters::Mortar<2> &mortar_parameters);

template std::tuple<std::vector<double>, double>
compute_interface_dimensions_circular<3>(
  const Triangulation<3>      &triangulation,
  const Mapping<3>            &mapping,
  const Parameters::Mortar<3> &mortar_parameters);

template std::pair<double, double>
compute_interface_dimensions_linear<2>(
  const Triangulation<2>      &triangulation,
  const Mapping<2>            &mapping,
  const Parameters::Mortar<2> &mortar_parameters);

template std::pair<double, double>
compute_interface_dimensions_linear<3>(
  const Triangulation<3>      &triangulation,
  const Mapping<3>            &mapping,
  const Parameters::Mortar<3> &mortar_parameters);

template Quadrature<2>
construct_quadrature(const Quadrature<2>         &quadrature,
                     const Parameters::Mortar<2> &mortar_parameters);

template Quadrature<3>
construct_quadrature(const Quadrature<3>         &quadrature,
                     const Parameters::Mortar<3> &mortar_parameters);

template std::vector<double>
compute_stage_heights<2>(const Triangulation<2>      &triangulation,
                         const Parameters::Mortar<2> &mortar_parameters);

template std::vector<double>
compute_stage_heights<3>(const Triangulation<3>      &triangulation,
                         const Parameters::Mortar<3> &mortar_parameters);

template void
mortar_workload_imbalance(const Triangulation<2>      &triangulation,
                          const Parameters::Mortar<2> &mortar_parameters,
                          const ConditionalOStream    &pcout);

template void
mortar_workload_imbalance(const Triangulation<3>      &triangulation,
                          const Parameters::Mortar<3> &mortar_parameters,
                          const ConditionalOStream    &pcout);
