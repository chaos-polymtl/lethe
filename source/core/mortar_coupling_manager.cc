// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/config.h>

#include <core/mortar_coupling_manager.h>

#include <deal.II/base/mpi_noncontiguous_partitioner.templates.h>

#include <deal.II/fe/fe_nothing.h>


/*-------------- MortarManagerBase -------------------------------*/

template <int dim>
bool
MortarManagerBase<dim>::is_mesh_aligned() const
{
  AssertThrow(dim != 1, ExcInternalError());

  const double tolerance = 1e-8;
  const double delta_0   = 2 * numbers::PI / n_subdivisions[0];

  return std::abs(rotation_angle / delta_0 -
                  std::round(rotation_angle / delta_0)) < tolerance;
}

template <int dim>
unsigned int
MortarManagerBase<dim>::get_n_total_mortars() const
{
  if constexpr (dim == 1)
    return 1;

  unsigned int n_total_subdivisions = n_subdivisions[0];

  // In 3D, besides the subdivisons in the plane perpendicular to the rotation
  // axis, we also need to account for the number of subdivisons along the
  // rotation axis direction
  if constexpr (dim == 3)
    n_total_subdivisions *= n_subdivisions[1];

  if (this->is_mesh_aligned()) // aligned
    {
      return n_total_subdivisions;
    }
  else // inside/outside
    {
      return 2 * n_total_subdivisions;
    }
}

template <int dim>
unsigned int
MortarManagerBase<dim>::get_n_mortars() const
{
  if constexpr (dim == 1)
    return 1;

  if (this->is_mesh_aligned()) // aligned
    {
      return 1;
    }
  else // inside/outside
    {
      return 2;
    }
}

template <int dim>
std::vector<unsigned int>
MortarManagerBase<dim>::get_mortar_indices(const Point<dim> &face_center) const
{
  if constexpr (dim == 1)
    return std::vector<unsigned int>{0};

  // Mesh alignment type and cell indexes
  const auto [type, id_in_plane, id_out_plane] = get_config(face_center);

  if (type == 0) // aligned
    {
      std::vector<unsigned int> indices;

      const unsigned int index = id_in_plane;

      AssertIndexRange(index, n_subdivisions[0]);

      indices.emplace_back(index + n_subdivisions[0] * id_out_plane);

      return indices;
    }
  else if (type == 1) // inside
    {
      std::vector<unsigned int> indices;

      for (unsigned int q = 0; q < 2; ++q)
        {
          const unsigned int index =
            (id_in_plane * 2 + 1 + q) % (n_subdivisions[0] * 2);

          AssertIndexRange(index, n_subdivisions[0] * 2);

          indices.emplace_back(index + n_subdivisions[0] * id_out_plane);
        }

      return indices;
    }
  else // outside
    {
      std::vector<unsigned int> indices;

      for (unsigned int q = 0; q < 2; ++q)
        {
          const unsigned int index = id_in_plane * 2 + q;

          AssertIndexRange(index, n_subdivisions[0] * 2);

          indices.emplace_back(index + n_subdivisions[0] * id_out_plane);
        }

      return indices;
    }
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
MortarManagerBase<dim>::get_n_points() const
{
  if constexpr (dim == 1)
    return 1;

  return get_n_mortars() * n_quadrature_points;
}

template <int dim>
std::vector<Point<dim>>
MortarManagerBase<dim>::get_points(const Point<dim> &face_center) const
{
  if constexpr (dim == 1)
    return std::vector<Point<dim>>{face_center};

  // Mesh alignment type and cell index
  const auto [type, id_in_plane, id_out_plane] = get_config(face_center);
  // Angle variation within each cell
  const double delta_0 = 2 * numbers::PI / n_subdivisions[0];
  double       delta_1 = 1.0;

  if constexpr (dim == 3)
    delta_1 = radius[1] / n_subdivisions[1];

  if (type == 0) // aligned
    {
      std::vector<Point<dim>> points;

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        {
          const auto x =
            from_1D((id_in_plane + quadrature.point(q)[0]) * delta_0);

          if constexpr (dim == 3)
            points.emplace_back(
              x[0], x[1], (id_out_plane + quadrature.point(q)[1]) * delta_1);
          else
            points.emplace_back(x);
        }

      return points;
    }
  else // Point at the inner boundary lies somewhere in the face of the outer
       // boundary cell
    {
      // rad_0: first cell vertex (fixed)
      // rad_1: shifted vertex
      // rad_2: last cell vertex (fixed)
      double rad_0, rad_1, rad_2;
      // Minimum rotation angle
      double rot_min =
        rotation_angle - std::floor(rotation_angle / delta_0) * delta_0;

      if (type == 2) // outside
        {
          rad_0 = id_in_plane * delta_0;
          rad_1 = id_in_plane * delta_0 + rot_min;
          rad_2 = (id_in_plane + 1) * delta_0;
        }
      else // inside
        {
          rad_0 = id_in_plane * delta_0 + rot_min;
          rad_1 = (id_in_plane + 1) * delta_0;
          rad_2 = (id_in_plane + 1) * delta_0 + rot_min;
        }

      std::vector<Point<dim>> points;

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        {
          const auto x =
            from_1D(rad_0 + quadrature.point(q)[0] * (rad_1 - rad_0));

          if constexpr (dim == 3)
            points.emplace_back(
              x[0], x[1], (id_out_plane + quadrature.point(q)[1]) * delta_1);
          else
            points.emplace_back(x);
        }

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        {
          const auto x =
            from_1D(rad_1 + quadrature.point(q)[0] * (rad_2 - rad_1));

          if constexpr (dim == 3)
            points.emplace_back(
              x[0], x[1], (id_out_plane + quadrature.point(q)[1]) * delta_1);
          else
            points.emplace_back(x);
        }

      return points;
    }
}

template <int dim>
std::vector<Point<std::max(1, dim - 1)>>
MortarManagerBase<dim>::get_points_ref(const Point<dim> &face_center) const
{
  if (dim == 1)
    return std::vector<Point<std::max(1, dim - 1)>>{
      Point<std::max(1, dim - 1)>()};

  const auto [type, _, __] = get_config(face_center);

  const double delta_0 = 2 * numbers::PI / n_subdivisions[0];

  if (type == 0) // aligned
    {
      std::vector<Point<std::max(1, dim - 1)>> points;

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        points.emplace_back(quadrature.point(q));

      return points;
    }
  else // inside/outside
    {
      double rad_0, rad_1, rad_2;

      double rot_min =
        (rotation_angle - std::floor(rotation_angle / delta_0) * delta_0) /
        delta_0;

      if (type == 2) // outside
        {
          rad_0 = 0.0;
          rad_1 = rot_min;
          rad_2 = 1.0;
        }
      else // inside
        {
          rad_0 = 0.0;
          rad_1 = 1.0 - rot_min;
          rad_2 = 1.0;
        }

      std::vector<Point<std::max(1, dim - 1)>> points;

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        {
          const double x = rad_0 + quadrature.point(q)[0] * (rad_1 - rad_0);

          if (dim == 2)
            points.emplace_back(x);
          else
            points.emplace_back(x, quadrature.point(q)[1]);
        }

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        {
          const double x = rad_1 + quadrature.point(q)[0] * (rad_2 - rad_1);

          if (dim == 2)
            points.emplace_back(x);
          else
            points.emplace_back(x, quadrature.point(q)[1]);
        }

      return points;
    }
}

template <int dim>
std::vector<double>
MortarManagerBase<dim>::get_weights(const Point<dim> &face_center) const
{
  if (dim == 1)
    return std::vector<double>{1.0};

  // Mesh alignment type and cell index
  const auto [type, id_in_plane, _] = get_config(face_center);
  // Angle variation within each cell
  const double delta_0 = 2 * numbers::PI / n_subdivisions[0];
  double       delta_1 = 1.0;

  if (dim == 3)
    delta_1 = radius[1] / n_subdivisions[1];

  if (type == 0) // aligned
    {
      std::vector<double> weights;

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        weights.emplace_back(radius[0] * quadrature.weight(q) * delta_0 *
                             delta_1);

      return weights;
    }
  else // inside/outside
    {
      double rad_0, rad_1, rad_2;

      double rot_min =
        rotation_angle - std::floor(rotation_angle / delta_0) * delta_0;

      if (type == 2) // outside
        {
          rad_0 = id_in_plane * delta_0;
          rad_1 = id_in_plane * delta_0 + rot_min;
          rad_2 = (id_in_plane + 1) * delta_0;
        }
      else // inside
        {
          rad_0 = id_in_plane * delta_0 + rot_min;
          rad_1 = (id_in_plane + 1) * delta_0;
          rad_2 = (id_in_plane + 1) * delta_0 + rot_min;
        }

      std::vector<double> weights;

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        weights.emplace_back(radius[0] * quadrature.weight(q) *
                             (rad_1 - rad_0) * delta_1);

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        weights.emplace_back(radius[0] * quadrature.weight(q) *
                             (rad_2 - rad_1) * delta_1);

      return weights;
    }
}

template <int dim>
std::vector<Tensor<1, dim, double>>
MortarManagerBase<dim>::get_normals(const Point<dim> &face_center) const
{
  // Coordinates of cell quadrature points
  const auto points = get_points(face_center);

  std::vector<Tensor<1, dim, double>> result;

  result.reserve(points.size());

  for (const auto &point : points)
    result.emplace_back(get_normal(point));

  return result;
}

template <int dim>
std::tuple<unsigned int, unsigned int, unsigned int>
MortarManagerBase<dim>::get_config(const Point<dim> &face_center) const
{
  const auto angle_cell_center = to_1D(face_center);

  // Alignment tolerance
  const double tolerance = 1e-8;
  // Angular variation in each cell
  const double delta_0 = 2 * numbers::PI / n_subdivisions[0];
  // Minimum rotation angle
  const double rot_min =
    rotation_angle - std::floor(rotation_angle / delta_0) * delta_0;

  AssertThrow(rot_min <= delta_0, ExcInternalError());

  // Point position in the cell
  const double segment = (angle_cell_center - delta_0 / 2) / delta_0;
  // Point position after rotation
  const double segment_rot =
    (angle_cell_center - delta_0 / 2 - rot_min) / delta_0;
  // Cell index in the direction of the rotation axis
  unsigned int id_out_plane = 0;

  if constexpr (dim == 3)
    {
      const double delta_1 = radius[1] / n_subdivisions[1];
      id_out_plane         = static_cast<unsigned int>(
        std::round((face_center[2] - delta_1 / 2) / delta_1));
    }

  if (this->is_mesh_aligned())
    {
      // Case 1: mesh is aligned
      return {0, static_cast<unsigned int>(std::round(segment)), id_out_plane};
    }
  else
    {
      // Case 2: mesh is not aligned
      if (std::abs(segment - std::round(segment)) < tolerance)
        // outer (fixed) domain
        return {2,
                static_cast<unsigned int>(std::round(segment)),
                id_out_plane};
      else
        // inner (rotated) domain
        return {1,
                // static_cast<unsigned int>(std::round(segment_rot)) %
                //   (2 * n_subdivisions[0]),
                // id_out_plane};
                (static_cast<unsigned int>(std::round(segment_rot)) +
                 2 * n_subdivisions[0]) %
                  (2 * n_subdivisions[0]),
                id_out_plane};
    }
}


/*-------------- Auxiliary Functions -------------------------------*/

#if DEAL_II_VERSION_GTE(9, 7, 0)
template <int dim>
std::tuple<unsigned int, double, double>
compute_n_subdivisions_and_radius(
  const Triangulation<dim>      &triangulation,
  const Mapping<dim>            &mapping,
  const Parameters::Mortar<dim> &mortar_parameters)
{
  // Number of subdivisions per process
  unsigned int n_subdivisions_local = 0;
  // Number of vertices at the boundary per process
  unsigned int n_vertices_local = 0;
  // Tolerance for rotor radius computation
  const double tolerance = 1e-8;
  // Min and max values for rotor radius computation
  double radius_min = std::numeric_limits<double>::max();
  double radius_max = 0;
  // Minimum rotation angle in initial mesh configuration
  double pre_rotation_min = std::numeric_limits<double>::max();

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

                      const auto vertices = mapping.get_vertices(cell, face_no);

                      for (unsigned int vertex_index = 0;
                           vertex_index < face->n_vertices();
                           vertex_index++)
                        {
                          n_vertices_local++;
                          const auto v = vertices[vertex_index];
                          double     radius_current =
                            v.distance(mortar_parameters.center_of_rotation);

                          // In 3D, the interface radius to be computed is
                          // actually with respect to the rotation axis. For
                          // this computation, we use the relation:
                          // radius_current = ||VC x d||/||d||, where VC is the
                          // current vertex minus the center of rotation, and d
                          // is the rotation axis
                          if constexpr (dim == 3)
                            {
                              const auto aux = cross_product_3d(
                                (v - mortar_parameters.center_of_rotation),
                                mortar_parameters.rotation_axis);
                              const double num = aux.norm();
                              const double den =
                                mortar_parameters.rotation_axis.norm();
                              radius_current = num / den;
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
                          n_vertices_local++;
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

  // Total number of faces
  const unsigned int n_subdivisions =
    Utilities::MPI::sum(n_subdivisions_local,
                        triangulation.get_mpi_communicator());

  // Min and max values over all processes
  radius_min =
    Utilities::MPI::min(radius_min, triangulation.get_mpi_communicator());
  radius_max =
    Utilities::MPI::max(radius_max, triangulation.get_mpi_communicator());

  pre_rotation_min =
    Utilities::MPI::min(pre_rotation_min, triangulation.get_mpi_communicator());

  AssertThrow(
    std::abs(radius_max - radius_min) < tolerance,
    ExcMessage(
      "The computed radius of the rotor mesh has a variation greater than "
      "the tolerance across the rotor domain, meaning that the prescribed "
      "center of rotation and the rotor geometry are not in accordance."));

  // Final radius value
  const double radius = radius_min;

  return {n_subdivisions, radius, pre_rotation_min};
}
#else
template <int dim>
std::tuple<unsigned int, double, double>
compute_n_subdivisions_and_radius(const Triangulation<dim> &,
                                  const Mapping<dim> &,
                                  const Parameters::Mortar<dim> &)
{
  AssertThrow(false,
              ExcMessage(
                "The mortar coupling requires deal.II 9.7 or more recent."));

  return {1, 1.0, 0.0};
}
#endif

template <int dim>
Quadrature<dim>
construct_quadrature(const Quadrature<dim>         &quadrature,
                     const Parameters::Mortar<dim> &mortar_parameters)
{
  const double oversampling_factor = mortar_parameters.oversampling_factor;

  for (unsigned int i = 1; i <= 10; ++i)
    if (quadrature == QGauss<dim>(i))
      return QGauss<dim>(i * oversampling_factor);

  AssertThrow(false, ExcNotImplemented());

  return quadrature;
}


/*-------------- MortarManagerCircle -------------------------------*/

template <int dim>
Point<dim>
MortarManagerCircle<dim>::from_1D(const double radiant) const
{
  return radius_to_point<dim>(this->radius[0], radiant + pre_rotation_angle);
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



/*-------------- CouplingOperator -------------------------------*/
#if DEAL_II_VERSION_GTE(9, 7, 0)
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
  , constraints(constraints)
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

#  ifdef DEBUG
  std::vector<double> vec_local_cells(n_sub_cells * 2, 0.0);
  std::vector<double> vec_ghost_cells(n_sub_cells * 2, 0.0);
#  endif

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto face_no : cell->face_indices())
        if ((cell->face(face_no)->boundary_id() == bid_m) ||
            (cell->face(face_no)->boundary_id() == bid_p))
          {
            const auto face = cell->face(face_no);

            // Indices of mortars on face of cell.
            const auto indices =
              mortar_manager->get_mortar_indices(get_face_center(cell, face));

            const auto local_dofs = this->get_dof_indices(cell);

            // Loop over the mortar indices at the rotor-stator
            // interface. The logic of local (rotor) and ghost (stator) is the
            // same as in the previous loop.
            for (unsigned int ii = 0; ii < indices.size(); ++ii)
              {
                unsigned int i = indices[ii];
                unsigned int id_local, id_ghost;

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

#  ifdef DEBUG
                vec_local_cells[id_local] += 1.0;
                vec_ghost_cells[id_ghost] += 1.0;
#  endif
                is_local_cell.emplace_back(id_local);
                is_ghost_cell.emplace_back(id_ghost);

                for (const auto i : local_dofs)
                  dof_indices.emplace_back(i);
              }

            // Weights of quadrature points
            const auto weights =
              mortar_manager->get_weights(get_face_center(cell, face));
            data.all_weights.insert(data.all_weights.end(),
                                    weights.begin(),
                                    weights.end());

            // Normals of quadrature points
            if constexpr (dim == 3)
              {
                const auto points =
                  mortar_manager->get_points(get_face_center(cell, face));
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
                    for (unsigned int i = 0, j = 0; i < dim; ++i)
                      if ((face_no / 2) != i)
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
                auto normals =
                  mortar_manager->get_normals(get_face_center(cell, face));
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
                      get_face_center(cell, face));

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

            // Penalty parmeter
            const Number penalty_parameter = compute_penalty_parameter(cell);

            // Store penalty parameter for all quadrature points
            for (unsigned int i = 0; i < mortar_manager->get_n_points(); ++i)
              data.all_penalty_parameter.emplace_back(penalty_parameter);
          }

#  ifdef DEBUG
  Utilities::MPI::sum(vec_local_cells,
                      dof_handler.get_mpi_communicator(),
                      vec_local_cells);
  Utilities::MPI::sum(vec_ghost_cells,
                      dof_handler.get_mpi_communicator(),
                      vec_ghost_cells);

  std::set<unsigned int> vec_local_cells_0;
  std::set<unsigned int> vec_local_cells_2;

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
#  endif

  // Setup communication
  partitioner.reinit(is_local_cell,
                     is_ghost_cell,
                     dof_handler.get_mpi_communicator());

  // Finalized penalty parameters
  const unsigned n_q_points =
    mortar_manager->get_n_points() / mortar_manager->get_n_mortars();
  std::vector<Number> all_penalty_parameter_ghost(
    data.all_penalty_parameter.size());
  partitioner.template export_to_ghosted_array<Number, 0>(
    data.all_penalty_parameter, all_penalty_parameter_ghost, n_q_points);
  for (unsigned int i = 0; i < data.all_penalty_parameter.size(); ++i)
    data.all_penalty_parameter[i] =
      std::max(data.all_penalty_parameter[i], all_penalty_parameter_ghost[i]);

  // Finialize DoF indices and update constraints
  dof_indices_ghost.resize(dof_indices.size());
  partitioner.template export_to_ghosted_array<types::global_dof_index, 0>(
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
  }
}
#else
template <int dim, typename Number>
CouplingOperator<dim, Number>::CouplingOperator(
  const Mapping<dim>                                        &mapping,
  const DoFHandler<dim>                                     &dof_handler,
  const AffineConstraints<Number>                           &constraints,
  const std::shared_ptr<CouplingEvaluationBase<dim, Number>> evaluator,
  const std::shared_ptr<MortarManagerBase<dim>>              mortar_manager,
  const unsigned int                                         bid_m,
  const unsigned int                                         bid_p,
  const double)
  : mapping(mapping)
  , dof_handler(dof_handler)
  , constraints(constraints)
  , bid_m(bid_m)
  , bid_p(bid_p)
  , evaluator(evaluator)
  , mortar_manager(mortar_manager)
{
  AssertThrow(false,
              ExcMessage(
                "The mortar coupling requires deal.II 9.7 or more recent."));
}
#endif

template <int dim, typename Number>
const AffineConstraints<Number> &
CouplingOperator<dim, Number>::get_affine_constraints() const
{
  return constraints_extended;
}

template <int dim, typename Number>
Number
CouplingOperator<dim, Number>::compute_penalty_factor(const unsigned int degree,
                                                      const Number factor) const
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

#if DEAL_II_VERSION_GTE(9, 7, 0)
template <int dim, typename Number>
template <typename VectorType>
void
CouplingOperator<dim, Number>::vmult_add(VectorType       &dst,
                                         const VectorType &src) const
{
  // 1) Evaluate
  unsigned int ptr_q = 0;

  Vector<Number> buffer;

  std::vector<Number> all_values_local(data.all_normals.size() * q_data_size);
  std::vector<Number> all_values_ghost(data.all_normals.size() * q_data_size);

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
          {
            // Quadrature points at the current face. Note: we process
            // all mortars of the face together here.
            const unsigned int n_q_points = mortar_manager->get_n_points();

            evaluator->local_reinit(
              cell,
              ArrayView<const Point<dim, Number>>(all_points_ref.data() + ptr_q,
                                                  n_q_points));

            buffer.reinit(n_dofs_per_cell);

            const auto local_dofs = this->get_dof_indices(cell);

            for (unsigned int i = 0; i < local_dofs.size(); ++i)
              buffer[i] = src[local_dofs[i]];

            evaluator->local_evaluate(data,
                                      buffer,
                                      ptr_q,
                                      1,
                                      all_values_local.data() +
                                        ptr_q * q_data_size);

            ptr_q += n_q_points;
          }

  // 2) Communicate
  const unsigned n_q_points =
    mortar_manager->get_n_points() / mortar_manager->get_n_mortars();
  partitioner.template export_to_ghosted_array<Number, 0>(
    ArrayView<const Number>(reinterpret_cast<Number *>(all_values_local.data()),
                            all_values_local.size()),
    ArrayView<Number>(reinterpret_cast<Number *>(all_values_ghost.data()),
                      all_values_ghost.size()),
    n_q_points * q_data_size);

  // 3) Integrate
  ptr_q = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
          {
            // Quadrature points at the current face. Note: we process
            // all mortars of the face together here.
            const unsigned int n_q_points = mortar_manager->get_n_points();

            evaluator->local_reinit(
              cell,
              ArrayView<const Point<dim, Number>>(all_points_ref.data() + ptr_q,
                                                  n_q_points));

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
            constraints.distribute_local_to_global(buffer, local_dofs, dst);

            ptr_q += n_q_points;
          }
  dst.compress(VectorOperation::add);
}
#else
template <int dim, typename Number>
template <typename VectorType>
void
CouplingOperator<dim, Number>::vmult_add(VectorType &, const VectorType &) const
{}
#endif

template <int dim, typename Number>
template <typename VectorType>
void
CouplingOperator<dim, Number>::add_diagonal_entries(VectorType &diagonal) const
{
  unsigned int ptr_q = 0;

  Vector<Number>      buffer, diagonal_local;
  std::vector<Number> all_values_local, all_values_ghost;

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
          {
            // Quadrature points at the current face. Note: we process
            // all mortars of the face together here.
            const unsigned int n_q_points = mortar_manager->get_n_points();

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
            constraints.distribute_local_to_global(diagonal_local,
                                                   local_dofs,
                                                   diagonal);

            ptr_q += n_q_points;
          }

  diagonal.compress(VectorOperation::add);
}

template <int dim, typename Number>
void
CouplingOperator<dim, Number>::add_sparsity_pattern_entries(
  SparsityPatternBase &dsp) const
{
  const auto constraints = &constraints_extended;

  for (unsigned int i = 0; i < dof_indices.size(); i += n_dofs_per_cell)
    {
      std::vector<types::global_dof_index> a(dof_indices.begin() + i,
                                             dof_indices.begin() + i +
                                               n_dofs_per_cell);
      std::vector<types::global_dof_index> b(dof_indices_ghost.begin() + i,
                                             dof_indices_ghost.begin() + i +
                                               n_dofs_per_cell);

      constraints->add_entries_local_to_global(a, b, dsp);
      constraints->add_entries_local_to_global(b, a, dsp);
    }
}

#if DEAL_II_VERSION_GTE(9, 7, 0)
template <int dim, typename Number>
void
CouplingOperator<dim, Number>::add_system_matrix_entries(
  TrilinosWrappers::SparseMatrix &system_matrix) const
{
  const auto constraints = &constraints_extended;

  std::vector<Number> all_values_local(data.all_normals.size() *
                                       n_dofs_per_cell * q_data_size);
  std::vector<Number> all_values_ghost(data.all_normals.size() *
                                       n_dofs_per_cell * q_data_size);

  unsigned int ptr_q = 0;

  Vector<Number> buffer;

  // 1) Evaluate
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
          {
            // Quadrature points at the current face. Note: we process
            // all mortars of the face together here.
            const unsigned int n_q_points = mortar_manager->get_n_points();

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
          }

  const unsigned n_q_points =
    mortar_manager->get_n_points() / mortar_manager->get_n_mortars();

  // 2) Communicate: export data from local to ghost side
  partitioner.template export_to_ghosted_array<Number, 0>(
    ArrayView<const Number>(reinterpret_cast<Number *>(all_values_local.data()),
                            all_values_local.size()),
    ArrayView<Number>(reinterpret_cast<Number *>(all_values_ghost.data()),
                      all_values_ghost.size()),
    n_dofs_per_cell * n_q_points * q_data_size);


  ptr_q                 = 0;
  unsigned int ptr_dofs = 0;

  // 3) Integrate
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
          {
            // Number of mortars attached to the current cell (i.e. 1 for
            // aligned rotor-stator meshes and 2 for non-aligned case)
            const unsigned int n_mortars = mortar_manager->get_n_mortars();

            // Loop over mortars
            for (unsigned int m = 0; m < n_mortars; ++m)
              {
                evaluator->local_reinit(cell,
                                        ArrayView<const Point<dim, Number>>(
                                          all_points_ref.data() + ptr_q,
                                          n_q_points));

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
                        constraints->distribute_local_to_global(
                          cell_matrix, local_dof_indices, system_matrix);
                      }
                    else // ghost cell -> local-ghost block
                      {
                        std::vector<types::global_dof_index>
                          local_dof_indices_ghost(dof_indices_ghost.begin() +
                                                    ptr_dofs,
                                                  dof_indices_ghost.begin() +
                                                    ptr_dofs + n_dofs_per_cell);

                        constraints->distribute_local_to_global(
                          cell_matrix,
                          local_dof_indices,
                          local_dof_indices_ghost,
                          system_matrix);
                      }
                  }

                ptr_dofs += n_dofs_per_cell;

                ptr_q += n_q_points;
              }
          }

  AssertDimension(ptr_q, data.all_normals.size());
  AssertDimension(ptr_dofs, dof_indices.size());
}
#else
template <int dim, typename Number>
void
CouplingOperator<dim, Number>::add_system_matrix_entries(
  TrilinosWrappers::SparseMatrix &) const
{
  AssertThrow(false,
              ExcMessage(
                "The mortar coupling requires deal.II 9.7 or more recent."));
}
#endif


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
CouplingOperator<1, double>::add_diagonal_entries(
  LinearAlgebra::distributed::Vector<double> &) const;
template void
CouplingOperator<2, double>::add_diagonal_entries(
  LinearAlgebra::distributed::Vector<double> &) const;
template void
CouplingOperator<3, double>::add_diagonal_entries(
  LinearAlgebra::distributed::Vector<double> &) const;

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

template std::tuple<unsigned int, double, double>
compute_n_subdivisions_and_radius<2>(
  const Triangulation<2>      &triangulation,
  const Mapping<2>            &mapping,
  const Parameters::Mortar<2> &mortar_parameters);

template std::tuple<unsigned int, double, double>
compute_n_subdivisions_and_radius<3>(
  const Triangulation<3>      &triangulation,
  const Mapping<3>            &mapping,
  const Parameters::Mortar<3> &mortar_parameters);

template Quadrature<2>
construct_quadrature(const Quadrature<2>         &quadrature,
                     const Parameters::Mortar<2> &mortar_parameters);

template Quadrature<3>
construct_quadrature(const Quadrature<3>         &quadrature,
                     const Parameters::Mortar<3> &mortar_parameters);
