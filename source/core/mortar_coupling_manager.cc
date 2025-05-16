// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/config.h>

#if DEAL_II_VERSION_GTE(9, 7, 0)

#  include <core/mortar_coupling_manager.h>


/*-------------- MortarManager -------------------------------*/

template <int dim>
MortarManager<dim>::MortarManager(const unsigned int n_subdivisions,
                                  const unsigned int n_quadrature_points,
                                  const double       radius,
                                  const double       rotation_angle)
  : n_subdivisions(n_subdivisions)
  , n_quadrature_points(n_quadrature_points)
  , radius(radius)
  , rotation_angle(rotation_angle)
  , quadrature(n_quadrature_points)
{}

template <int dim>
bool
MortarManager<dim>::is_mesh_aligned() const
{
  const double tolerance = 1e-8;
  const double delta     = 2 * numbers::PI / n_subdivisions;

  return std::abs(rotation_angle / delta - std::round(rotation_angle / delta)) <
         tolerance;
}

template <int dim>
unsigned int
MortarManager<dim>::get_n_points() const
{
  if (this->is_mesh_aligned()) // aligned
    {
      return n_subdivisions * n_quadrature_points;
    }
  else // inside/outside
    {
      return 2 * n_subdivisions * n_quadrature_points;
    }
}

template <int dim>
unsigned int
MortarManager<dim>::get_n_points(const double &angle_cell_center) const
{
  (void)angle_cell_center;

  if (this->is_mesh_aligned()) // aligned
    {
      return n_quadrature_points;
    }
  else // inside/outside
    {
      return 2 * n_quadrature_points;
    }
}

template <int dim>
std::vector<unsigned int>
MortarManager<dim>::get_indices(const double &angle_cell_center) const
{
  // Mesh alignment type and cell index
  const auto [type, id] = get_config(angle_cell_center);

  if (type == 0) // aligned
    {
      std::vector<unsigned int> indices;

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        {
          const unsigned int index = id * n_quadrature_points + q;

          AssertIndexRange(index, get_n_points());

          indices.emplace_back(index);
        }

      return indices;
    }
  else if (type == 1) // inside
    {
      std::vector<unsigned int> indices;

      for (unsigned int q = 0; q < n_quadrature_points * 2; ++q)
        {
          const unsigned int index =
            (id * n_quadrature_points * 2 + n_quadrature_points + q) %
            get_n_points();

          AssertIndexRange(index, get_n_points());

          indices.emplace_back(index);
        }

      return indices;
    }
  else // outside
    {
      std::vector<unsigned int> indices;

      for (unsigned int q = 0; q < n_quadrature_points * 2; ++q)
        {
          const unsigned int index = id * n_quadrature_points * 2 + q;

          AssertIndexRange(index, get_n_points());

          indices.emplace_back(index);
        }

      return indices;
    }
}

template <int dim>
std::vector<Point<dim>>
MortarManager<dim>::get_points(const double rad) const
{
  // Mesh alignment type and cell index
  const auto [type, id] = get_config(rad);
  // Angle variation within each cell
  const double delta = 2 * numbers::PI / n_subdivisions;

  if (type == 0) // aligned
    {
      std::vector<Point<dim>> points;

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        points.emplace_back(
          radius_to_point<dim>(radius, (id + quadrature.point(q)[0]) * delta));

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
        rotation_angle - std::floor(rotation_angle / delta) * delta;

      if (type == 2) // outside
        {
          rad_0 = id * delta;
          rad_1 = id * delta + rot_min;
          rad_2 = (id + 1) * delta;
        }
      else // inside
        {
          rad_0 = id * delta + rot_min;
          rad_1 = (id + 1) * delta;
          rad_2 = (id + 1) * delta + rot_min;
        }

      std::vector<Point<dim>> points;

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        points.emplace_back(radius_to_point<dim>(
          radius, rad_0 + quadrature.point(q)[0] * (rad_1 - rad_0)));

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        points.emplace_back(radius_to_point<dim>(
          radius, rad_1 + quadrature.point(q)[0] * (rad_2 - rad_1)));

      return points;
    }
}

template <int dim>
std::vector<Point<1>>
MortarManager<dim>::get_points_ref(const double angle_cell_center) const
{
  const auto [type, id] = get_config(angle_cell_center);

  const double delta = 2 * numbers::PI / n_subdivisions;

  if (type == 0) // aligned
    {
      std::vector<Point<1>> points;

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        points.emplace_back(quadrature.point(q));

      return points;
    }
  else // inside/outside
    {
      double rad_0, rad_1, rad_2;

      double rot_min =
        (rotation_angle - std::floor(rotation_angle / delta) * delta) / delta;

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

      std::vector<Point<1>> points;

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        points.emplace_back(rad_0 + quadrature.point(q)[0] * (rad_1 - rad_0));

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        points.emplace_back(rad_1 + quadrature.point(q)[0] * (rad_2 - rad_1));

      return points;
    }
}

template <int dim>
std::vector<double>
MortarManager<dim>::get_weights(const double &angle_cell_center) const
{
  // Mesh alignment type and cell index
  const auto [type, id] = get_config(angle_cell_center);
  // Angle variation within each cell
  const double delta = 2 * numbers::PI / n_subdivisions;

  if (type == 0) // aligned
    {
      std::vector<double> points;

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        points.emplace_back(radius * quadrature.weight(q) * delta);

      return points;
    }
  else // inside/outside
    {
      double rad_0, rad_1, rad_2;

      double rot_min =
        rotation_angle - std::floor(rotation_angle / delta) * delta;

      if (type == 2) // outside
        {
          rad_0 = id * delta;
          rad_1 = id * delta + rot_min;
          rad_2 = (id + 1) * delta;
        }
      else // inside
        {
          rad_0 = id * delta + rot_min;
          rad_1 = (id + 1) * delta;
          rad_2 = (id + 1) * delta + rot_min;
        }

      std::vector<double> points;

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        points.emplace_back(radius * quadrature.weight(q) * (rad_1 - rad_0));

      for (unsigned int q = 0; q < n_quadrature_points; ++q)
        points.emplace_back(radius * quadrature.weight(q) * (rad_2 - rad_1));

      return points;
    }
}

template <int dim>
std::vector<Tensor<1, dim, double>>
MortarManager<dim>::get_normals(const double &angle_cell_center) const
{
  // Coordinates of cell quadrature points
  const auto points = get_points(angle_cell_center);

  std::vector<Tensor<1, dim, double>> result;

  for (const auto &point : points)
    result.emplace_back(point / point.norm());

  return result;
}

template <int dim>
std::pair<unsigned int, unsigned int>
MortarManager<dim>::get_config(const double &rad) const
{
  // Aalignment tolerance
  const double tolerance = 1e-8;
  // Angular variation in each cell
  const double delta = 2 * numbers::PI / n_subdivisions;
  // Minimum rotation angle
  double rot_min = rotation_angle - std::floor(rotation_angle / delta) * delta;
  // Point position in the cell
  const double segment = (rad - delta / 2) / delta;
  // Point position after rotation
  const double segment_rot = (rad - delta / 2 - rot_min) / delta;

  if (this->is_mesh_aligned())
    {
      // Case 1: mesh is aligned
      return {0, std::round(segment)};
    }
  else
    {
      // Case 2: mesh is not aligned
      if (std::abs(segment - std::round(segment)) < tolerance)
        // outer (fixed) domain
        return {2, std::round(segment)};
      else
        // inner (rotated) domain
        return {1,
                static_cast<unsigned int>(std::round(segment_rot)) %
                  (2 * n_subdivisions)};
    }
}



/*-------------- CouplingOperator -------------------------------*/

template <int dim, typename Number>
CouplingOperator<dim, Number>::CouplingOperator(
  const Mapping<dim>                                        &mapping,
  const DoFHandler<dim>                                     &dof_handler,
  const AffineConstraints<Number>                           &constraints,
  const Quadrature<dim>                                      quadrature,
  const std::shared_ptr<CouplingEvaluationBase<dim, Number>> evaluator,
  const unsigned int                                         n_subdivisions,
  const double                                               radius,
  const double                                               rotation_angle,
  const unsigned int                                         bid_rotor,
  const unsigned int                                         bid_stator,
  const double                                               sip_factor)
  : mapping(mapping)
  , dof_handler(dof_handler)
  , constraints(constraints)
  , quadrature(quadrature)
  , bid_rotor(bid_rotor)
  , bid_stator(bid_stator)
  , evaluator(evaluator)
{
  this->N                    = evaluator->data_size();
  this->relevant_dof_indices = evaluator->get_relevant_dof_indices();
  this->n_dofs_per_cell      = this->relevant_dof_indices.size();

  data.penalty_factor =
    compute_penalty_factor(dof_handler.get_fe().degree, sip_factor);

  // Create manager at the quadrature point level
  mortar_manager_q = std::make_shared<MortarManager<dim>>(
    n_subdivisions,
    quadrature.get_tensor_basis()[0].size(),
    radius,
    rotation_angle);

  // Create manager at the cell level
  mortar_manager_cell = std::make_shared<MortarManager<dim>>(n_subdivisions,
                                                             1,
                                                             radius,
                                                             rotation_angle);

  // Number of quadrature points
  const unsigned int n_points = mortar_manager_q->get_n_points();
  // Number of cells
  const unsigned int n_sub_cells = mortar_manager_cell->get_n_points();

  std::vector<types::global_dof_index> is_local;
  std::vector<types::global_dof_index> is_ghost;
  std::vector<types::global_dof_index> is_local_cell;
  std::vector<types::global_dof_index> is_ghost_cell;

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto face_no : cell->face_indices())
        if ((cell->face(face_no)->boundary_id() == bid_rotor) ||
            (cell->face(face_no)->boundary_id() == bid_stator))
          {
            const auto face = cell->face(face_no);

            // Indices of quadrature points
            const auto indices_q =
              mortar_manager_q->get_indices(get_angle_cell_center(cell, face));

            /* Loop over the quadrature points, storing indices for both sides.
             * We assume that the rotor side is the 'local' reference, and the
             * stator side is the 'ghost reference. */
            for (unsigned int ii = 0; ii < indices_q.size(); ++ii)
              {
                unsigned int i = indices_q[ii];
                unsigned int id_local, id_ghost;

                if (face->boundary_id() == bid_rotor)
                  {
                    id_local = i;
                    id_ghost = i + n_points;
                  }
                else if (face->boundary_id() == bid_stator)
                  {
                    id_local = i + n_points;
                    id_ghost = i;
                  }

                is_local.emplace_back(id_local);
                is_ghost.emplace_back(id_ghost);
              }

            // Indices of cells/DoFs on them
            const auto indices = mortar_manager_cell->get_indices(
              get_angle_cell_center(cell, face));

            const auto local_dofs = this->get_dof_indices(cell);

            /* Loop over the DoFs indices of the cells at the rotor-stator
             * interface. The logic of local (rotor) and ghost (stator) is the
             * same as in the previous loop. */
            for (unsigned int ii = 0; ii < indices.size(); ++ii)
              {
                unsigned int i = indices[ii];
                unsigned int id_local, id_ghost;

                if (face->boundary_id() == bid_rotor)
                  {
                    id_local = i;
                    id_ghost = i + n_sub_cells;
                  }
                else if (face->boundary_id() == bid_stator)
                  {
                    id_local = i + n_sub_cells;
                    id_ghost = i;
                  }

                is_local_cell.emplace_back(id_local);
                is_ghost_cell.emplace_back(id_ghost);

                for (const auto i : local_dofs)
                  dof_indices.emplace_back(i);
              }

            // Weights of quadrature points
            const auto weights =
              mortar_manager_q->get_weights(get_angle_cell_center(cell, face));
            data.all_weights.insert(data.all_weights.end(),
                                    weights.begin(),
                                    weights.end());

            // Normals of quadrature points
            if (false)
              {
                const auto points = mortar_manager_q->get_points(
                  get_angle_cell_center(cell, face));
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
                    if ((face_no / 2) == 0)
                      {
                        quad.emplace_back(p[1]);
                      }
                    else if ((face_no / 2) == 1)
                      {
                        quad.emplace_back(p[0]);
                      }
                    else
                      {
                        AssertThrow(false, ExcInternalError());
                      }
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
                auto normals = mortar_manager_q->get_normals(
                  get_angle_cell_center(cell, face));
                if (face->boundary_id() == bid_stator)
                  for (auto &normal : normals)
                    normal *= -1.0;
                data.all_normals.insert(data.all_normals.end(),
                                        normals.begin(),
                                        normals.end());

                auto points = mortar_manager_q->get_points_ref(
                  get_angle_cell_center(cell, face));

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

            // Penalty parmeter
            const Number penalty_parameter = compute_penalty_parameter(cell);

            // Store penalty parameter for all quadrature points
            for (unsigned int i = 0; i < mortar_manager_q->get_n_points(
                                           get_angle_cell_center(cell, face));
                 ++i)
              data.all_penalty_parameter.emplace_back(penalty_parameter);
          }

  // Setup communication
  partitioner.reinit(is_local, is_ghost, dof_handler.get_mpi_communicator());
  partitioner_cell.reinit(is_local_cell,
                          is_ghost_cell,
                          dof_handler.get_mpi_communicator());

  // Finalized penalty parameters
  std::vector<Number> all_penalty_parameter_ghost(
    data.all_penalty_parameter.size());
  partitioner.template export_to_ghosted_array<Number, 1>(
    data.all_penalty_parameter, all_penalty_parameter_ghost);
  for (unsigned int i = 0; i < data.all_penalty_parameter.size(); ++i)
    data.all_penalty_parameter[i] =
      std::min(data.all_penalty_parameter[i], all_penalty_parameter_ghost[i]);

  // Finialize DoF indices and update constraints
  dof_indices_ghost.resize(dof_indices.size());
  partitioner_cell.template export_to_ghosted_array<types::global_dof_index, 0>(
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
      const Number factor = 0.5; /*Assuming that we are within the domain*/
      for (unsigned int q = 0; q < face_quadrature.size(); ++q)
        surface_area += fe_face_values.JxW(q) * factor;
    }

  return surface_area / volume;
}

template <int dim, typename Number>
double
CouplingOperator<dim, Number>::get_angle_cell_center(
  const typename Triangulation<dim>::cell_iterator &cell,
  const typename Triangulation<dim>::face_iterator &face) const
{
  if (false)
    return point_to_angle(face->center());
  else
    return point_to_angle(mapping.transform_unit_to_real_cell(
      cell,
      MappingQ1<dim>().transform_real_to_unit_cell(cell, face->center())));
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
void
CouplingOperator<dim, Number>::vmult_add(VectorType       &dst,
                                         const VectorType &src) const
{
  // 1) Evaluate
  unsigned int ptr_q = 0;

  Vector<Number> buffer;

  std::vector<Number> all_value_m(data.all_normals.size() * N);
  std::vector<Number> all_value_p(data.all_normals.size() * N);

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_rotor) ||
            (face->boundary_id() == bid_stator))
          {
            /* Number of quadrature points at the cell(rotor)/cell(stator)
             * interaction. For non-aligned meshes, this value indicates the
             * number of quadrature points at both rotor and stator cells. */
            const unsigned int n_q_points =
              mortar_manager_q->get_n_points(get_angle_cell_center(cell, face));

            evaluator->local_reinit(
              cell,
              ArrayView<const Point<dim, Number>>(all_points_ref.data() + ptr_q,
                                                  n_q_points));

            buffer.reinit(n_dofs_per_cell);

            const auto local_dofs = this->get_dof_indices(cell);

            for (unsigned int i = 0; i < local_dofs.size(); ++i)
              buffer[i] = src[local_dofs[i]];

            evaluator->local_evaluate(
              data, buffer, ptr_q, 1, all_value_m.data() + ptr_q * N);

            ptr_q += n_q_points;
          }

  // 2) Communicate
  partitioner.template export_to_ghosted_array<Number, 0>(
    ArrayView<const Number>(reinterpret_cast<Number *>(all_value_m.data()),
                            all_value_m.size()),
    ArrayView<Number>(reinterpret_cast<Number *>(all_value_p.data()),
                      all_value_p.size()),
    N);

  // 3) Integrate
  ptr_q = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_rotor) ||
            (face->boundary_id() == bid_stator))
          {
            // Quadrature points at the cell(rotor)/cell(stator) interaction
            const unsigned int n_q_points =
              mortar_manager_q->get_n_points(get_angle_cell_center(cell, face));

            evaluator->local_reinit(
              cell,
              ArrayView<const Point<dim, Number>>(all_points_ref.data() + ptr_q,
                                                  n_q_points));

            buffer.reinit(n_dofs_per_cell);
            evaluator->local_integrate(data,
                                       buffer,
                                       ptr_q,
                                       1,
                                       all_value_m.data() + ptr_q * N,
                                       all_value_p.data() + ptr_q * N);

            const auto local_dofs = this->get_dof_indices(cell);
            constraints.distribute_local_to_global(buffer, local_dofs, dst);

            ptr_q += n_q_points;
          }

  dst.compress(VectorOperation::add);
}

template <int dim, typename Number>
void
CouplingOperator<dim, Number>::add_diagonal_entries(VectorType &diagonal) const
{
  unsigned int ptr_q = 0;

  Vector<Number>      buffer, diagonal_local;
  std::vector<Number> all_value_m, all_value_p;

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_rotor) ||
            (face->boundary_id() == bid_stator))
          {
            // Quadrature points at the cell(rotor)/cell(stator) interaction
            const unsigned int n_q_points =
              mortar_manager_q->get_n_points(get_angle_cell_center(cell, face));

            evaluator->local_reinit(
              cell,
              ArrayView<const Point<dim, Number>>(all_points_ref.data() + ptr_q,
                                                  n_q_points));

            buffer.reinit(n_dofs_per_cell);
            diagonal_local.reinit(n_dofs_per_cell);
            all_value_m.resize(n_q_points * N);
            all_value_p.resize(n_q_points * N);

            for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                  buffer[j] = static_cast<Number>(i == j);

                evaluator->local_evaluate(
                  data, buffer, ptr_q, 1, all_value_m.data());

                buffer.reinit(n_dofs_per_cell);
                evaluator->local_integrate(data,
                                           buffer,
                                           ptr_q,
                                           1,
                                           all_value_m.data(),
                                           all_value_p.data());

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

template <int dim, typename Number>
void
CouplingOperator<dim, Number>::add_system_matrix_entries(
  TrilinosWrappers::SparseMatrix &system_matrix) const
{
  const auto constraints = &constraints_extended;

  std::vector<Number> all_value_m(data.all_normals.size() * n_dofs_per_cell *
                                  N);
  std::vector<Number> all_value_p(data.all_normals.size() * n_dofs_per_cell *
                                  N);

  unsigned int ptr_q = 0;

  Vector<Number> buffer;

  // 1) Evaluate
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_rotor) ||
            (face->boundary_id() == bid_stator))
          {
            // Quadrature points at the cell(rotor)/cell(stator) interaction
            const unsigned int n_q_points =
              mortar_manager_q->get_n_points(get_angle_cell_center(cell, face));

            evaluator->local_reinit(
              cell,
              ArrayView<const Point<dim, Number>>(all_points_ref.data() + ptr_q,
                                                  n_q_points));

            buffer.reinit(n_dofs_per_cell);

            for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                  buffer[j] = static_cast<Number>(i == j);

                evaluator->local_evaluate(data,
                                          buffer,
                                          ptr_q,
                                          n_dofs_per_cell,
                                          all_value_m.data() +
                                            (ptr_q * n_dofs_per_cell + i) * N);
              }

            ptr_q += n_q_points;
          }

  const unsigned n_q_points =
    Utilities::pow(quadrature.get_tensor_basis()[0].size(), dim - 1);

  // 2) Communicate
  partitioner_cell.template export_to_ghosted_array<Number, 0>(
    ArrayView<const Number>(reinterpret_cast<Number *>(all_value_m.data()),
                            all_value_m.size()),
    ArrayView<Number>(reinterpret_cast<Number *>(all_value_p.data()),
                      all_value_p.size()),
    n_dofs_per_cell * n_q_points * N);


  ptr_q                 = 0;
  unsigned int ptr_dofs = 0;

  // 3) Integrate
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_rotor) ||
            (face->boundary_id() == bid_stator))
          {
            const unsigned int n_sub_cells = mortar_manager_cell->get_n_points(
              get_angle_cell_center(cell, face));

            for (unsigned int sc = 0; sc < n_sub_cells; ++sc)
              {
                const unsigned int n_q_points =
                  mortar_manager_q->get_n_points(
                    get_angle_cell_center(cell, face)) /
                  n_sub_cells;

                evaluator->local_reinit(cell,
                                        ArrayView<const Point<dim, Number>>(
                                          all_points_ref.data() + ptr_q,
                                          n_q_points));

                for (unsigned int bb = 0; bb < 2; ++bb)
                  {
                    FullMatrix<Number> cell_matrix(n_dofs_per_cell,
                                                   n_dofs_per_cell);

                    for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                      {
                        buffer.reinit(n_dofs_per_cell);
                        if (bb == 0)
                          evaluator->local_integrate(
                            data,
                            buffer,
                            ptr_q,
                            n_dofs_per_cell,
                            all_value_m.data() +
                              (ptr_q * n_dofs_per_cell + i) * N,
                            nullptr);
                        else
                          evaluator->local_integrate(
                            data,
                            buffer,
                            ptr_q,
                            n_dofs_per_cell,
                            nullptr,
                            all_value_p.data() +
                              (ptr_q * n_dofs_per_cell + i) * N);

                        for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                          cell_matrix[j][i] = buffer[j];
                      }


                    std::vector<types::global_dof_index> a(
                      dof_indices.begin() + ptr_dofs,
                      dof_indices.begin() + ptr_dofs + n_dofs_per_cell);

                    if (bb == 0)
                      {
                        constraints->distribute_local_to_global(cell_matrix,
                                                                a,
                                                                system_matrix);
                      }
                    else
                      {
                        std::vector<types::global_dof_index> b(
                          dof_indices_ghost.begin() + ptr_dofs,
                          dof_indices_ghost.begin() + ptr_dofs +
                            n_dofs_per_cell);

                        constraints->distribute_local_to_global(cell_matrix,
                                                                a,
                                                                b,
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



/*-------------- CouplingEvaluationBaseSIPG -------------------------------*/



/**
 * @brief Compute the number of subdivisions at the rotor-stator interface and the rotor radius
 * @param[in] dof_handler DoFHandler associated to the triangulation
 * @param[in] mortar_parameters The information about the mortar method
 * control, including the rotor mesh parameters
 *
 * @return n_subdivisions Number of cells at the interface between inner
 * and outer domains
 * @return radius Radius at the interface between inner and outer domains
 */
template <int dim>
std::pair<unsigned int, double>
compute_n_subdivisions_and_radius(
  const DoFHandler<dim>         &dof_handler,
  const Parameters::Mortar<dim> &mortar_parameters)
{
  // Number of subdivisions per process
  unsigned int n_subdivisions_local = 0;
  // Number of vertices at the boundary per process
  unsigned int n_vertices_local = 0;
  // Tolerance for rotor radius computation
  const double tolerance = 1e-8;
  // Min and max values for rotor radius computation
  double radius_min = 1e12;
  double radius_max = 1e-12;

  // Check number of faces and vertices at the rotor-stator interface
  for (const auto &cell :
       dof_handler.get_triangulation().active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (const auto &face : cell->face_iterators())
            {
              if (face->at_boundary())
                {
                  if (face->boundary_id() ==
                      mortar_parameters.rotor_boundary_id)
                    {
                      n_subdivisions_local++;
                      for (unsigned int vertex_index = 0;
                           vertex_index < face->n_vertices();
                           vertex_index++)
                        {
                          n_vertices_local++;
                          auto   v = face->vertex(vertex_index);
                          double radius_current =
                            v.distance(mortar_parameters.center_of_rotation);
                          radius_min = std::min(radius_min, radius_current);
                          radius_max = std::max(radius_max, radius_current);
                        }
                    }
                }
            }
        }
    }

  // Total number of faces
  const unsigned int n_subdivisions =
    Utilities::MPI::sum(n_subdivisions_local,
                        dof_handler.get_mpi_communicator());

  // Min and max values over all processes
  radius_min =
    Utilities::MPI::min(radius_min, dof_handler.get_mpi_communicator());
  radius_max =
    Utilities::MPI::max(radius_max, dof_handler.get_mpi_communicator());

  AssertThrow(
    std::abs(radius_max - radius_min) < tolerance,
    ExcMessage(
      "The computed radius of the rotor mesh has a variation greater than "
      "the tolerance across the rotor domain, meaning that the prescribed "
      "center of rotation and the rotor geometry are not in accordance."));

  // Final radius value
  const double radius = radius_min;

  return {n_subdivisions, radius};
}

/**
 * @brief Construct oversampled quadrature
 *
 * @param[in] quadrature Quadrature for local cell operations
 * @param[in] mortar_parameters The information about the mortar method
 * control, including the rotor mesh parameters
 *
 * @return Quadrature oversampled
 */
template <int dim>
static Quadrature<dim>
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

template <int dim, int n_components, typename Number>
CouplingEvaluationBaseSIPG<dim, n_components, Number>::
  CouplingEvaluationBaseSIPG(const Mapping<dim>    &mapping,
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
CouplingEvaluationBaseSIPG<dim, n_components, Number>::data_size() const
{
  return n_components * 2;
}

template <int dim, int n_components, typename Number>
const std::vector<unsigned int> &
CouplingEvaluationBaseSIPG<dim, n_components, Number>::
  get_relevant_dof_indices() const
{
  return relevant_dof_indices;
}

template <int dim, int n_components, typename Number>
void
CouplingEvaluationBaseSIPG<dim, n_components, Number>::local_reinit(
  const typename Triangulation<dim>::cell_iterator &cell,
  const ArrayView<const Point<dim, Number>>        &points) const
{
  this->phi_m.reinit(cell, points);
}

template <int dim, int n_components, typename Number>
void
CouplingEvaluationBaseSIPG<dim, n_components, Number>::local_evaluate(
  const CouplingEvaluationData<dim, Number> &data,
  const Vector<Number>                      &buffer,
  const unsigned int                         ptr_q,
  const unsigned int                         q_stride,
  Number                                    *all_value_m) const
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

      // Store data in buffer
      BufferRW<Number> buffer_m(all_value_m, q * 2 * n_components * q_stride);

      buffer_m.write(value_m);
      buffer_m.write(gradient_m);
    }
}

template <int dim, int n_components, typename Number>
void
CouplingEvaluationBaseSIPG<dim, n_components, Number>::local_integrate(
  const CouplingEvaluationData<dim, Number> &data,
  Vector<Number>                            &buffer,
  const unsigned int                         ptr_q,
  const unsigned int                         q_stride,
  Number                                    *all_value_m,
  Number                                    *all_value_p) const
{
  for (const auto q : this->phi_m.quadrature_point_indices())
    {
      const unsigned int q_index = ptr_q + q;

      BufferRW<Number> buffer_m(all_value_m, q * 2 * n_components * q_stride);
      BufferRW<Number> buffer_p(all_value_p, q * 2 * n_components * q_stride);

      const auto value_m           = buffer_m.template read<value_type>();
      const auto value_p           = buffer_p.template read<value_type>();
      const auto normal_gradient_m = buffer_m.template read<value_type>();
      const auto normal_gradient_p = buffer_p.template read<value_type>();

      const auto JxW               = data.all_weights[q_index];
      const auto penalty_parameter = data.all_penalty_parameter[q_index];
      const auto normal            = data.all_normals[q_index];

      const auto value_jump = (value_m - value_p);
      const auto gradient_normal_avg =
        (normal_gradient_m - normal_gradient_p) * 0.5;

      const double sigma = penalty_parameter * data.penalty_factor;

      // - (n avg(∇v), jump(u))
      this->phi_m.submit_gradient(outer(-value_jump, normal) * 0.5 * JxW, q);

      // + (jump(v), σ jump(u) - avg(∇u) n)
      this->phi_m.submit_value((value_jump * sigma - gradient_normal_avg) * JxW,
                               q);
    }

  this->phi_m.test_and_sum(buffer,
                           EvaluationFlags::values |
                             EvaluationFlags::gradients);
}


/*-------------- Explicit Instantiations -------------------------------*/
template class MortarManager<2>;
template class MortarManager<3>;

template class CouplingOperator<2, double>;
template class CouplingOperator<3, double>;

template class CouplingEvaluationBaseSIPG<2, 1, double>;
template class CouplingEvaluationBaseSIPG<2, 2, double>;
template class CouplingEvaluationBaseSIPG<2, 3, double>;
template class CouplingEvaluationBaseSIPG<3, 1, double>;
template class CouplingEvaluationBaseSIPG<3, 3, double>;
template class CouplingEvaluationBaseSIPG<3, 4, double>;

#endif
