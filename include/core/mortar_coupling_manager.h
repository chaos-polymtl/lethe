// SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_core_mortar_coupling_manager_h
#define lethe_core_mortar_coupling_manager_h

#include <core/utilities.h>

#include <deal.II/base/mpi_noncontiguous_partitioner.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.templates.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

using namespace dealii;

/**
 * @brief Base class for the mortar manager
 * @param n_subdivisions Number of cells at the interface between inner and outer domains
 * @param n_quadrature_points Number of quadrature points per cell
 * @param radius Radius at the interface between inner and outer domains
 * @param rotate_pi Rotation angle for the inner domain
 */
template <int dim>
class MortarManager
{
public:
  MortarManager(const unsigned int n_subdivisions,
                const unsigned int n_quadrature_points,
                const double       radius,
                const double       rotate_pi)
    : n_subdivisions(n_subdivisions)
    , n_quadrature_points(n_quadrature_points)
    , radius(radius)
    , rotate_pi(rotate_pi)
    , quadrature(n_quadrature_points)
  {}

  /**
   * @brief Verify if cells of the inner and outer domains are aligned
   */
  bool
  is_mesh_aligned() const
  {
    const double tolerance = 1e-8;
    const double delta     = 2 * numbers::PI / n_subdivisions;

    return std::abs(rotate_pi / delta - std::round(rotate_pi / delta)) <
           tolerance;
  }

  /**
   * @brief Returns the total number of quadrature points at the inner/outer boundary interface
   */
  unsigned int
  get_n_points() const
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

  /**
   * @brief Returns the number of quadrature points at each inner/outer cell matching pair
   * @param[in] rad Angular coordinate of cell center
   */
  unsigned int
  get_n_points(const double &rad) const
  {
    (void)rad;

    if (this->is_mesh_aligned()) // aligned
      {
        return n_quadrature_points;
      }
    else // inside/outside
      {
        return 2 * n_quadrature_points;
      }
  }

  /**
   * @brief Returns the indices of all quadrature points at both sides of the interface
   * @param[in] rad Angular coordinate of cell center
   */
  std::vector<unsigned int>
  get_indices(const double &rad) const
  {
    // mesh alignment type and cell index
    const auto [type, id] = get_config(rad);

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

  /**
   * @brief Returns the coordinates of the quadrature points at both sides of the inerface
   * @param[in] rad Angular coordinate of cell center
   * @param[out] points Coordinate of quadrature points of the cell
   */
  std::vector<Point<dim>>
  get_points(const double rad) const
  {
    // mesh alignment type and cell index
    const auto [type, id] = get_config(rad);
    // angle variation within each cell
    const double delta = 2 * numbers::PI / n_subdivisions;

    if (type == 0) // aligned
      {
        std::vector<Point<dim>> points;

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          points.emplace_back(
            rad_to_point<dim>(radius, (id + quadrature.point(q)[0]) * delta));

        return points;
      }
    else // point at the inner boundary lies somewhere in the face of the outer boundary cell
      {
        // rad_0: first cell vertex (fixed)
        // rad_1: shifted vertex
        // rad_2: last cell vertex (fixed)
        double rad_0, rad_1, rad_2;
        // minimum rotation angle
        double rot_min = rotate_pi - std::floor(rotate_pi / delta) * delta;

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
          points.emplace_back(rad_to_point<dim>(radius,
                                                rad_0 + quadrature.point(q)[0] *
                                                          (rad_1 - rad_0)));

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          points.emplace_back(rad_to_point<dim>(radius,
                                                rad_1 + quadrature.point(q)[0] *
                                                          (rad_2 - rad_1)));

        return points;
      }
  }

  std::vector<Point<1>>
  get_points_ref(const double rad) const
  {
    const auto [type, id] = get_config(rad);

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
          (rotate_pi - std::floor(rotate_pi / delta) * delta) / delta;

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

  /**
   * @brief Returns the weights of the quadrature points at both sides of the interface
   * @param[in] rad Angular coordinate of cell center
   * @param[out] points Angular weights of quadrature points of the cell 
   */
  std::vector<double>
  get_weights(const double &rad) const
  {
    // mesh alignment type and cell index
    const auto [type, id] = get_config(rad);
    // angle variation within each cell
    const double delta = 2 * numbers::PI / n_subdivisions;

    if (type == 0) // aligned
      {
        std::vector<double> points;

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          points.emplace_back(quadrature.weight(q) * delta);

        return points;
      }
    else // inside/outside
      {
        double rad_0, rad_1, rad_2;

        double rot_min = rotate_pi - std::floor(rotate_pi / delta) * delta;

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
          points.emplace_back(quadrature.weight(q) * (rad_1 - rad_0));

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          points.emplace_back(quadrature.weight(q) * (rad_2 - rad_1));

        return points;
      }
  }

  /**
   * @brief Returns the normal vector for the quadrature points
   * @param[in] rad Angular coordinate of cell center
   * @param[out] result Normal vectors of the cell quadrature points
   */
  std::vector<Tensor<1, dim, double>>
  get_normals(const double &rad) const
  {
    // Coordinates of cell quadrature points
    const auto points = get_points(rad);

    std::vector<Tensor<1, dim, double>> result;

    for (const auto &point : points)
      result.emplace_back(point / point.norm());

    return result;
  }

private:
  /**
   * @brief Returns the mesh alignement type and cell index
   * @param[in] rad Angular coordinate of cell center
   * @param[out] type Cell configuration type at the interface
   *                  type = 0: mesh aligned
   *                  type = 1: mesh not aligned, inner domain (allows rotation)
   *                  type = 2: mesh not aligned, outer domain (fixed)
   * @param[out] id Index of the cell in which lies the rotated cell center
   */
  std::pair<unsigned int, unsigned int>
  get_config(const double &rad) const
  {
    // alignment tolerance
    const double tolerance = 1e-8;
    // angular variation in each cell
    const double delta = 2 * numbers::PI / n_subdivisions;
    // minimum rotation angle
    double rot_min = rotate_pi - std::floor(rotate_pi / delta) * delta;
    // point position in the cell
    const double segment = (rad - delta / 2) / delta;
    // point position after rotation
    const double segment_rot = (rad - delta / 2 - rot_min) / delta;

    if (this->is_mesh_aligned())
      {
        // case 1: mesh is aligned
        return {0, std::round(segment)};
      }
    else
      {
        // case 2: mesh is not aligned
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

  const unsigned int n_subdivisions;
  const unsigned int n_quadrature_points;
  const double       radius;
  const double       rotate_pi;
  QGauss<1>          quadrature;
};



template <int dim, typename Number>
Number
contract(const Tensor<1, dim, Number> &grad,
         const Tensor<1, dim, Number> &normal)
{
  return grad * normal;
}

template <int dim, typename Number>
Tensor<1, dim, Number>
contract(const Tensor<2, dim, Number> &grad,
         const Tensor<1, dim, Number> &normal)
{
  return grad * normal;
}

template <int n_components, int dim, typename Number>
Tensor<1, n_components, Number>
contract(const Tensor<1, n_components, Tensor<1, dim, Number>> &grad,
         const Tensor<1, dim, Number>                          &normal)
{
  Tensor<1, n_components, Number> result;

  for (int r = 0; r < n_components; ++r)
    result[r] = grad[r] * normal;

  return result;
}

template <int dim, typename Number>
Tensor<1, dim, Number>
outer(const Number &value, const Tensor<1, dim, Number> &normal)
{
  return value * normal;
}

template <int dim, typename Number>
Tensor<2, dim, Number>
outer(const Tensor<1, dim, Number> &value, const Tensor<1, dim, Number> &normal)
{
  Tensor<2, dim, Number> result;

  for (unsigned int c = 0; c < dim; ++c)
    result[c] = value[c] * normal;

  return result;
}

template <int n_components, int dim, typename Number>
Tensor<1, n_components, Tensor<1, dim, Number>>
outer(const Tensor<1, n_components, Number> &value,
      const Tensor<1, dim, Number>          &normal)
{
  Tensor<1, n_components, Tensor<1, dim, Number>> result;

  for (unsigned int c = 0; c < n_components; ++c)
    result[c] = value[c] * normal;

  return result;
}

template <int dim, int dim_, typename Number>
inline DEAL_II_ALWAYS_INLINE void
symm_scalar_product_add(Tensor<1, dim_, Tensor<1, dim, Number>> &v_gradient,
                        const Tensor<2, dim, Number>            &u_gradient,
                        const Number                            &factor)
{
  for (unsigned int d = 0; d < dim; ++d)
    v_gradient[d][d] += u_gradient[d][d] * factor;

  for (unsigned int e = 0; e < dim; ++e)
    for (unsigned int d = e + 1; d < dim; ++d)
      {
        const auto tmp = (u_gradient[d][e] + u_gradient[e][d]) * (factor * 0.5);
        v_gradient[d][e] += tmp;
        v_gradient[e][d] += tmp;
      }
}

template <int dim, int n_components, typename Number>
class CouplingOperator
{
public:
  using FEPointIntegrator = FEPointEvaluation<n_components, dim, dim, Number>;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  CouplingOperator(const Mapping<dim>              &mapping,
                   const DoFHandler<dim>           &dof_handler,
                   const AffineConstraints<Number> &constraints,
                   const Quadrature<dim>            quadrature,
                   const unsigned int               n_subdivisions,
                   const double                     radius,
                   const double                     rotate_pi,
                   const unsigned int               bid_0,
                   const unsigned int               bid_1,
                   const double                     sip_factor = 1.0);

  const AffineConstraints<Number> &
  get_affine_constraints() const;

  void
  vmult_add(VectorType &dst, const VectorType &src) const;

  void
  add_diagonal_entries(VectorType &diagonal) const;

  void
  add_sparsity_pattern_entries(TrilinosWrappers::SparsityPattern &dsp) const;

  void
  add_system_matrix_entries(
    TrilinosWrappers::SparseMatrix &system_matrix) const;

private:
  Number
  compute_penalty_factor(const unsigned int degree, const Number factor) const;

  Number
  compute_penalty_parameter(
    const typename Triangulation<dim>::cell_iterator &cell) const;

  double
  get_rad(const typename Triangulation<dim>::cell_iterator &cell,
          const typename Triangulation<dim>::face_iterator &face) const;

  const Mapping<dim>              &mapping;
  const DoFHandler<dim>           &dof_handler;
  const AffineConstraints<Number> &constraints;
  const Quadrature<dim>            quadrature;

  std::vector<std::tuple<std::vector<double>,
                         typename Triangulation<dim>::active_cell_iterator,
                         std::vector<Point<dim>>,
                         typename Triangulation<dim>::active_cell_iterator,
                         std::vector<Point<dim>>,
                         std::vector<Tensor<1, dim, Number>>,
                         Number>>
    all_intersections;

  Number penalty_factor;

  std::shared_ptr<MortarManager<dim>> mortar_manager_q;
  std::shared_ptr<MortarManager<dim>> mortar_manager_cell;

  unsigned int bid_0;
  unsigned int bid_1;

  Utilities::MPI::NoncontiguousPartitioner partitioner;
  Utilities::MPI::NoncontiguousPartitioner partitioner_cell;

  std::vector<types::global_dof_index> dof_indices;
  std::vector<types::global_dof_index> dof_indices_ghost;

  std::vector<Number>                 all_penalty_parameter;
  std::vector<Number>                 all_weights;
  std::vector<Point<dim, Number>>     all_points_ref;
  std::vector<Tensor<1, dim, Number>> all_normals;

  AffineConstraints<Number> constraints_extended;
};

template <int dim, int n_components, typename Number>
CouplingOperator<dim, n_components, Number>::CouplingOperator(
  const Mapping<dim>              &mapping,
  const DoFHandler<dim>           &dof_handler,
  const AffineConstraints<Number> &constraints,
  const Quadrature<dim>            quadrature,
  const unsigned int               n_subdivisions,
  const double                     radius,
  const double                     rotate_pi,
  const unsigned int               bid_0,
  const unsigned int               bid_1,
  const double                     sip_factor)
  : mapping(mapping)
  , dof_handler(dof_handler)
  , constraints(constraints)
  , quadrature(quadrature)
{
  this->bid_0 = bid_0;
  this->bid_1 = bid_1;

  penalty_factor =
    compute_penalty_factor(dof_handler.get_fe().degree, sip_factor);

  mortar_manager_q = std::make_shared<MortarManager<dim>>(
    n_subdivisions, quadrature.get_tensor_basis()[0].size(), radius, rotate_pi);

  mortar_manager_cell =
    std::make_shared<MortarManager<dim>>(n_subdivisions, 1, radius, rotate_pi);

  const unsigned int n_points    = mortar_manager_q->get_n_points();
  const unsigned int n_sub_cells = mortar_manager_cell->get_n_points();

  std::vector<types::global_dof_index> is_local;
  std::vector<types::global_dof_index> is_ghost;
  std::vector<types::global_dof_index> is_local_cell;
  std::vector<types::global_dof_index> is_ghost_cell;

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto face_no : cell->face_indices())
        if ((cell->face(face_no)->boundary_id() == bid_0) ||
            (cell->face(face_no)->boundary_id() == bid_1))
          {
            const auto face = cell->face(face_no);

            // indices
            const auto indices_q =
              mortar_manager_q->get_indices(get_rad(cell, face));
            for (unsigned int ii = 0; ii < indices_q.size(); ++ii)
              {
                unsigned int i = indices_q[ii];
                unsigned int id_local, id_ghost;

                if (face->boundary_id() == bid_0)
                  {
                    id_local = i;
                    id_ghost = i + n_points;
                  }
                else if (face->boundary_id() == bid_1)
                  {
                    id_local = i + n_points;
                    id_ghost = i;
                  }

                is_local.emplace_back(id_local);
                is_ghost.emplace_back(id_ghost);
              }

            // indices of cells/DoFs on them
            const auto indices =
              mortar_manager_cell->get_indices(get_rad(cell, face));

            std::vector<types::global_dof_index> local_dofs(
              dof_handler.get_fe().n_dofs_per_cell());
            cell->get_dof_indices(local_dofs);

            for (unsigned int ii = 0; ii < indices.size(); ++ii)
              {
                unsigned int i = indices[ii];
                unsigned int id_local, id_ghost;

                if (face->boundary_id() == bid_0)
                  {
                    id_local = i;
                    id_ghost = i + n_sub_cells;
                  }
                else if (face->boundary_id() == bid_1)
                  {
                    id_local = i + n_sub_cells;
                    id_ghost = i;
                  }

                is_local_cell.emplace_back(id_local);
                is_ghost_cell.emplace_back(id_ghost);

                for (const auto i : local_dofs)
                  dof_indices.emplace_back(i);
              }

            // weights
            const auto weights =
              mortar_manager_q->get_weights(get_rad(cell, face));
            all_weights.insert(all_weights.end(),
                               weights.begin(),
                               weights.end());

            // normals
            auto normals = mortar_manager_q->get_normals(get_rad(cell, face));
            if (face->boundary_id() == bid_0)
              for (auto &normal : normals)
                normal *= -1.0;
            all_normals.insert(all_normals.end(),
                               normals.begin(),
                               normals.end());

            // points (also convert real to unit coordinates)
            if (false)
              {
                const auto points =
                  mortar_manager_q->get_points(get_rad(cell, face));
                std::vector<Point<dim, Number>> points_ref(points.size());
                mapping.transform_points_real_to_unit_cell(cell,
                                                           points,
                                                           points_ref);
                all_points_ref.insert(all_points_ref.end(),
                                      points_ref.begin(),
                                      points_ref.end());
              }
            else
              {
                auto points =
                  mortar_manager_q->get_points_ref(get_rad(cell, face));

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

            // penalty parmeter
            const Number penalty_parameter = compute_penalty_parameter(cell);

            for (unsigned int i = 0;
                 i < mortar_manager_q->get_n_points(get_rad(cell, face));
                 ++i)
              all_penalty_parameter.emplace_back(penalty_parameter);
          }

  // setup communication
  partitioner.reinit(is_local, is_ghost, dof_handler.get_mpi_communicator());
  partitioner_cell.reinit(is_local_cell,
                          is_ghost_cell,
                          dof_handler.get_mpi_communicator());

  // finalized penalty parameters
  std::vector<Number> all_penalty_parameter_ghost(all_penalty_parameter.size());
  partitioner.template export_to_ghosted_array<Number, 1>(
    all_penalty_parameter, all_penalty_parameter_ghost);
  for (unsigned int i = 0; i < all_penalty_parameter.size(); ++i)
    all_penalty_parameter[i] =
      std::min(all_penalty_parameter[i], all_penalty_parameter_ghost[i]);

  // finialize DoF indices
  dof_indices_ghost.resize(dof_indices.size());
  partitioner_cell.template export_to_ghosted_array<types::global_dof_index, 0>(
    dof_indices, dof_indices_ghost, dof_handler.get_fe().n_dofs_per_cell());

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

template <int dim, int n_components, typename Number>
const AffineConstraints<Number> &
CouplingOperator<dim, n_components, Number>::get_affine_constraints() const
{
  return constraints_extended;
}

template <int dim, int n_components, typename Number>
Number
CouplingOperator<dim, n_components, Number>::compute_penalty_factor(
  const unsigned int degree,
  const Number       factor) const
{
  return factor * (degree + 1.0) * (degree + 1.0);
}

template <int dim, int n_components, typename Number>
Number
CouplingOperator<dim, n_components, Number>::compute_penalty_parameter(
  const typename Triangulation<dim>::cell_iterator &cell) const
{
  const auto        &fe     = dof_handler.get_fe();
  const unsigned int degree = fe.degree;

  dealii::QGauss<dim>   quadrature(degree + 1);
  dealii::FEValues<dim> fe_values(mapping,
                                  fe,
                                  quadrature,
                                  dealii::update_JxW_values);

  dealii::QGauss<dim - 1>   face_quadrature(degree + 1);
  dealii::FEFaceValues<dim> fe_face_values(mapping,
                                           fe,
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
        (cell->at_boundary(f) && !cell->has_periodic_neighbor(f)) ? 1. : 0.5;
      for (unsigned int q = 0; q < face_quadrature.size(); ++q)
        surface_area += fe_face_values.JxW(q) * factor;
    }

  return surface_area / volume;
}

template <int dim, int n_components, typename Number>
double
CouplingOperator<dim, n_components, Number>::get_rad(
  const typename Triangulation<dim>::cell_iterator &cell,
  const typename Triangulation<dim>::face_iterator &face) const
{
  if (false)
    return point_to_rad(face->center());
  else
    return point_to_rad(mapping.transform_unit_to_real_cell(
      cell,
      MappingQ1<dim>().transform_real_to_unit_cell(cell, face->center())));
}

template <int dim, int n_components, typename Number>
void
CouplingOperator<dim, n_components, Number>::vmult_add(
  VectorType       &dst,
  const VectorType &src) const
{
  // 1) evaluate
  unsigned int ptr_q = 0;

  FEPointIntegrator phi_m(mapping,
                          dof_handler.get_fe(),
                          update_values | update_gradients);

  Vector<Number> buffer;

  std::vector<typename FEPointIntegrator::value_type> all_value_m(
    all_normals.size() * 2);
  std::vector<typename FEPointIntegrator::value_type> all_value_p(
    all_normals.size() * 2);

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_0) || (face->boundary_id() == bid_1))
          {
            const unsigned int n_q_points =
              mortar_manager_q->get_n_points(get_rad(cell, face));

            phi_m.reinit(cell,
                         ArrayView<const Point<dim, Number>>(
                           all_points_ref.data() + ptr_q, n_q_points));

            buffer.reinit(cell->get_fe().n_dofs_per_cell());
            cell->get_dof_values(src, buffer);

            phi_m.evaluate(buffer,
                           EvaluationFlags::values |
                             EvaluationFlags::gradients);

            for (const auto q : phi_m.quadrature_point_indices())
              {
                const unsigned int q_index = ptr_q + q;

                AssertIndexRange(q_index, all_normals.size());
                AssertIndexRange(q_index * 2 + 1, all_value_m.size());

                const auto normal     = all_normals[q_index];
                const auto value_m    = phi_m.get_value(q);
                const auto gradient_m = contract(phi_m.get_gradient(q), normal);

                all_value_m[q_index * 2 + 0] = value_m;
                all_value_m[q_index * 2 + 1] = gradient_m;
              }

            ptr_q += n_q_points;
          }

  // 2) communicate
  partitioner.template export_to_ghosted_array<Number, 2u * n_components>(
    ArrayView<const Number>(reinterpret_cast<Number *>(all_value_m.data()),
                            all_value_m.size() * n_components),
    ArrayView<Number>(reinterpret_cast<Number *>(all_value_p.data()),
                      all_value_p.size() * n_components));

  // 3) integrate
  ptr_q = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_0) || (face->boundary_id() == bid_1))
          {
            const unsigned int n_q_points =
              mortar_manager_q->get_n_points(get_rad(cell, face));

            phi_m.reinit(cell,
                         ArrayView<const Point<dim, Number>>(
                           all_points_ref.data() + ptr_q, n_q_points));

            for (const auto q : phi_m.quadrature_point_indices())
              {
                const unsigned int q_index = ptr_q + q;

                const auto value_m           = all_value_m[q_index * 2 + 0];
                const auto value_p           = all_value_p[q_index * 2 + 0];
                const auto normal_gradient_m = all_value_m[q_index * 2 + 1];
                const auto normal_gradient_p = all_value_p[q_index * 2 + 1];
                const auto JxW               = all_weights[q_index];
                const auto penalty_parameter = all_penalty_parameter[q_index];
                const auto normal            = all_normals[q_index];

                const auto jump_value = (value_m - value_p) * 0.5 * JxW;
                const auto avg_gradient =
                  (normal_gradient_m + normal_gradient_p) * 0.5 * JxW;

                const double sigma = penalty_parameter * penalty_factor;

                phi_m.submit_gradient(-outer(jump_value, normal), q);
                phi_m.submit_value(jump_value * sigma * 2.0 - avg_gradient, q);
              }

            buffer.reinit(cell->get_fe().n_dofs_per_cell());
            phi_m.test_and_sum(buffer,
                               EvaluationFlags::values |
                                 EvaluationFlags::gradients);

            std::vector<types::global_dof_index> local_dofs(
              dof_handler.get_fe().n_dofs_per_cell());
            cell->get_dof_indices(local_dofs);
            constraints.distribute_local_to_global(buffer, local_dofs, dst);

            ptr_q += n_q_points;
          }

  dst.compress(VectorOperation::add);
}

template <int dim, int n_components, typename Number>
void
CouplingOperator<dim, n_components, Number>::add_diagonal_entries(
  VectorType &diagonal) const
{
  unsigned int ptr_q = 0;

  FEPointIntegrator phi_m(mapping,
                          dof_handler.get_fe(),
                          update_values | update_gradients);

  Vector<Number>                                      buffer, diagonal_local;
  std::vector<typename FEPointIntegrator::value_type> all_value_m, all_value_p;

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_0) || (face->boundary_id() == bid_1))
          {
            const unsigned int n_dofs_per_cell =
              cell->get_fe().n_dofs_per_cell();
            const unsigned int n_q_points =
              mortar_manager_q->get_n_points(get_rad(cell, face));

            phi_m.reinit(cell,
                         ArrayView<const Point<dim, Number>>(
                           all_points_ref.data() + ptr_q, n_q_points));

            buffer.reinit(n_dofs_per_cell);
            diagonal_local.reinit(n_dofs_per_cell);
            all_value_m.resize(n_q_points * 2);
            all_value_p.resize(n_q_points * 2);

            for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                  buffer[j] = static_cast<Number>(i == j);

                phi_m.evaluate(buffer,
                               EvaluationFlags::values |
                                 EvaluationFlags::gradients);

                for (const auto q : phi_m.quadrature_point_indices())
                  {
                    const unsigned int q_index = ptr_q + q;

                    AssertIndexRange(q_index, all_normals.size());

                    const auto normal  = all_normals[q_index];
                    const auto value_m = phi_m.get_value(q);
                    const auto gradient_m =
                      contract(phi_m.get_gradient(q), normal);

                    all_value_m[q * 2 + 0] = value_m;
                    all_value_m[q * 2 + 1] = gradient_m;
                  }

                for (const auto q : phi_m.quadrature_point_indices())
                  {
                    const unsigned int q_index = ptr_q + q;

                    const auto value_m           = all_value_m[q * 2 + 0];
                    const auto value_p           = all_value_p[q * 2 + 0];
                    const auto normal_gradient_m = all_value_m[q * 2 + 1];
                    const auto normal_gradient_p = all_value_p[q * 2 + 1];

                    const auto JxW = all_weights[q_index];
                    const auto penalty_parameter =
                      all_penalty_parameter[q_index];
                    const auto normal = all_normals[q_index];

                    const auto jump_value = (value_m - value_p) * 0.5 * JxW;
                    const auto avg_gradient =
                      (normal_gradient_m + normal_gradient_p) * 0.5 * JxW;

                    const double sigma = penalty_parameter * penalty_factor;

                    phi_m.submit_gradient(-outer(jump_value, normal), q);
                    phi_m.submit_value(jump_value * sigma * 2.0 - avg_gradient,
                                       q);
                  }

                phi_m.test_and_sum(buffer,
                                   EvaluationFlags::values |
                                     EvaluationFlags::gradients);

                diagonal_local[i] = buffer[i];
              }

            std::vector<types::global_dof_index> local_dofs(
              dof_handler.get_fe().n_dofs_per_cell());
            cell->get_dof_indices(local_dofs);
            constraints.distribute_local_to_global(diagonal_local,
                                                   local_dofs,
                                                   diagonal);

            ptr_q += n_q_points;
          }

  diagonal.compress(VectorOperation::add);
}

template <int dim, int n_components, typename Number>
void
CouplingOperator<dim, n_components, Number>::add_sparsity_pattern_entries(
  TrilinosWrappers::SparsityPattern &dsp) const
{
  const auto constraints = &constraints_extended;

  const unsigned int n_dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell();

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

template <int dim, int n_components, typename Number>
void
CouplingOperator<dim, n_components, Number>::add_system_matrix_entries(
  TrilinosWrappers::SparseMatrix &system_matrix) const
{
  const auto constraints = &constraints_extended;

  const unsigned int n_dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell();

  std::vector<typename FEPointIntegrator::value_type> all_value_m(
    all_normals.size() * n_dofs_per_cell * 2);
  std::vector<typename FEPointIntegrator::value_type> all_value_p(
    all_normals.size() * n_dofs_per_cell * 2);

  unsigned int ptr_q = 0;

  FEPointIntegrator phi_m(mapping,
                          dof_handler.get_fe(),
                          update_values | update_gradients);
  Vector<Number>    buffer;

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_0) || (face->boundary_id() == bid_1))
          {
            const unsigned int n_q_points =
              mortar_manager_q->get_n_points(get_rad(cell, face));

            phi_m.reinit(cell,
                         ArrayView<const Point<dim, Number>>(
                           all_points_ref.data() + ptr_q, n_q_points));

            buffer.reinit(n_dofs_per_cell);

            for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                  buffer[j] = static_cast<Number>(i == j);

                phi_m.evaluate(buffer,
                               EvaluationFlags::values |
                                 EvaluationFlags::gradients);

                for (const auto q : phi_m.quadrature_point_indices())
                  {
                    const unsigned int q_index = ptr_q + q;

                    AssertIndexRange(q_index, all_normals.size());
                    AssertIndexRange(q_index * 2 + 1, all_value_m.size());

                    const auto normal  = all_normals[q_index];
                    const auto value_m = phi_m.get_value(q);
                    const auto gradient_m =
                      contract(phi_m.get_gradient(q), normal);

                    all_value_m[(q_index * n_dofs_per_cell + i) * 2 + 0] =
                      value_m;
                    all_value_m[(q_index * n_dofs_per_cell + i) * 2 + 1] =
                      gradient_m;
                  }
              }

            ptr_q += n_q_points;
          }

  const unsigned n_q_points =
    Utilities::pow(quadrature.get_tensor_basis()[0].size(), dim - 1);

  partitioner_cell.template export_to_ghosted_array<Number, 0>(
    ArrayView<const Number>(reinterpret_cast<Number *>(all_value_m.data()),
                            all_value_m.size() * n_components),
    ArrayView<Number>(reinterpret_cast<Number *>(all_value_p.data()),
                      all_value_p.size() * n_components),
    dof_handler.get_fe().n_dofs_per_cell() * n_q_points * 2 * n_components);


  ptr_q                 = 0;
  unsigned int ptr_dofs = 0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == bid_0) || (face->boundary_id() == bid_1))
          {
            const unsigned int n_sub_cells =
              mortar_manager_cell->get_n_points(get_rad(cell, face));

            for (unsigned int sc = 0; sc < n_sub_cells; ++sc)
              {
                const unsigned int n_q_points =
                  mortar_manager_q->get_n_points(get_rad(cell, face)) /
                  n_sub_cells;

                phi_m.reinit(cell,
                             ArrayView<const Point<dim, Number>>(
                               all_points_ref.data() + ptr_q, n_q_points));

                for (unsigned int bb = 0; bb < 2; ++bb)
                  {
                    FullMatrix<Number> cell_matrix(n_dofs_per_cell,
                                                   n_dofs_per_cell);

                    for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                      {
                        for (const auto q : phi_m.quadrature_point_indices())
                          {
                            const unsigned int q_index = ptr_q + q;

                            const auto value_m =
                              (bb == 0) ?
                                all_value_m[(q_index * n_dofs_per_cell + i) *
                                              2 +
                                            0] :
                                typename FEPointIntegrator::value_type();
                            const auto value_p =
                              (bb == 1) ?
                                all_value_p[(q_index * n_dofs_per_cell + i) *
                                              2 +
                                            0] :
                                typename FEPointIntegrator::value_type();
                            const auto normal_gradient_m =
                              (bb == 0) ?
                                all_value_m[(q_index * n_dofs_per_cell + i) *
                                              2 +
                                            1] :
                                typename FEPointIntegrator::value_type();
                            const auto normal_gradient_p =
                              (bb == 1) ?
                                all_value_p[(q_index * n_dofs_per_cell + i) *
                                              2 +
                                            1] :
                                typename FEPointIntegrator::value_type();
                            const auto JxW = all_weights[q_index];
                            const auto penalty_parameter =
                              all_penalty_parameter[q_index];
                            const auto normal = all_normals[q_index];

                            const auto jump_value =
                              (value_m - value_p) * 0.5 * JxW;
                            const auto avg_gradient =
                              (normal_gradient_m + normal_gradient_p) * 0.5 *
                              JxW;

                            const double sigma =
                              penalty_parameter * penalty_factor;

                            phi_m.submit_gradient(-outer(jump_value, normal),
                                                  q);
                            phi_m.submit_value(jump_value * sigma * 2.0 -
                                                 avg_gradient,
                                               q);
                          }

                        buffer.reinit(cell->get_fe().n_dofs_per_cell());
                        phi_m.test_and_sum(buffer,
                                           EvaluationFlags::values |
                                             EvaluationFlags::gradients);

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

  AssertDimension(ptr_q, all_normals.size());
  AssertDimension(ptr_dofs, dof_indices.size());
}

#endif
