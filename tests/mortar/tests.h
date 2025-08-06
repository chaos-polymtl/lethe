// SPDX-FileCopyrightText: Copyright (c) 2021 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

// Lethe
#include <core/mortar_coupling_manager.h>

using namespace dealii;

/**
 * @brief Create grid with hyper_shell and hyper_shell geometries
 * Only the merged triangulation (and resulting grid) are stored
 */
template <int dim>
void
hyper_shell_with_hyper_shell(const double        radius,
                             Triangulation<dim> &tria,
                             const double        tolerance = 0.0)
{
  const double r1_i = radius * 0.25;
  const double r1_o = radius * 0.5;
  const double r2_i = radius * 0.5;
  const double r2_o = radius * 1.0;

  // inner domain triangulation
  Triangulation<2> circle_one;
  GridGenerator::hyper_shell(circle_one, Point<2>(), r1_i, r1_o, 6, true);
  // outer domain triangulation
  Triangulation<2> circle_two;
  GridGenerator::hyper_shell(circle_two, Point<2>(), r2_i, r2_o, 6, true);

  // shift boundary id of circle two
  for (const auto &face : circle_two.active_face_iterators())
    if (face->at_boundary())
      face->set_boundary_id(face->boundary_id() + 2);

  // create unique triangulation
  Triangulation<2> temp;
  GridGenerator::merge_triangulations(
    circle_one, circle_two, temp, tolerance, true, true);
  temp.set_manifold(0, SphericalManifold<2>(Point<2>()));

  if constexpr (dim == 3)
    GridGenerator::extrude_triangulation(temp, 3, radius, tria, true);
  if constexpr (dim == 2)
    tria.copy_triangulation(temp);

  // store manifolds in merged triangulation
  if (dim == 2)
    tria.set_manifold(0, SphericalManifold<dim>(Point<dim>()));
  else
    tria.set_manifold(0, CylindricalManifold<dim>(2));
}

/**
 * @brief Create grid with hyper_ball_balanced and hyper_cube_with_cylindrical_hole geometries
 * Only the merged triangulation (and resulting grid) are stored
 */
template <int dim>
void
hyper_cube_with_cylindrical_hole(const double        radius,
                                 const double        outer_radius,
                                 const double        rotate,
                                 Triangulation<dim> &tria)
{
  Triangulation<2> tria_0, tria_1;

  // inner domain triangulation
  GridGenerator::hyper_ball_balanced(tria_0, {}, radius);
  GridTools::rotate(rotate, tria_0);

  // outer domain triangulation
  GridGenerator::hyper_cube_with_cylindrical_hole(
    tria_1, radius, outer_radius, 5.0, 1.0, true);

  // shift boundary IDs # in outer grid
  for (const auto &face : tria_1.active_face_iterators())
    if (face->at_boundary())
      {
        face->set_boundary_id(face->boundary_id() + 1);
        face->set_manifold_id(face->manifold_id() + 2);
      }

  // create unique triangulation
  Triangulation<2> temp;
  GridGenerator::merge_triangulations(tria_0, tria_1, temp, 0, true, true);
  temp.set_manifold(0, SphericalManifold<2>(Point<2>()));
  temp.set_manifold(1, FlatManifold<2>());
  temp.set_manifold(2, SphericalManifold<2>(Point<2>()));

  if constexpr (dim == 3)
    GridGenerator::extrude_triangulation(temp, 2, radius, tria, true);
  if constexpr (dim == 2)
    tria.copy_triangulation(temp);

  // store manifolds in merged triangulation
  if (dim == 2)
    tria.set_manifold(0, SphericalManifold<dim>(Point<dim>()));
  else
    tria.set_manifold(0, CylindricalManifold<dim>(2));
  tria.set_manifold(1, FlatManifold<dim>());
  if (dim == 2)
    tria.set_manifold(2, SphericalManifold<dim>(Point<dim>()));
  else
    tria.set_manifold(2, CylindricalManifold<dim>(2));
}

/**
 * @brief Create grid with hyper_ball_balanced and hyper_cube_with_cylindrical_hole geometries
 * Only the merged triangulation (and resulting grid) are stored
 * Consider a tolerance of 1e-9 when merging boundary points
 */
template <int dim>
void
hyper_cube_with_cylindrical_hole_with_tolerance(const double radius,
                                                const double outer_radius,
                                                const double rotate,
                                                Triangulation<dim> &tria)
{
  Triangulation<dim> tria_0, tria_1;

  // inner domain triangulation
  GridGenerator::hyper_ball_balanced(tria_0, {}, radius);
  GridTools::rotate(rotate, tria_0);

  // outer domain triangulation
  GridGenerator::hyper_cube_with_cylindrical_hole(
    tria_1, radius, outer_radius, 5.0, 1.0, true);

  // shift boundary IDs # in outer grid
  for (const auto &face : tria_1.active_face_iterators())
    if (face->at_boundary())
      {
        face->set_boundary_id(face->boundary_id() + 1);
        face->set_manifold_id(face->manifold_id() + 2);
      }

  // crete unique triangulation
  GridGenerator::merge_triangulations(tria_0, tria_1, tria, 1e-9, true, false);
  // store manifolds in merged triangulaiton
  tria.set_manifold(0, SphericalManifold<dim>(Point<dim>()));
  tria.set_manifold(1, FlatManifold<dim>());
  tria.set_manifold(2, SphericalManifold<dim>(Point<dim>()));
}

/**
 * @brief TODO
 */
void
split_hyper_cube(Triangulation<2> &tria,
                 const double      left,
                 const double      right,
                 const double      mid,
                 const double      tolerance = 0.0)
{
  Triangulation<2> tria_0, tria_1;

  // inner domain triangulation
  GridGenerator::subdivided_hyper_rectangle(tria_0,
                                            std::vector<unsigned int>{1, 1},
                                            Point<2>(left, left),
                                            Point<2>(mid, right),
                                            true);

  // outer domain triangulation
  GridGenerator::subdivided_hyper_rectangle(tria_1,
                                            std::vector<unsigned int>{1, 1},
                                            Point<2>(mid, left),
                                            Point<2>(right, right),
                                            true);

  // shift boundary IDs # in outer grid
  for (const auto &face : tria_1.active_face_iterators())
    if (face->at_boundary())
      face->set_boundary_id(face->boundary_id() + 4);

  // create unique triangulation
  GridGenerator::merge_triangulations(
    tria_0, tria_1, tria, tolerance, true, true);
}

/**
 * @brief TODO
 */
void
split_hyper_cube(Triangulation<1> &tria,
                 const double      left,
                 const double      right,
                 const double      mid,
                 const double      tolerance = 0.0)
{
  Triangulation<1> tria_0, tria_1;

  // inner domain triangulation
  GridGenerator::subdivided_hyper_rectangle(
    tria_0, std::vector<unsigned int>{1}, Point<1>(left), Point<1>(mid), true);

  // outer domain triangulation
  GridGenerator::subdivided_hyper_rectangle(
    tria_1, std::vector<unsigned int>{1}, Point<1>(mid), Point<1>(right), true);

  // shift boundary IDs # in outer grid
  for (const auto &cell : tria_1.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary())
        face->set_boundary_id(face->boundary_id() + 2);

  // create unique triangulation
  GridGenerator::merge_triangulations(
    tria_0, tria_1, tria, tolerance, true, true);
}

/**
 * @brief TODO
 */
template <int dim>
void
split_hyper_cube(Triangulation<dim> &tria,
                 const double        left      = 0.0,
                 const double        right     = 1.0,
                 const double        tolerance = 0.0)
{
  split_hyper_cube(tria, left, right, (left + right) / 2.0, tolerance);
}

template <int dim>
Quadrature<dim>
construct_quadrature(const Quadrature<dim> &quad)
{
  const double oversampling_factor = 2.0; // make parameter

  for (unsigned int i = 1; i <= 10; ++i)
    if (quad == QGauss<dim>(i))
      return QGauss<dim>(i * oversampling_factor);

  AssertThrow(false, ExcNotImplemented());

  return quad;
}

template <int dim, typename Number>
class CouplingEvaluationStokes : public CouplingEvaluationBase<dim, Number>
{
public:
  using FEPointIntegratorU = FEPointEvaluation<dim, dim, dim, Number>;
  using FEPointIntegratorP = FEPointEvaluation<1, dim, dim, Number>;

  using u_value_type = typename FEPointIntegratorU::value_type;

  CouplingEvaluationStokes(const Mapping<dim>    &mapping,
                           const DoFHandler<dim> &dof_handler,
                           const bool weak_pressure_gradient_term   = true,
                           const bool weak_velocity_divergence_term = true)
    : fe_sub_u(dof_handler.get_fe().base_element(
                 dof_handler.get_fe().component_to_base_index(0).first),
               dim)
    , fe_sub_p(dof_handler.get_fe().base_element(
                 dof_handler.get_fe().component_to_base_index(dim).first),
               1)
    , phi_u_m(mapping, fe_sub_u, update_values | update_gradients)
    , phi_p_m(mapping, fe_sub_p, update_values)
    , weak_pressure_gradient_term(weak_pressure_gradient_term)
    , weak_velocity_divergence_term(weak_velocity_divergence_term)
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

  unsigned int
  data_size() const override
  {
    return 4 * dim;
  }

  const std::vector<unsigned int> &
  get_relevant_dof_indices() const override
  {
    return relevant_dof_indices;
  }

  void
  local_reinit(const typename Triangulation<dim>::cell_iterator &cell,
               const ArrayView<const Point<dim, Number>> &points) const override
  {
    this->phi_u_m.reinit(cell, points);
    this->phi_p_m.reinit(cell, points);
  }

  void
  local_evaluate(const CouplingEvaluationData<dim, Number> &data,
                 const Vector<Number>                      &buffer,
                 const unsigned int                         ptr_q,
                 const unsigned int                         q_stride,
                 Number *all_value_m) const override
  {
    AssertDimension(buffer.size(),
                    fe_sub_u.n_dofs_per_cell() + fe_sub_p.n_dofs_per_cell());

    ArrayView<const Number> buffer_u(buffer.data() + 0,
                                     fe_sub_u.n_dofs_per_cell());
    ArrayView<const Number> buffer_p(buffer.data() + fe_sub_u.n_dofs_per_cell(),
                                     fe_sub_p.n_dofs_per_cell());

    this->phi_u_m.evaluate(buffer_u,
                           EvaluationFlags::values |
                             EvaluationFlags::gradients);
    this->phi_p_m.evaluate(buffer_p, EvaluationFlags::values);

    for (const auto q : this->phi_u_m.quadrature_point_indices())
      {
        const unsigned int q_index = ptr_q + q;

        const auto normal = data.all_normals[q_index];

        const auto value_m    = this->phi_u_m.get_value(q);
        const auto gradient_m = contract(this->phi_u_m.get_gradient(q), normal);
        const auto p_value_m  = this->phi_p_m.get_value(q) * normal;

        BufferRW<Number> buffer_m(all_value_m, q * 4 * dim * q_stride);

        buffer_m.write(value_m);
        buffer_m.write(gradient_m);
        buffer_m.write(p_value_m);
      }
  }

  void
  local_integrate(const CouplingEvaluationData<dim, Number> &data,
                  Vector<Number>                            &buffer,
                  const unsigned int                         ptr_q,
                  const unsigned int                         q_stride,
                  Number                                    *all_value_m,
                  Number *all_value_p) const override
  {
    for (const auto q : this->phi_u_m.quadrature_point_indices())
      {
        const unsigned int q_index = ptr_q + q;

        BufferRW<Number> buffer_m(all_value_m, q * 4 * dim * q_stride);
        BufferRW<Number> buffer_p(all_value_p, q * 4 * dim * q_stride);

        const auto value_m           = buffer_m.template read<u_value_type>();
        const auto value_p           = buffer_p.template read<u_value_type>();
        const auto normal_gradient_m = buffer_m.template read<u_value_type>();
        const auto normal_gradient_p = buffer_p.template read<u_value_type>();
        const auto normal_p_value_m  = buffer_m.template read<u_value_type>();
        const auto normal_p_value_p  = buffer_p.template read<u_value_type>();

        const auto JxW               = data.all_weights[q_index];
        const auto penalty_parameter = data.all_penalty_parameter[q_index];
        const auto normal            = data.all_normals[q_index];

        const auto u_value_avg  = (value_m + value_p) * 0.5;
        const auto u_value_jump = value_m - value_p;
        const auto u_gradient_avg =
          (normal_gradient_m - normal_gradient_p) * 0.5;
        const auto p_value_avg = (normal_p_value_m - normal_p_value_p) * 0.5;

        typename FEPointIntegratorU::value_type u_normal_gradient_avg_result =
          {};
        typename FEPointIntegratorU::value_type u_value_jump_result = {};
        typename FEPointIntegratorP::value_type p_value_jump_result = {};

        const double sigma = penalty_parameter * data.penalty_factor;

        if (true /*Laplace term*/)
          {
            // - (n avg(∇v), jump(u))
            u_normal_gradient_avg_result -= u_value_jump;

            // - (jump(v), avg(∇u) n)
            u_value_jump_result -= u_gradient_avg;

            // + (jump(v), σ jump(u))
            u_value_jump_result += sigma * u_value_jump;
          }

        if (weak_pressure_gradient_term)
          {
            // + (jump(v), avg(p) n)
            u_value_jump_result += p_value_avg;
          }
        else
          {
            // nothing to do
          }

        if (weak_velocity_divergence_term)
          {
            // + (jump(q), avg(u) n)
            if constexpr (dim == 1)
              p_value_jump_result += u_value_avg * normal[0];
            else
              p_value_jump_result += u_value_avg * normal;
          }
        else
          {
            // - (avg(q), jump(u) n)
            if constexpr (dim == 1)
              p_value_jump_result -= 0.5 * u_value_jump * normal[0];
            else
              p_value_jump_result -= 0.5 * u_value_jump * normal;
          }

        phi_u_m.submit_gradient(outer(u_normal_gradient_avg_result, normal) *
                                  0.5 * JxW,
                                q);
        phi_u_m.submit_value(u_value_jump_result * JxW, q);
        phi_p_m.submit_value(p_value_jump_result * JxW, q);
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

  const FESystem<dim>        fe_sub_u;
  const FESystem<dim>        fe_sub_p;
  mutable FEPointIntegratorU phi_u_m;
  mutable FEPointIntegratorP phi_p_m;

  const bool weak_pressure_gradient_term;
  const bool weak_velocity_divergence_term;

  std::vector<unsigned int> relevant_dof_indices;
};

/**
 * @brief Base class for the operator base
 */
template <int dim,
          int n_components,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class OperatorBase : public Subscriptor
{
public:
  using FECellIntegrator =
    FEEvaluation<dim, -1, 0, n_components, Number, VectorizedArrayType>;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  OperatorBase() = default;

  OperatorBase(const Mapping<dim>              &mapping,
               const DoFHandler<dim>           &dof_handler,
               const AffineConstraints<Number> &constraints,
               const Quadrature<dim>           &quadrature)
  {
    reinit(mapping, dof_handler, constraints, quadrature);
  }

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<dim>           &quadrature)
  {
    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
    data.mapping_update_flags =
      update_quadrature_points | update_gradients | update_values;

    matrix_free.reinit(mapping, dof_handler, constraints, quadrature, data);

    valid_system = false;
  }

  /**
   * @brief Create coupling operator
   */
  void
  add_coupling(const std::shared_ptr<MortarManagerBase<dim>> mortar_manager,
               const unsigned int                            bid_0,
               const unsigned int                            bid_1,
               const double                                  sip_factor = 1.0)
  {
    const std::shared_ptr<CouplingEvaluationBase<dim, Number>>
      coupling_evaluator =
        std::make_shared<CouplingEvaluationSIPG<dim, n_components, Number>>(
          *matrix_free.get_mapping_info().mapping,
          matrix_free.get_dof_handler());

    coupling_operator = std::make_shared<CouplingOperator<dim, Number>>(
      *matrix_free.get_mapping_info().mapping,
      matrix_free.get_dof_handler(),
      matrix_free.get_affine_constraints(),
      coupling_evaluator,
      mortar_manager,
      bid_0,
      bid_1,
      sip_factor);
  }

  virtual types::global_dof_index
  m() const
  {
    if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
      return this->matrix_free.get_dof_handler().n_dofs(
        this->matrix_free.get_mg_level());
    else
      return this->matrix_free.get_dof_handler().n_dofs();
  }

  Number
  el(unsigned int, unsigned int) const
  {
    DEAL_II_NOT_IMPLEMENTED();
    return 0;
  }

  void
  initialize_dof_vector(VectorType &dst) const
  {
    matrix_free.initialize_dof_vector(dst);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    src.update_ghost_values();

    matrix_free.cell_loop(
      &OperatorBase<dim, n_components, Number, VectorizedArrayType>::
        do_vmult_cell,
      this,
      dst,
      src,
      true);

    // apply coupling terms
    if (coupling_operator)
      {
        // apply constraints
        // TODO: only apply relevant constraints
        const auto &constraints = coupling_operator->get_affine_constraints();
        constraints.distribute(const_cast<VectorType &>(src));
        src.update_ghost_values();

        coupling_operator->vmult_add(dst, src);

        constraints.set_zero(const_cast<VectorType &>(src));
        constraints.set_zero(dst);
      }

    src.zero_out_ghost_values();
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    vmult(dst, src);
  }

  void
  compute_inverse_diagonal(VectorType &diagonal) const
  {
    matrix_free.initialize_dof_vector(diagonal);

    MatrixFreeTools::compute_diagonal(
      matrix_free,
      diagonal,
      &OperatorBase<dim, n_components, Number, VectorizedArrayType>::
        do_vmult_cell_single,
      this);

    // add coupling terms
    if (coupling_operator)
      coupling_operator->add_diagonal_entries(diagonal);

    for (auto &i : diagonal)
      i = (i != 0.0) ? (1.0 / i) : 1.0;
  }

  const TrilinosWrappers::SparseMatrix &
  get_system_matrix() const
  {
    initialize_system_matrix();

    return system_matrix;
  }

  void
  initialize_system_matrix() const
  {
    const auto &dof_handler = matrix_free.get_dof_handler();

    const auto constraints = (coupling_operator != nullptr) ?
                               (&coupling_operator->get_affine_constraints()) :
                               (&matrix_free.get_affine_constraints());

    if (system_matrix.m() == 0 || system_matrix.n() == 0)
      {
        system_matrix.clear();

        TrilinosWrappers::SparsityPattern dsp;

        dsp.reinit(dof_handler.locally_owned_dofs(),
                   dof_handler.get_mpi_communicator());

        DoFTools::make_sparsity_pattern(dof_handler, dsp, *constraints);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_sparsity_pattern_entries(dsp);

        dsp.compress();

        system_matrix.reinit(dsp);
      }

    if (this->valid_system == false)
      {
        system_matrix = 0.0;

        MatrixFreeTools::compute_matrix(
          matrix_free,
          *constraints,
          system_matrix,
          &OperatorBase<dim, n_components, Number, VectorizedArrayType>::
            do_vmult_cell_single,
          this);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_system_matrix_entries(system_matrix);

        system_matrix.compress(VectorOperation::add);

        this->valid_system = true;
      }
  }


protected:
  void
  do_vmult_cell(const MatrixFree<dim, Number>               &data,
                VectorType                                  &dst,
                const VectorType                            &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FECellIntegrator phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);

        do_vmult_cell_single(phi);
        phi.distribute_local_to_global(dst);
      }
  }

  virtual void
  do_vmult_cell_single(FECellIntegrator &phi) const = 0;

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  mutable TrilinosWrappers::SparseMatrix       system_matrix;
  mutable bool                                 valid_system;

  std::shared_ptr<CouplingOperator<dim, Number>> coupling_operator;

private:
};


template <int dim,
          int n_components,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class PoissonOperator
  : public OperatorBase<dim, n_components, Number, VectorizedArrayType>
{
public:
  PoissonOperator() = default;

  PoissonOperator(const MappingQ<dim>             &mapping,
                  const DoFHandler<dim>           &dof_handler,
                  const AffineConstraints<Number> &constraints,
                  const Quadrature<dim>           &quadrature)
    : OperatorBase<dim, n_components, Number, VectorizedArrayType>(mapping,
                                                                   dof_handler,
                                                                   constraints,
                                                                   quadrature)
  {}

  void
  do_vmult_cell_single(
    typename OperatorBase<dim, n_components, Number, VectorizedArrayType>::
      FECellIntegrator &phi) const override
  {
    phi.evaluate(EvaluationFlags::gradients);

    for (unsigned int q = 0; q < phi.n_q_points; ++q)
      phi.submit_gradient(phi.get_gradient(q), q);

    phi.integrate(EvaluationFlags::gradients);
  }
};


template <int dim,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class StokesOperator
  : public OperatorBase<dim, dim + 1, Number, VectorizedArrayType>
{
public:
  using BaseClass = OperatorBase<dim, dim + 1, Number, VectorizedArrayType>;
  using FECellIntegrator = typename BaseClass::FECellIntegrator;

  StokesOperator(const MappingQ<dim>             &mapping,
                 const DoFHandler<dim>           &dof_handler,
                 const AffineConstraints<Number> &constraints,
                 const Quadrature<dim>           &quadrature,
                 const double                     delta_1_scaling)
    : BaseClass(mapping, dof_handler, constraints, quadrature)
    , delta_1_scaling(delta_1_scaling)
  {}

  void
  do_vmult_cell_single(FECellIntegrator &phi) const override
  {
    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    const auto cell = phi.get_current_cell_index();

    VectorizedArrayType delta_1;
    for (unsigned int v = 0;
         v < this->matrix_free.n_active_entries_per_cell_batch(cell);
         ++v)
      delta_1[v] =
        delta_1_scaling *
        this->matrix_free.get_cell_iterator(cell, v)->minimum_vertex_distance();

    for (unsigned int q = 0; q < phi.n_q_points; ++q)
      {
        typename FECellIntegrator::value_type    value_result    = {};
        typename FECellIntegrator::gradient_type gradient_result = {};

        const auto value    = phi.get_value(q);
        const auto gradient = phi.get_gradient(q);

        const VectorizedArray<Number>                 p_value = value[dim];
        const Tensor<1, dim, VectorizedArray<Number>> p_gradient =
          gradient[dim];

        Tensor<2, dim, VectorizedArray<Number>> u_gradient;

        for (unsigned int d = 0; d < dim; ++d)
          u_gradient[d] = gradient[d];

        // a)     (ε(v), 2νε(u))
        if (true)
          symm_scalar_product_add(gradient_result,
                                  u_gradient,
                                  VectorizedArrayType(2.0));
        else
          for (unsigned int d = 0; d < dim; ++d)
            gradient_result[d] += u_gradient[d];

        // b)   - (div(v), p)
        for (unsigned int d = 0; d < dim; ++d)
          gradient_result[d][d] -= p_value;

        // c)     (q, div(u))
        for (unsigned int d = 0; d < dim; ++d)
          value_result[dim] -= u_gradient[d][d];

        // d) δ_1 (∇q, ∇p)
        gradient_result[dim] = -delta_1 * p_gradient;

        phi.submit_value(value_result, q);
        phi.submit_gradient(gradient_result, q);
      }

    phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  }

private:
  const double delta_1_scaling;
};


/**
 * @brief Base class for the operator base
 */
template <int dim,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class GeneralStokesOperator : public Subscriptor
{
public:
  using FECellIntegratorU =
    FEEvaluation<dim, -1, 0, dim, Number, VectorizedArrayType>;
  using FECellIntegratorP =
    FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;


  GeneralStokesOperator(const MappingQ<dim>             &mapping,
                        const DoFHandler<dim>           &dof_handler,
                        const AffineConstraints<Number> &constraints,
                        const Quadrature<dim>           &quadrature,
                        const double                     delta_1_scaling,
                        const bool weak_velocity_divergence_term = false)
    : delta_1_scaling(delta_1_scaling)
    , weak_velocity_divergence_term(weak_velocity_divergence_term)
  {
    reinit(mapping, dof_handler, constraints, quadrature);
  }

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<dim>           &quadrature)
  {
    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
    data.mapping_update_flags =
      update_quadrature_points | update_gradients | update_values;

    matrix_free.reinit(mapping, dof_handler, constraints, quadrature, data);

    valid_system = false;
  }


  /**
   * @brief Create coupling operator
   *
   * @param[in] n_subdivisions Number of cells at the interface between inner
   * and outer domains
   * @param[in] radius Radius at the interface between inner and outer domains
   * @param[in] rotate_pi Rotation angle for the inner domain
   * @param[in] bid_0 Boundary ID of inner domain (rotor)
   * @param[in] bid_1 Boundary ID of outer domain (stator)
   * @param[in] sip_factor Penalty factor (akin to symmetric interior penalty
   * factor in SIPG)
   */
  void
  add_coupling(const std::shared_ptr<MortarManagerBase<dim>> mortar_manager,
               const unsigned int                            bid_0,
               const unsigned int                            bid_1,
               const double                                  sip_factor = 1.0)
  {
    const bool is_p_disc = matrix_free.get_dof_handler()
                             .get_fe()
                             .base_element(matrix_free.get_dof_handler()
                                             .get_fe()
                                             .component_to_base_index(dim)
                                             .first)
                             .n_dofs_per_vertex() == 0;

    const std::shared_ptr<CouplingEvaluationBase<dim, Number>>
      coupling_evaluator =
        std::make_shared<CouplingEvaluationStokes<dim, Number>>(
          *matrix_free.get_mapping_info().mapping,
          matrix_free.get_dof_handler(),
          !is_p_disc,
          weak_velocity_divergence_term);

    coupling_operator = std::make_shared<CouplingOperator<dim, Number>>(
      *matrix_free.get_mapping_info().mapping,
      matrix_free.get_dof_handler(),
      matrix_free.get_affine_constraints(),
      coupling_evaluator,
      mortar_manager,
      bid_0,
      bid_1,
      sip_factor);
  }

  virtual types::global_dof_index
  m() const
  {
    if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
      return this->matrix_free.get_dof_handler().n_dofs(
        this->matrix_free.get_mg_level());
    else
      return this->matrix_free.get_dof_handler().n_dofs();
  }

  Number
  el(unsigned int, unsigned int) const
  {
    DEAL_II_NOT_IMPLEMENTED();
    return 0;
  }

  void
  initialize_dof_vector(VectorType &dst) const
  {
    matrix_free.initialize_dof_vector(dst);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    src.update_ghost_values();

    matrix_free.cell_loop(
      &GeneralStokesOperator<dim, Number, VectorizedArrayType>::do_vmult_cell,
      this,
      dst,
      src,
      true);

    // apply coupling terms
    if (coupling_operator)
      coupling_operator->vmult_add(dst, src);

    src.zero_out_ghost_values();
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    vmult(dst, src);
  }

  void
  compute_inverse_diagonal(VectorType &diagonal) const
  {
    matrix_free.initialize_dof_vector(diagonal);

    MatrixFreeTools::internal::
      ComputeMatrixScratchData<dim, VectorizedArray<Number>, false>
        data_cell;

    data_cell.dof_numbers               = {0, 0};
    data_cell.quad_numbers              = {0, 0};
    data_cell.n_components              = {dim, 1};
    data_cell.first_selected_components = {0, dim};
    data_cell.batch_type                = {0, 0};

    data_cell
      .op_create = [&](const std::pair<unsigned int, unsigned int> &range) {
      std::vector<
        std::unique_ptr<FEEvaluationData<dim, VectorizedArray<Number>, false>>>
        phi;

      phi.emplace_back(
        std::make_unique<FECellIntegratorU>(matrix_free, range, 0, 0, 0));

      phi.emplace_back(
        std::make_unique<FECellIntegratorP>(matrix_free, range, 0, 0, dim));

      return phi;
    };

    data_cell.op_reinit = [](auto &phi, const unsigned batch) {
      static_cast<FECellIntegratorU &>(*phi[0]).reinit(batch);
      static_cast<FECellIntegratorP &>(*phi[1]).reinit(batch);
    };

    data_cell.op_compute = [&](auto &phi) {
      auto &phi_0 = static_cast<FECellIntegratorU &>(*phi[0]);
      auto &phi_1 = static_cast<FECellIntegratorP &>(*phi[1]);

      do_vmult_cell_single(phi_0, phi_1);
    };

    std::vector<VectorType *> diagonal_global_components(1);
    diagonal_global_components[0] = &diagonal;

    MatrixFreeTools::internal::compute_diagonal(
      matrix_free, data_cell, {}, {}, diagonal, diagonal_global_components);

    // add coupling terms
    if (coupling_operator)
      coupling_operator->add_diagonal_entries(diagonal);

    for (auto &i : diagonal)
      i = (i != 0.0) ? (1.0 / i) : 1.0;
  }

  const TrilinosWrappers::SparseMatrix &
  get_system_matrix() const
  {
    initialize_system_matrix();

    return system_matrix;
  }

  void
  initialize_system_matrix() const
  {
    const auto &dof_handler = matrix_free.get_dof_handler();

    auto constraints = &matrix_free.get_affine_constraints();

    AffineConstraints<Number> affine_constraints_tmp;

    if (coupling_operator)
      {
        affine_constraints_tmp.copy_from(
          coupling_operator->get_affine_constraints());

        affine_constraints_tmp.close();
      }

    if (system_matrix.m() == 0 || system_matrix.n() == 0)
      {
        system_matrix.clear();

        TrilinosWrappers::SparsityPattern dsp;

        dsp.reinit(dof_handler.locally_owned_dofs(),
                   dof_handler.get_mpi_communicator());

        DoFTools::make_sparsity_pattern(dof_handler, dsp, *constraints);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_sparsity_pattern_entries(dsp);

        dsp.compress();

        system_matrix.reinit(dsp);
      }

    if (this->valid_system == false)
      {
        system_matrix = 0.0;

        MatrixFreeTools::internal::
          ComputeMatrixScratchData<dim, VectorizedArray<Number>, false>
            data_cell;

        data_cell.dof_numbers               = {0, 0};
        data_cell.quad_numbers              = {0, 0};
        data_cell.n_components              = {dim, 1};
        data_cell.first_selected_components = {0, dim};
        data_cell.batch_type                = {0, 0};

        data_cell.op_create =
          [&](const std::pair<unsigned int, unsigned int> &range) {
            std::vector<std::unique_ptr<
              FEEvaluationData<dim, VectorizedArray<Number>, false>>>
              phi;

            phi.emplace_back(
              std::make_unique<FECellIntegratorU>(matrix_free, range, 0, 0, 0));

            phi.emplace_back(std::make_unique<FECellIntegratorP>(
              matrix_free, range, 0, 0, dim));

            return phi;
          };

        data_cell.op_reinit = [](auto &phi, const unsigned batch) {
          static_cast<FECellIntegratorU &>(*phi[0]).reinit(batch);
          static_cast<FECellIntegratorP &>(*phi[1]).reinit(batch);
        };

        data_cell.op_compute = [&](auto &phi) {
          auto &phi_0 = static_cast<FECellIntegratorU &>(*phi[0]);
          auto &phi_1 = static_cast<FECellIntegratorP &>(*phi[1]);

          do_vmult_cell_single(phi_0, phi_1);
        };

        MatrixFreeTools::internal::compute_matrix(
          matrix_free, *constraints, data_cell, {}, {}, system_matrix);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_system_matrix_entries(system_matrix);

        system_matrix.compress(VectorOperation::add);

        this->valid_system = true;
      }
  }


private:
  void
  do_vmult_cell(const MatrixFree<dim, Number>               &matrix_free,
                VectorType                                  &dst,
                const VectorType                            &src,
                const std::pair<unsigned int, unsigned int> &range) const
  {
    FECellIntegratorU integrator_u(matrix_free, range, 0, 0, 0);
    FECellIntegratorP integrator_p(matrix_free, range, 0, 0, dim);
    for (unsigned cell = range.first; cell < range.second; ++cell)
      {
        integrator_u.reinit(cell);
        integrator_p.reinit(cell);

        integrator_u.read_dof_values(src);
        integrator_p.read_dof_values(src);
        do_vmult_cell_single(integrator_u, integrator_p);
        integrator_u.distribute_local_to_global(dst);
        integrator_p.distribute_local_to_global(dst);
      }
  }

  void
  do_vmult_cell_single(FECellIntegratorU &phi_u, FECellIntegratorP &phi_p) const
  {
    phi_u.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    const auto cell = phi_u.get_current_cell_index();

    VectorizedArrayType delta_1;
    for (unsigned int v = 0;
         v < this->matrix_free.n_active_entries_per_cell_batch(cell);
         ++v)
      delta_1[v] =
        delta_1_scaling *
        this->matrix_free.get_cell_iterator(cell, v)->minimum_vertex_distance();

    for (unsigned int q = 0; q < phi_u.n_q_points; ++q)
      {
        typename FECellIntegratorP::value_type    p_value_result    = {};
        typename FECellIntegratorP::gradient_type p_gradient_result = {};
        typename FECellIntegratorU::value_type    u_value_result    = {};
        typename FECellIntegratorU::gradient_type u_gradient_result = {};

        const auto p_value    = phi_p.get_value(q);
        const auto p_gradient = phi_p.get_gradient(q);

        const auto u_value    = phi_u.get_value(q);
        const auto u_gradient = phi_u.get_gradient(q);

        // (ε(v), 2νε(u))
        if (false)
          {
            symm_scalar_product_add(u_gradient_result,
                                    u_gradient,
                                    VectorizedArrayType(2.0));
          }
        else
          {
            u_gradient_result = u_gradient;
          }

        // - (div(v), p)
        if constexpr (dim == 1)
          u_gradient_result[0] -= p_value;
        else
          for (unsigned int d = 0; d < dim; ++d)
            u_gradient_result[d][d] -= p_value;

        if (weak_velocity_divergence_term)
          {
            // - (∇q, u)
            if constexpr (dim == 1)
              p_gradient_result[0] -= u_value;
            else
              p_gradient_result -= u_value;
          }
        else
          {
            // + (q, div(u))
            if constexpr (dim == 1)
              p_value_result += u_gradient[0];
            else
              for (unsigned int d = 0; d < dim; ++d)
                p_value_result += u_gradient[d][d];
          }

        // δ_1 (∇q, ∇p)
        if (delta_1_scaling != 0.0)
          p_gradient_result += delta_1 * p_gradient;

        phi_p.submit_value(p_value_result, q);
        phi_p.submit_gradient(p_gradient_result, q);

        phi_u.submit_value(u_value_result, q);
        phi_u.submit_gradient(u_gradient_result, q);
      }

    phi_u.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  }

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  mutable TrilinosWrappers::SparseMatrix       system_matrix;
  mutable bool                                 valid_system;

  std::shared_ptr<CouplingOperator<dim, Number>> coupling_operator;

  const double delta_1_scaling;
  const bool   weak_velocity_divergence_term;
};



/**
 * @brief Base class for the operator base
 */
template <int dim,
          int n_components,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class PoissonOperatorDG : public Subscriptor
{
public:
  using FECellIntegrator =
    FEEvaluation<dim, -1, 0, n_components, Number, VectorizedArrayType>;
  using FEFaceIntegrator =
    FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  PoissonOperatorDG() = default;

  PoissonOperatorDG(const Mapping<dim>              &mapping,
                    const DoFHandler<dim>           &dof_handler,
                    const AffineConstraints<Number> &constraints,
                    const Quadrature<dim>           &quadrature)
  {
    reinit(mapping, dof_handler, constraints, quadrature);
  }

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<dim>           &quadrature)
  {
    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
    data.mapping_update_flags =
      update_quadrature_points | update_gradients | update_JxW_values;
    data.mapping_update_flags_inner_faces =
      update_values | update_gradients | update_JxW_values;
    data.mapping_update_flags_boundary_faces =
      update_values | update_gradients | update_JxW_values;

    matrix_free.reinit(mapping, dof_handler, constraints, quadrature, data);

    valid_system = false;

    compute_penalty_parameters();

    panalty_factor = compute_pentaly_factor(dof_handler.get_fe().degree, 1.0);
  }

  /**
   * @brief Create coupling operator
   */
  void
  add_coupling(const unsigned int n_subdivisions,
               const double       radius,
               const double       rotate_pi,
               const unsigned int bid_0,
               const unsigned int bid_1,
               const double       sip_factor = 1.0)
  {
    const std::shared_ptr<MortarManagerBase<dim>> mortar_manager =
      std::make_shared<MortarManagerCircle<dim>>(n_subdivisions,
                                                 radius,
                                                 matrix_free.get_quadrature(),
                                                 rotate_pi);

    const std::shared_ptr<CouplingEvaluationBase<dim, Number>>
      coupling_evaluator =
        std::make_shared<CouplingEvaluationSIPG<dim, n_components, Number>>(
          *matrix_free.get_mapping_info().mapping,
          matrix_free.get_dof_handler());

    coupling_operator = std::make_shared<CouplingOperator<dim, Number>>(
      *matrix_free.get_mapping_info().mapping,
      matrix_free.get_dof_handler(),
      matrix_free.get_affine_constraints(),
      coupling_evaluator,
      mortar_manager,
      bid_0,
      bid_1,
      sip_factor);
  }

  virtual types::global_dof_index
  m() const
  {
    if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
      return this->matrix_free.get_dof_handler().n_dofs(
        this->matrix_free.get_mg_level());
    else
      return this->matrix_free.get_dof_handler().n_dofs();
  }

  Number
  el(unsigned int, unsigned int) const
  {
    DEAL_II_NOT_IMPLEMENTED();
    return 0;
  }

  void
  initialize_dof_vector(VectorType &dst) const
  {
    matrix_free.initialize_dof_vector(dst);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    src.update_ghost_values();

    matrix_free.loop(
      &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
        do_vmult_cell,
      &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
        local_apply_face,
      &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
        local_apply_boundary,
      this,
      dst,
      src,
      true);

    // apply coupling terms
    if (coupling_operator)
      coupling_operator->vmult_add(dst, src);

    src.zero_out_ghost_values();
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    vmult(dst, src);
  }

  void
  compute_inverse_diagonal(VectorType &diagonal) const
  {
    matrix_free.initialize_dof_vector(diagonal);

    MatrixFreeTools::compute_diagonal(
      matrix_free,
      diagonal,
      &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
        do_vmult_cell_single,
      &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
        local_apply_face_cell,
      &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
        local_apply_boundary_cell,
      this);

    // add coupling terms
    if (coupling_operator)
      coupling_operator->add_diagonal_entries(diagonal);

    for (auto &i : diagonal)
      i = (i != 0.0) ? (1.0 / i) : 1.0;
  }

  const TrilinosWrappers::SparseMatrix &
  get_system_matrix() const
  {
    initialize_system_matrix();

    return system_matrix;
  }

  void
  initialize_system_matrix() const
  {
    const auto &dof_handler = matrix_free.get_dof_handler();

    const auto constraints = (coupling_operator != nullptr) ?
                               (&coupling_operator->get_affine_constraints()) :
                               (&matrix_free.get_affine_constraints());

    if (system_matrix.m() == 0 || system_matrix.n() == 0)
      {
        system_matrix.clear();

        TrilinosWrappers::SparsityPattern dsp;

        dsp.reinit(dof_handler.locally_owned_dofs(),
                   dof_handler.get_mpi_communicator());

        DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, *constraints);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_sparsity_pattern_entries(dsp);

        dsp.compress();

        system_matrix.reinit(dsp);
      }

    if (this->valid_system == false)
      {
        system_matrix = 0.0;

        MatrixFreeTools::compute_matrix(
          matrix_free,
          *constraints,
          system_matrix,
          &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
            do_vmult_cell_single,
          &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
            local_apply_face_cell,
          &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
            local_apply_boundary_cell,
          this);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_system_matrix_entries(system_matrix);

        system_matrix.compress(VectorOperation::add);

        this->valid_system = true;
      }
  }


protected:
  void
  do_vmult_cell(const MatrixFree<dim, Number>               &data,
                VectorType                                  &dst,
                const VectorType                            &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FECellIntegrator phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);

        do_vmult_cell_single(phi);
        phi.distribute_local_to_global(dst);
      }
  }

  void
  local_apply_face(
    const MatrixFree<dim, Number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceIntegrator phi_m(data, true);
    FEFaceIntegrator phi_p(data, false);

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi_m.reinit(face);
        phi_p.reinit(face);
        phi_m.read_dof_values(src);
        phi_p.read_dof_values(src);
        local_apply_face_cell(phi_m, phi_p);
        phi_m.distribute_local_to_global(dst);
        phi_p.distribute_local_to_global(dst);
      }
  }

  void
  local_apply_boundary(
    const MatrixFree<dim, Number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceIntegrator phi(data, true);
    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi.reinit(face);
        phi.read_dof_values(src);
        local_apply_boundary_cell(phi);
        phi.distribute_local_to_global(dst);
      }
  }

  void
  do_vmult_cell_single(FECellIntegrator &phi) const
  {
    phi.evaluate(EvaluationFlags::gradients);

    for (unsigned int q = 0; q < phi.n_q_points; ++q)
      phi.submit_gradient(phi.get_gradient(q), q);

    phi.integrate(EvaluationFlags::gradients);
  }

  void
  local_apply_face_cell(FEFaceIntegrator &phi_m, FEFaceIntegrator &phi_p) const
  {
    phi_m.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    const auto sigma = std::max(phi_m.read_cell_data(penalty_parameters),
                                phi_p.read_cell_data(penalty_parameters)) *
                       panalty_factor;

    for (const auto q : phi_m.quadrature_point_indices())
      {
        const auto average_value =
          (phi_m.get_value(q) - phi_p.get_value(q)) * 0.5;
        const auto average_valgrad =
          average_value * 2. * sigma -
          (phi_m.get_normal_derivative(q) + phi_p.get_normal_derivative(q)) *
            0.5;

        phi_m.submit_normal_derivative(-average_value, q);
        phi_p.submit_normal_derivative(-average_value, q);
        phi_m.submit_value(average_valgrad, q);
        phi_p.submit_value(-average_valgrad, q);
      }

    phi_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  }

  void
  local_apply_boundary_cell(FEFaceIntegrator &phi) const
  {
    if ((phi.boundary_id() == 0 || phi.boundary_id() == 5))
      {
        const VectorizedArrayType zero = 0.0;

        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          phi.begin_dof_values()[i] = zero;
        return;
      }

    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    const auto sigma = phi.read_cell_data(penalty_parameters) * panalty_factor;

    for (const auto q : phi.quadrature_point_indices())
      {
        const auto average_value = phi.get_value(q);
        const auto average_valgrad =
          average_value * sigma * 2.0 - phi.get_normal_derivative(q);

        phi.submit_normal_derivative(-average_value, q);
        phi.submit_value(average_valgrad, q);
      }

    phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  }

  void
  compute_penalty_parameters()
  {
    const unsigned int n_cells =
      matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
    penalty_parameters.resize(n_cells);

    const auto &mapping = *matrix_free.get_mapping_info().mapping;
    const auto &fe      = matrix_free.get_dof_handler().get_fe();

    FEValues<dim> fe_values(mapping,
                            fe,
                            matrix_free.get_quadrature(),
                            update_JxW_values);

    FEFaceValues<dim> fe_face_values(mapping,
                                     fe,
                                     matrix_free.get_face_quadrature(),
                                     update_JxW_values);

    for (unsigned int cell = 0; cell < n_cells; ++cell)
      for (unsigned int v = 0;
           v < matrix_free.n_active_entries_per_cell_batch(cell);
           ++v)
        {
          const auto dealii_cell = matrix_free.get_cell_iterator(cell, v);
          fe_values.reinit(dealii_cell);

          // compute cell volume
          Number volume = 0.0;
          for (const auto q : fe_values.quadrature_point_indices())
            volume += fe_values.JxW(q);

          // compute surface area
          Number surface_area = 0.0;
          for (const auto f : dealii_cell->face_indices())
            {
              fe_face_values.reinit(dealii_cell, f);

              const Number factor = (dealii_cell->at_boundary(f) &&
                                     !dealii_cell->has_periodic_neighbor(f)) ?
                                      1. :
                                      0.5;

              for (const auto q : fe_face_values.quadrature_point_indices())
                surface_area += fe_face_values.JxW(q) * factor;
            }

          penalty_parameters[cell][v] = surface_area / volume;
        }
  }

  Number
  compute_pentaly_factor(const unsigned int degree, const Number factor) const
  {
    return factor * (degree + 1.0) * (degree + 1.0);
  }

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  mutable TrilinosWrappers::SparseMatrix       system_matrix;
  mutable bool                                 valid_system;

  AlignedVector<VectorizedArrayType> penalty_parameters;
  VectorizedArrayType                panalty_factor;

  std::shared_ptr<CouplingOperator<dim, Number>> coupling_operator;
};


/**
 * @brief Base class for the operator base
 */
template <int dim,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class GeneralStokesOperatorDG : public Subscriptor
{
public:
  using FECellIntegratorU =
    FEEvaluation<dim, -1, 0, dim, Number, VectorizedArrayType>;
  using FEFaceIntegratorU =
    FEFaceEvaluation<dim, -1, 0, dim, Number, VectorizedArrayType>;
  using FECellIntegratorP =
    FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;
  using FEFaceIntegratorP =
    FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;


  GeneralStokesOperatorDG(const MappingQ<dim>             &mapping,
                          const DoFHandler<dim>           &dof_handler,
                          const AffineConstraints<Number> &constraints,
                          const Quadrature<dim>           &quadrature,
                          const double                     sip_factor = 1.0,
                          const bool   weak_pressure_gradient_term    = true,
                          const bool   weak_velocity_divergence_term  = true,
                          const double delta_1_scaling                = 0.0)
    : weak_pressure_gradient_term(weak_pressure_gradient_term)
    , weak_velocity_divergence_term(weak_velocity_divergence_term)
    , delta_1_scaling(delta_1_scaling)
  {
    reinit(mapping, dof_handler, constraints, quadrature, sip_factor);
  }

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<dim>           &quadrature,
         const double                     sip_factor = 1.0)
  {
    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
    data.mapping_update_flags =
      update_quadrature_points | update_gradients | update_values;
    data.mapping_update_flags_inner_faces =
      update_values | update_gradients | update_JxW_values;
    data.mapping_update_flags_boundary_faces =
      update_values | update_gradients | update_JxW_values;

    matrix_free.reinit(mapping, dof_handler, constraints, quadrature, data);

    valid_system = false;

    this->sip_factor = sip_factor;

    compute_penalty_parameters();

    panalty_factor =
      compute_pentaly_factor(dof_handler.get_fe().degree, sip_factor);
  }

  /**
   * @brief Create coupling operator
   */
  void
  add_coupling(const std::shared_ptr<MortarManagerBase<dim>> mortar_manager,
               const unsigned int                            bid_0,
               const unsigned int                            bid_1)
  {
    const std::shared_ptr<CouplingEvaluationBase<dim, Number>>
      coupling_evaluator =
        std::make_shared<CouplingEvaluationStokes<dim, Number>>(
          *matrix_free.get_mapping_info().mapping,
          matrix_free.get_dof_handler(),
          weak_pressure_gradient_term,
          weak_velocity_divergence_term);

    coupling_operator = std::make_shared<CouplingOperator<dim, Number>>(
      *matrix_free.get_mapping_info().mapping,
      matrix_free.get_dof_handler(),
      matrix_free.get_affine_constraints(),
      coupling_evaluator,
      mortar_manager,
      bid_0,
      bid_1,
      sip_factor);

    coupling_bids.insert(bid_0);
    coupling_bids.insert(bid_1);

    compute_penalty_parameters();
  }

  virtual types::global_dof_index
  m() const
  {
    if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
      return this->matrix_free.get_dof_handler().n_dofs(
        this->matrix_free.get_mg_level());
    else
      return this->matrix_free.get_dof_handler().n_dofs();
  }

  Number
  el(unsigned int, unsigned int) const
  {
    DEAL_II_NOT_IMPLEMENTED();
    return 0;
  }

  void
  initialize_dof_vector(VectorType &dst) const
  {
    matrix_free.initialize_dof_vector(dst);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    src.update_ghost_values();

    matrix_free.loop(
      &GeneralStokesOperatorDG<dim, Number, VectorizedArrayType>::do_vmult_cell,
      &GeneralStokesOperatorDG<dim, Number, VectorizedArrayType>::
        local_apply_face,
      &GeneralStokesOperatorDG<dim, Number, VectorizedArrayType>::
        local_apply_boundary,
      this,
      dst,
      src,
      true);

    // apply coupling terms
    if (coupling_operator)
      {
        // apply constraints
        // TODO: only apply relevant constraints
        const auto &constraints = coupling_operator->get_affine_constraints();
        constraints.distribute(const_cast<VectorType &>(src));
        src.update_ghost_values();

        coupling_operator->vmult_add(dst, src);

        constraints.set_zero(const_cast<VectorType &>(src));
        constraints.set_zero(dst);
      }

    src.zero_out_ghost_values();
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    vmult(dst, src);
  }

  void
  compute_inverse_diagonal(VectorType &diagonal) const
  {
    (void)diagonal;

    AssertThrow(false, ExcNotImplemented());
  }

  const TrilinosWrappers::SparseMatrix &
  get_system_matrix() const
  {
    initialize_system_matrix();

    return system_matrix;
  }

  void
  initialize_system_matrix() const
  {
    const auto &dof_handler = matrix_free.get_dof_handler();

    auto constraints = &matrix_free.get_affine_constraints();

    AffineConstraints<Number> affine_constraints_tmp;

    if (coupling_operator)
      {
        affine_constraints_tmp.copy_from(
          coupling_operator->get_affine_constraints());
        affine_constraints_tmp.close();
      }

    if (system_matrix.m() == 0 || system_matrix.n() == 0)
      {
        system_matrix.clear();

        TrilinosWrappers::SparsityPattern dsp;

        dsp.reinit(dof_handler.locally_owned_dofs(),
                   dof_handler.get_mpi_communicator());

        DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, *constraints);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_sparsity_pattern_entries(dsp);

        dsp.compress();

        system_matrix.reinit(dsp);
      }

    if (this->valid_system == false)
      {
        system_matrix = 0.0;

        MatrixFreeTools::internal::
          ComputeMatrixScratchData<dim, VectorizedArray<Number>, false>
            data_cell;

        data_cell.dof_numbers               = {0, 0};
        data_cell.quad_numbers              = {0, 0};
        data_cell.n_components              = {dim, 1};
        data_cell.first_selected_components = {0, dim};
        data_cell.batch_type                = {0, 0};

        data_cell.op_create =
          [&](const std::pair<unsigned int, unsigned int> &range) {
            std::vector<std::unique_ptr<
              FEEvaluationData<dim, VectorizedArray<Number>, false>>>
              phi;

            phi.emplace_back(
              std::make_unique<FECellIntegratorU>(matrix_free, range, 0, 0, 0));

            phi.emplace_back(std::make_unique<FECellIntegratorP>(
              matrix_free, range, 0, 0, dim));

            return phi;
          };

        data_cell.op_reinit = [](auto &phi, const unsigned batch) {
          static_cast<FECellIntegratorU &>(*phi[0]).reinit(batch);
          static_cast<FECellIntegratorP &>(*phi[1]).reinit(batch);
        };

        data_cell.op_compute = [&](auto &phi) {
          auto &phi_0 = static_cast<FECellIntegratorU &>(*phi[0]);
          auto &phi_1 = static_cast<FECellIntegratorP &>(*phi[1]);

          do_vmult_cell_single(phi_0, phi_1);
        };

        MatrixFreeTools::internal::
          ComputeMatrixScratchData<dim, VectorizedArray<Number>, true>
            data_face;

        data_face.dof_numbers               = {0, 0, 0, 0};
        data_face.quad_numbers              = {0, 0, 0, 0};
        data_face.n_components              = {dim, dim, 1, 1};
        data_face.first_selected_components = {0, 0, dim, dim};
        data_face.batch_type                = {1, 2, 1, 2};

        data_face.op_create =
          [&](const std::pair<unsigned int, unsigned int> &range) {
            std::vector<std::unique_ptr<
              FEEvaluationData<dim, VectorizedArray<Number>, true>>>
              phi;

            phi.emplace_back(std::make_unique<FEFaceIntegratorU>(
              matrix_free, range, true, 0, 0, 0));

            phi.emplace_back(std::make_unique<FEFaceIntegratorU>(
              matrix_free, range, false, 0, 0, 0));

            phi.emplace_back(std::make_unique<FEFaceIntegratorP>(
              matrix_free, range, true, 0, 0, dim));

            phi.emplace_back(std::make_unique<FEFaceIntegratorP>(
              matrix_free, range, false, 0, 0, dim));

            return phi;
          };

        data_face.op_reinit = [](auto &phi, const unsigned batch) {
          static_cast<FEFaceIntegratorU &>(*phi[0]).reinit(batch);
          static_cast<FEFaceIntegratorU &>(*phi[1]).reinit(batch);
          static_cast<FEFaceIntegratorP &>(*phi[2]).reinit(batch);
          static_cast<FEFaceIntegratorP &>(*phi[3]).reinit(batch);
        };

        data_face.op_compute = [&](auto &phi) {
          auto &phi_0 = static_cast<FEFaceIntegratorU &>(*phi[0]);
          auto &phi_1 = static_cast<FEFaceIntegratorU &>(*phi[1]);
          auto &phi_2 = static_cast<FEFaceIntegratorP &>(*phi[2]);
          auto &phi_3 = static_cast<FEFaceIntegratorP &>(*phi[3]);

          local_apply_face_cell(phi_0, phi_1, phi_2, phi_3);
        };

        MatrixFreeTools::internal::
          ComputeMatrixScratchData<dim, VectorizedArray<Number>, true>
            data_boundary;

        data_boundary.dof_numbers               = {0, 0};
        data_boundary.quad_numbers              = {0, 0};
        data_boundary.n_components              = {dim, 1};
        data_boundary.first_selected_components = {0, dim};
        data_boundary.batch_type                = {1, 1};

        data_boundary.op_create =
          [&](const std::pair<unsigned int, unsigned int> &range) {
            std::vector<std::unique_ptr<
              FEEvaluationData<dim, VectorizedArray<Number>, true>>>
              phi;

            phi.emplace_back(std::make_unique<FEFaceIntegratorU>(
              matrix_free, range, true, 0, 0, 0));
            phi.emplace_back(std::make_unique<FEFaceIntegratorP>(
              matrix_free, range, true, 0, 0, dim));

            return phi;
          };

        data_boundary.op_reinit = [](auto &phi, const unsigned batch) {
          static_cast<FEFaceIntegratorU &>(*phi[0]).reinit(batch);
          static_cast<FEFaceIntegratorP &>(*phi[1]).reinit(batch);
        };

        data_boundary.op_compute = [&](auto &phi) {
          auto &phi_0 = static_cast<FEFaceIntegratorU &>(*phi[0]);
          auto &phi_1 = static_cast<FEFaceIntegratorP &>(*phi[1]);

          local_apply_boundary_cell(phi_0, phi_1);
        };

        MatrixFreeTools::internal::compute_matrix(matrix_free,
                                                  *constraints,
                                                  data_cell,
                                                  data_face,
                                                  data_boundary,
                                                  system_matrix);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_system_matrix_entries(system_matrix);

        system_matrix.compress(VectorOperation::add);

        this->valid_system = true;
      }
  }


private:
  void
  do_vmult_cell(const MatrixFree<dim, Number>               &matrix_free,
                VectorType                                  &dst,
                const VectorType                            &src,
                const std::pair<unsigned int, unsigned int> &range) const
  {
    FECellIntegratorU integrator_u(matrix_free, range, 0, 0, 0);
    FECellIntegratorP integrator_p(matrix_free, range, 0, 0, dim);
    for (unsigned cell = range.first; cell < range.second; ++cell)
      {
        integrator_u.reinit(cell);
        integrator_p.reinit(cell);

        integrator_u.read_dof_values(src);
        integrator_p.read_dof_values(src);
        do_vmult_cell_single(integrator_u, integrator_p);
        integrator_u.distribute_local_to_global(dst);
        integrator_p.distribute_local_to_global(dst);
      }
  }

  void
  local_apply_face(
    const MatrixFree<dim, Number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceIntegratorU phi_u_m(data, face_range, true, 0, 0, 0);
    FEFaceIntegratorU phi_u_p(data, face_range, false, 0, 0, 0);
    FEFaceIntegratorP phi_p_m(data, face_range, true, 0, 0, dim);
    FEFaceIntegratorP phi_p_p(data, face_range, false, 0, 0, dim);

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi_u_m.reinit(face);
        phi_u_p.reinit(face);
        phi_p_m.reinit(face);
        phi_p_p.reinit(face);
        phi_u_m.read_dof_values(src);
        phi_u_p.read_dof_values(src);
        phi_p_m.read_dof_values(src);
        phi_p_p.read_dof_values(src);
        local_apply_face_cell(phi_u_m, phi_u_p, phi_p_m, phi_p_p);
        phi_u_m.distribute_local_to_global(dst);
        phi_u_p.distribute_local_to_global(dst);
        phi_p_m.distribute_local_to_global(dst);
        phi_p_p.distribute_local_to_global(dst);
      }
  }

  void
  local_apply_boundary(
    const MatrixFree<dim, Number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceIntegratorU phi_u(data, face_range, true, 0, 0, 0);
    FEFaceIntegratorP phi_p(data, face_range, true, 0, 0, dim);
    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi_u.reinit(face);
        phi_p.reinit(face);
        phi_u.read_dof_values(src);
        phi_p.read_dof_values(src);
        local_apply_boundary_cell(phi_u, phi_p);
        phi_u.distribute_local_to_global(dst);
        phi_p.distribute_local_to_global(dst);
      }
  }

  void
  do_vmult_cell_single(FECellIntegratorU &phi_u, FECellIntegratorP &phi_p) const
  {
    phi_u.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    const auto cell = phi_u.get_current_cell_index();

    VectorizedArrayType delta_1;
    for (unsigned int v = 0;
         v < this->matrix_free.n_active_entries_per_cell_batch(cell);
         ++v)
      delta_1[v] =
        delta_1_scaling *
        this->matrix_free.get_cell_iterator(cell, v)->minimum_vertex_distance();

    for (unsigned int q = 0; q < phi_u.n_q_points; ++q)
      {
        typename FECellIntegratorP::value_type    p_value_result    = {};
        typename FECellIntegratorP::gradient_type p_gradient_result = {};
        typename FECellIntegratorU::value_type    u_value_result    = {};
        typename FECellIntegratorU::gradient_type u_gradient_result = {};

        const auto p_value    = phi_p.get_value(q);
        const auto p_gradient = phi_p.get_gradient(q);
        const auto u_value    = phi_u.get_value(q);
        const auto u_gradient = phi_u.get_gradient(q);

        if (true /*Laplace term*/)
          {
            // (∇v, ∇u)
            u_gradient_result += u_gradient;
          }

        if (weak_pressure_gradient_term)
          {
            // - (div(v), p)
            for (unsigned int d = 0; d < dim; ++d)
              u_gradient_result[d][d] -= p_value;
          }
        else
          {
            // + (v, ∇p)
            u_value_result += p_gradient;
          }

        if (weak_velocity_divergence_term)
          {
            // - (∇q, u)
            p_gradient_result -= u_value * vel_div_sign;
          }
        else
          {
            // + (q, div(u))
            for (unsigned int d = 0; d < dim; ++d)
              p_value_result += u_gradient[d][d] * vel_div_sign;
          }

        // δ_1 (∇q, ∇p)
        if (delta_1_scaling != 0.0)
          p_gradient_result += delta_1 * p_gradient;

        phi_p.submit_value(p_value_result, q);
        phi_p.submit_gradient(p_gradient_result, q);

        phi_u.submit_value(u_value_result, q);
        phi_u.submit_gradient(u_gradient_result, q);
      }

    phi_u.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  }

  void
  local_apply_face_cell(FEFaceIntegratorU &phi_u_m,
                        FEFaceIntegratorU &phi_u_p,
                        FEFaceIntegratorP &phi_p_m,
                        FEFaceIntegratorP &phi_p_p) const
  {
    phi_u_m.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_u_p.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p_m.evaluate(EvaluationFlags::values);
    phi_p_p.evaluate(EvaluationFlags::values);

    const auto sigma = std::max(phi_u_m.read_cell_data(penalty_parameters),
                                phi_u_p.read_cell_data(penalty_parameters)) *
                       panalty_factor;

    VectorizedArrayType mask = 1.0;

    const unsigned int face = phi_u_m.get_current_cell_index();
    for (unsigned int v = 0;
         v < matrix_free.n_active_entries_per_face_batch(face);
         ++v)
      if (matrix_free.get_face_iterator(face, v, true).first->material_id() !=
          matrix_free.get_face_iterator(face, v, false).first->material_id())
        mask[v] = 0.0;

    for (const auto q : phi_u_m.quadrature_point_indices())
      {
        const auto u_value_avg =
          (phi_u_m.get_value(q) + phi_u_p.get_value(q)) * 0.5;
        const auto u_value_jump = phi_u_m.get_value(q) - phi_u_p.get_value(q);
        const auto u_gradient_avg =
          (phi_u_m.get_gradient(q) + phi_u_p.get_gradient(q)) * 0.5;
        const auto p_value_avg =
          (phi_p_m.get_value(q) + phi_p_p.get_value(q)) * 0.5;
        const auto normal = phi_u_m.normal_vector(q);

        typename FECellIntegratorU::value_type u_normal_gradient_avg_result =
          {};
        typename FECellIntegratorU::value_type u_value_jump_result = {};
        typename FECellIntegratorP::value_type p_value_jump_result = {};

        if (true /*Laplace term*/)
          {
            // - (n avg(∇v), jump(u))
            u_normal_gradient_avg_result -= u_value_jump;

            // - (jump(v), avg(∇u) n)
            u_value_jump_result -= u_gradient_avg * normal;

            // + (jump(v), σ jump(u))
            u_value_jump_result += sigma * u_value_jump;
          }

        if (weak_pressure_gradient_term)
          {
            // + (jump(v), avg(p) n)
            u_value_jump_result += p_value_avg * normal;
          }
        else
          {
            // nothing to do
          }

        if (weak_velocity_divergence_term)
          {
            // + (jump(q), avg(u) n)
            p_value_jump_result += u_value_avg * normal * vel_div_sign;
          }
        else
          {
            // - (avg(q), jump(u) n)
            p_value_jump_result -= 0.5 * u_value_jump * normal;
          }

        phi_u_m.submit_normal_derivative(u_normal_gradient_avg_result * 0.5, q);
        phi_u_p.submit_normal_derivative(u_normal_gradient_avg_result * 0.5, q);
        phi_u_m.submit_value(u_value_jump_result, q);
        phi_u_p.submit_value(-u_value_jump_result, q);
        phi_p_m.submit_value(p_value_jump_result, q);
        phi_p_p.submit_value(-p_value_jump_result, q);
      }

    phi_u_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_u_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p_m.integrate(EvaluationFlags::values);
    phi_p_p.integrate(EvaluationFlags::values);
  }

  void
  local_apply_boundary_cell(FEFaceIntegratorU &phi_u_m,
                            FEFaceIntegratorP &phi_p_m) const
  {
    if (coupling_bids.find(phi_u_m.boundary_id()) != coupling_bids.end())
      {
        for (unsigned int i = 0; i < phi_u_m.dofs_per_cell; ++i)
          phi_u_m.begin_dof_values()[i] = 0.0;
        for (unsigned int i = 0; i < phi_p_m.dofs_per_cell; ++i)
          phi_p_m.begin_dof_values()[i] = 0.0;
        return;
      }

    phi_u_m.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p_m.evaluate(EvaluationFlags::values);

    const auto sigma =
      phi_u_m.read_cell_data(penalty_parameters) * panalty_factor;

    for (const auto q : phi_u_m.quadrature_point_indices())
      {
        const auto u_value_avg =
          (phi_u_m.get_value(q) - phi_u_m.get_value(q)) * 0.5;
        const auto u_value_jump = phi_u_m.get_value(q) + phi_u_m.get_value(q);
        const auto u_gradient_avg =
          (phi_u_m.get_gradient(q) + phi_u_m.get_gradient(q)) * 0.5;
        const auto p_value_avg =
          (phi_p_m.get_value(q) + phi_p_m.get_value(q)) * 0.5;
        const auto normal = phi_u_m.normal_vector(q);

        typename FECellIntegratorU::value_type u_normal_gradient_avg_result =
          {};
        typename FECellIntegratorU::value_type u_value_jump_result = {};
        typename FECellIntegratorP::value_type p_value_jump_result = {};

        if (true /*Laplace term*/)
          {
            // - (n avg(∇v), jump(u))
            u_normal_gradient_avg_result -= u_value_jump;

            // - (jump(v), avg(∇u) n)
            u_value_jump_result -= u_gradient_avg * normal;

            // + (jump(v), σ jump(u))
            u_value_jump_result += sigma * u_value_jump;
          }

        if (weak_pressure_gradient_term)
          {
            // + (jump(v), avg(p) n)
            u_value_jump_result += p_value_avg * normal;
          }
        else
          {
            // nothing to do
          }

        if (weak_velocity_divergence_term)
          {
            // + (jump(q), avg(u) n)
            p_value_jump_result += u_value_avg * normal * vel_div_sign;
          }
        else
          {
            // nothing to do
          }

        phi_u_m.submit_normal_derivative(u_normal_gradient_avg_result * 0.5, q);
        phi_u_m.submit_value(u_value_jump_result, q);
        phi_p_m.submit_value(p_value_jump_result, q);
      }

    phi_u_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p_m.integrate(EvaluationFlags::values);
  }

  void
  compute_penalty_parameters()
  {
    const unsigned int n_cells =
      matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
    penalty_parameters.resize(n_cells);

    const auto &mapping = *matrix_free.get_mapping_info().mapping;
    const auto &fe      = matrix_free.get_dof_handler().get_fe();

    FEValues<dim> fe_values(mapping,
                            fe,
                            matrix_free.get_quadrature(),
                            update_JxW_values);

    FEFaceValues<dim> fe_face_values(mapping,
                                     fe,
                                     matrix_free.get_face_quadrature(),
                                     update_JxW_values);

    for (unsigned int cell = 0; cell < n_cells; ++cell)
      for (unsigned int v = 0;
           v < matrix_free.n_active_entries_per_cell_batch(cell);
           ++v)
        {
          const auto dealii_cell = matrix_free.get_cell_iterator(cell, v);
          fe_values.reinit(dealii_cell);

          // compute cell volume
          Number volume = 0.0;
          for (const auto q : fe_values.quadrature_point_indices())
            volume += fe_values.JxW(q);

          // compute surface area
          Number surface_area = 0.0;
          for (const auto f : dealii_cell->face_indices())
            {
              fe_face_values.reinit(dealii_cell, f);

              const Number factor =
                (dealii_cell->at_boundary(f) &&
                 !dealii_cell->has_periodic_neighbor(f) &&
                 (coupling_bids.find(dealii_cell->face(f)->boundary_id()) ==
                  coupling_bids.end())) ?
                  1. :
                  0.5;

              for (const auto q : fe_face_values.quadrature_point_indices())
                surface_area += fe_face_values.JxW(q) * factor;
            }

          penalty_parameters[cell][v] = surface_area / volume;
        }
  }

  Number
  compute_pentaly_factor(const unsigned int degree, const Number factor) const
  {
    return factor * (degree + 1.0) * (degree + 1.0);
  }

  const bool   weak_pressure_gradient_term;
  const bool   weak_velocity_divergence_term;
  const double delta_1_scaling;

  const double vel_div_sign = +1.0;

  mutable double sip_factor;

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  mutable TrilinosWrappers::SparseMatrix       system_matrix;
  mutable bool                                 valid_system;

  AlignedVector<VectorizedArrayType> penalty_parameters;
  VectorizedArrayType                panalty_factor;

  std::shared_ptr<CouplingOperator<dim, Number>> coupling_operator;

  std::set<unsigned int> coupling_bids;
};


template <int dim>
class MyMortarManagerCircle : public MortarManagerBase<dim>
{
public:
  template <int dim2>
  MyMortarManagerCircle(const std::vector<unsigned int> &n_subdivisions,
                        const std::vector<double>       &radius,
                        const Quadrature<dim2>          &quadrature,
                        const double                     rotation_angle);

protected:
  Point<dim>
  from_1D(const double rad) const override;

  double
  to_1D(const Point<dim> &point) const override;

  Tensor<1, dim, double>
  get_normal(const Point<dim> &point) const override;
};


template <int dim>
template <int dim2>
MyMortarManagerCircle<dim>::MyMortarManagerCircle(
  const std::vector<unsigned int> &n_subdivisions,
  const std::vector<double>       &radius,
  const Quadrature<dim2>          &quadrature,
  const double                     rotation_angle)
  : MortarManagerBase<dim>(n_subdivisions, radius, quadrature, rotation_angle)
{}

template <int dim>
Point<dim>
MyMortarManagerCircle<dim>::from_1D(const double radiant) const
{
  return radius_to_point<dim>(this->radius[0], radiant);
}

template <int dim>
double
MyMortarManagerCircle<dim>::to_1D(const Point<dim> &point) const
{
  return point_to_angle(point, Point<dim>());
}

template <int dim>
Tensor<1, dim, double>
MyMortarManagerCircle<dim>::get_normal(const Point<dim> &point_in) const
{
  Point<dim> point = point_in;

  if (dim == 3)
    point[2] = 0.0;

  return point / point.norm();
}
