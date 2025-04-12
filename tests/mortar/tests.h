// SPDX-FileCopyrightText: Copyright (c) 2021 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

// Lethe
#include <core/mortar_coupling_manager.h>

using namespace dealii;

/**
 * @brief Create grid with hyper_shell and hyper_shell geometries
 * Only the merged triangulation (and resulting grid) are stored
 */
template <int dim>
void
hyper_shell_with_hyper_shell(const double radius, Triangulation<dim> &tria)
{
  const double r1_i = radius * 0.25;
  const double r1_o = radius * 0.5;
  const double r2_i = radius * 0.5;
  const double r2_o = radius * 1.0;

  // inner domain triangulation
  Triangulation<dim> circle_one;
  GridGenerator::hyper_shell(circle_one,
                             (dim == 2) ? Point<dim>(0, 0) :
                                          Point<dim>(0, 0, 0),
                             r1_i,
                             r1_o,
                             6,
                             true);
  // outer domain triangulation
  Triangulation<dim> circle_two;
  GridGenerator::hyper_shell(circle_two,
                             (dim == 2) ? Point<dim>(0, 0) :
                                          Point<dim>(0, 0, 0),
                             r2_i,
                             r2_o,
                             6,
                             true);

  // shift boundary id of circle two
  for (const auto &face : circle_two.active_face_iterators())
    if (face->at_boundary())
      face->set_boundary_id(face->boundary_id() + 2);

  // create unique triangulation
  GridGenerator::merge_triangulations(
    circle_one, circle_two, tria, 0, true, true);
  // store manifolds in merged triangulation
  tria.set_manifold(0,
                    SphericalManifold<dim>((dim == 2) ? Point<dim>(0, 0) :
                                                        Point<dim>(0, 0, 0)));
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

  // create unique triangulation
  GridGenerator::merge_triangulations(tria_0, tria_1, tria, 0, true, true);
  // store manifolds in merged triangulation
  tria.set_manifold(0, tria_0.get_manifold(0));
  tria.set_manifold(1, tria_0.get_manifold(1));
  tria.set_manifold(2, tria_1.get_manifold(0));
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
  tria.set_manifold(0, tria_0.get_manifold(0));
  tria.set_manifold(1, tria_0.get_manifold(1));
  tria.set_manifold(2, tria_1.get_manifold(0));
}

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
  add_coupling(const unsigned int n_subdivisions,
               const double       radius,
               const double       rotate_pi,
               const unsigned int bid_0,
               const unsigned int bid_1,
               const double       sip_factor = 1.0)
  {
    coupling_operator =
      std::make_shared<CouplingOperator<dim, n_components, Number>>(
        *matrix_free.get_mapping_info().mapping,
        matrix_free.get_dof_handler(),
        matrix_free.get_affine_constraints(),
        matrix_free.get_quadrature(),
        n_subdivisions,
        radius,
        rotate_pi,
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
                   dof_handler.get_communicator());

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

  std::shared_ptr<CouplingOperator<dim, n_components, Number>>
    coupling_operator;
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
        0.01 *
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
          value_result[dim] += u_gradient[d][d];

        // d) δ_1 (∇q, ∇p)
        gradient_result[dim] = delta_1 * p_gradient;

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
  using FECellIntegrator =
    FEEvaluation<dim, -1, 0, dim + 1, Number, VectorizedArrayType>;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;


  GeneralStokesOperator(const MappingQ<dim>             &mapping,
                        const DoFHandler<dim>           &dof_handler,
                        const AffineConstraints<Number> &constraints,
                        const Quadrature<dim>           &quadrature,
                        const double                     delta_1_scaling)
    : delta_1_scaling(delta_1_scaling)
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
  add_coupling(const unsigned int n_subdivisions,
               const double       radius,
               const double       rotate_pi,
               const unsigned int bid_0,
               const unsigned int bid_1,
               const double       sip_factor = 1.0)
  {
    coupling_operator_v = std::make_shared<CouplingOperator<dim, dim, Number>>(
      *matrix_free.get_mapping_info().mapping,
      matrix_free.get_dof_handler(),
      matrix_free.get_affine_constraints(),
      matrix_free.get_quadrature(),
      n_subdivisions,
      radius,
      rotate_pi,
      bid_0,
      bid_1,
      sip_factor,
      0);
    coupling_operator_p = std::make_shared<CouplingOperator<dim, 1, Number>>(
      *matrix_free.get_mapping_info().mapping,
      matrix_free.get_dof_handler(),
      matrix_free.get_affine_constraints(),
      matrix_free.get_quadrature(),
      n_subdivisions,
      radius,
      rotate_pi,
      bid_0,
      bid_1,
      sip_factor,
      dim);
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
    if (coupling_operator_v)
      coupling_operator_v->vmult_add(dst, src);

    if (coupling_operator_p)
      coupling_operator_p->vmult_add(dst, src);

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
      &GeneralStokesOperator<dim, Number, VectorizedArrayType>::
        do_vmult_cell_single,
      this);

    // add coupling terms
    if (coupling_operator_v)
      coupling_operator_v->add_diagonal_entries(diagonal);

    if (coupling_operator_p)
      coupling_operator_p->add_diagonal_entries(diagonal);

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

    if (coupling_operator_v)
      {
        affine_constraints_tmp.copy_from(
          coupling_operator_v->get_affine_constraints());

        if (coupling_operator_p)
          affine_constraints_tmp.merge(
            coupling_operator_p->get_affine_constraints(),
            AffineConstraints<Number>::MergeConflictBehavior::left_object_wins,
            true);

        affine_constraints_tmp.close();
      }

    if (system_matrix.m() == 0 || system_matrix.n() == 0)
      {
        system_matrix.clear();

        TrilinosWrappers::SparsityPattern dsp;

        dsp.reinit(dof_handler.locally_owned_dofs(),
                   dof_handler.get_communicator());

        DoFTools::make_sparsity_pattern(dof_handler, dsp, *constraints);

        // apply coupling terms
        if (coupling_operator_v)
          coupling_operator_v->add_sparsity_pattern_entries(dsp);

        if (coupling_operator_p)
          coupling_operator_p->add_sparsity_pattern_entries(dsp);

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
          &GeneralStokesOperator<dim, Number, VectorizedArrayType>::
            do_vmult_cell_single,
          this);

        // apply coupling terms
        if (coupling_operator_v)
          coupling_operator_v->add_system_matrix_entries(system_matrix);

        if (coupling_operator_p)
          coupling_operator_p->add_system_matrix_entries(system_matrix);

        system_matrix.compress(VectorOperation::add);

        this->valid_system = true;
      }
  }


private:
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
  do_vmult_cell_single(FECellIntegrator &phi) const
  {
    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    const auto cell = phi.get_current_cell_index();

    VectorizedArrayType delta_1;
    for (unsigned int v = 0;
         v < this->matrix_free.n_active_entries_per_cell_batch(cell);
         ++v)
      delta_1[v] =
        0.01 *
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
          value_result[dim] += u_gradient[d][d];

        // d) δ_1 (∇q, ∇p)
        gradient_result[dim] = delta_1 * p_gradient;

        phi.submit_value(value_result, q);
        phi.submit_gradient(gradient_result, q);
      }

    phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  }

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  mutable TrilinosWrappers::SparseMatrix       system_matrix;
  mutable bool                                 valid_system;

  std::shared_ptr<CouplingOperator<dim, dim, Number>> coupling_operator_v;
  std::shared_ptr<CouplingOperator<dim, 1, Number>>   coupling_operator_p;

  const double delta_1_scaling;
};
