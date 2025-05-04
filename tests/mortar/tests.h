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
        construct_quadrature(matrix_free.get_quadrature()),
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

private:
  static Quadrature<dim>
  construct_quadrature(const Quadrature<dim> &quad)
  {
    const double oversamling_factor = 2.0; // make parameter

    for (unsigned int i = 1; i <= 10; ++i)
      if (quad == QGauss<dim>(i))
        return QGauss<dim>(i * oversamling_factor);

    AssertThrow(false, ExcNotImplemented());

    return quad;
  }
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

    const bool is_p_disc = matrix_free.get_dof_handler()
                             .get_fe()
                             .base_element(matrix_free.get_dof_handler()
                                             .get_fe()
                                             .component_to_base_index(dim)
                                             .first)
                             .n_dofs_per_vertex() == 0;

    if (is_p_disc == false)
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

        const auto u_gradient = phi_u.get_gradient(q);

        // a)     (ε(v), 2νε(u))
        if (true)
          {
            symm_scalar_product_add(u_gradient_result,
                                    u_gradient,
                                    VectorizedArrayType(2.0));
          }
        else
          {
            u_gradient_result = u_gradient;
          }

        // b)   - (div(v), p)
        for (unsigned int d = 0; d < dim; ++d)
          u_gradient_result[d][d] -= p_value;

        // c)     (q, div(u))
        for (unsigned int d = 0; d < dim; ++d)
          p_value_result -= u_gradient[d][d];

        // d) δ_1 (∇q, ∇p)
        if (delta_1_scaling != 0.0)
          p_gradient_result = -delta_1 * p_gradient;

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

  std::shared_ptr<CouplingOperator<dim, dim, Number>> coupling_operator_v;
  std::shared_ptr<CouplingOperator<dim, 1, Number>>   coupling_operator_p;

  const double delta_1_scaling;
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

    comute_penalty_parameters();

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
                   dof_handler.get_communicator());

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
  comute_penalty_parameters()
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

  std::shared_ptr<CouplingOperator<dim, n_components, Number>>
    coupling_operator;
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
    data.mapping_update_flags_inner_faces =
      update_values | update_gradients | update_JxW_values;
    data.mapping_update_flags_boundary_faces =
      update_values | update_gradients | update_JxW_values;

    matrix_free.reinit(mapping, dof_handler, constraints, quadrature, data);

    valid_system = false;

    comute_penalty_parameters();

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

    coupling_bids.insert(bid_0);
    coupling_bids.insert(bid_1);
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
    if (coupling_operator_v)
      coupling_operator_v->vmult_add(dst, src);

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

    if (coupling_operator_v)
      {
        affine_constraints_tmp.copy_from(
          coupling_operator_v->get_affine_constraints());

        affine_constraints_tmp.close();
      }

    if (system_matrix.m() == 0 || system_matrix.n() == 0)
      {
        system_matrix.clear();

        TrilinosWrappers::SparsityPattern dsp;

        dsp.reinit(dof_handler.locally_owned_dofs(),
                   dof_handler.get_communicator());

        DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, *constraints);

        // apply coupling terms
        if (coupling_operator_v)
          coupling_operator_v->add_sparsity_pattern_entries(dsp);

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
        if (coupling_operator_v)
          coupling_operator_v->add_system_matrix_entries(system_matrix);

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

    for (unsigned int q = 0; q < phi_u.n_q_points; ++q)
      {
        typename FECellIntegratorP::value_type    p_value_result    = {};
        typename FECellIntegratorP::gradient_type p_gradient_result = {};
        typename FECellIntegratorU::value_type    u_value_result    = {};
        typename FECellIntegratorU::gradient_type u_gradient_result = {};

        const auto p_value    = phi_p.get_value(q);
        const auto u_value    = phi_u.get_value(q);
        const auto u_gradient = phi_u.get_gradient(q);

        // a)     (∇v, ∇u)
        u_gradient_result = u_gradient;

        // b)   - (div(v), p)
        for (unsigned int d = 0; d < dim; ++d)
          u_gradient_result[d][d] -= p_value;

        // c)   - (∇q, u)
        p_gradient_result = -u_value;

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

    for (const auto q : phi_u_m.quadrature_point_indices())
      {
        const auto average_value =
          (phi_u_m.get_value(q) - phi_u_p.get_value(q)) * 0.5;
        const auto average_normal_value_u =
          average_value * phi_u_m.normal_vector(q);
        const auto average_value_p =
          (phi_p_m.get_value(q) - phi_p_p.get_value(q)) * 0.5;
        const auto avg_normal_gradient = (phi_u_m.get_normal_derivative(q) +
                                          phi_u_p.get_normal_derivative(q)) *
                                         0.5;

        const auto average_valgrad = average_value * 2. * sigma -
                                     avg_normal_gradient +
                                     average_value_p * phi_u_m.normal_vector(q);

        phi_u_m.submit_normal_derivative(-average_value, q);
        phi_u_p.submit_normal_derivative(-average_value, q);
        phi_u_m.submit_value(average_valgrad, q);
        phi_u_p.submit_value(-average_valgrad, q);
        phi_p_m.submit_value(average_normal_value_u, q);
        phi_p_p.submit_value(-average_normal_value_u, q);
      }

    phi_u_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_u_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p_m.integrate(EvaluationFlags::values);
    phi_p_p.integrate(EvaluationFlags::values);
  }

  void
  local_apply_boundary_cell(FEFaceIntegratorU &phi_u,
                            FEFaceIntegratorP &phi_p) const
  {
    if (coupling_bids.find(phi_u.boundary_id()) != coupling_bids.end())
      {
        for (unsigned int i = 0; i < phi_u.dofs_per_cell; ++i)
          phi_u.begin_dof_values()[i] = 0.0;
        for (unsigned int i = 0; i < phi_p.dofs_per_cell; ++i)
          phi_p.begin_dof_values()[i] = 0.0;
        return;
      }

    phi_u.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.evaluate(EvaluationFlags::values);

    const auto sigma =
      phi_u.read_cell_data(penalty_parameters) * panalty_factor;

    for (const auto q : phi_u.quadrature_point_indices())
      {
        const auto average_value = phi_u.get_value(q) * 0.0;
        const auto average_normal_value_u =
          average_value * phi_u.normal_vector(q);
        const auto average_value_p = phi_p.get_value(q);

        const auto average_valgrad = average_value * sigma * 2.0 -
                                     phi_u.get_normal_derivative(q) +
                                     average_value_p * phi_u.normal_vector(q);

        phi_u.submit_normal_derivative(-average_value, q);
        phi_u.submit_value(average_valgrad, q);
        phi_p.submit_value(average_normal_value_u, q);
      }

    phi_u.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.integrate(EvaluationFlags::values);
  }

  void
  comute_penalty_parameters()
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

  std::shared_ptr<CouplingOperator<dim, dim, Number>> coupling_operator_v;

  std::set<unsigned int> coupling_bids;
};
