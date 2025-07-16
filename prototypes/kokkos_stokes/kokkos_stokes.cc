// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Solve stabilized Stokes with equal-order elements:
// (ε(v), 2νε(u)) - (div(v), p) + (q, div(u)) - (v, f) + δ_1 (∇q, ∇p - f) = 0
//                                                       +-----stab-----+

#include <deal.II/base/convergence_table.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_base.templates.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/multigrid.templates.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

static unsigned int counter = 0;



template <int dim, typename Number, typename MemorySpace>
class MGTransferMFWrapper
  : public MGTransferBase<
      LinearAlgebra::distributed::Vector<Number, MemorySpace>>
{
public:
  using VectorType = LinearAlgebra::distributed::Vector<Number, MemorySpace>;
  using VectorTypeHost =
    LinearAlgebra::distributed::Vector<Number, dealii::MemorySpace::Host>;

  MGTransferMFWrapper(
    const MGLevelObject<MGTwoLevelTransfer<dim, VectorTypeHost>> &mg_transfers,
    const std::function<void(const unsigned int, VectorTypeHost &)>
      &initialize_dof_vector)
    : transfer(mg_transfers, initialize_dof_vector)
  {}

  template <typename Number2>
  void
  copy_to_mg(
    const DoFHandler<dim>                                          &dof_handler,
    MGLevelObject<VectorType>                                      &dst,
    const LinearAlgebra::distributed::Vector<Number2, MemorySpace> &src) const
  {
    MGLevelObject<VectorTypeHost> dst_host(dst.min_level(), dst.max_level());
    LinearAlgebra::distributed::Vector<Number2, dealii::MemorySpace::Host>
      src_host;

    copy_to_host(src_host, src);
    for (unsigned int l = dst.min_level(); l < dst.max_level(); ++l)
      copy_to_host(dst_host[l], dst[l]);

    transfer.copy_to_mg(dof_handler, dst_host, src_host);

    for (unsigned int l = dst.min_level(); l <= dst.max_level(); ++l)
      copy_from_host(dst[l], dst_host[l]);
  }

  template <typename Number2>
  void
  copy_from_mg(const DoFHandler<dim> &dof_handler,
               LinearAlgebra::distributed::Vector<Number2, MemorySpace> &dst,
               const MGLevelObject<VectorType> &src) const
  {
    LinearAlgebra::distributed::Vector<Number2, dealii::MemorySpace::Host>
                                  dst_host;
    MGLevelObject<VectorTypeHost> src_host(src.min_level(), src.max_level());

    copy_to_host(dst_host, dst);
    for (unsigned int l = src.min_level(); l <= src.max_level(); ++l)
      copy_to_host(src_host[l], src[l]);

    transfer.copy_from_mg(dof_handler, dst_host, src_host);

    copy_from_host(dst, dst_host);
  }

  void
  prolongate(const unsigned int to_level,
             VectorType        &dst,
             const VectorType  &src) const override
  {
    VectorTypeHost dst_host;
    VectorTypeHost src_host;

    copy_to_host(dst_host, dst);
    copy_to_host(src_host, src);

    transfer.prolongate(to_level, dst_host, src_host);

    copy_from_host(dst, dst_host);
  }

  void
  restrict_and_add(const unsigned int from_level,
                   VectorType        &dst,
                   const VectorType  &src) const override
  {
    VectorTypeHost dst_host;
    VectorTypeHost src_host;

    copy_to_host(dst_host, dst);
    copy_to_host(src_host, src);

    transfer.restrict_and_add(from_level, dst_host, src_host);

    copy_from_host(dst, dst_host);
  }

private:
  const MGTransferMF<dim, Number, dealii::MemorySpace::Host> transfer;

  template <typename Number2>
  void
  copy_to_host(
    LinearAlgebra::distributed::Vector<Number2, dealii::MemorySpace::Host> &dst,
    const LinearAlgebra::distributed::Vector<Number2, MemorySpace> &src) const
  {
    LinearAlgebra::ReadWriteVector<Number2> rw_vector(
      src.get_partitioner()->locally_owned_range());
    rw_vector.import_elements(src, VectorOperation::insert);

    dst.reinit(src.get_partitioner());
    dst.import_elements(rw_vector, VectorOperation::insert);
  }

  template <typename Number2>
  void
  copy_from_host(
    LinearAlgebra::distributed::Vector<Number2, MemorySpace> &dst,
    const LinearAlgebra::distributed::Vector<Number2, dealii::MemorySpace::Host>
      &src) const
  {
    LinearAlgebra::ReadWriteVector<Number2> rw_vector(
      src.get_partitioner()->locally_owned_range());
    rw_vector.import_elements(src, VectorOperation::insert);

    if (dst.size() == 0)
      dst.reinit(src.get_partitioner());
    dst.import_elements(rw_vector, VectorOperation::insert);
  }
};

template <int dim, int fe_degree, typename Number, typename MemorySpace>
class StokesOperator;

template <int dim, int fe_degree, typename Number>
class StokesOperator<dim, fe_degree, Number, MemorySpace::Host>
  : public EnableObserverPointer
{
  using FECellIntegrator =
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, Number>;

public:
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>;

  StokesOperator() = default;

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<1>             &quadrature,
         const LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
           &delta_1_dealii)
  {
    typename MatrixFree<dim, Number>::AdditionalData additional_data;
    additional_data.mapping_update_flags =
      update_JxW_values | update_values | update_gradients;

    matrix_free.reinit(
      mapping, dof_handler, constraints, quadrature, additional_data);

    delta_1.resize(matrix_free.n_cell_batches());

    for (unsigned int cell = 0; cell < matrix_free.n_cell_batches(); ++cell)
      for (unsigned int v = 0;
           v < matrix_free.n_active_entries_per_cell_batch(cell);
           ++v)
        delta_1[cell][v] =
          delta_1_dealii[matrix_free.get_cell_iterator(cell, v)
                           ->active_cell_index()];
  }

  void
  initialize_dof_vector(VectorType &vec) const
  {
    matrix_free.initialize_dof_vector(vec);
  }

  types::global_dof_index
  m() const
  {
    return matrix_free.get_dof_handler().n_dofs();
  }

  Number
  el(unsigned int, unsigned int) const
  {
    DEAL_II_NOT_IMPLEMENTED();
    return 0;
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    matrix_free.cell_loop(&StokesOperator::local_apply, this, dst, src, true);
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    AssertThrow(false, ExcNotImplemented());

    (void)dst;
    (void)src;
  }

  void
  compute_inverse_diagonal(VectorType &diagonal_global) const
  {
    this->initialize_dof_vector(diagonal_global);

    MatrixFreeTools::compute_diagonal(matrix_free,
                                      diagonal_global,
                                      &StokesOperator::do_vmult_cell_single,
                                      this);

    for (auto &e : diagonal_global)
      e = (e == 0.0) ? 1.0 : (1.0 / e);
  }

private:
  void
  local_apply(const MatrixFree<dim, Number>               &data,
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
    const auto delta_1 = this->delta_1[phi.get_current_cell_index()];

    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
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

  MatrixFree<dim, Number> matrix_free;

  std::vector<VectorizedArray<Number>> delta_1;
};



template <int dim, int fe_degree, typename Number>
class StokesOperatorQuad
{
  using FECellIntegrator =
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, Number>;

public:
  DEAL_II_HOST_DEVICE
  StokesOperatorQuad(const Number *delta_1)
    : delta_1(delta_1)
  {}

  DEAL_II_HOST_DEVICE
  StokesOperatorQuad(
    const typename Portable::MatrixFree<dim, Number>::Data *data,
    const Number                                           *delta_1)
    : data(data)
    , delta_1(delta_1)
  {}

  DEAL_II_HOST_DEVICE void
  operator()(FECellIntegrator *phi, const int q_point) const
  {
    const int  cell_index = phi->get_current_cell_index();
    const auto delta_1    = this->delta_1[cell_index];

    typename FECellIntegrator::value_type    value_result    = {};
    typename FECellIntegrator::gradient_type gradient_result = {};

    const auto value    = phi->get_value(q_point);
    const auto gradient = phi->get_gradient(q_point);

    const Number                 p_value    = value[dim];
    const Tensor<1, dim, Number> p_gradient = gradient[dim];

    Tensor<2, dim, Number> u_gradient;

    for (unsigned int d = 0; d < dim; ++d)
      u_gradient[d] = gradient[d];

    // a)     (ε(v), 2νε(u))
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

    phi->submit_value(value_result, q_point);
    phi->submit_gradient(gradient_result, q_point);
  }

  DEAL_II_HOST_DEVICE void
  set_matrix_free_data(
    const typename Portable::MatrixFree<dim, Number>::Data &data)
  {
    this->data = &data;
  }

private:
  const typename Portable::MatrixFree<dim, Number>::Data *data;
  const Number                                           *delta_1;
  static const unsigned int                               n_q_points =
    dealii::Utilities::pow(fe_degree + 1, dim);
};

template <int dim, int fe_degree, typename Number>
class StokesOperatorLocal
{
public:
  StokesOperatorLocal(const Number *delta_1)
    : delta_1(delta_1)
  {}

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceVector<Number>                   &src,
             Portable::DeviceVector<Number>                         &dst) const
  {
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, Number> phi(
      data);
    phi.read_dof_values(src);
    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi.apply_for_each_quad_point(
      StokesOperatorQuad<dim, fe_degree, Number>(data, delta_1));
    phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi.distribute_local_to_global(dst);
  }

  static const unsigned int n_local_dofs =
    Utilities::pow(fe_degree + 1, dim) * (dim + 1);
  static const unsigned int n_q_points = Utilities::pow(fe_degree + 1, dim);

private:
  const Number *delta_1;
};

template <int dim, int fe_degree, typename Number>
class StokesOperator<dim, fe_degree, Number, MemorySpace::Default>
  : public EnableObserverPointer
{
public:
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;

  StokesOperator() = default;

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<1>             &quadrature,
         const LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
           &delta_1_host_dealii)
  {
    typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
    additional_data.mapping_update_flags =
      update_JxW_values | update_values | update_gradients;

    matrix_free.reinit(
      mapping, dof_handler, constraints, quadrature, additional_data);

    LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> delta_1_host(
      dof_handler.get_triangulation().n_active_cells());

    unsigned int c = 0;
    for (const auto &color : matrix_free.get_colored_graph())
      for (const auto &cell : color)
        delta_1_host[c++] = delta_1_host_dealii[cell->active_cell_index()];

    delta_1.reinit(dof_handler.get_triangulation().n_active_cells());

    LinearAlgebra::ReadWriteVector<Number> rw_vector(
      dof_handler.get_triangulation().n_active_cells());
    rw_vector.import_elements(delta_1_host, VectorOperation::insert);
    delta_1.import_elements(rw_vector, VectorOperation::insert);
  }

  types::global_dof_index
  m() const
  {
    return matrix_free.get_dof_handler().n_dofs();
  }

  Number
  el(unsigned int, unsigned int) const
  {
    DEAL_II_NOT_IMPLEMENTED();
    return 0;
  }

  template <typename VectorType2>
  void
  initialize_dof_vector(VectorType2 &vec) const
  {
    matrix_free.initialize_dof_vector(vec);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = 0.0; // TODO: annoying
    StokesOperatorLocal<dim, fe_degree, Number> local_operator(
      delta_1.get_values());
    matrix_free.cell_loop(local_operator, src, dst);
    matrix_free.copy_constrained_values(src, dst); // TODO: annoying
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    AssertThrow(false, ExcNotImplemented());

    (void)dst;
    (void)src;
  }

  void
  compute_inverse_diagonal(VectorType &diagonal_global) const
  {
    matrix_free.initialize_dof_vector(diagonal_global);
    StokesOperatorQuad<dim, fe_degree, Number> stokes_operator_quad(
      delta_1.get_values());
    MatrixFreeTools::
      compute_diagonal<dim, fe_degree, fe_degree + 1, dim + 1, Number>(
        matrix_free,
        diagonal_global,
        stokes_operator_quad,
        EvaluationFlags::values | EvaluationFlags::gradients,
        EvaluationFlags::values | EvaluationFlags::gradients);

    Number *diagonal_global_ptr = diagonal_global.get_values();

    Kokkos::parallel_for(
      "lethe::invert_vector",
      Kokkos::RangePolicy<MemorySpace::Default::kokkos_space::execution_space>(
        0, diagonal_global.locally_owned_size()),
      KOKKOS_LAMBDA(int i) {
        diagonal_global_ptr[i] = 1.0 / diagonal_global_ptr[i];
      });
  }

private:
  Portable::MatrixFree<dim, Number> matrix_free;

  VectorType delta_1;
};



template <int dim, typename T>
class RHSFunction : public Function<dim, T>
{
public:
  RHSFunction()
    : Function<dim, T>(dim + 1)
  {}

  virtual T
  value(const Point<dim, T> &p, const unsigned int c = 0) const override
  {
    const double a = numbers::PI;
    const double x = p[0];
    const double y = p[1];

    if (c == 0)
      return 2 * a * a *
               (std::sin(a * x) * std::sin(a * x) -
                std::cos(a * x) * std::cos(a * x)) *
               std::sin(a * y) * std::cos(a * y) +
             4 * a * a * std::sin(a * x) * std::sin(a * x) * std::sin(a * y) *
               std::cos(a * y) +
             a * std::sin(a * y) * std::cos(a * x);
    else if (c == 1)
      return -2 * a * a *
               (std::sin(a * y) * std::sin(a * y) -
                std::cos(a * y) * std::cos(a * y)) *
               std::sin(a * x) * std::cos(a * x) -
             4 * a * a * std::sin(a * x) * std::sin(a * y) * std::sin(a * y) *
               std::cos(a * x) +
             a * std::sin(a * x) * std::cos(a * y);
    else if (c == 2)
      return 0.0;

    AssertThrow(false, ExcNotImplemented());

    return 0.0;
  }
};



template <int dim, typename T>
class ExactSolution : public Function<dim, T>
{
public:
  ExactSolution()
    : Function<dim, T>(dim + 1)
  {}

  virtual T
  value(const Point<dim, T> &p, const unsigned int c = 0) const override
  {
    const double a = numbers::PI;
    const double x = p[0];
    const double y = p[1];

    if (c == 0)
      return std::sin(a * x) * std::sin(a * x) * std::cos(a * y) *
             std::sin(a * y);
    else if (c == 1)
      return -std::cos(a * x) * std::sin(a * x) * std::sin(a * y) *
             std::sin(a * y);
    else if (c == 2)
      return std::sin(a * x) * std::sin(a * y);

    AssertThrow(false, ExcNotImplemented());

    return 0.0;
  }
};



template <unsigned int dim, const int degree, typename MemorySpace>
void
run(const unsigned int n_refinements, ConvergenceTable &table)
{
  const MPI_Comm comm = MPI_COMM_WORLD;

  using Number     = double;
  using VectorType = LinearAlgebra::distributed::Vector<Number, MemorySpace>;
  using VectorTypeHost = LinearAlgebra::distributed::Vector<Number>;

  const bool   use_multigrid   = true;
  const double delta_1_scaling = 0.1;

  if (std::is_same_v<MemorySpace, dealii::MemorySpace::Host>)
    table.add_value("version", "host");
  else
    table.add_value("version", "default");

  table.add_value("fe_degree", degree);
  table.add_value("n_refinements", n_refinements);

  parallel::distributed::Triangulation<dim> tria(comm);

  GridGenerator::hyper_cube(tria, -1.0, +1.0);
  tria.refine_global(n_refinements);

  const MappingQ1<dim> mapping;
  const FE_Q<dim>      fe_q(degree);
  const FESystem<dim>  fe(fe_q, dim + 1);
  const QGauss<dim>    quadrature(degree + 1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  table.add_value("n_dofs", dof_handler.n_dofs());

  AffineConstraints<Number> constraints;
  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();

  VectorTypeHost all_delta_1(tria.n_active_cells());

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      all_delta_1[cell->active_cell_index()] =
        delta_1_scaling * cell->minimum_vertex_distance();

  StokesOperator<dim, degree, Number, MemorySpace> stokes_operator;

  stokes_operator.reinit(mapping,
                         dof_handler,
                         constraints,
                         quadrature.get_tensor_basis()[0],
                         all_delta_1);

  VectorType src, dst;

  stokes_operator.initialize_dof_vector(src);
  stokes_operator.initialize_dof_vector(dst);

  {
    VectorTypeHost src_host(src.get_partitioner());

    RHSFunction<dim, Number> rhs_func;

    VectorTools::create_right_hand_side<dim, dim>(
      mapping, dof_handler, quadrature, rhs_func, src_host, constraints);

    FEValues<dim>              fe_values(mapping,
                            fe,
                            quadrature,
                            update_values | update_gradients |
                              update_JxW_values | update_quadrature_points);
    FEValuesViews::Vector<dim> velocities(fe_values, 0);
    FEValuesViews::Scalar<dim> pressure(fe_values, dim);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          const double delta_1 = all_delta_1[cell->active_cell_index()];

          Vector<double>                       rhs_local(fe.n_dofs_per_cell());
          std::vector<types::global_dof_index> indices(fe.n_dofs_per_cell());

          cell->get_dof_indices(indices);

          for (const unsigned int q : fe_values.quadrature_point_indices())
            {
              const auto JxW   = fe_values.JxW(q);
              const auto point = fe_values.quadrature_point(q);

              Tensor<1, dim> source;
              for (unsigned int d = 0; d < dim; ++d)
                source[d] = rhs_func.value(point, d);

              for (const unsigned int i : fe_values.dof_indices())
                rhs_local(i) +=
                  delta_1 * source * pressure.gradient(i, q) * JxW;
            }

          constraints.distribute_local_to_global(rhs_local, indices, src_host);
        }

    src_host.compress(VectorOperation::add);

    LinearAlgebra::ReadWriteVector<Number> rw_vector(
      src.get_partitioner()->locally_owned_range());
    rw_vector.import_elements(src_host, VectorOperation::insert);
    src.import_elements(rw_vector, VectorOperation::insert);

    dst = 0.0;
  }

  ReductionControl        solver_control;
  SolverGMRES<VectorType> solver(solver_control);

  if (!use_multigrid)
    {
      DiagonalMatrix<VectorType> preconditioner;
      stokes_operator.compute_inverse_diagonal(preconditioner.get_vector());

      solver.solve(stokes_operator, dst, src, preconditioner);
    }
  else
    {
      using LevelMatrixType = StokesOperator<dim, degree, Number, MemorySpace>;
      using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
      using SmootherType =
        PreconditionRelaxation<LevelMatrixType, SmootherPreconditionerType>;
      using MGTransferType = typename std::conditional<
        std::is_same_v<MemorySpace, dealii::MemorySpace::Host>,
        MGTransferMF<dim, Number>,
        MGTransferMFWrapper<dim, Number, MemorySpace>>::type;

      const auto coarse_grid_triangulations =
        MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
          dof_handler.get_triangulation());

      const unsigned int min_level = 0;
      const unsigned int max_level = coarse_grid_triangulations.size() - 1;

      MGLevelObject<DoFHandler<dim>> mg_dof_handlers(min_level, max_level);
      MGLevelObject<AffineConstraints<Number>> mg_constraints(min_level,
                                                              max_level);
      MGLevelObject<LevelMatrixType> mg_matrices(min_level, max_level);

      MGLevelObject<MGTwoLevelTransfer<dim, VectorTypeHost>> mg_transfers(
        min_level, max_level);

      // level operators
      for (unsigned int level = min_level; level <= max_level; ++level)
        {
          const auto &tria        = *coarse_grid_triangulations[level];
          auto       &dof_handler = mg_dof_handlers[level];
          auto       &constraint  = mg_constraints[level];

          dof_handler.reinit(tria);
          dof_handler.distribute_dofs(fe);

          constraint.reinit(dof_handler.locally_owned_dofs(),
                            DoFTools::extract_locally_relevant_dofs(
                              dof_handler));

          DoFTools::make_zero_boundary_constraints(dof_handler, constraint);
          constraint.close();

          VectorTypeHost all_delta_1(tria.n_active_cells());

          for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->is_locally_owned())
              all_delta_1[cell->active_cell_index()] =
                delta_1_scaling * cell->minimum_vertex_distance();

          mg_matrices[level].reinit(mapping,
                                    dof_handler,
                                    constraint,
                                    quadrature.get_tensor_basis()[0],
                                    all_delta_1);
        }

      mg::Matrix<VectorType> mg_matrix(mg_matrices);

      // transfer operator
      for (unsigned int level = min_level; level < max_level; ++level)
        mg_transfers[level + 1].reinit(mg_dof_handlers[level + 1],
                                       mg_dof_handlers[level],
                                       mg_constraints[level + 1],
                                       mg_constraints[level]);

      MGTransferType mg_transfer(mg_transfers, [&](const auto l, auto &vec) {
        mg_matrices[l].initialize_dof_vector(vec);
      });

      // smoother
      MGLevelObject<typename SmootherType::AdditionalData> smoother_data(
        min_level, max_level);

      for (unsigned int level = min_level; level <= max_level; ++level)
        {
          smoother_data[level].preconditioner =
            std::make_shared<SmootherPreconditionerType>();
          mg_matrices[level].compute_inverse_diagonal(
            smoother_data[level].preconditioner->get_vector());
          smoother_data[level].smoothing_range     = 20;
          smoother_data[level].n_iterations        = 5;
          smoother_data[level].eig_cg_n_iterations = 20;
          smoother_data[level].eigenvalue_algorithm =
            SmootherType::AdditionalData::EigenvalueAlgorithm::power_iteration;
          smoother_data[level].constraints.copy_from(mg_constraints[level]);
        }

      MGSmootherPrecondition<LevelMatrixType, SmootherType, VectorType>
        mg_smoother;
      mg_smoother.initialize(mg_matrices, smoother_data);

      for (unsigned int level = min_level; level <= max_level; ++level)
        {
          VectorType vec;
          mg_matrices[level].initialize_dof_vector(vec);
          mg_smoother.smoothers[level].estimate_eigenvalues(vec);
        }

      // coarse-grid solver
      MGCoarseGridApplySmoother<VectorType> mg_coarse;
      mg_coarse.initialize(mg_smoother);

      // put everything together
      Multigrid<VectorType> mg(
        mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);

      PreconditionMG<dim, VectorType, MGTransferType> preconditioner(
        dof_handler, mg, mg_transfer);

      // solve
      solver.solve(stokes_operator, dst, src, preconditioner);
    }

  table.add_value("n_iterations", solver_control.last_step());

  {
    VectorTypeHost dst_host(dst.get_partitioner());

    LinearAlgebra::ReadWriteVector<Number> rw_vector(
      src.get_partitioner()->locally_owned_range());
    rw_vector.import_elements(dst, VectorOperation::insert);
    dst_host.import_elements(rw_vector, VectorOperation::insert);

    ExactSolution<dim, Number> exact_solution;

    const ComponentSelectFunction<dim> u_mask(std::make_pair(0, dim), dim + 1);
    const ComponentSelectFunction<dim> p_mask(dim, dim + 1);

    Vector<Number> cell_wise_error;

    dst_host.update_ghost_values();
    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      dst_host,
                                      exact_solution,
                                      cell_wise_error,
                                      quadrature,
                                      VectorTools::NormType::L2_norm,
                                      &u_mask);

    const auto error_u =
      VectorTools::compute_global_error(tria,
                                        cell_wise_error,
                                        VectorTools::NormType::L2_norm);

    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      dst_host,
                                      exact_solution,
                                      cell_wise_error,
                                      quadrature,
                                      VectorTools::NormType::L2_norm,
                                      &p_mask);

    const auto error_p =
      VectorTools::compute_global_error(tria,
                                        cell_wise_error,
                                        VectorTools::NormType::L2_norm);

    table.add_value("error_u", error_u);
    table.set_scientific("error_u", true);
    table.add_value("error_p", error_p);
    table.set_scientific("error_p", true);
  }

  {
    VectorTypeHost dst_host(dst.get_partitioner());

    LinearAlgebra::ReadWriteVector<Number> rw_vector(
      src.get_partitioner()->locally_owned_range());
    rw_vector.import_elements(dst, VectorOperation::insert);
    dst_host.import_elements(rw_vector, VectorOperation::insert);

    std::string file_name = "solution_" + std::to_string(counter++) + ".vtu";

    DataOut<dim> data_out;

    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    std::vector<std::string> labels(dim + 1, "u");
    labels[dim] = "p";

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim + 1, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation[dim] =
      DataComponentInterpretation::component_is_scalar;

    data_out.add_data_vector(dof_handler,
                             dst_host,
                             labels,
                             data_component_interpretation);
    data_out.build_patches(mapping,
                           degree + 1,
                           DataOut<dim>::CurvedCellRegion::curved_inner_cells);
    data_out.write_vtu_in_parallel(file_name, MPI_COMM_WORLD);
  }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  const unsigned int dim           = 2;
  const unsigned int fe_degree     = 3;
  unsigned int       n_refinements = 3;

  ConvergenceTable table;

  run<dim, fe_degree, MemorySpace::Host>(n_refinements, table);
  run<dim, fe_degree, MemorySpace::Default>(n_refinements, table);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    table.write_text(std::cout);
}
