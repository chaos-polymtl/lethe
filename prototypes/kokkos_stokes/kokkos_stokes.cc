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
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

static unsigned int counter = 0;

template <int dim, int fe_degree, typename Number, typename MemorySpace>
class StokesOperator;

template <int dim, int fe_degree, typename Number>
class StokesOperator<dim, fe_degree, Number, MemorySpace::Host>
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

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    matrix_free.cell_loop(&StokesOperator::local_apply, this, dst, src, true);
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
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;

  StokesOperatorQuad(const VectorType &delta_1)
    : delta_1(delta_1.get_values())
  {}

  StokesOperatorQuad(
    const typename Portable::MatrixFree<dim, Number>::Data *data,
    int                                                     cell,
    const VectorType                                       &delta_1)
    : data(data)
    , cell(cell)
    , delta_1(delta_1.get_values())
  {}

  DEAL_II_HOST_DEVICE void
  operator()(FECellIntegrator *phi, const int q_point) const
  {
    const auto delta_1 = this->delta_1[data->local_q_point_id(cell, 1, 0)];

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

  DEAL_II_HOST_DEVICE
  void
  set_cell(int cell)
  {
    this->cell = cell;
  }

  DEAL_II_HOST_DEVICE void
  set_matrix_free_data(
    const typename Portable::MatrixFree<dim, Number>::Data &data)
  {
    this->data = &data;
  }

private:
  const typename Portable::MatrixFree<dim, Number>::Data *data;
  int                                                     cell;
  const Number                                           *delta_1;
};

template <int dim, int fe_degree, typename Number>
class StokesOperatorLocal
{
public:
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;

  StokesOperatorLocal(const VectorType &delta_1)
    : delta_1(delta_1)
  {}

  DEAL_II_HOST_DEVICE void
  operator()(const unsigned int                                      cell,
             const typename Portable::MatrixFree<dim, Number>::Data *gpu_data,
             Portable::SharedData<dim, Number> *shared_data,
             const Number                      *src,
             Number                            *dst) const
  {
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, Number> phi(
      /*cell,*/ gpu_data, shared_data);
    phi.read_dof_values(src);
    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi.apply_for_each_quad_point(
      StokesOperatorQuad<dim, fe_degree, Number>(gpu_data, cell, delta_1));
    phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi.distribute_local_to_global(dst);
  }

  static const unsigned int n_local_dofs =
    Utilities::pow(fe_degree + 1, dim) * (dim + 1);
  static const unsigned int n_q_points = Utilities::pow(fe_degree + 1, dim);

private:
  const VectorType &delta_1;
};

template <int dim, int fe_degree, typename Number>
class StokesOperator<dim, fe_degree, Number, MemorySpace::Default>
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
    rw_vector.import(delta_1_host, VectorOperation::insert);
    delta_1.import(rw_vector, VectorOperation::insert);
  }

  void
  initialize_dof_vector(VectorType &vec) const
  {
    matrix_free.initialize_dof_vector(vec);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = 0.0; // TODO: annoying
    StokesOperatorLocal<dim, fe_degree, Number> local_operator(delta_1);
    matrix_free.cell_loop(local_operator, src, dst);
    matrix_free.copy_constrained_values(src, dst); // TODO: annoying
  }

  void
  compute_inverse_diagonal(VectorType &diagonal_global) const
  {
    matrix_free.initialize_dof_vector(diagonal_global);
    StokesOperatorQuad<dim, fe_degree, Number> laplace_operator_quad(delta_1);
    MatrixFreeTools::
      compute_diagonal<dim, fe_degree, fe_degree + 1, dim + 1, Number>(
        matrix_free,
        diagonal_global,
        laplace_operator_quad,
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

  LinearAlgebra::distributed::Vector<Number> all_delta_1(tria.n_active_cells());

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      all_delta_1[cell->active_cell_index()] =
        delta_1_scaling * cell->minimum_vertex_distance();

  StokesOperator<dim, degree, Number, MemorySpace> laplace_operator;

  laplace_operator.reinit(mapping,
                          dof_handler,
                          constraints,
                          quadrature.get_tensor_basis()[0],
                          all_delta_1);

  VectorType src, dst;

  laplace_operator.initialize_dof_vector(src);
  laplace_operator.initialize_dof_vector(dst);

  {
    LinearAlgebra::distributed::Vector<Number> src_host(src.get_partitioner());

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

    LinearAlgebra::ReadWriteVector<Number> rw_vector(
      src.get_partitioner()->locally_owned_range());
    rw_vector.import(src_host, VectorOperation::insert);
    src.import(rw_vector, VectorOperation::insert);

    dst = 0.0;
  }

  // PreconditionIdentity preconditioner;
  DiagonalMatrix<VectorType> preconditioner;
  laplace_operator.compute_inverse_diagonal(preconditioner.get_vector());

  ReductionControl        solver_control;
  SolverGMRES<VectorType> solver(solver_control);
  solver.solve(laplace_operator, dst, src, preconditioner);

  table.add_value("n_iterations", solver_control.last_step());

  {
    LinearAlgebra::distributed::Vector<Number> dst_host(dst.get_partitioner());

    LinearAlgebra::ReadWriteVector<Number> rw_vector(
      src.get_partitioner()->locally_owned_range());
    rw_vector.import(dst, VectorOperation::insert);
    dst_host.import(rw_vector, VectorOperation::insert);

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
    LinearAlgebra::distributed::Vector<Number> dst_host(dst.get_partitioner());

    LinearAlgebra::ReadWriteVector<Number> rw_vector(
      src.get_partitioner()->locally_owned_range());
    rw_vector.import(dst, VectorOperation::insert);
    dst_host.import(rw_vector, VectorOperation::insert);

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

  table.write_text(std::cout);
}
