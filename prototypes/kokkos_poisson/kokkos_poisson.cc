/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2024 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * A simple matrix-free solver for solving a vector Laplace problem
 * using Kokkos via Portable::MatrixFree and Porable::FEEvaluation.
 *
 * ---------------------------------------------------------------------*/

#include <deal.II/base/convergence_table.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/repartitioning_policy_tools.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/matrix_free/cuda_fe_evaluation.h>
#include <deal.II/matrix_free/cuda_matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

static unsigned int counter = 0;

template <int dim,
          int fe_degree,
          int n_components,
          typename Number,
          typename MemorySpace>
class LaplaceOperator;

template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperator<dim, fe_degree, n_components, Number, MemorySpace::Host>
{
public:
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>;

  LaplaceOperator() = default;

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<1>             &quadrature)
  {
    typename MatrixFree<dim, Number>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_gradients;

    matrix_free.reinit(
      mapping, dof_handler, constraints, quadrature, additional_data);
  }

  void
  initialize_dof_vector(VectorType &vec) const
  {
    matrix_free.initialize_dof_vector(vec);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    matrix_free.cell_loop(&LaplaceOperator::local_apply, this, dst, src, true);
  }

private:
  void
  local_apply(const MatrixFree<dim, Number>               &data,
              VectorType                                  &dst,
              const VectorType                            &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, Number> phi(data);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);

        phi.read_dof_values_plain(src);
        phi.evaluate(EvaluationFlags::gradients);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate(EvaluationFlags::gradients);
        phi.distribute_local_to_global(dst);
      }
  }

  MatrixFree<dim, Number> matrix_free;
};



template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperatorQuad
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(CUDAWrappers::
               FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, Number>
                      *fe_eval,
             const int q_point) const
  {
    fe_eval->submit_gradient(fe_eval->get_gradient(q_point), q_point);
  }
};

template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperatorLocal
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(
    const unsigned int                                          cell,
    const typename CUDAWrappers::MatrixFree<dim, Number>::Data *gpu_data,
    CUDAWrappers::SharedData<dim, Number>                      *shared_data,
    const Number                                               *src,
    Number                                                     *dst) const
  {
    (void)cell; // TODO?

    CUDAWrappers::
      FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, Number>
        fe_eval(
          /*cell,*/ gpu_data, shared_data);
    fe_eval.read_dof_values(src);
    fe_eval.evaluate(false, true);
    fe_eval.apply_for_each_quad_point(
      LaplaceOperatorQuad<dim, fe_degree, n_components, Number>());
    fe_eval.integrate(false, true);
    fe_eval.distribute_local_to_global(dst);
  }
  static const unsigned int n_dofs_1d    = fe_degree + 1;
  static const unsigned int n_local_dofs = Utilities::pow(fe_degree + 1, dim);
  static const unsigned int n_q_points   = Utilities::pow(fe_degree + 1, dim);
};

template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperator<dim,
                      fe_degree,
                      n_components,
                      Number,
                      MemorySpace::Default>
{
public:
  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;

  LaplaceOperator() = default;

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<1>             &quadrature)
  {
    typename CUDAWrappers::MatrixFree<dim, Number>::AdditionalData
      additional_data;
    additional_data.mapping_update_flags = update_JxW_values | update_gradients;

    matrix_free.reinit(
      mapping, dof_handler, constraints, quadrature, additional_data);
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
    LaplaceOperatorLocal<dim, fe_degree, n_components, Number> local_operator;
    matrix_free.cell_loop(local_operator, src, dst);
    matrix_free.copy_constrained_values(src, dst); // TODO: annoying
  }

private:
  CUDAWrappers::MatrixFree<dim, Number> matrix_free;
};



template <int dim, typename T>
class AnalyticalFunction : public Function<dim, T>
{
public:
  AnalyticalFunction(const unsigned int n_components)
    : Function<dim, T>(n_components)
  {}

  virtual T
  value(const Point<dim, T> &p, const unsigned int component = 0) const override
  {
    double temp = 0.0;

    for (unsigned int d = 0; d < dim; ++d)
      temp += std::sin(p[d]);

    return temp * (1.0 + component);
  }
};



template <unsigned int dim,
          const int    degree,
          int          n_components,
          typename MemorySpace>
void
run(const unsigned int n_refinements, ConvergenceTable &table)
{
  const MPI_Comm comm = MPI_COMM_WORLD;

  using Number     = double;
  using VectorType = LinearAlgebra::distributed::Vector<Number, MemorySpace>;

  parallel::distributed::Triangulation<dim> tria(comm);

  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_refinements);

  const MappingQ1<dim> mapping;
  const FE_Q<dim>      fe_q(degree);
  const FESystem<dim>  fe(fe_q, n_components);
  const QGauss<dim>    quadrature(degree + 1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<Number> constraints;
  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  constraints.close();

  LaplaceOperator<dim, degree, n_components, Number, MemorySpace>
    laplace_operator;

  laplace_operator.reinit(mapping,
                          dof_handler,
                          constraints,
                          quadrature.get_tensor_basis()[0]);

  VectorType src, dst;

  laplace_operator.initialize_dof_vector(src);
  laplace_operator.initialize_dof_vector(dst);

  {
    LinearAlgebra::distributed::Vector<Number> src_host(src.get_partitioner());

    VectorTools::create_right_hand_side<dim, dim>(
      mapping,
      dof_handler,
      quadrature,
      AnalyticalFunction<dim, Number>(n_components),
      src_host,
      constraints);

    LinearAlgebra::ReadWriteVector<Number> rw_vector(
      src.get_partitioner()->locally_owned_range());
    rw_vector.import(src_host, VectorOperation::insert);
    src.import(rw_vector, VectorOperation::insert);

    dst = 0.0;
  }

  PreconditionIdentity preconditioner;

  ReductionControl     solver_control;
  SolverCG<VectorType> solver(solver_control);
  solver.solve(laplace_operator, dst, src, preconditioner);

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

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(dst_host, "solution");
    data_out.build_patches(mapping,
                           degree + 1,
                           DataOut<dim>::CurvedCellRegion::curved_inner_cells);
    data_out.write_vtu_in_parallel(file_name, MPI_COMM_WORLD);
  }

  table.add_value("fe_degree", degree);
  table.add_value("n_refinements", n_refinements);
  table.add_value("n_components", n_components);
  table.add_value("n_dofs", dof_handler.n_dofs());

  if (std::is_same_v<MemorySpace, dealii::MemorySpace::Host>)
    table.add_value("version", "host");
  else
    table.add_value("version", "default");

  table.add_value("norm", dst.l2_norm());
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  const unsigned int dim           = 2;
  const unsigned int fe_degree     = 3;
  unsigned int       n_refinements = 3;

  ConvergenceTable table;

  run<dim, fe_degree, 1, MemorySpace::Host>(n_refinements, table);
  run<dim, fe_degree, 1, MemorySpace::Default>(n_refinements, table);
  run<dim, fe_degree, dim, MemorySpace::Host>(n_refinements, table);
  run<dim, fe_degree, dim, MemorySpace::Default>(n_refinements, table);
  run<dim, fe_degree, dim + 1, MemorySpace::Host>(n_refinements, table);
  run<dim, fe_degree, dim + 1, MemorySpace::Default>(n_refinements, table);

  table.write_text(std::cout);
}
