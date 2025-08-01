// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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

template <int dim,
          int fe_degree,
          int n_components,
          typename Number,
          typename MemorySpace>
class LaplaceOperator;

template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperator<dim, fe_degree, n_components, Number, MemorySpace::Host>
  : public EnableObserverPointer
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
  initialize_dof_vector(VectorType &vec) const
  {
    matrix_free.initialize_dof_vector(vec);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    matrix_free.cell_loop(&LaplaceOperator::local_apply, this, dst, src, true);
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
                                      &LaplaceOperator::do_vmult_cell_single,
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
    FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, Number> phi(data);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);

        phi.read_dof_values_plain(src);
        do_vmult_cell_single(phi);
        phi.distribute_local_to_global(dst);
      }
  }

  void
  do_vmult_cell_single(
    FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, Number> &phi)
    const
  {
    phi.evaluate(EvaluationFlags::gradients);
    for (unsigned int q = 0; q < phi.n_q_points; ++q)
      phi.submit_gradient(phi.get_gradient(q), q);
    phi.integrate(EvaluationFlags::gradients);
  }

  MatrixFree<dim, Number> matrix_free;
};



template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperatorQuad
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, Number>
             *phi,
    const int q_point) const
  {
    phi->submit_gradient(phi->get_gradient(q_point), q_point);
  }

  DEAL_II_HOST_DEVICE
  void
  set_cell(int)
  {}

  DEAL_II_HOST_DEVICE void
  set_matrix_free_data(const typename Portable::MatrixFree<dim, Number>::Data &)
  {}
};

template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperatorLocal
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceVector<Number>                   &src,
             Portable::DeviceVector<Number>                         &dst) const
  {
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, Number>
      phi(data);
    phi.read_dof_values(src);
    phi.evaluate(EvaluationFlags::gradients);
    phi.apply_for_each_quad_point(
      LaplaceOperatorQuad<dim, fe_degree, n_components, Number>());
    phi.integrate(EvaluationFlags::gradients);
    phi.distribute_local_to_global(dst);
  }

  static const unsigned int n_local_dofs =
    Utilities::pow(fe_degree + 1, dim) * n_components;
  static const unsigned int n_q_points = Utilities::pow(fe_degree + 1, dim);
};

template <int dim, int fe_degree, int n_components, typename Number>
class LaplaceOperator<dim,
                      fe_degree,
                      n_components,
                      Number,
                      MemorySpace::Default> : public EnableObserverPointer
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
    typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_JxW_values | update_gradients;

    matrix_free.reinit(
      mapping, dof_handler, constraints, quadrature, additional_data);
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
    LaplaceOperatorLocal<dim, fe_degree, n_components, Number> local_operator;
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
    LaplaceOperatorQuad<dim, fe_degree, n_components, Number>
      laplace_operator_quad;
    MatrixFreeTools::
      compute_diagonal<dim, fe_degree, fe_degree + 1, n_components, Number>(
        matrix_free,
        diagonal_global,
        laplace_operator_quad,
        EvaluationFlags::gradients,
        EvaluationFlags::gradients);

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
  using VectorTypeHost = LinearAlgebra::distributed::Vector<Number>;

  const bool use_multigrid = true;

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
    VectorTypeHost src_host(src.get_partitioner());

    VectorTools::create_right_hand_side<dim, dim>(
      mapping,
      dof_handler,
      quadrature,
      AnalyticalFunction<dim, Number>(n_components),
      src_host,
      constraints);

    LinearAlgebra::ReadWriteVector<Number> rw_vector(
      src.get_partitioner()->locally_owned_range());
    rw_vector.import_elements(src_host, VectorOperation::insert);
    src.import_elements(rw_vector, VectorOperation::insert);

    dst = 0.0;
  }

  ReductionControl     solver_control;
  SolverCG<VectorType> solver(solver_control);

  if (!use_multigrid)
    {
      DiagonalMatrix<VectorType> preconditioner;
      laplace_operator.compute_inverse_diagonal(preconditioner.get_vector());

      solver.solve(laplace_operator, dst, src, preconditioner);
    }
  else
    {
      using LevelMatrixType =
        LaplaceOperator<dim, degree, n_components, Number, MemorySpace>;
      using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
      using SmootherType               = PreconditionChebyshev<LevelMatrixType,
                                                 VectorType,
                                                 SmootherPreconditionerType>;
      using MGTransferType             = typename std::conditional<
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
          auto &dof_handler = mg_dof_handlers[level];
          auto &constraint  = mg_constraints[level];

          dof_handler.reinit(*coarse_grid_triangulations[level]);
          dof_handler.distribute_dofs(fe);

          constraint.reinit(dof_handler.locally_owned_dofs(),
                            DoFTools::extract_locally_relevant_dofs(
                              dof_handler));

          DoFTools::make_zero_boundary_constraints(dof_handler, constraint);
          constraint.close();

          mg_matrices[level].reinit(mapping,
                                    dof_handler,
                                    constraint,
                                    quadrature.get_tensor_basis()[0]);
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
          smoother_data[level].degree              = 5;
          smoother_data[level].eig_cg_n_iterations = 20;
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
      solver.solve(laplace_operator, dst, src, preconditioner);
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
  table.add_value("n_iterations", solver_control.last_step());

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

#if DEAL_II_VERSION_GTE(9, 7, 0)
  run<dim, fe_degree, dim, MemorySpace::Host>(n_refinements, table);
  run<dim, fe_degree, dim, MemorySpace::Default>(n_refinements, table);
  run<dim, fe_degree, dim + 1, MemorySpace::Host>(n_refinements, table);
  run<dim, fe_degree, dim + 1, MemorySpace::Default>(n_refinements, table);
#endif

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    table.write_text(std::cout);
}
