// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief PoissonOperator: test in a multigrid setting.
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../mortar_functions/generate_grids.h"
#include "../mortar_functions/mortar_operator.h"

using namespace dealii;


/**
 * @brief Parameters for the geometric multigrid
 */
struct GMGParameters
{
  struct CoarseSolverParameters
  {
    std::string  type            = "cg_with_amg"; // "cg";
    unsigned int maxiter         = 10000;
    double       abstol          = 1e-20;
    double       reltol          = 1e-4;
    unsigned int smoother_sweeps = 1;
    unsigned int n_cycles        = 1;
    std::string  smoother_type   = "ILU";
  };

  struct SmootherParameters
  {
    std::string  type                = "chebyshev";
    double       smoothing_range     = 20;
    unsigned int degree              = 5;
    unsigned int eig_cg_n_iterations = 20;
  };

  SmootherParameters     smoother;
  CoarseSolverParameters coarse_solver;

  unsigned int maxiter = 10000;
  double       abstol  = 1e-20;
  double       reltol  = 1e-4;
};

template <typename VectorType,
          int dim,
          typename SystemMatrixType,
          typename LevelMatrixType,
          typename MGTransferType>
static void
mg_solve(SolverControl                        &solver_control,
         VectorType                           &dst,
         const VectorType                     &src,
         const GMGParameters                  &mg_data,
         const DoFHandler<dim>                &dof,
         const SystemMatrixType               &fine_matrix,
         const MGLevelObject<LevelMatrixType> &mg_matrices,
         const MGTransferType                 &mg_transfer)
{
  AssertThrow(mg_data.smoother.type == "chebyshev", ExcNotImplemented());

  const unsigned int min_level = mg_matrices.min_level();
  const unsigned int max_level = mg_matrices.max_level();

  using Number                     = typename VectorType::value_type;
  using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
  using SmootherType               = PreconditionChebyshev<LevelMatrixType,
                                             VectorType,
                                             SmootherPreconditionerType>;
  using PreconditionerType = PreconditionMG<dim, VectorType, MGTransferType>;

  // Initialize level operators
  mg::Matrix<VectorType> mg_matrix(mg_matrices);

  // Initialize smoothers
  MGLevelObject<typename SmootherType::AdditionalData> smoother_data(min_level,
                                                                     max_level);

  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      smoother_data[level].preconditioner =
        std::make_shared<SmootherPreconditionerType>();
      mg_matrices[level].compute_inverse_diagonal(
        smoother_data[level].preconditioner->get_vector());
      smoother_data[level].smoothing_range = mg_data.smoother.smoothing_range;
      smoother_data[level].degree          = mg_data.smoother.degree;
      smoother_data[level].eig_cg_n_iterations =
        mg_data.smoother.eig_cg_n_iterations;
    }

  MGSmootherPrecondition<LevelMatrixType, SmootherType, VectorType> mg_smoother;
  mg_smoother.initialize(mg_matrices, smoother_data);

  // Initialize coarse-grid solver
  ReductionControl coarse_grid_solver_control(mg_data.coarse_solver.maxiter,
                                              mg_data.coarse_solver.abstol,
                                              mg_data.coarse_solver.reltol,
                                              false,
                                              false);
  SolverGMRES<VectorType> coarse_grid_solver(coarse_grid_solver_control);

  PreconditionIdentity precondition_identity;
  PreconditionChebyshev<LevelMatrixType, VectorType, DiagonalMatrix<VectorType>>
    precondition_chebyshev;

#ifdef DEAL_II_WITH_TRILINOS
  TrilinosWrappers::PreconditionAMG precondition_amg;
#endif

  std::unique_ptr<MGCoarseGridBase<VectorType>> mg_coarse;

  if (mg_data.coarse_solver.type == "cg_with_amg")
    {
      // CG with AMG as preconditioner

#ifdef DEAL_II_WITH_TRILINOS
      TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
      amg_data.smoother_sweeps = mg_data.coarse_solver.smoother_sweeps;
      amg_data.n_cycles        = mg_data.coarse_solver.n_cycles;
      amg_data.smoother_type   = mg_data.coarse_solver.smoother_type.c_str();

      // CG with AMG as preconditioner
      precondition_amg.initialize(mg_matrices[min_level].get_system_matrix(),
                                  amg_data);

      mg_coarse = std::make_unique<
        MGCoarseGridIterativeSolver<VectorType,
                                    SolverGMRES<VectorType>,
                                    LevelMatrixType,
                                    decltype(precondition_amg)>>(
        coarse_grid_solver, mg_matrices[min_level], precondition_amg);
#else
      AssertThrow(false, ExcNotImplemented());
#endif
    }
  else
    {
      AssertThrow(false, ExcNotImplemented());
    }

  // Create multigrid object
  Multigrid<VectorType> mg(mg_matrix,
                           *mg_coarse,
                           mg_transfer,
                           mg_smoother,
                           mg_smoother,
                           min_level,
                           max_level);

  // Convert it to a preconditioner
  PreconditionerType preconditioner(dof, mg, mg_transfer);

  // Finally, solve
  SolverGMRES<VectorType>(solver_control)
    .solve(fine_matrix, dst, src, preconditioner);
}


template <int dim, typename Number = double>
void
test(const unsigned int n_refinements, const unsigned int fe_degree_fine)
{
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  const unsigned int min_level = 0;
  const unsigned int max_level = n_refinements;

  const double radius = 1.0;
  const double rotate = 360. / (4 * Utilities::pow(2, n_refinements + 1)) * 0.9;
  const double rotate_pi = 2 * numbers::PI * rotate / 360.0;

  const MPI_Comm comm = MPI_COMM_WORLD;

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(comm) == 0);

  MGLevelObject<std::unique_ptr<MappingQ<dim>>> mappings(min_level, max_level);
  MGLevelObject<std::unique_ptr<Quadrature<dim>>> quads(min_level, max_level);
  MGLevelObject<parallel::distributed::Triangulation<dim>> triangulations(
    min_level, max_level, comm);
  MGLevelObject<DoFHandler<dim>>           dof_handlers(min_level, max_level);
  MGLevelObject<AffineConstraints<Number>> constraints(min_level, max_level);
  MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers(min_level,
                                                               max_level);
  MGLevelObject<PoissonOperator<dim, 1, Number>>     operators(min_level,
                                                           max_level);

  std::unique_ptr<Mapping<dim>> mapping_;

  // set up levels
  for (auto l = min_level; l <= max_level; ++l)
    {
      auto &tria        = triangulations[l];
      auto &dof_handler = dof_handlers[l];
      auto &constraint  = constraints[l];
      auto &op          = operators[l];

      std::unique_ptr<FiniteElement<dim>> fe;

      fe          = std::make_unique<FE_Q<dim>>(fe_degree_fine);
      quads[l]    = std::make_unique<QGauss<dim>>(fe_degree_fine + 1);
      mappings[l] = std::make_unique<MappingQ<dim>>(fe_degree_fine);

      if (l == max_level)
        mapping_ = mappings[l]->clone();

      // set up triangulation
      hyper_cube_with_cylindrical_hole(radius, 2.0, rotate_pi, tria);
      tria.refine_global(l);

      // set up dofhandler
      dof_handler.reinit(tria);
      dof_handler.distribute_dofs(*fe);

      // set up constraints
      const IndexSet locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);
      constraint.reinit(dof_handler.locally_owned_dofs(),
                        locally_relevant_dofs);
      DoFTools::make_zero_boundary_constraints(dof_handler, 1, constraint);
      DoFTools::make_zero_boundary_constraints(dof_handler, 2, constraint);
      DoFTools::make_zero_boundary_constraints(dof_handler, 3, constraint);
      DoFTools::make_zero_boundary_constraints(dof_handler, 4, constraint);
      constraint.close();

      // set up operator
      op.reinit(*mappings[l], dof_handler, constraint, *quads[l]);

      const auto mortar_manager =
        std::make_shared<MortarManagerCircle<dim>>(4 * Utilities::pow(2, l + 1),
                                                   radius,
                                                   construct_quadrature(
                                                     *quads[l]),
                                                   rotate_pi);

      op.add_coupling(mortar_manager, 0, 5);
    }

  GMGParameters mg_data; // TODO

  ReductionControl solver_control(
    mg_data.maxiter, mg_data.abstol, mg_data.reltol, false, false);

  if (true)
    {
      // set up transfer operator
      for (unsigned int l = min_level; l < max_level; ++l)
        transfers[l + 1].reinit(dof_handlers[l + 1],
                                dof_handlers[l],
                                constraints[l + 1],
                                constraints[l]);

      MGTransferGlobalCoarsening<dim, VectorType> transfer(
        transfers, [&](const auto l, auto &vec) {
          operators[l].initialize_dof_vector(vec);
        });

      VectorType dst, src;
      operators[max_level].initialize_dof_vector(dst);
      operators[max_level].initialize_dof_vector(src);
      VectorTools::create_right_hand_side(*mapping_,
                                          dof_handlers[max_level],
                                          *quads[max_level],
                                          Functions::ConstantFunction<dim>(1.0),
                                          src,
                                          constraints[max_level]);

      mg_solve(solver_control,
               dst,
               src,
               mg_data,
               dof_handlers[max_level],
               operators[max_level],
               operators,
               transfer);

      pcout << dim << ' ' << fe_degree_fine << ' ' << n_refinements << ' '
            << solver_control.last_step() << std::endl;

      static unsigned int counter = 0;

      MGLevelObject<VectorType> results(min_level, max_level);

      transfer.interpolate_to_mg(dof_handlers[max_level], results, dst);

      for (unsigned int l = min_level; l <= max_level; ++l)
        {
          DataOut<dim> data_out;

          DataOutBase::VtkFlags flags;
          flags.write_higher_order_cells = true;
          data_out.set_flags(flags);

          data_out.attach_dof_handler(dof_handlers[l]);
          data_out.add_data_vector(
            results[l],
            "solution",
            DataOut_DoFData<dim, dim>::DataVectorType::type_dof_data);

          Vector<double> ranks(triangulations[l].n_active_cells());
          ranks = Utilities::MPI::this_mpi_process(comm);
          data_out.add_data_vector(ranks, "ranks");

          data_out.build_patches(*mapping_, fe_degree_fine);

          std::string file_name = "test." + std::to_string(dim) + "." +
                                  std::to_string(counter) + "." +
                                  std::to_string(l) + ".vtu";
          data_out.write_vtu_in_parallel(file_name, comm);
        }

      counter++;
    }
  else
    {
      static unsigned int counter = 0;

      MGLevelObject<VectorType> results(min_level, max_level);

      for (unsigned int l = min_level; l <= max_level; ++l)
        {
          VectorType dst, src;
          operators[l].initialize_dof_vector(dst);
          operators[l].initialize_dof_vector(src);
          VectorTools::create_right_hand_side(*mapping_,
                                              dof_handlers[l],
                                              *quads[l],
                                              Functions::ConstantFunction<dim>(
                                                1.0),
                                              src,
                                              constraints[l]);

          TrilinosWrappers::SolverDirect preconditioner;
          preconditioner.initialize(operators[l].get_system_matrix());
          SolverGMRES<VectorType> solver(solver_control);
          solver.solve(operators[l], dst, src, preconditioner);

          pcout << dim << ' ' << fe_degree_fine << ' ' << n_refinements << ' '
                << solver_control.last_step() << std::endl;


          DataOut<dim> data_out;

          DataOutBase::VtkFlags flags;
          flags.write_higher_order_cells = true;
          data_out.set_flags(flags);

          data_out.attach_dof_handler(dof_handlers[l]);
          data_out.add_data_vector(
            dst,
            "solution",
            DataOut_DoFData<dim, dim>::DataVectorType::type_dof_data);
          data_out.build_patches(*mapping_, fe_degree_fine);

          std::ofstream output("test." + std::to_string(dim) + "." +
                               std::to_string(counter) + "." +
                               std::to_string(l) + ".vtk");
          data_out.write_vtk(output);
        }

      counter++;
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  try
    {
      test<2>(3, 3);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
