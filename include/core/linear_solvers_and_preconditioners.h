/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
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
 * ---------------------------------------------------------------------
 */

#ifndef lethe_linear_solvers_and_preconditioners_h
#define lethe_linear_solvers_and_preconditioners_h

#include <solvers/mf_navier_stokes_operators.h>
#include <solvers/simulation_parameters.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/matrix_free/operators.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>


template <int dim>
void
setup_ilu_preconditioner(TrilinosWrappers::PreconditionILU &ilu_preconditioner,
                         SimulationParameters<dim> &        parameters,
                         const TrilinosWrappers::SparseMatrix &system_matrix,
                         const int &current_preconditioner_fill_level)
{
  const double ilu_atol = parameters.linear_solver.ilu_precond_atol;
  const double ilu_rtol = parameters.linear_solver.ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    current_preconditioner_fill_level, ilu_atol, ilu_rtol, 0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner.initialize(system_matrix, preconditionerOptions);
}

template <int dim>
void
setup_amg_preconditioner(TrilinosWrappers::PreconditionAMG &amg_preconditioner,
                         SimulationParameters<dim> &        parameters,
                         const TrilinosWrappers::SparseMatrix &system_matrix,
                         const int &current_preconditioner_fill_level,
                         const DoFHandler<dim> &dof_handler)
{
  // Constant modes for velocity
  std::vector<std::vector<bool>> constant_modes;

  // Constant modes include pressure since everything is in the same matrix
  ComponentMask components(dim + 1, true);
  DoFTools::extract_constant_modes(dof_handler, components, constant_modes);

  TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
  amg_data.constant_modes = constant_modes;

  const bool elliptic              = false;
  bool       higher_order_elements = false;
  if (parameters.fem_parameters.velocity_order > 1)
    higher_order_elements = true;
  const unsigned int n_cycles = parameters.linear_solver.amg_n_cycles;
  const bool         w_cycle  = parameters.linear_solver.amg_w_cycles;
  const double       aggregation_threshold =
    parameters.linear_solver.amg_aggregation_threshold;
  const unsigned int smoother_sweeps =
    parameters.linear_solver.amg_smoother_sweeps;
  const unsigned int smoother_overlap =
    parameters.linear_solver.amg_smoother_overlap;
  const bool                                        output_details = false;
  const char *                                      smoother_type  = "ILU";
  const char *                                      coarse_type    = "ILU";
  TrilinosWrappers::PreconditionAMG::AdditionalData preconditionerOptions(
    elliptic,
    higher_order_elements,
    n_cycles,
    w_cycle,
    aggregation_threshold,
    constant_modes,
    smoother_sweeps,
    smoother_overlap,
    output_details,
    smoother_type,
    coarse_type);

  Teuchos::ParameterList              parameter_ml;
  std::unique_ptr<Epetra_MultiVector> distributed_constant_modes;
  preconditionerOptions.set_parameters(parameter_ml,
                                       distributed_constant_modes,
                                       system_matrix);
  const double ilu_fill = current_preconditioner_fill_level;
  const double ilu_atol = parameters.linear_solver.amg_precond_ilu_atol;
  const double ilu_rtol = parameters.linear_solver.amg_precond_ilu_rtol;
  parameter_ml.set("smoother: ifpack level-of-fill", ilu_fill);
  parameter_ml.set("smoother: ifpack absolute threshold", ilu_atol);
  parameter_ml.set("smoother: ifpack relative threshold", ilu_rtol);

  parameter_ml.set("coarse: ifpack level-of-fill", ilu_fill);
  parameter_ml.set("coarse: ifpack absolute threshold", ilu_atol);
  parameter_ml.set("coarse: ifpack relative threshold", ilu_rtol);
  amg_preconditioner = std::make_shared<TrilinosWrappers::PreconditionAMG>();
  amg_preconditioner.initialize(system_matrix, parameter_ml);
}

// Multigrid global coarsening algorithm and parameters
struct MultigridParameters
{
  struct
  {
    std::string  type            = "gmres_with_amg";
    unsigned int maxiter         = 2000;
    double       abstol          = 1e-14;
    double       reltol          = 1e-4;
    unsigned int smoother_sweeps = 1;
    unsigned int n_cycles        = 1;
    std::string  smoother_type   = "ILU";
  } coarse_solver;

  struct
  {
    std::string  type         = "relaxation";
    unsigned int n_iterations = 10;
    double       relaxation   = 0.5;
  } smoother;
};

// setup_ls_multigrid_preconditioner(
//   this->ls_multigrid_preconditioner,
//   this->simulation_parameters,
//   *this->system_operator,
//   *this->mapping,
//   this->dof_handler,
//   *this->cell_quadrature,
//   this->present_solution,
//   this->forcing_function,
//   this->simulation_parameters.physical_properties_manager
//     .get_viscosity_scale());

template <int dim,
          typename VectorType,
          typename OperatorType,
          typename LSTransferType>
void
setup_ls_multigrid_preconditioner(
  std::shared_ptr<PreconditionMG<dim, VectorType, LSTransferType>>
                                   preconditioner,
  const SimulationParameters<dim> &simulation_parameters,
  const OperatorType &             system_matrix,
  const Mapping<dim> &             mapping,
  const DoFHandler<dim> &          dof_handler,
  const Quadrature<dim> &          quadrature,
  const VectorType &               solution,
  const Function<dim> *            forcing_function,
  const double                     viscosity)
{
  // Only needed to deduce the OperatorType
  (void)system_matrix;
  (void)simulation_parameters;

  using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
  using SmootherType =
    PreconditionRelaxation<OperatorType, SmootherPreconditionerType>;
  using PreconditionerType =
    PreconditionMG<dim, VectorType, MGTransferMatrixFree<dim, double>>;

  MGLevelObject<std::shared_ptr<NavierStokesSUPGPSPGOperator<dim, double>>>
                                                     operators;
  std::shared_ptr<MGTransferMatrixFree<dim, double>> mg_transfer =
    std::make_shared<MGTransferMatrixFree<dim, double>>();
  MGLevelObject<VectorType>                    mg_solution;
  MGLevelObject<std::shared_ptr<OperatorType>> mg_interface_in;
  MGLevelObject<std::shared_ptr<OperatorType>> mg_interface_out;
  MGLevelObject<AffineConstraints<double>>     level_constraints;
  MGConstrainedDoFs                            mg_constrained_dofs;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_operators;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_interface_in;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_interface_out;

  const unsigned int n_h_levels =
    dof_handler.get_triangulation().n_global_levels();

  unsigned int minlevel = 0;
  unsigned int maxlevel = n_h_levels - 1;

  operators.resize(0, n_h_levels - 1);
  mg_solution.resize(0, n_h_levels - 1);
  level_constraints.resize(0, n_h_levels - 1);
  ls_mg_interface_in.resize(0, n_h_levels - 1);
  ls_mg_interface_out.resize(0, n_h_levels - 1);
  ls_mg_operators.resize(0, n_h_levels - 1);

  std::set<types::boundary_id> dirichlet_boundary_ids = {0, 1, 2, 3, 4, 5};

  mg_constrained_dofs.initialize(dof_handler);
  mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                     dirichlet_boundary_ids);

  std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>> partitioners(
    dof_handler.get_triangulation().n_global_levels());

  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    {
      const IndexSet relevant_dofs =
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);

      level_constraints[level].reinit(relevant_dofs);
      level_constraints[level].add_lines(
        mg_constrained_dofs.get_boundary_indices(level));
      level_constraints[level].close();

      operators[level] =
        std::make_unique<NavierStokesSUPGPSPGOperator<dim, double>>();

      operators[level]->reinit(mapping,
                               dof_handler,
                               level_constraints[level],
                               quadrature,
                               forcing_function,
                               viscosity,
                               level);

      operators[level]->initialize_dof_vector(mg_solution[level]);

      ls_mg_operators[level].initialize(*operators[level]);
      ls_mg_interface_in[level].initialize(*operators[level]);
      ls_mg_interface_out[level].initialize(*operators[level]);

      partitioners[level] = operators[level]->get_vector_partitioner();
    }

  mg_transfer->initialize_constraints(mg_constrained_dofs);
  mg_transfer->build(dof_handler, partitioners);
  mg_transfer->interpolate_to_mg(dof_handler, mg_solution, solution);

  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    {
      mg_solution[level].update_ghost_values();
      operators[level]->evaluate_non_linear_term(mg_solution[level]);
    }

  ConditionalOStream pcout(std::cout,
                           (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                            0));

  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    pcout << "   MG Level " << level << ": " << dof_handler.n_dofs(level)
          << " DoFs, " << dof_handler.get_triangulation().n_cells(level)
          << " cells" << std::endl;

  mg::Matrix<VectorType> mg_matrix(ls_mg_operators);

  MGLevelObject<typename SmootherType::AdditionalData> smoother_data(minlevel,
                                                                     maxlevel);

  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    {
      smoother_data[level].preconditioner =
        std::make_shared<DiagonalMatrix<VectorType>>();
      operators[level]->compute_inverse_diagonal(
        smoother_data[level].preconditioner->get_vector());
      smoother_data[level].n_iterations = 10;
      smoother_data[level].relaxation   = 0.5;
    }

  MGSmootherPrecondition<OperatorType, SmootherType, VectorType> mg_smoother;
  mg_smoother.initialize(operators, smoother_data);

  ReductionControl coarse_grid_solver_control(2000, 1e-14, 1e-4, false, false);
  SolverGMRES<VectorType> coarse_grid_solver(coarse_grid_solver_control);

  std::shared_ptr<MGCoarseGridBase<VectorType>> mg_coarse;

  TrilinosWrappers::PreconditionAMG                 precondition_amg;
  TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
  amg_data.smoother_sweeps = 1;
  amg_data.n_cycles        = 1;
  amg_data.smoother_type   = "ILU";

  precondition_amg.initialize(operators[minlevel]->get_system_matrix(),
                              amg_data);

  mg_coarse =
    std::make_unique<MGCoarseGridIterativeSolver<VectorType,
                                                 SolverGMRES<VectorType>,
                                                 OperatorType,
                                                 decltype(precondition_amg)>>(
      coarse_grid_solver, *operators[minlevel], precondition_amg);

  mg::Matrix<VectorType> mg_interface_matrix_in(ls_mg_interface_in);
  mg::Matrix<VectorType> mg_interface_matrix_out(ls_mg_interface_out);

  std::shared_ptr<Multigrid<VectorType>> mg =
    std::make_shared<Multigrid<VectorType>>(
      mg_matrix, *mg_coarse, *mg_transfer, mg_smoother, mg_smoother);

  if (dof_handler.get_triangulation().has_hanging_nodes())
    mg->set_edge_matrices(mg_interface_matrix_in, mg_interface_matrix_out);


  preconditioner =
    std::make_shared<PreconditionMG<dim, VectorType, LSTransferType>>(
      dof_handler, *mg, *mg_transfer);
}

#endif
