// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/vof_linear_subequations_solver.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_solver.h>


template <int dim>
void
VOFLinearSubequationsSolver<dim>::setup_dofs()
{
  // Get MPI communicator
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  // Distribute and renumber DoFs
  this->dof_handler->distribute_dofs(*this->fe);
  DoFRenumbering::Cuthill_McKee(*this->dof_handler);

  // Get locally owned and relevant DoFs
  this->locally_owned_dofs = this->dof_handler->locally_owned_dofs();
  this->locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(*this->dof_handler);

  // Constraints
  this->constraints.clear();
  this->constraints.reinit(this->locally_owned_dofs,
                           this->locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(*this->dof_handler,
                                          this->constraints);

  // Add periodic boundary conditions for the linear subequation
  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions_vof.type)
    {
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            *this->dof_handler,
            id,
            this->simulation_parameters.boundary_conditions_vof
              .periodic_neighbor_id.at(id),
            this->simulation_parameters.boundary_conditions_vof
              .periodic_direction.at(id),
            this->constraints);
        }
    }

  this->constraints.close();

  // Sparsity pattern
  DynamicSparsityPattern dsp(this->locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(*this->dof_handler,
                                  dsp,
                                  this->constraints,
                                  false);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             this->locally_owned_dofs,
                                             mpi_communicator,
                                             this->locally_relevant_dofs);

  // Reinitialize system matrix and right-hand side (rhs)
  this->system_matrix.reinit(this->locally_owned_dofs,
                             this->locally_owned_dofs,
                             dsp,
                             mpi_communicator);
  this->system_rhs.reinit(this->locally_owned_dofs, mpi_communicator);

  // Reinitialize solution vectors
  this->present_solution->reinit(this->locally_owned_dofs,
                                 this->locally_relevant_dofs,
                                 mpi_communicator);
  this->evaluation_point = *this->present_solution;


  if (this->subequation_verbosity != Parameters::Verbosity::quiet)
    {
      std::string subequation_string =
        this->subequations_interface.get_subequation_string(
          this->subequation_id);
      this->pcout << "   Number of " << subequation_string
                  << " degrees of freedom: " << this->dof_handler->n_dofs()
                  << std::endl;
    }

  // Provide DoFHandler and solutions to the subequations interface
  this->subequations_interface.set_dof_handler(this->subequation_id,
                                               this->dof_handler);
  this->subequations_interface.set_solution(this->subequation_id,
                                            this->present_solution);
}


template <int dim>
void
VOFLinearSubequationsSolver<dim>::solve_linear_system_and_update_solution()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  const AffineConstraints<double> &constraints_used = this->constraints;

  const bool verbose(
    this->subequation_verbosity != Parameters::Verbosity::quiet &&
    this->linear_solver_verbosity != Parameters::Verbosity::quiet);

  if (verbose)
    {
      std::string subequation_string =
        this->subequations_interface.get_subequation_string(
          this->subequation_id);

      this->pcout << "  -Solving " << subequation_string << ":" << std::endl;
    }

  // Set tolerance
  const double normalize_metric =
    simulation_parameters.non_linear_solver.at(PhysicsID::VOF)
        .normalize_residual_by_volume ?
      GridTools::volume(*this->triangulation, *this->mapping) :
      1.0;
  const double linear_solver_tolerance =
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .minimum_residual /
    normalize_metric;

  const double non_normalized_linear_solver_tolerance =
    linear_solver_tolerance * normalize_metric;

  // Solution vector
  GlobalVectorType completely_distributed_solution(this->locally_owned_dofs,
                                                   mpi_communicator);

  if (verbose)
    {
      this->pcout << "    -Tolerance of iterative solver is: "
                  << linear_solver_tolerance << std::endl;
    }

  // ILU preconditioner
  const unsigned int ilu_fill =
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  this->ilu_preconditioner =
    std::make_shared<TrilinosWrappers::PreconditionILU>();

  this->ilu_preconditioner->initialize(this->system_matrix,
                                       preconditionerOptions);

  // CG solver
  SolverControl solver_control(
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF).max_iterations,
    non_normalized_linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverCG solver(solver_control);

  solver.solve(this->system_matrix,
               completely_distributed_solution,
               this->system_rhs,
               *this->ilu_preconditioner);

  if (verbose)
    {
      this->pcout << "    -Iterative solver took " << solver_control.last_step()
                  << " steps to reach a residual norm of "
                  << solver_control.last_value() / normalize_metric
                  << std::endl;
    }

  // Update constraints vector
  constraints_used.distribute(completely_distributed_solution);

  // Update solution vectors
  *this->present_solution = completely_distributed_solution;
  this->evaluation_point  = *this->present_solution;
}


template <int dim>
void
VOFLinearSubequationsSolver<dim>::solve()
{
  check_dependencies_validity();
  assemble_system_matrix_and_rhs();
  solve_linear_system_and_update_solution();
  this->subequations_interface.set_solution_valid(this->subequation_id);
}


template class VOFLinearSubequationsSolver<2>;
template class VOFLinearSubequationsSolver<3>;
