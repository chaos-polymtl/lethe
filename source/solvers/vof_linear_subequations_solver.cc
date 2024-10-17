// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/vof_linear_subequations_solver.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_solver.h>

template <int dim, typename ScratchDataType>
void
VOFLinearSubequationsSolver<dim, ScratchDataType>::setup_dofs()
{
  // Get MPI communicator
  auto mpi_communicator = this->triangulation->get_communicator();

  // Distribute and renumber DoFs
  this->dof_handler.distribute_dofs(*this->fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  // Get locally owned and relevant DoFs
  this->locally_owned_dofs = this->dof_handler.locally_owned_dofs();
  this->locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(this->dof_handler);

  // Constraints
  this->constraints.clear();
  this->constraints.reinit(this->locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(this->dof_handler, this->constraints);
  this->constraints.close();

  // Sparsity pattern
  DynamicSparsityPattern dsp(this->locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(this->dof_handler,
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
  this->present_solution.reinit(this->locally_owned_dofs,
                                this->locally_relevant_dofs,
                                mpi_communicator);
  this->evaluation_point = this->present_solution;


  if (this->subequation_verbosity != Parameters::Verbosity::quiet)
    {
      std::string subequation_string =
        this->subequations->get_subequation_string(this->subequation_id);
      this->pcout << "   Number of " << subequation_string
                  << " degrees of freedom: " << this->dof_handler.n_dofs()
                  << std::endl;
    }

  // Provide DoFHandler and solutions to the subequations interface
  this->subequations->set_dof_handler(
    VOFSubequationsID::phase_gradient_projection, &this->dof_handler);
  this->subequations->set_solution(VOFSubequationsID::phase_gradient_projection,
                                   &this->present_solution);
}

template <int dim, typename ScratchDataType>
void
VOFLinearSubequationsSolver<dim, ScratchDataType>::assemble_system_matrix()
{
  this->system_matrix = 0;
  this->assembler =
    this->subequations->template assembler_cast<ScratchDataType>(
      this->subequation_id,
      this->simulation_parameters.multiphysics.vof_parameters);

  const DoFHandler<dim> *dof_handler_vof =
    this->multiphysics->get_dof_handler(PhysicsID::VOF);

  auto scratch_data_parent = this->subequations->scratch_data_cast(
    VOFSubequationsID::phase_gradient_projection,
    *this->fe,
    *this->cell_quadrature,
    *this->mapping,
    dof_handler_vof->get_fe());
  ScratchDataType *scratch_data =
    dynamic_cast<ScratchDataType *>(scratch_data_parent.get());

  WorkStream::run(
    this->dof_handler.begin_active(),
    this->dof_handler.end(),
    *this,
    &VOFLinearSubequationsSolver::assemble_local_system_matrix,
    &VOFLinearSubequationsSolver::copy_local_matrix_to_global_matrix,
    *scratch_data,
    StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                              this->cell_quadrature->size()));

  this->system_matrix.compress(VectorOperation::add);
}

template <int dim, typename ScratchDataType>
void
VOFLinearSubequationsSolver<dim, ScratchDataType>::
  copy_local_matrix_to_global_matrix(const StabilizedMethodsCopyData &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->constraints;
  constraints_used.distribute_local_to_global(copy_data.local_matrix,
                                              copy_data.local_dof_indices,
                                              this->system_matrix);
}

template <int dim, typename ScratchDataType>
void
VOFLinearSubequationsSolver<dim, ScratchDataType>::assemble_system_rhs()
{
  this->system_rhs = 0;

  const DoFHandler<dim> *dof_handler_vof =
    this->multiphysics->get_dof_handler(PhysicsID::VOF);

  auto scratch_data_parent = this->subequations->scratch_data_cast(
    VOFSubequationsID::phase_gradient_projection,
    *this->fe,
    *this->cell_quadrature,
    *this->mapping,
    dof_handler_vof->get_fe());

  ScratchDataType *scratch_data =
    dynamic_cast<ScratchDataType *>(scratch_data_parent.get());

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &VOFLinearSubequationsSolver::assemble_local_system_rhs,
                  &VOFLinearSubequationsSolver::copy_local_rhs_to_global_rhs,
                  *scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  this->system_rhs.compress(VectorOperation::add);
}

template <int dim, typename ScratchDataType>
void
VOFLinearSubequationsSolver<dim, ScratchDataType>::copy_local_rhs_to_global_rhs(
  const StabilizedMethodsCopyData &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->constraints;
  constraints_used.distribute_local_to_global(copy_data.local_rhs,
                                              copy_data.local_dof_indices,
                                              this->system_rhs);
}

template <int dim, typename ScratchDataType>
void
VOFLinearSubequationsSolver<dim, ScratchDataType>::
  solve_linear_system_and_update_solution(const bool &is_post_mesh_adaptation)
{
  auto mpi_communicator = this->triangulation->get_communicator();

  const AffineConstraints<double> &constraints_used = this->constraints;

  const bool verbose(
    this->subequation_verbosity != Parameters::Verbosity::quiet &&
    this->linear_solver_verbosity != Parameters::Verbosity::quiet &&
    !is_post_mesh_adaptation);

  if (verbose)
    {
      std::string subequation_string =
        this->subequations->get_subequation_string(this->subequation_id);

      this->pcout << "  -Solving " << subequation_string << ":" << std::endl;
    }

  // Set tolerance
  const double linear_solver_tolerance =
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .minimum_residual;

  // Solution vector
  GlobalVectorType completely_distributed_solution(this->locally_owned_dofs,
                                                   mpi_communicator);
  completely_distributed_solution = this->present_solution;

  if (verbose)
    {
      this->pcout << "    -Tolerance of iterative solver is: "
                  << linear_solver_tolerance << std::endl;
    }

  // ILU preconditioner
  const double ilu_fill =
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
    linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverCG solver(solver_control);

  solver.solve(this->system_matrix,
               completely_distributed_solution,
               this->system_rhs,
               *this->ilu_preconditioner);

  if (verbose)
    {
      this->pcout << "    -Iterative solver took: "
                  << solver_control.last_step()
                  << " steps to reach a residual norm of "
                  << solver_control.last_value() << std::endl;
    }

  // Update constraints vector
  constraints_used.distribute(completely_distributed_solution);

  // Update solution vectors
  this->present_solution = completely_distributed_solution;
  this->evaluation_point = this->present_solution;
}

template <int dim, typename ScratchDataType>
void
VOFLinearSubequationsSolver<dim, ScratchDataType>::solve(
  const bool &is_post_mesh_adaptation)
{
  assemble_system_matrix();
  assemble_system_rhs();
  solve_linear_system_and_update_solution(is_post_mesh_adaptation);
}

template class VOFLinearSubequationsSolver<
  2,
  VOFPhaseGradientProjectionScratchData<2>>;
template class VOFLinearSubequationsSolver<
  3,
  VOFPhaseGradientProjectionScratchData<3>>;
