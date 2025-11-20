// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/solutions_output.h>

#include <solvers/vof_algebraic_interface_reinitialization.h>
#include <solvers/vof_assemblers.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::setup_dofs()
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
  define_non_zero_constraints();
  define_zero_constraints();

  // Sparsity pattern
  DynamicSparsityPattern dsp(this->locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(*this->dof_handler,
                                  dsp,
                                  this->nonzero_constraints,
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
  this->previous_solution.reinit(this->locally_owned_dofs,
                                 this->locally_relevant_dofs,
                                 mpi_communicator);
  this->newton_update.reinit(this->locally_owned_dofs, mpi_communicator);
  this->local_evaluation_point.reinit(this->locally_owned_dofs,
                                      mpi_communicator);
  this->previous_local_evaluation_point.reinit(this->locally_owned_dofs,
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
VOFAlgebraicInterfaceReinitialization<dim>::define_zero_constraints()
{
  this->zero_constraints.clear();
  this->zero_constraints.reinit(this->locally_owned_dofs,
                                this->locally_relevant_dofs);

  DoFTools::make_hanging_node_constraints(*this->dof_handler,
                                          this->zero_constraints);

  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            *this->dof_handler,
            id,
            this->simulation_parameters.boundary_conditions.periodic_neighbor_id
              .at(id),
            this->simulation_parameters.boundary_conditions.periodic_direction
              .at(id),
            this->zero_constraints);
        }
    }
  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions_vof.type)
    {
      if (type == BoundaryConditions::BoundaryType::vof_dirichlet)
        {
          VectorTools::interpolate_boundary_values(
            *this->dof_handler,
            id,
            Functions::ZeroFunction<dim>(),
            this->zero_constraints);
        }
    }
  this->zero_constraints.close();
}


template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::define_non_zero_constraints()
{
  this->nonzero_constraints.clear();
  this->nonzero_constraints.reinit(this->locally_owned_dofs,
                                   this->locally_relevant_dofs);

  DoFTools::make_hanging_node_constraints(*this->dof_handler,
                                          this->nonzero_constraints);

  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            *this->dof_handler,
            id,
            this->simulation_parameters.boundary_conditions.periodic_neighbor_id
              .at(id),
            this->simulation_parameters.boundary_conditions.periodic_direction
              .at(id),
            this->nonzero_constraints);
        }
    }
  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions_vof.type)
    {
      if (type == BoundaryConditions::BoundaryType::vof_dirichlet)
        {
          VectorTools::interpolate_boundary_values(
            *this->dof_handler,
            id,
            *this->simulation_parameters.boundary_conditions_vof.phase_fraction
               .at(id),
            this->nonzero_constraints);
        }
    }
  this->nonzero_constraints.close();
}


template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::set_initial_conditions()
{
  // Get VOF DoFHandler
  const DoFHandler<dim> &dof_handler_vof =
    this->subequations_interface.get_vof_dof_handler();

  // Interpolate VOF solution to algebraic interface reinitialization
  FETools::interpolate(dof_handler_vof,
                       this->subequations_interface.get_vof_solution(),
                       *this->dof_handler,
                       this->nonzero_constraints,
                       this->newton_update);

  *this->present_solution = this->newton_update;
  this->previous_solution =
    *this->present_solution; // We only have 1 previous solution (bdf1)
  this->previous_local_evaluation_point =
    this->previous_solution; // For steady-state criterion evaluation

  // Initial condition
  if (this->simulation_parameters.multiphysics.vof_parameters
        .regularization_method.algebraic_interface_reinitialization
        .output_reinitialization_steps)
    {
      this->pvdhandler.times_and_names.clear();
      write_output_results(0);
    }
}


template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::assemble_system_matrix()
{
  // Reinitialize system matrix
  this->system_matrix = 0;

  // Get projected VOF phase gradient and curvature DoFHandler
  const DoFHandler<dim> &dof_handler_vof_phase_fraction_gradient =
    this->subequations_interface.get_dof_handler(
      VOFSubequationsID::phase_gradient_projection);
  const DoFHandler<dim> &dof_handler_curvature =
    this->subequations_interface.get_dof_handler(
      VOFSubequationsID::curvature_projection);

  // Initialize FEValues for algebraic interface reinitialization, VOF
  // phase fraction gradient and curvature projections
  FEValues<dim> fe_values_algebraic_reinitialization(*this->mapping,
                                                     *this->fe,
                                                     *this->cell_quadrature,
                                                     update_values |
                                                       update_gradients |
                                                       update_JxW_values);
  FEValues<dim> fe_values_phase_gradient_projection(
    *this->mapping,
    dof_handler_vof_phase_fraction_gradient.get_fe(),
    *this->cell_quadrature,
    update_values);
  FEValues<dim> fe_values_curvature_projection(*this->mapping,
                                               dof_handler_curvature.get_fe(),
                                               *this->cell_quadrature,
                                               update_values);

  // Initialize size of arrays
  const unsigned int n_q_points =
    fe_values_algebraic_reinitialization.get_quadrature().size();
  const unsigned int n_dofs_per_cell =
    fe_values_algebraic_reinitialization.get_fe().n_dofs_per_cell();

  // Initialize local matrix
  FullMatrix<double> local_matrix(n_dofs_per_cell, n_dofs_per_cell);

  // Initialize local dof indices array
  std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);

  // Extractor for phase fraction gradient vector
  FEValuesExtractors::Vector phase_fraction_gradients(0);

  // Initialize phase fraction, projected phase gradient and projected curvature
  // solution arrays
  std::vector<double>         present_phase_fraction_values(n_q_points);
  std::vector<Tensor<1, dim>> present_phase_gradient_projection_values(
    n_q_points);
  std::vector<double>         present_curvature_values(n_q_points);
  std::vector<Tensor<1, dim>> present_reinitialized_phase_gradient_values(
    n_q_points);

  // Initialize shape function arrays
  std::vector<double>         phi(n_dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_phi(n_dofs_per_cell);

  // Compute diffusivity coefficient with the smallest cell size
  const double h_min = identify_minimum_cell_size(
    *this->mapping,
    this->subequations_interface.get_dof_handler(
      VOFSubequationsID::algebraic_interface_reinitialization),
    *this->cell_quadrature,
    this->triangulation->get_mpi_communicator());
  const double diffusivity_coefficient = compute_diffusivity(h_min);

  // BDF coefficients for pseudo time-stepping
  const Parameters::SimulationControl::TimeSteppingMethod method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;
  const Vector<double> bdf_coefficient_vector =
    calculate_bdf_coefficients(method, this->time_step_vector);

  // Assemble system matrix
  for (const auto &cell : this->dof_handler->active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Reinitialize local matrix
          local_matrix = 0;

          // Get phase gradient projection and curvature cell iterators
          typename DoFHandler<dim>::active_cell_iterator
            vof_phase_gradient_projection_cell(
              &(*this->triangulation),
              cell->level(),
              cell->index(),
              &dof_handler_vof_phase_fraction_gradient);
          typename DoFHandler<dim>::active_cell_iterator
            curvature_projection_cell(&(*this->triangulation),
                                      cell->level(),
                                      cell->index(),
                                      &dof_handler_curvature);

          // Reinitialize FEValues with corresponding cell
          fe_values_algebraic_reinitialization.reinit(cell);
          fe_values_phase_gradient_projection.reinit(
            vof_phase_gradient_projection_cell);
          fe_values_curvature_projection.reinit(curvature_projection_cell);

          // Get vector of JxW
          std::vector<double> JxW_vec =
            fe_values_algebraic_reinitialization.get_JxW_values();

          // Get get present phase fraction, projected VOF phase gradient and
          // curvature values
          fe_values_algebraic_reinitialization.get_function_values(
            this->evaluation_point, present_phase_fraction_values);
          fe_values_phase_gradient_projection[phase_fraction_gradients]
            .get_function_values(
              this->subequations_interface.get_solution(
                VOFSubequationsID::phase_gradient_projection),
              present_phase_gradient_projection_values);
          fe_values_curvature_projection.get_function_values(
            this->subequations_interface.get_solution(
              VOFSubequationsID::curvature_projection),
            present_curvature_values);
          fe_values_algebraic_reinitialization.get_function_gradients(
            this->evaluation_point,
            present_reinitialized_phase_gradient_values);

          // Set tolerance to avoid division by zero
          const double tolerance = 1e-12;

          // Loop over quadrature points
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              // Shape functions
              for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
                {
                  phi[k] =
                    fe_values_algebraic_reinitialization.shape_value(k, q);
                  grad_phi[k] =
                    fe_values_algebraic_reinitialization.shape_grad(k, q);
                }

              // Extract phase fraction value, projected phase fraction gradient
              // and curvature
              const double phase_fraction = present_phase_fraction_values[q];
              const Tensor<1, dim> projected_vof_phase_gradient =
                present_phase_gradient_projection_values[q];
              const double         curvature = present_curvature_values[q];
              const Tensor<1, dim> reinitialized_phase_gradient =
                present_reinitialized_phase_gradient_values[q];

              // Compute normal vector with projected phase gradient
              const double projected_vof_phase_gradient_norm =
                projected_vof_phase_gradient.norm() + tolerance;
              const Tensor<1, dim> interface_normal =
                (projected_vof_phase_gradient /
                 projected_vof_phase_gradient_norm);

              // Local matrix assembly
              for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                    {
                      local_matrix(i, j) +=
                        (
                          // Time-stepping term
                          phi[i] * phi[j] * bdf_coefficient_vector[0] +
                          // Compressive term
                          phi[i] *
                            (scalar_product(grad_phi[j], interface_normal) -
                             2 * phase_fraction *
                               scalar_product(grad_phi[j], interface_normal) -
                             2 * phi[j] *
                               scalar_product(reinitialized_phase_gradient,
                                              interface_normal) -
                             (phi[j] - 2 * phase_fraction * phi[j]) *
                               curvature) +
                          // Diffusive term
                          scalar_product(grad_phi[i],
                                         (diffusivity_coefficient *
                                          scalar_product(grad_phi[j],
                                                         interface_normal) *
                                          interface_normal))) *
                        JxW_vec[q];
                    }
                }
            }
          // Distribute the local contributions to the global system
          cell->get_dof_indices(local_dof_indices);
          this->zero_constraints.distribute_local_to_global(
            local_matrix, local_dof_indices, this->system_matrix);
        }
    }
  this->system_matrix.compress(VectorOperation::add);
}


template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::assemble_system_rhs()
{
  // Reinitialize system right-hand side (rhs)
  this->system_rhs = 0;

  // Get projected phase gradient, curvature and reinitialized phase gradient
  // DoFHandler
  const DoFHandler<dim> &dof_handler_vof_phase_fraction_gradient =
    this->subequations_interface.get_dof_handler(
      VOFSubequationsID::phase_gradient_projection);
  const DoFHandler<dim> &dof_handler_curvature =
    this->subequations_interface.get_dof_handler(
      VOFSubequationsID::curvature_projection);

  // Initialize FEValues for algebraic interface reinitialization, VOF phase
  // fraction gradient and curvature projections
  FEValues<dim> fe_values_algebraic_reinitialization(*this->mapping,
                                                     *this->fe,
                                                     *this->cell_quadrature,
                                                     update_values |
                                                       update_gradients |
                                                       update_JxW_values);
  FEValues<dim> fe_values_vof_phase_gradient_projection(
    *this->mapping,
    dof_handler_vof_phase_fraction_gradient.get_fe(),
    *this->cell_quadrature,
    update_values);
  FEValues<dim> fe_values_curvature_projection(*this->mapping,
                                               dof_handler_curvature.get_fe(),
                                               *this->cell_quadrature,
                                               update_values);

  // Initialize size of arrays
  const unsigned int n_q_points =
    fe_values_algebraic_reinitialization.get_quadrature().size();
  const unsigned int n_dofs_per_cell =
    fe_values_algebraic_reinitialization.get_fe().n_dofs_per_cell();

  //  Initialize local rhs
  Vector<double> local_rhs(n_dofs_per_cell);

  // Initialize local dof indices array
  std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);

  // Extractor for phase fraction gradient vector
  FEValuesExtractors::Vector phase_fraction_gradients(0);

  // Initialize phase fraction and projected phase fraction gradient and
  // curvature solution arrays
  std::vector<double>         present_phase_fraction_values(n_q_points);
  std::vector<double>         previous_phase_fraction_values(n_q_points);
  std::vector<Tensor<1, dim>> present_vof_phase_gradient_projection_values(
    n_q_points);
  std::vector<double>         present_curvature_values(n_q_points);
  std::vector<Tensor<1, dim>> present_reinitialized_phase_gradient_values(
    n_q_points);

  // Initialize shape function arrays
  std::vector<double>         phi(n_dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_phi(n_dofs_per_cell);

  // Compute diffusivity coefficient with the smallest cell size
  const double h_min = identify_minimum_cell_size(
    *this->mapping,
    this->subequations_interface.get_dof_handler(
      VOFSubequationsID::algebraic_interface_reinitialization),
    *this->cell_quadrature,
    this->triangulation->get_mpi_communicator());
  const double diffusivity_coefficient = compute_diffusivity(h_min);

  // BDF coefficients for pseudo time-stepping
  const Parameters::SimulationControl::TimeSteppingMethod method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;
  const Vector<double> bdf_coefficient_vector =
    calculate_bdf_coefficients(method, this->time_step_vector);

  for (const auto &cell : this->dof_handler->active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Reinitialize local rhs
          local_rhs = 0;

          // Get phase gradient projection cell iterator
          typename DoFHandler<dim>::active_cell_iterator
            vof_phase_gradient_projection_cell(
              &(*this->triangulation),
              cell->level(),
              cell->index(),
              &dof_handler_vof_phase_fraction_gradient);
          typename DoFHandler<dim>::active_cell_iterator
            curvature_projection_cell(&(*this->triangulation),
                                      cell->level(),
                                      cell->index(),
                                      &dof_handler_curvature);

          // Reinitialize FEValues with corresponding cells
          fe_values_algebraic_reinitialization.reinit(cell);
          fe_values_vof_phase_gradient_projection.reinit(
            vof_phase_gradient_projection_cell);
          fe_values_curvature_projection.reinit(curvature_projection_cell);

          // Get vector of Jacobi determinant times the quadrature weights
          std::vector<double> JxW_vec =
            fe_values_algebraic_reinitialization.get_JxW_values();

          // Get VOF phase fraction and projected phase gradient values
          fe_values_algebraic_reinitialization.get_function_values(
            this->evaluation_point, present_phase_fraction_values);
          fe_values_algebraic_reinitialization.get_function_values(
            this->previous_solution, previous_phase_fraction_values);
          fe_values_vof_phase_gradient_projection[phase_fraction_gradients]
            .get_function_values(
              this->subequations_interface.get_solution(
                VOFSubequationsID::phase_gradient_projection),
              present_vof_phase_gradient_projection_values);
          fe_values_curvature_projection.get_function_values(
            this->subequations_interface.get_solution(
              VOFSubequationsID::curvature_projection),
            present_curvature_values);
          fe_values_algebraic_reinitialization.get_function_gradients(
            this->evaluation_point,
            present_reinitialized_phase_gradient_values);

          // Set tolerance to avoid division by zero
          const double tolerance = 1e-12;

          // Loop over quadrature points
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              // Shape functions
              for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
                {
                  phi[k] =
                    fe_values_algebraic_reinitialization.shape_value(k, q);
                  grad_phi[k] =
                    fe_values_algebraic_reinitialization.shape_grad(k, q);
                }

              // Extract present and previous phase fraction values and
              // projected phase gradient values
              std::vector<double> phase_fraction_values(2);
              phase_fraction_values[0]    = present_phase_fraction_values[q];
              const double phase_fraction = phase_fraction_values[0];
              phase_fraction_values[1]    = previous_phase_fraction_values[q];
              const Tensor<1, dim> projected_vof_phase_gradient =
                present_vof_phase_gradient_projection_values[q];
              const double         curvature = present_curvature_values[q];
              const Tensor<1, dim> reinitialized_phase_gradient =
                present_reinitialized_phase_gradient_values[q];

              // Compute normal vector with projected phase gradient
              const double projected_vof_phase_gradient_norm =
                projected_vof_phase_gradient.norm() + tolerance;
              const Tensor<1, dim> interface_normal =
                (projected_vof_phase_gradient /
                 projected_vof_phase_gradient_norm);

              // Local rhs assembly
              for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                {
                  // Time-stepping term (bdf1)
                  for (unsigned int p = 0; p < 2; ++p)
                    {
                      local_rhs(i) -= bdf_coefficient_vector[p] *
                                      (phase_fraction_values[p] * phi[i]) *
                                      JxW_vec[q];
                    }

                  local_rhs(i) -=
                    (
                      // Compressive term
                      phi[i] *
                        (scalar_product(reinitialized_phase_gradient,
                                        interface_normal) -
                         2 * phase_fraction *
                           scalar_product(reinitialized_phase_gradient,
                                          interface_normal) -
                         phase_fraction * (1 - phase_fraction) * curvature) +
                      // Diffusive term
                      (scalar_product(
                        grad_phi[i],
                        (diffusivity_coefficient *
                         scalar_product(reinitialized_phase_gradient,
                                        interface_normal) *
                         interface_normal)))) *
                    JxW_vec[q];
                }
            }
          // Distribute the local contributions to the global system
          cell->get_dof_indices(local_dof_indices);
          this->zero_constraints.distribute_local_to_global(local_rhs,
                                                            local_dof_indices,
                                                            this->system_rhs);
        }
    }
  this->system_rhs.compress(VectorOperation::add);
}


template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::solve_linear_system()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  const AffineConstraints<double> &constraints_used = this->zero_constraints;

  const bool verbose(
    this->simulation_parameters.vof_subequations_linear_solvers
      .at(VOFSubequationsID::algebraic_interface_reinitialization)
      .verbosity != Parameters::Verbosity::quiet);

  // Get residual conditions
  const double absolute_residual =
    this->simulation_parameters.vof_subequations_linear_solvers
      .at(VOFSubequationsID::algebraic_interface_reinitialization)
      .minimum_residual;
  const double relative_residual =
    this->simulation_parameters.vof_subequations_linear_solvers
      .at(VOFSubequationsID::algebraic_interface_reinitialization)
      .relative_residual;

  // Set linear solver tolerance
  const double rescale_metric   = this->get_residual_rescale_metric();
  const double current_residual = this->system_rhs.l2_norm() / rescale_metric;
  const double linear_solver_tolerance =
    std::max(relative_residual * current_residual, absolute_residual);
  const double non_rescaled_linear_solver_tolerance =
    linear_solver_tolerance * rescale_metric;

  if (verbose)
    {
      this->pcout << "  -Tolerance of iterative solver is: "
                  << linear_solver_tolerance << std::endl;
    }

  // ILU preconditioner
  const unsigned int ilu_fill =
    this->simulation_parameters.vof_subequations_linear_solvers
      .at(VOFSubequationsID::algebraic_interface_reinitialization)
      .ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.vof_subequations_linear_solvers
      .at(VOFSubequationsID::algebraic_interface_reinitialization)
      .ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.vof_subequations_linear_solvers
      .at(VOFSubequationsID::algebraic_interface_reinitialization)
      .ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  TrilinosWrappers::PreconditionILU ilu_preconditioner;

  ilu_preconditioner.initialize(this->system_matrix, preconditionerOptions);

  GlobalVectorType completely_distributed_solution(this->locally_owned_dofs,
                                                   mpi_communicator);

  SolverControl solver_control(
    this->simulation_parameters.vof_subequations_linear_solvers
      .at(VOFSubequationsID::algebraic_interface_reinitialization)
      .max_iterations,
    non_rescaled_linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    false,
    this->simulation_parameters.vof_subequations_linear_solvers
      .at(VOFSubequationsID::algebraic_interface_reinitialization)
      .max_krylov_vectors);

  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);

  solver.solve(this->system_matrix,
               completely_distributed_solution,
               this->system_rhs,
               ilu_preconditioner);

  if (verbose)
    {
      this->pcout << "  -Iterative solver took: " << solver_control.last_step()
                  << " steps to reach a residual norm of "
                  << solver_control.last_value() / rescale_metric << std::endl;
    }

  // Update constraints and newton vectors
  constraints_used.distribute(completely_distributed_solution);
  this->newton_update = completely_distributed_solution;
}


template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::solve()
{
  // Check if all required solutions are valid
  check_dependencies_validity();

  // Compute time-step
  this->current_time_step = compute_time_step();

  double current_time_step_inv = 1. / this->current_time_step;
  this->time_step_vector[0]    = this->current_time_step;

  // Set initial conditions from the VOF solution (step 0)
  set_initial_conditions();

  // Reinitialization steps counter
  unsigned int step = 1;

  // Define verbosity
  const bool verbose(
    this->subequation_verbosity != Parameters::Verbosity::quiet ||
    this->simulation_parameters.vof_subequations_non_linear_solvers
        .at(VOFSubequationsID::algebraic_interface_reinitialization)
        .verbosity != Parameters::Verbosity::quiet);
  std::string subequation_string =
    this->subequations_interface.get_subequation_string(this->subequation_id);

  if (verbose)
    {
      announce_string(this->pcout, "VOF algebraic interface reinitialization");

      this->pcout << "Current algebraic reinitialization time-step: "
                  << this->current_time_step << std::endl;
    }

  // For debugging purposes
  if (this->simulation_parameters.multiphysics.vof_parameters
        .regularization_method.algebraic_interface_reinitialization
        .output_reinitialization_steps)
    write_output_results(step);

  // Iterate until a stop criterion is met
  while (continue_iterating(current_time_step_inv, step))
    {
      if (verbose)
        {
          this->pcout << "-Solving " << subequation_string << ", step " << step
                      << std::endl;
        }

      // Update previous solution
      this->previous_solution               = *this->present_solution;
      this->previous_local_evaluation_point = this->local_evaluation_point;

      // Update non-zero constraints
      define_non_zero_constraints();

      // Solve non-linear equation
      this->solve_non_linear_system();

      // For debugging purposes
      if (this->simulation_parameters.multiphysics.vof_parameters
            .regularization_method.algebraic_interface_reinitialization
            .output_reinitialization_steps)
        write_output_results(step);

      // Increment step number
      step++;
    }

  if (verbose)
    this->pcout << "The solver took: " << step - 1 << " reinitialization steps"
                << std::endl;

  this->subequations_interface.set_solution_valid(this->subequation_id);
}

template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::write_output_results(
  const unsigned int step)
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();
  const std::string folder =
    this->simulation_parameters.simulation_control.output_folder +
    "/algebraic-reinitialization-steps-output/";

  const std::string file_name =
    "algebraic-reinitialization-" +
    this->simulation_parameters.simulation_control.output_name + "-it" +
    Utilities::int_to_string(this->simulation_control->get_iteration_number(),
                             5);

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      1, DataComponentInterpretation::component_is_scalar);

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    vector_data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);

  DataOut<dim> data_out;

  // Attach solution data to DataOut object
  data_out.attach_dof_handler(*this->dof_handler);
  data_out.add_data_vector(*this->present_solution, "reinit_phase_fraction");
  data_out.add_data_vector(this->subequations_interface.get_vof_dof_handler(),
                           this->subequations_interface.get_vof_solution(),
                           "vof_phase_fraction",
                           data_component_interpretation);
  std::vector<std::string> vof_gradient_solution_names(dim,
                                                       "vof_phase_gradient");
  data_out.add_data_vector(this->subequations_interface.get_dof_handler(
                             VOFSubequationsID::phase_gradient_projection),
                           this->subequations_interface.get_solution(
                             VOFSubequationsID::phase_gradient_projection),
                           vof_gradient_solution_names,
                           vector_data_component_interpretation);
  data_out.add_data_vector(this->subequations_interface.get_dof_handler(
                             VOFSubequationsID::curvature_projection),
                           this->subequations_interface.get_solution(
                             VOFSubequationsID::curvature_projection),
                           "vof_curvature",
                           data_component_interpretation);

  data_out.build_patches(*this->mapping,
                         1, // subdivision
                         DataOut<dim>::curved_inner_cells);

  write_vtu_and_pvd<dim>(this->pvdhandler,
                         data_out,
                         folder,
                         file_name,
                         step, // time,
                         step,
                         1, // group_files,
                         mpi_communicator);
}


template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::check_dependencies_validity()
{
  AssertThrow(this->subequations_interface.get_solution_validity(
                VOFSubequationsID::phase_gradient_projection),
              PhaseGradientProjectionIsInvalid(
                this->subequations_interface.get_subequation_string(
                  VOFSubequationsID::phase_gradient_projection),
                this->subequations_interface.get_subequation_string(
                  this->subequation_id)));
  AssertThrow(this->subequations_interface.get_solution_validity(
                VOFSubequationsID::curvature_projection),
              CurvatureProjectionIsInvalid(
                this->subequations_interface.get_subequation_string(
                  VOFSubequationsID::curvature_projection),
                this->subequations_interface.get_subequation_string(
                  this->subequation_id)));
}

template class VOFAlgebraicInterfaceReinitialization<2>;
template class VOFAlgebraicInterfaceReinitialization<3>;
