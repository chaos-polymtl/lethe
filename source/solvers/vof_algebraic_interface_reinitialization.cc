// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/vof_algebraic_interface_reinitialization.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_solver.h>

template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::setup_dofs()
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
  define_non_zero_constraints();
  define_zero_constraints();

  // Sparsity pattern
  DynamicSparsityPattern dsp(this->locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(this->dof_handler,
                                  dsp,
                                  this->nonzero_constraints,
                                  true);
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
  this->previous_solution.reinit(this->locally_owned_dofs,
                                 this->locally_relevant_dofs,
                                 mpi_communicator);
  this->newton_update.reinit(this->locally_owned_dofs, mpi_communicator);
  this->local_evaluation_point.reinit(
    this->locally_owned_dofs,
    mpi_communicator);
  this->evaluation_point = this->present_solution;

  if (this->subequation_verbosity != Parameters::Verbosity::quiet)
    {
      std::string subequation_string =
        this->subequations_interface->get_subequation_string(
          this->subequation_id);
      this->pcout << "   Number of " << subequation_string
                  << " degrees of freedom: " << this->dof_handler.n_dofs()
                  << std::endl;
    }

  // Provide DoFHandler and solutions to the subequations interface
  this->subequations_interface->set_dof_handler(this->subequation_id,
                                                &this->dof_handler);
  this->subequations_interface->set_solution(this->subequation_id,
                                             &this->present_solution);
}


template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::define_zero_constraints()
{
  this->zero_constraints.clear();
  this->zero_constraints.reinit(this->locally_relevant_dofs);

  DoFTools::make_hanging_node_constraints(this->dof_handler,
                                          this->zero_constraints);

  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            this->dof_handler,
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
            this->dof_handler,
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
  this->nonzero_constraints.reinit(this->locally_relevant_dofs);

  DoFTools::make_hanging_node_constraints(this->dof_handler,
                                          this->nonzero_constraints);

  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            this->dof_handler,
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
            this->dof_handler,
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
  const DoFHandler<dim> *dof_handler_vof =
    this->multiphysics_interface->get_dof_handler(PhysicsID::VOF);

  // Interpolate VOF solution to algebraic interface reinitialization
  VectorTools::interpolate_to_different_mesh(
    *dof_handler_vof,
    *this->multiphysics_interface->get_solution(PhysicsID::VOF),
    this->dof_handler,
    this->local_evaluation_point);

  this->nonzero_constraints.distribute(this->local_evaluation_point);
  this->present_solution = this->local_evaluation_point;
  this->previous_solution =
    this->present_solution; // We only have 1 previous solution (bdf1)
}


template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::assemble_system_matrix()
{
  // Reinitialize system matrix
  this->system_matrix = 0;

  // Get projected phase gradient DoFHandler
  const DoFHandler<dim> *dof_handler_phase_fraction_gradient =
    this->subequations_interface->get_dof_handler(
      VOFSubequationsID::phase_gradient_projection);

  // Initialize FEValues for interface algebraic reinitialization and VOF
  FEValues<dim> fe_values_algebraic_reinitialization(*this->mapping,
                                                     *this->fe,
                                                     *this->cell_quadrature,
                                                     update_values |
                                                       update_gradients |
                                                       update_JxW_values);
  FEValues<dim> fe_values_phase_gradient_projection(
    *this->mapping,
    dof_handler_phase_fraction_gradient->get_fe(),
    *this->cell_quadrature,
    update_values);

  // Initialize size of arrays
  const double n_q_points =
    fe_values_algebraic_reinitialization.get_quadrature().size();
  const double n_dofs_per_cell =
    fe_values_algebraic_reinitialization.get_fe().n_dofs_per_cell();

  // Initialize local matrix
  FullMatrix<double> local_matrix(n_dofs_per_cell, n_dofs_per_cell);

  // Initialize local dof indices array
  std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);

  // Extractor for phase fraction gradient vector
  FEValuesExtractors::Vector phase_fraction_gradients(0);

  // Initialize phase fraction solution array
  std::vector<double>         present_phase_fraction_values(n_q_points);
  std::vector<Tensor<1, dim>> present_phase_gradient_projection_values(
    n_q_points);

  // Initialize shape function arrays
  std::vector<double>         phi(n_dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_phi(n_dofs_per_cell);

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Reinitialize local matrix
          local_matrix = 0;

          // Get phase gradient projection cell iterator
          typename DoFHandler<dim>::active_cell_iterator
            phase_gradient_projection_cell(&(*this->triangulation),
                                           cell->level(),
                                           cell->index(),
                                           dof_handler_phase_fraction_gradient);

          // Reinitialize FEValues with corresponding cells
          fe_values_algebraic_reinitialization.reinit(cell);
          fe_values_phase_gradient_projection.reinit(
            phase_gradient_projection_cell);

          // Get vector of Jacobi determinant times the quadrature weights
          std::vector<double> JxW_vec =
            fe_values_algebraic_reinitialization.get_JxW_values();

          // Compute cell size
          auto &fe_algebraic_reinitialization =
            fe_values_algebraic_reinitialization.get_fe();
          const double h =
            compute_cell_diameter<dim>(compute_cell_measure_with_JxW(JxW_vec),
                                       fe_algebraic_reinitialization.degree);

          // Compute diffusivity coefficient
          const double diffusivity_coefficient =
            2 * std::pow(h, 0.9); // TODO AMISHGA change power with parameter

          // Get get present phase fraction and projected phase gradient values
          fe_values_algebraic_reinitialization.get_function_values(
            this->evaluation_point, present_phase_fraction_values);
          fe_values_phase_gradient_projection[phase_fraction_gradients]
            .get_function_values(
              *this->subequations_interface->get_solution(
                VOFSubequationsID::phase_gradient_projection),
              present_phase_gradient_projection_values);

          // Set tolerance to avoid division by zero
          const double tolerance = 1e-12;

          // BDF coefficients for pseudo time-stepping
          const Parameters::SimulationControl::TimeSteppingMethod method =
            Parameters::SimulationControl::TimeSteppingMethod::bdf1;
          const Vector<double> bdf_coefficient_vector =
            calculate_bdf_coefficients(
              method, this->simulation_control.get_time_steps_vector());

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

              // Extract phase fraction value
              const double phase_fraction = present_phase_fraction_values[q];

              // Compute normal vector with projected phase gradient
              const Tensor<1, dim> projected_vof_phase_gradient =
                present_phase_gradient_projection_values[q];
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
                          scalar_product(grad_phi[i],
                                         // Compressive term
                                         (phi[j] * (1 - 2 * phase_fraction) *
                                            interface_normal
                                          // Diffusive term
                                          - diffusivity_coefficient *
                                              scalar_product(grad_phi[j],
                                                             interface_normal) *
                                              interface_normal))) *
                        JxW_vec[q];
                    }
                }
            }
          // Distribute the local contributions to the global system
          cell->get_dof_indices(local_dof_indices);
          this->nonzero_constraints.distribute_local_to_global(
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

  // Get projected phase gradient DoFHandler
  const DoFHandler<dim> *dof_handler_phase_fraction_gradient =
    this->subequations_interface->get_dof_handler(
      VOFSubequationsID::phase_gradient_projection);

  // Initialize FEValues for interface algebraic reinitialization and phase
  // fraction gradient projeciton
  FEValues<dim> fe_values_algebraic_reinitialization(*this->mapping,
                                                     *this->fe,
                                                     *this->cell_quadrature,
                                                     update_values |
                                                       update_gradients |
                                                       update_JxW_values);
  FEValues<dim> fe_values_phase_gradient_projection(
    *this->mapping,
    dof_handler_phase_fraction_gradient->get_fe(),
    *this->cell_quadrature,
    update_values);

  // Initialize size of arrays
  const double n_q_points =
    fe_values_algebraic_reinitialization.get_quadrature().size();
  const double n_dofs_per_cell =
    fe_values_algebraic_reinitialization.get_fe().n_dofs_per_cell();

  //  Initialize local rhs
  Vector<double> local_rhs(n_dofs_per_cell);

  // Initialize local dof indices array
  std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);

  // Extractor for phase fraction gradient vector
  FEValuesExtractors::Vector phase_fraction_gradients(0);

  // Initialize phase fraction and projected phase fraction gradient solution
  // arrays
  std::vector<double>         present_phase_fraction_values(n_q_points);
  std::vector<double>         previous_phase_fraction_values(n_q_points);
  std::vector<Tensor<1, dim>> present_phase_gradient_projection_values(
    n_q_points);

  // Initialize shape function arrays
  std::vector<double>         phi(n_dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_phi(n_dofs_per_cell);

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Reinitialize local rhs
          local_rhs = 0;

          // Get phase gradient projection cell iterator
          typename DoFHandler<dim>::active_cell_iterator
            phase_gradient_projection_cell(&(*this->triangulation),
                                           cell->level(),
                                           cell->index(),
                                           dof_handler_phase_fraction_gradient);

          // Reinitialize FEValues with corresponding cells
          fe_values_algebraic_reinitialization.reinit(cell);
          fe_values_phase_gradient_projection.reinit(
            phase_gradient_projection_cell);

          // Get vector of Jacobi determinant times the quadrature weights
          std::vector<double> JxW_vec =
            fe_values_algebraic_reinitialization.get_JxW_values();


          // Compute cell size
          auto &fe_algebraic_reinitialization =
            fe_values_algebraic_reinitialization.get_fe();
          const double h =
            compute_cell_diameter<dim>(compute_cell_measure_with_JxW(JxW_vec),
                                       fe_algebraic_reinitialization.degree);

          // Compute diffusivity coefficient
          const double diffusivity_coefficient =
            2 * std::pow(h, 0.9); // TODO AMISHGA change power with parameter

          // Get VOF phase fraction and projected phase gradient values
          fe_values_algebraic_reinitialization.get_function_values(
            this->evaluation_point, present_phase_fraction_values);
          fe_values_algebraic_reinitialization.get_function_values(
            this->previous_solution, previous_phase_fraction_values);
          fe_values_phase_gradient_projection[phase_fraction_gradients]
            .get_function_values(
              *this->subequations_interface->get_solution(
                VOFSubequationsID::phase_gradient_projection),
              present_phase_gradient_projection_values);

          // Set tolerance to avoid division by zero
          const double tolerance = 1e-12;

          // BDF coefficients for pseudo time-stepping
          const Parameters::SimulationControl::TimeSteppingMethod method =
            Parameters::SimulationControl::TimeSteppingMethod::bdf1;
          const Vector<double> bdf_coefficient_vector =
            calculate_bdf_coefficients(
              method, this->simulation_control.get_time_steps_vector());

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
              phase_fraction_values[0] = present_phase_fraction_values[q];
              phase_fraction_values[1] = previous_phase_fraction_values[q];
              const Tensor<1, dim> projected_vof_phase_gradient =
                present_phase_gradient_projection_values[q];

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
                    (scalar_product(
                      grad_phi[i],
                      (
                        // Compressive term
                        (phase_fraction_values[0] *
                         (1 - phase_fraction_values[0])) *
                          interface_normal
                        // Diffusive term
                        - diffusivity_coefficient *
                            scalar_product(projected_vof_phase_gradient,
                                           interface_normal) *
                            interface_normal))) *
                    JxW_vec[q];
                }
            }
          // Distribute the local contributions to the global system
          cell->get_dof_indices(local_dof_indices);
          this->nonzero_constraints.distribute_local_to_global(
            local_rhs, local_dof_indices, this->system_rhs);
        }
    }
  this->system_rhs.compress(VectorOperation::add);
}


template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::solve_linear_system(
  const bool initial_step,
  const bool /*renewed_matrix*/)
{
  auto mpi_communicator = this->triangulation->get_communicator();

  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;

  const bool verbose(
    this->subequation_verbosity != Parameters::Verbosity::quiet &&
    this->linear_solver_verbosity != Parameters::Verbosity::quiet);

  if (verbose)
    {
      std::string subequation_string =
        this->subequations_interface->get_subequation_string(
          this->subequation_id);

      this->pcout << "  -Solving " << subequation_string << ":" << std::endl;
    }

  // Get residual conditions
  const double absolute_residual =
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .minimum_residual;
  const double relative_residual =
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .relative_residual;

  // Set tolerance
  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);

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

  GlobalVectorType completely_distributed_solution(this->locally_owned_dofs,
                                                   mpi_communicator);

  SolverControl solver_control(
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF).max_iterations,
    linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    false,
    this->simulation_parameters.linear_solver.at(PhysicsID::VOF)
      .max_krylov_vectors);

  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);

  solver.solve(this->system_matrix,
               completely_distributed_solution,
               this->system_rhs,
               *this->ilu_preconditioner);

  if (verbose)
    {
      this->pcout << "  -Iterative solver took: " << solver_control.last_step()
                  << " steps to reach a residual norm of "
                  << solver_control.last_value() << std::endl;
    }

  // Update constraints and newton vectors
  constraints_used.distribute(completely_distributed_solution);
  this->newton_update = completely_distributed_solution;
}


template <int dim>
void
VOFAlgebraicInterfaceReinitialization<dim>::solve(
  const bool & /*is_post_mesh_adaptation = false*/)
{
  // Get time-step
  double current_time_step_inv =
    1. / this->simulation_control.get_time_steps_vector()[0];

  // TODO AMISHGA ERASE
  this->pcout << "Time-step: " << this->simulation_control.get_time_steps_vector()[0] << std::endl;

  // Set tolerance for the stop criterion of the pseudo-time-stepping scheme
  double tolerance = 1e-6;

  // Set initial conditions from the VOF solution
  set_initial_conditions();

  // Solve a first time-step
  this->solve_non_linear_system(false);

  // Number of iterations counter
  unsigned int it = 1;

  // Iterate until the stop criterion is met
  while (continue_iterating(current_time_step_inv, tolerance))
    {
      // Update previous solution
      this->previous_solution = this->present_solution;

      // Update non-zero constraints
      define_non_zero_constraints(); // TODO AMISHGA check if necessary

      // Solve non-linear equation
      this->solve_non_linear_system(false);

      it += 1;
    }

  if (this->subequation_verbosity != Parameters::Verbosity::quiet)
    this->pcout << "The solver took: " << it << " reinitialization steps"
                << std::endl;
}


template class VOFAlgebraicInterfaceReinitialization<2>;
template class VOFAlgebraicInterfaceReinitialization<3>;
