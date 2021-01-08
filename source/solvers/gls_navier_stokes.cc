/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#include "solvers/gls_navier_stokes.h"

#include "core/bdf.h"
#include "core/grids.h"
#include "core/manifolds.h"
#include "core/sdirk.h"
#include "core/time_integration_utilities.h"

// Constructor for class GLSNavierStokesSolver
template <int dim>
GLSNavierStokesSolver<dim>::GLSNavierStokesSolver(
  NavierStokesSolverParameters<dim> &p_nsparam)
  : NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>(p_nsparam)
{}

template <int dim>
GLSNavierStokesSolver<dim>::~GLSNavierStokesSolver()
{
  this->dof_handler.clear();
}



template <int dim>
void
GLSNavierStokesSolver<dim>::set_solution_vector(double value)
{
  auto &present_solution = this->get_present_solution();
  present_solution       = value;
}

template <int dim>
void
GLSNavierStokesSolver<dim>::setup_dofs_cfd()
{
  TimerOutput::Scope t(this->computing_timer, "setup_dofs");

  // Clear the preconditioner before the matrix they are associated with is
  // cleared
  amg_preconditioner.reset();
  ilu_preconditioner.reset();

  // Now reset system matrix
  system_matrix.clear();

  this->dof_handler.distribute_dofs(this->fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  this->locally_owned_dofs = this->dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                          this->locally_relevant_dofs);

  const MappingQ<dim>        mapping(this->velocity_fem_degree,
                              this->nsparam.fem_parameters.qmapping_all);
  FEValuesExtractors::Vector velocities(0);

  // Non-zero constraints
  auto &nonzero_constraints =
    this->get_nonzero_constraints(Parameters::Multiphysics::ID::fluid);
  {
    nonzero_constraints.clear();

    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            nonzero_constraints);
    for (unsigned int i_bc = 0; i_bc < this->nsparam.boundary_conditions.size;
         ++i_bc)
      {
        if (this->nsparam.boundary_conditions.type[i_bc] ==
            BoundaryConditions::BoundaryType::noslip)
          {
            VectorTools::interpolate_boundary_values(
              mapping,
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              dealii::Functions::ZeroFunction<dim>(dim + 1),
              nonzero_constraints,
              this->fe.component_mask(velocities));
          }
        else if (this->nsparam.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert(
              this->nsparam.boundary_conditions.id[i_bc]);
            VectorTools::compute_no_normal_flux_constraints(
              this->dof_handler,
              0,
              no_normal_flux_boundaries,
              nonzero_constraints);
          }
        else if (this->nsparam.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::function)
          {
            VectorTools::interpolate_boundary_values(
              mapping,
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              NavierStokesFunctionDefined<dim>(
                &this->nsparam.boundary_conditions.bcFunctions[i_bc].u,
                &this->nsparam.boundary_conditions.bcFunctions[i_bc].v,
                &this->nsparam.boundary_conditions.bcFunctions[i_bc].w),
              nonzero_constraints,
              this->fe.component_mask(velocities));
          }

        else if (this->nsparam.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              this->nsparam.boundary_conditions.periodic_id[i_bc],
              this->nsparam.boundary_conditions.periodic_direction[i_bc],
              nonzero_constraints);
          }
      }
  }
  nonzero_constraints.close();

  {
    this->zero_constraints.clear();
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            this->zero_constraints);

    for (unsigned int i_bc = 0; i_bc < this->nsparam.boundary_conditions.size;
         ++i_bc)
      {
        if (this->nsparam.boundary_conditions.type[i_bc] ==
            BoundaryConditions::BoundaryType::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert(
              this->nsparam.boundary_conditions.id[i_bc]);
            VectorTools::compute_no_normal_flux_constraints(
              this->dof_handler,
              0,
              no_normal_flux_boundaries,
              this->zero_constraints);
          }
        else if (this->nsparam.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              this->nsparam.boundary_conditions.periodic_id[i_bc],
              this->nsparam.boundary_conditions.periodic_direction[i_bc],
              this->zero_constraints);
          }
        else // if(nsparam.boundaryConditions.boundaries[i_bc].type==Parameters::noslip
             // || Parameters::function)
          {
            VectorTools::interpolate_boundary_values(
              mapping,
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              dealii::Functions::ZeroFunction<dim>(dim + 1),
              this->zero_constraints,
              this->fe.component_mask(velocities));
          }
      }
  }
  this->zero_constraints.close();

  TrilinosWrappers::MPI::Vector &present_solution =
    this->get_present_solution();
  present_solution.reinit(this->locally_owned_dofs,
                          this->locally_relevant_dofs,
                          this->mpi_communicator);
  this->solution_m1.reinit(this->locally_owned_dofs,
                           this->locally_relevant_dofs,
                           this->mpi_communicator);
  this->solution_m2.reinit(this->locally_owned_dofs,
                           this->locally_relevant_dofs,
                           this->mpi_communicator);
  this->solution_m3.reinit(this->locally_owned_dofs,
                           this->locally_relevant_dofs,
                           this->mpi_communicator);

  TrilinosWrappers::MPI::Vector &newton_update = this->get_newton_update();
  newton_update.reinit(this->locally_owned_dofs, this->mpi_communicator);
  TrilinosWrappers::MPI::Vector &system_rhs = this->get_system_rhs();
  system_rhs.reinit(this->locally_owned_dofs, this->mpi_communicator);
  TrilinosWrappers::MPI::Vector &local_evaluation_point =
    this->get_local_evaluation_point();
  local_evaluation_point.reinit(this->locally_owned_dofs,
                                this->mpi_communicator);

  DynamicSparsityPattern dsp(this->locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(this->dof_handler,
                                  dsp,
                                  nonzero_constraints,
                                  false);
  SparsityTools::distribute_sparsity_pattern(
    dsp,
    this->dof_handler.locally_owned_dofs(),
    this->mpi_communicator,
    this->locally_relevant_dofs);
  system_matrix.reinit(this->locally_owned_dofs,
                       this->locally_owned_dofs,
                       dsp,
                       this->mpi_communicator);

  if (this->nsparam.post_processing.calculate_average_velocities)
    {
      AssertThrow(this->nsparam.mesh_adaptation.type ==
                    Parameters::MeshAdaptation::Type::none,
                  ExcMessage(
                    "Time-averaging velocities and calculating reynolds "
                    "stresses are currently unavailable for mesh "
                    "adaptation."));

      this->average_velocities.initialize_vectors(this->locally_owned_dofs,
                                                  this->locally_relevant_dofs,
                                                  this->fe.n_dofs_per_vertex(),
                                                  this->mpi_communicator);

      if (this->nsparam.restart_parameters.checkpoint)
        {
          this->average_velocities.initialize_checkpoint_vectors(
            this->locally_owned_dofs,
            this->locally_relevant_dofs,
            this->mpi_communicator);
        }
    }

  double global_volume = GridTools::volume(*this->triangulation);

  this->pcout << "   Number of active cells:       "
              << this->triangulation->n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: "
              << this->dof_handler.n_dofs() << std::endl;
  this->pcout << "   Volume of triangulation:      " << global_volume
              << std::endl;
}

template <int dim>
template <bool                                              assemble_matrix,
          Parameters::SimulationControl::TimeSteppingMethod scheme,
          Parameters::VelocitySource::VelocitySourceType    velocity_source>
void
GLSNavierStokesSolver<dim>::assembleGLS()
{
  if (assemble_matrix)
    system_matrix = 0;
  auto &system_rhs = this->get_system_rhs();
  system_rhs       = 0;

  double         viscosity = this->nsparam.physical_properties.viscosity;
  Function<dim> *l_forcing_function = this->forcing_function;

  QGauss<dim>         quadrature_formula(this->number_quadrature_points);
  const MappingQ<dim> mapping(this->velocity_fem_degree,
                              this->nsparam.fem_parameters.qmapping_all);
  FEValues<dim>       fe_values(mapping,
                          this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients |
                            update_hessians);
  const unsigned int  dofs_per_cell = this->fe.dofs_per_cell;
  const unsigned int  n_q_points    = quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  FullMatrix<double>               local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>                   local_rhs(dofs_per_cell);
  std::vector<Vector<double>> rhs_force(n_q_points, Vector<double>(dim + 1));
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<Tensor<1, dim>>          present_velocity_values(n_q_points);
  std::vector<Tensor<2, dim>>          present_velocity_gradients(n_q_points);
  std::vector<double>                  present_pressure_values(n_q_points);
  std::vector<Tensor<1, dim>>          present_pressure_gradients(n_q_points);
  std::vector<Tensor<1, dim>>          present_velocity_laplacians(n_q_points);
  std::vector<Tensor<2, dim>>          present_velocity_hess(n_q_points);

  Tensor<1, dim> force;
  Tensor<1, dim> beta_force = this->beta;

  // Velocity dependent source term
  //----------------------------------
  // Angular velocity of the rotating frame. This is always a 3D vector even in
  // 2D.
  Tensor<1, dim> omega_vector;

  double omega_z  = this->nsparam.velocitySource.omega_z;
  omega_vector[0] = this->nsparam.velocitySource.omega_x;
  omega_vector[1] = this->nsparam.velocitySource.omega_y;
  if (dim == 3)
    omega_vector[2] = this->nsparam.velocitySource.omega_z;

  std::vector<double>         div_phi_u(dofs_per_cell);
  std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
  std::vector<Tensor<3, dim>> hess_phi_u(dofs_per_cell);
  std::vector<Tensor<1, dim>> laplacian_phi_u(dofs_per_cell);
  std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
  std::vector<double>         phi_p(dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_phi_p(dofs_per_cell);

  // Values at previous time step for transient schemes
  std::vector<Tensor<1, dim>> p1_velocity_values(n_q_points);
  std::vector<Tensor<1, dim>> p2_velocity_values(n_q_points);
  std::vector<Tensor<1, dim>> p3_velocity_values(n_q_points);

  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Time steps and inverse time steps which is used for numerous calculations
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Vector for the BDF coefficients
  // The coefficients are stored in the following fashion :
  // 0 - n+1
  // 1 - n
  // 2 - n-1
  // 3 - n-2
  Vector<double> bdf_coefs;

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
      scheme == Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
    bdf_coefs = bdf_coefficients(1, time_steps_vector);

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    bdf_coefs = bdf_coefficients(2, time_steps_vector);

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    bdf_coefs = bdf_coefficients(3, time_steps_vector);

  // Matrix of coefficients for the SDIRK methods
  // The lines store the information required for each step
  // Column 0 always refer to outcome of the step that is being calculated
  // Column 1 always refer to step n
  // Column 2+ refer to intermediary steps
  FullMatrix<double> sdirk_coefs;
  if (is_sdirk2(scheme))
    sdirk_coefs = sdirk_coefficients(2, dt);

  if (is_sdirk3(scheme))
    sdirk_coefs = sdirk_coefficients(3, dt);

  // Element size
  double h;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          auto &evaluation_point = this->get_evaluation_point();
          fe_values.reinit(cell);

          if (dim == 2)
            h = std::sqrt(4. * cell->measure() / M_PI) /
                this->velocity_fem_degree;
          else if (dim == 3)
            h = pow(6 * cell->measure() / M_PI, 1. / 3.) /
                this->velocity_fem_degree;

          local_matrix = 0;
          local_rhs    = 0;

          // Gather velocity (values, gradient and laplacian)
          fe_values[velocities].get_function_values(evaluation_point,
                                                    present_velocity_values);
          fe_values[velocities].get_function_gradients(
            evaluation_point, present_velocity_gradients);
          fe_values[velocities].get_function_laplacians(
            evaluation_point, present_velocity_laplacians);

          // Gather pressure (values, gradient)
          fe_values[pressure].get_function_values(evaluation_point,
                                                  present_pressure_values);
          fe_values[pressure].get_function_gradients(
            evaluation_point, present_pressure_gradients);

          std::vector<Point<dim>> quadrature_points =
            fe_values.get_quadrature_points();

          // Calculate forcing term if there is a forcing function
          if (l_forcing_function)
            l_forcing_function->vector_value_list(quadrature_points, rhs_force);

          // Gather the previous time steps depending on the number of stages
          // of the time integration scheme
          if (scheme !=
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            fe_values[velocities].get_function_values(this->solution_m1,
                                                      p1_velocity_values);

          if (time_stepping_method_has_two_stages(scheme))
            fe_values[velocities].get_function_values(this->solution_m2,
                                                      p2_velocity_values);

          if (time_stepping_method_has_three_stages(scheme))
            fe_values[velocities].get_function_values(this->solution_m3,
                                                      p3_velocity_values);

          // Loop over the quadrature points
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              // Gather into local variables the relevant fields
              const Tensor<1, dim> velocity = present_velocity_values[q];
              const Tensor<2, dim> velocity_gradient =
                present_velocity_gradients[q];
              const double present_velocity_divergence =
                trace(velocity_gradient);
              const Tensor<1, dim> p1_velocity = p1_velocity_values[q];
              const Tensor<1, dim> p2_velocity = p2_velocity_values[q];
              const Tensor<1, dim> p3_velocity = p3_velocity_values[q];
              const double current_pressure    = present_pressure_values[q];



              // Calculation of the magnitude of the velocity for the
              // stabilization parameter
              const double u_mag =
                std::max(velocity.norm(), 1e-12 * GLS_u_scale);

              // Store JxW in local variable for faster access;
              const double JxW = fe_values.JxW(q);

              // Calculation of the GLS stabilization parameter. The
              // stabilization parameter used is different if the simulation is
              // steady or unsteady. In the unsteady case it includes the value
              // of the time-step
              const double tau =
                is_steady(scheme) ?
                  1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                                 9 * std::pow(4 * viscosity / (h * h), 2)) :
                  1. /
                    std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                              9 * std::pow(4 * viscosity / (h * h), 2));

              // Gather the shape functions, their gradient and their laplacian
              // for the velocity and the pressure
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  div_phi_u[k]  = fe_values[velocities].divergence(k, q);
                  grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                  phi_u[k]      = fe_values[velocities].value(k, q);
                  hess_phi_u[k] = fe_values[velocities].hessian(k, q);
                  phi_p[k]      = fe_values[pressure].value(k, q);
                  grad_phi_p[k] = fe_values[pressure].gradient(k, q);

                  for (int d = 0; d < dim; ++d)
                    laplacian_phi_u[k][d] = trace(hess_phi_u[k][d]);
                }

              // Establish the force vector
              for (int i = 0; i < dim; ++i)
                {
                  const unsigned int component_i =
                    this->fe.system_to_component_index(i).first;
                  force[i] = rhs_force[q](component_i);
                }
              // Correct force to include the dynamic forcing term for flow
              // control
              force = force + beta_force;

              // Calculate the strong residual for GLS stabilization
              auto strong_residual =
                velocity_gradient * velocity + present_pressure_gradients[q] -
                viscosity * present_velocity_laplacians[q] - force;

              if (velocity_source ==
                  Parameters::VelocitySource::VelocitySourceType::srf)
                {
                  if (dim == 2)
                    {
                      strong_residual +=
                        2 * omega_z * (-1.) * cross_product_2d(velocity);
                      auto centrifugal =
                        omega_z * (-1.) *
                        cross_product_2d(
                          omega_z * (-1.) *
                          cross_product_2d(quadrature_points[q]));
                      strong_residual += centrifugal;
                    }
                  else // dim == 3
                    {
                      strong_residual +=
                        2 * cross_product_3d(omega_vector, velocity);
                      strong_residual += cross_product_3d(
                        omega_vector,
                        cross_product_3d(omega_vector, quadrature_points[q]));
                    }
                }

              /* Adjust the strong residual in cases where the scheme is
               transient.
               The BDF schemes require values at previous time steps which are
               stored in the p1, p2 and p3 vectors. The SDIRK scheme require the
               values at the different stages, which are also stored in the same
               arrays.
               */

              if (scheme ==
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
                  scheme == Parameters::SimulationControl::TimeSteppingMethod::
                              steady_bdf)
                strong_residual += bdf_coefs[0] * velocity +
                                   bdf_coefs[1] * p1_velocity_values[q];

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                strong_residual += bdf_coefs[0] * velocity +
                                   bdf_coefs[1] * p1_velocity +
                                   bdf_coefs[2] * p2_velocity;

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                strong_residual +=
                  bdf_coefs[0] * velocity + bdf_coefs[1] * p1_velocity +
                  bdf_coefs[2] * p2_velocity + bdf_coefs[3] * p3_velocity;


              if (is_sdirk_step1(scheme))
                strong_residual += sdirk_coefs[0][0] * velocity +
                                   sdirk_coefs[0][1] * p1_velocity;

              if (is_sdirk_step2(scheme))
                {
                  strong_residual += sdirk_coefs[1][0] * velocity +
                                     sdirk_coefs[1][1] * p1_velocity +
                                     sdirk_coefs[1][2] * p2_velocity;
                }

              if (is_sdirk_step3(scheme))
                {
                  strong_residual += sdirk_coefs[2][0] * velocity +
                                     sdirk_coefs[2][1] * p1_velocity +
                                     sdirk_coefs[2][2] * p2_velocity +
                                     sdirk_coefs[2][3] * p3_velocity;
                }

              // Matrix assembly
              if (assemble_matrix)
                {
                  // We loop over the column first to prevent recalculation of
                  // the strong jacobian in the inner loop
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      const auto phi_u_j      = phi_u[j];
                      const auto grad_phi_u_j = grad_phi_u[j];
                      const auto phi_p_j      = phi_p[j];
                      const auto grad_phi_p_j = grad_phi_p[j];



                      auto strong_jac =
                        (velocity_gradient * phi_u_j + grad_phi_u_j * velocity +
                         grad_phi_p_j - viscosity * laplacian_phi_u[j]);

                      if (is_bdf(scheme))
                        strong_jac += phi_u_j * bdf_coefs[0];
                      if (is_sdirk(scheme))
                        strong_jac += phi_u_j * sdirk_coefs[0][0];

                      if (velocity_source ==
                          Parameters::VelocitySource::VelocitySourceType::srf)
                        {
                          if (dim == 2)
                            strong_jac +=
                              2 * omega_z * (-1.) * cross_product_2d(phi_u_j);
                          else if (dim == 3)
                            strong_jac +=
                              2 * cross_product_3d(omega_vector, phi_u_j);
                        }

                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                          const auto phi_u_i      = phi_u[i];
                          const auto grad_phi_u_i = grad_phi_u[i];
                          const auto phi_p_i      = phi_p[i];
                          const auto grad_phi_p_i = grad_phi_p[i];


                          local_matrix(i, j) +=
                            (
                              // Momentum terms
                              viscosity *
                                scalar_product(grad_phi_u_j, grad_phi_u_i) +
                              velocity_gradient * phi_u_j * phi_u_i +
                              grad_phi_u_j * velocity * phi_u_i -
                              div_phi_u[i] * phi_p_j +
                              // Continuity
                              phi_p_i * div_phi_u[j]) *
                            JxW;

                          // Mass matrix
                          if (is_bdf(scheme))
                            local_matrix(i, j) +=
                              phi_u_j * phi_u_i * bdf_coefs[0] * JxW;

                          if (is_sdirk(scheme))
                            local_matrix(i, j) +=
                              phi_u_j * phi_u_i * sdirk_coefs[0][0] * JxW;

                          // PSPG GLS term
                          local_matrix(i, j) +=
                            tau * (strong_jac * grad_phi_p_i) * JxW;

                          if (velocity_source == Parameters::VelocitySource::
                                                   VelocitySourceType::srf)
                            {
                              if (dim == 2)
                                local_matrix(i, j) +=
                                  2 * omega_z * (-1.) *
                                  cross_product_2d(phi_u_j) * phi_u_i * JxW;

                              else if (dim == 3)
                                local_matrix(i, j) +=
                                  2 * cross_product_3d(omega_vector, phi_u_j) *
                                  phi_u_i * JxW;
                            }


                          // PSPG TAU term is currently disabled because it does
                          // not alter the matrix sufficiently
                          // local_matrix(i, j) +=
                          //  -tau * tau * tau * 4 / h / h *
                          //  (velocity *phi_u_j) *
                          //  strong_residual * grad_phi_p_i *
                          //  fe_values.JxW(q);

                          // Jacobian is currently incomplete
                          if (SUPG)
                            {
                              local_matrix(i, j) +=
                                tau *
                                (strong_jac * (grad_phi_u_i * velocity) +
                                 strong_residual * (grad_phi_u_i * phi_u_j)) *
                                JxW;

                              // SUPG TAU term is currently disabled because it
                              // does not alter the matrix sufficiently
                              // local_matrix(i, j)
                              // +=
                              //   -strong_residual
                              //   * (grad_phi_u_i
                              //   *
                              //   velocity)
                              //   * tau * tau *
                              //   tau * 4 / h / h
                              //   *
                              //   (velocity
                              //   *phi_u_j) *
                              //   fe_values.JxW(q);
                            }
                        }
                    }
                }

              // Assembly of the right-hand side
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  const auto phi_u_i      = phi_u[i];
                  const auto grad_phi_u_i = grad_phi_u[i];
                  const auto phi_p_i      = phi_p[i];
                  const auto grad_phi_p_i = grad_phi_p[i];
                  const auto div_phi_u_i  = div_phi_u[i];


                  // Navier-Stokes Residual
                  local_rhs(i) +=
                    (
                      // Momentum
                      -viscosity *
                        scalar_product(velocity_gradient, grad_phi_u_i) -
                      velocity_gradient * velocity * phi_u_i +
                      current_pressure * div_phi_u_i + force * phi_u_i -
                      // Continuity
                      present_velocity_divergence * phi_p_i) *
                    JxW;

                  // Residual associated with BDF schemes
                  if (scheme == Parameters::SimulationControl::
                                  TimeSteppingMethod::bdf1 ||
                      scheme == Parameters::SimulationControl::
                                  TimeSteppingMethod::steady_bdf)
                    local_rhs(i) -=
                      bdf_coefs[0] * (velocity - p1_velocity) * phi_u_i * JxW;

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                    local_rhs(i) -= (bdf_coefs[0] * (velocity * phi_u_i) +
                                     bdf_coefs[1] * (p1_velocity * phi_u_i) +
                                     bdf_coefs[2] * (p2_velocity * phi_u_i)) *
                                    JxW;

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                    local_rhs(i) -= (bdf_coefs[0] * (velocity * phi_u_i) +
                                     bdf_coefs[1] * (p1_velocity * phi_u_i) +
                                     bdf_coefs[2] * (p2_velocity * phi_u_i) +
                                     bdf_coefs[3] * (p3_velocity * phi_u_i)) *
                                    JxW;

                  // Residuals associated with SDIRK schemes
                  if (is_sdirk_step1(scheme))
                    local_rhs(i) -=
                      (sdirk_coefs[0][0] * (velocity * phi_u_i) +
                       sdirk_coefs[0][1] * (p1_velocity * phi_u_i)) *
                      JxW;

                  if (is_sdirk_step2(scheme))
                    {
                      local_rhs(i) -=
                        (sdirk_coefs[1][0] * (velocity * phi_u_i) +
                         sdirk_coefs[1][1] * (p1_velocity * phi_u_i) +
                         sdirk_coefs[1][2] *
                           (p2_velocity_values[q] * phi_u_i)) *
                        JxW;
                    }

                  if (is_sdirk_step3(scheme))
                    {
                      local_rhs(i) -=
                        (sdirk_coefs[2][0] * (velocity * phi_u_i) +
                         sdirk_coefs[2][1] * (p1_velocity * phi_u_i) +
                         sdirk_coefs[2][2] * (p2_velocity * phi_u_i) +
                         sdirk_coefs[2][3] * (p3_velocity * phi_u_i)) *
                        JxW;
                    }

                  if (velocity_source ==
                      Parameters::VelocitySource::VelocitySourceType::srf)
                    {
                      if (dim == 2)
                        {
                          local_rhs(i) += -2 * omega_z * (-1.) *
                                          cross_product_2d(velocity) * phi_u_i *
                                          JxW;
                          auto centrifugal =
                            omega_z * (-1.) *
                            cross_product_2d(
                              omega_z * (-1.) *
                              cross_product_2d(quadrature_points[q]));
                          local_rhs(i) += -centrifugal * phi_u_i * JxW;
                        }
                      else if (dim == 3)
                        {
                          local_rhs(i) +=
                            -2 * cross_product_3d(omega_vector, velocity) *
                            phi_u_i * JxW;
                          local_rhs(i) +=
                            -cross_product_3d(
                              omega_vector,
                              cross_product_3d(omega_vector,
                                               quadrature_points[q])) *
                            phi_u_i * JxW;
                        }
                    }

                  // PSPG GLS term
                  local_rhs(i) += -tau * (strong_residual * grad_phi_p_i) * JxW;

                  // SUPG GLS term
                  if (SUPG)
                    {
                      local_rhs(i) +=
                        -tau * (strong_residual * (grad_phi_u_i * velocity)) *
                        JxW;
                    }
                }
            }

          cell->get_dof_indices(local_dof_indices);

          // The non-linear solver assumes that the nonzero constraints have
          // already been applied to the solution
          const AffineConstraints<double> &constraints_used =
            this->zero_constraints;
          // initial_step ? nonzero_constraints : zero_constraints;
          if (assemble_matrix)
            {
              constraints_used.distribute_local_to_global(local_matrix,
                                                          local_rhs,
                                                          local_dof_indices,
                                                          system_matrix,
                                                          system_rhs);
            }
          else
            {
              constraints_used.distribute_local_to_global(local_rhs,
                                                          local_dof_indices,
                                                          system_rhs);
            }
        }
    }
  if (assemble_matrix)
    system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

/**
 * Set the initial condition using a L2 or a viscous solver
 **/
template <int dim>
void
GLSNavierStokesSolver<dim>::set_initial_condition(
  Parameters::InitialConditionType initial_condition_type,
  bool                             restart)
{
  set_initial_condition_cfd(initial_condition_type, restart);
}

/**
 * Set the initial condition using a L2 or a viscous solver
 **/
template <int dim>
void
GLSNavierStokesSolver<dim>::set_initial_condition_cfd(
  Parameters::InitialConditionType initial_condition_type,
  bool                             restart)
{
  if (restart)
    {
      this->pcout << "************************" << std::endl;
      this->pcout << "---> Simulation Restart " << std::endl;
      this->pcout << "************************" << std::endl;
      this->read_checkpoint();
    }
  else if (initial_condition_type ==
           Parameters::InitialConditionType::L2projection)
    {
      assemble_L2_projection();
      solve_system_GMRES(true, 1e-15, 1e-15, true);
      auto &present_solution = this->get_present_solution();
      auto &newton_update    = this->get_newton_update();
      present_solution       = newton_update;
      this->finish_time_step();
      this->postprocess(true);
    }
  else if (initial_condition_type == Parameters::InitialConditionType::nodal)
    {
      this->set_nodal_values();
      this->finish_time_step();
      this->postprocess(true);
    }

  else if (initial_condition_type == Parameters::InitialConditionType::viscous)
    {
      this->set_nodal_values();
      double viscosity = this->nsparam.physical_properties.viscosity;
      this->nsparam.physical_properties.viscosity =
        this->nsparam.initial_condition->viscosity;
      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::steady, false, true);
      this->finish_time_step();
      this->postprocess(true);
      this->nsparam.physical_properties.viscosity = viscosity;
    }
  else
    {
      throw std::runtime_error("GLSNS - Initial condition could not be set");
    }
}

template <int dim>
void
GLSNavierStokesSolver<dim>::assemble_L2_projection()
{
  system_matrix    = 0;
  auto &system_rhs = this->get_system_rhs();
  system_rhs       = 0;
  QGauss<dim>         quadrature_formula(this->number_quadrature_points);
  const MappingQ<dim> mapping(this->velocity_fem_degree,
                              this->nsparam.fem_parameters.qmapping_all);
  FEValues<dim>       fe_values(mapping,
                          this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  const unsigned int  dofs_per_cell = this->fe.dofs_per_cell;
  const unsigned int  n_q_points    = quadrature_formula.size();
  FullMatrix<double>  local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>      local_rhs(dofs_per_cell);
  std::vector<Vector<double>>          initial_velocity(n_q_points,
                                                        Vector<double>(dim + 1));
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const FEValuesExtractors::Vector     velocities(0);
  const FEValuesExtractors::Scalar     pressure(dim);

  Tensor<1, dim> rhs_initial_velocity_pressure;
  double         rhs_initial_pressure;

  std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
  std::vector<double>         phi_p(dofs_per_cell);

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          local_matrix = 0;
          local_rhs    = 0;
          this->nsparam.initial_condition->uvwp.vector_value_list(
            fe_values.get_quadrature_points(), initial_velocity);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_p[k] = fe_values[pressure].value(k, q);
                  phi_u[k] = fe_values[velocities].value(k, q);
                }

              // Establish the rhs tensor operator
              for (int i = 0; i < dim; ++i)
                {
                  const unsigned int component_i =
                    this->fe.system_to_component_index(i).first;
                  rhs_initial_velocity_pressure[i] =
                    initial_velocity[q](component_i);
                }
              rhs_initial_pressure = initial_velocity[q](dim);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix(i, j) +=
                        (phi_u[j] * phi_u[i]) * fe_values.JxW(q);
                      local_matrix(i, j) +=
                        (phi_p[j] * phi_p[i]) * fe_values.JxW(q);
                    }
                  local_rhs(i) += (phi_u[i] * rhs_initial_velocity_pressure +
                                   phi_p[i] * rhs_initial_pressure) *
                                  fe_values.JxW(q);
                }
            }

          cell->get_dof_indices(local_dof_indices);
          auto &nonzero_constraints = this->get_nonzero_constraints();
          const AffineConstraints<double> &constraints_used =
            nonzero_constraints;
          constraints_used.distribute_local_to_global(local_matrix,
                                                      local_rhs,
                                                      local_dof_indices,
                                                      system_matrix,
                                                      system_rhs);
        }
    }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
GLSNavierStokesSolver<dim>::assemble_matrix_and_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_system");

  if (this->nsparam.velocitySource.type ==
      Parameters::VelocitySource::VelocitySourceType::none)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::steady_bdf,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else
        throw std::runtime_error(
          "The time stepping method provided is not supported by this solver");
    }

  else if (this->nsparam.velocitySource.type ==
           Parameters::VelocitySource::VelocitySourceType::srf)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
        assembleGLS<
          true,
          Parameters::SimulationControl::TimeSteppingMethod::steady_bdf,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else
        throw std::runtime_error(
          "The time stepping method provided is not supported by this solver");
    }


  if (this->simulation_control->is_first_assembly())
    {
      auto &system_rhs = this->get_system_rhs();
      this->simulation_control->provide_residual(system_rhs.l2_norm());
    }
}
template <int dim>
void
GLSNavierStokesSolver<dim>::assemble_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_rhs");

  if (this->nsparam.velocitySource.type ==
      Parameters::VelocitySource::VelocitySourceType::none)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::steady_bdf,
          Parameters::VelocitySource::VelocitySourceType::none>();
      else
        throw std::runtime_error(
          "The time stepping method provided is not supported by this solver");
    }
  if (this->nsparam.velocitySource.type ==
      Parameters::VelocitySource::VelocitySourceType::srf)
    {
      if (time_stepping_method ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::bdf3,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
        assembleGLS<
          false,
          Parameters::SimulationControl::TimeSteppingMethod::steady_bdf,
          Parameters::VelocitySource::VelocitySourceType::srf>();
      else
        throw std::runtime_error(
          "The time stepping method provided is not supported by this solver");
    }
}

template <int dim>
void
GLSNavierStokesSolver<dim>::solve_linear_system(const bool initial_step,
                                                const bool renewed_matrix)
{
  solve_linear_system_cfd(initial_step, renewed_matrix);
}


template <int dim>
void
GLSNavierStokesSolver<dim>::solve_linear_system_cfd(const bool initial_step,
                                                    const bool renewed_matrix)
{
  const double absolute_residual = this->nsparam.linear_solver.minimum_residual;
  const double relative_residual =
    this->nsparam.linear_solver.relative_residual;

  if (this->nsparam.linear_solver.solver ==
      Parameters::LinearSolver::SolverType::gmres)
    solve_system_GMRES(initial_step,
                       absolute_residual,
                       relative_residual,
                       renewed_matrix);
  else if (this->nsparam.linear_solver.solver ==
           Parameters::LinearSolver::SolverType::bicgstab)
    solve_system_BiCGStab(initial_step,
                          absolute_residual,
                          relative_residual,
                          renewed_matrix);
  else if (this->nsparam.linear_solver.solver ==
           Parameters::LinearSolver::SolverType::amg)
    solve_system_AMG(initial_step,
                     absolute_residual,
                     relative_residual,
                     renewed_matrix);
  else if (this->nsparam.linear_solver.solver ==
           Parameters::LinearSolver::SolverType::tfqmr)
    solve_system_TFQMR(initial_step,
                       absolute_residual,
                       relative_residual,
                       renewed_matrix);
  else if (this->nsparam.linear_solver.solver ==
           Parameters::LinearSolver::SolverType::direct)
    solve_system_direct(initial_step,
                        absolute_residual,
                        relative_residual,
                        renewed_matrix);
  else
    throw(std::runtime_error("This solver is not allowed"));
}

template <int dim>
void
GLSNavierStokesSolver<dim>::setup_ILU()
{
  TimerOutput::Scope t(this->computing_timer, "setup_ILU");

  const double ilu_fill = this->nsparam.linear_solver.ilu_precond_fill;
  const double ilu_atol = this->nsparam.linear_solver.ilu_precond_atol;
  const double ilu_rtol = this->nsparam.linear_solver.ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner->initialize(system_matrix, preconditionerOptions);
}

template <int dim>
void
GLSNavierStokesSolver<dim>::setup_AMG()
{
  TimerOutput::Scope t(this->computing_timer, "setup_AMG");

  std::vector<std::vector<bool>> constant_modes;
  // Constant modes include pressure since everything is in the same matrix
  std::vector<bool> velocity_components(dim + 1, true);
  velocity_components[dim] = true;
  DoFTools::extract_constant_modes(this->dof_handler,
                                   velocity_components,
                                   constant_modes);

  TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
  amg_data.constant_modes = constant_modes;

  const bool elliptic              = false;
  bool       higher_order_elements = false;
  if (this->velocity_fem_degree > 1)
    higher_order_elements = true;
  const unsigned int n_cycles = this->nsparam.linear_solver.amg_n_cycles;
  const bool         w_cycle  = this->nsparam.linear_solver.amg_w_cycles;
  const double       aggregation_threshold =
    this->nsparam.linear_solver.amg_aggregation_threshold;
  const unsigned int smoother_sweeps =
    this->nsparam.linear_solver.amg_smoother_sweeps;
  const unsigned int smoother_overlap =
    this->nsparam.linear_solver.amg_smoother_overlap;
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
  const double ilu_fill = this->nsparam.linear_solver.amg_precond_ilu_fill;
  const double ilu_atol = this->nsparam.linear_solver.amg_precond_ilu_atol;
  const double ilu_rtol = this->nsparam.linear_solver.amg_precond_ilu_rtol;
  parameter_ml.set("smoother: ifpack level-of-fill", ilu_fill);
  parameter_ml.set("smoother: ifpack absolute threshold", ilu_atol);
  parameter_ml.set("smoother: ifpack relative threshold", ilu_rtol);

  parameter_ml.set("coarse: ifpack level-of-fill", ilu_fill);
  parameter_ml.set("coarse: ifpack absolute threshold", ilu_atol);
  parameter_ml.set("coarse: ifpack relative threshold", ilu_rtol);
  amg_preconditioner = std::make_shared<TrilinosWrappers::PreconditionAMG>();
  amg_preconditioner->initialize(system_matrix, parameter_ml);
}

template <int dim>
void
GLSNavierStokesSolver<dim>::solve_system_GMRES(const bool   initial_step,
                                               const double absolute_residual,
                                               const double relative_residual,
                                               const bool   renewed_matrix)
{
  auto &system_rhs          = this->get_system_rhs();
  auto &nonzero_constraints = this->get_nonzero_constraints();

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    false, this->nsparam.linear_solver.max_krylov_vectors);


  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);

  if (renewed_matrix || !ilu_preconditioner)
    setup_ILU();

  {
    TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

    solver.solve(system_matrix,
                 completely_distributed_solution,
                 system_rhs,
                 *ilu_preconditioner);

    if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }
  }
  constraints_used.distribute(completely_distributed_solution);
  auto &newton_update = this->get_newton_update();
  newton_update       = completely_distributed_solution;
}

template <int dim>
void
GLSNavierStokesSolver<dim>::solve_system_BiCGStab(
  const bool   initial_step,
  const double absolute_residual,
  const double relative_residual,
  const bool   renewed_matrix)
{
  TimerOutput::Scope t(this->computing_timer, "solve");
  auto &             system_rhs          = this->get_system_rhs();
  auto &             nonzero_constraints = this->get_nonzero_constraints();

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);
  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);
  TrilinosWrappers::SolverBicgstab solver(solver_control);

  if (renewed_matrix || !ilu_preconditioner)
    setup_ILU();

  {
    TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

    solver.solve(system_matrix,
                 completely_distributed_solution,
                 system_rhs,
                 *ilu_preconditioner);

    if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }
    constraints_used.distribute(completely_distributed_solution);
    auto &newton_update = this->get_newton_update();
    newton_update       = completely_distributed_solution;
  }
}

template <int dim>
void
GLSNavierStokesSolver<dim>::solve_system_AMG(const bool   initial_step,
                                             const double absolute_residual,
                                             const double relative_residual,
                                             const bool   renewed_matrix)
{
  auto &system_rhs          = this->get_system_rhs();
  auto &nonzero_constraints = this->get_nonzero_constraints();

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;

  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);
  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    false, this->nsparam.linear_solver.max_krylov_vectors);

  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);

  if (renewed_matrix || !amg_preconditioner)
    setup_AMG();

  {
    TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

    solver.solve(system_matrix,
                 completely_distributed_solution,
                 system_rhs,
                 *amg_preconditioner);

    if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }

    constraints_used.distribute(completely_distributed_solution);

    auto &newton_update = this->get_newton_update();
    newton_update       = completely_distributed_solution;
  }
}


template <int dim>
void
GLSNavierStokesSolver<dim>::solve_system_direct(const bool   initial_step,
                                                const double absolute_residual,
                                                const double relative_residual,
                                                const bool /*renewed_matrix*/)
{
  auto &system_rhs          = this->get_system_rhs();
  auto &nonzero_constraints = this->get_nonzero_constraints();

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);
  TrilinosWrappers::SolverDirect solver(solver_control);

  solver.initialize(system_matrix);
  solver.solve(completely_distributed_solution, system_rhs);
  constraints_used.distribute(completely_distributed_solution);
  auto &newton_update = this->get_newton_update();
  newton_update       = completely_distributed_solution;
}


template <int dim>
void
GLSNavierStokesSolver<dim>::solve_system_TFQMR(const bool   initial_step,
                                               const double absolute_residual,
                                               const double relative_residual,
                                               const bool   renewed_matrix)
{
  auto &system_rhs          = this->get_system_rhs();
  auto &nonzero_constraints = this->get_nonzero_constraints();

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  TrilinosWrappers::SolverTFQMR solver(solver_control);

  if (renewed_matrix || !ilu_preconditioner)
    setup_ILU();

  {
    TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

    solver.solve(system_matrix,
                 completely_distributed_solution,
                 system_rhs,
                 *ilu_preconditioner);

    if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }
  }
  constraints_used.distribute(completely_distributed_solution);
  auto &newton_update = this->get_newton_update();
  newton_update       = completely_distributed_solution;
}

template <int dim>
void
GLSNavierStokesSolver<dim>::solve()
{
  read_mesh_and_manifolds(this->triangulation,
                          this->nsparam.mesh,
                          this->nsparam.manifolds_parameters,
                          this->nsparam.restart_parameters.restart,
                          this->nsparam.boundary_conditions);

  this->setup_dofs();
  this->set_initial_condition(this->nsparam.initial_condition->type,
                              this->nsparam.restart_parameters.restart);

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);
      this->dynamic_flow_control();

      if (this->simulation_control->is_at_start())
        this->first_iteration();
      else
        {
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
            refine_mesh();
          this->iterate();
        }
      this->postprocess(false);
      this->finish_time_step();
    }


  this->finish_simulation();
}


// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is
// valid before we actually compile the solver This greatly helps with debugging
template class GLSNavierStokesSolver<2>;
template class GLSNavierStokesSolver<3>;
