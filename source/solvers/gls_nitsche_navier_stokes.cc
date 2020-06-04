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
 * Author: Bruno Blais, Carole-Anne Daunais, Valérie Bibeau, Polytechnique Montreal, 2019-
 */

#include "solvers/gls_nitsche_navier_stokes.h"

#include "core/bdf.h"
#include "core/grids.h"
#include "core/manifolds.h"
#include "core/sdirk.h"
#include "core/time_integration_utilities.h"

// Constructor for class GLSNitscheNavierStokesSolver
template <int dim>
GLSNitscheNavierStokesSolver<dim>::GLSNitscheNavierStokesSolver(
  NavierStokesSolverParameters<dim> &p_nsparam,
  const unsigned int                 p_degreeVelocity,
  const unsigned int                 p_degreePressure)
  : NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>(
      p_nsparam,
      p_degreeVelocity,
      p_degreePressure)
{}

template <int dim>
GLSNitscheNavierStokesSolver<dim>::~GLSNitscheNavierStokesSolver()
{
  this->dof_handler.clear();
}



template <int dim>
void
GLSNitscheNavierStokesSolver<dim>::set_solution_vector(double value)
{
  this->present_solution = value;
}

template <int dim>
void
GLSNitscheNavierStokesSolver<dim>::setup_dofs()
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

  const MappingQ<dim>        mapping(this->degreeVelocity_,
                              this->nsparam.fem_parameters.qmapping_all);
  FEValuesExtractors::Vector velocities(0);

  // Non-zero constraints
  {
    this->nonzero_constraints.clear();

    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            this->nonzero_constraints);
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
              ZeroFunction<dim>(dim + 1),
              this->nonzero_constraints,
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
              this->nonzero_constraints);
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
              this->nonzero_constraints,
              this->fe.component_mask(velocities));
          }

        else if (this->nsparam.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints<DoFHandler<dim>>(
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              this->nsparam.boundary_conditions.periodic_id[i_bc],
              this->nsparam.boundary_conditions.periodic_direction[i_bc],
              this->nonzero_constraints);
          }
      }
  }
  this->nonzero_constraints.close();

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
            DoFTools::make_periodicity_constraints<DoFHandler<dim>>(
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
              ZeroFunction<dim>(dim + 1),
              this->zero_constraints,
              this->fe.component_mask(velocities));
          }
      }
  }
  this->zero_constraints.close();

  this->present_solution.reinit(this->locally_owned_dofs,
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

  this->newton_update.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->system_rhs.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->local_evaluation_point.reinit(this->locally_owned_dofs,
                                      this->mpi_communicator);

  DynamicSparsityPattern dsp(this->locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(this->dof_handler,
                                  dsp,
                                  this->nonzero_constraints,
                                  false);
#if !(DEAL_II_VERSION_GTE(9, 2, 0))
  SparsityTools::distribute_sparsity_pattern(
    dsp,
    this->dof_handler.compute_n_locally_owned_dofs_per_processor(),
    this->mpi_communicator,
    this->locally_relevant_dofs);
  system_matrix.reinit(this->locally_owned_dofs,
                       this->locally_owned_dofs,
                       dsp,
                       this->mpi_communicator);
#else
  SparsityTools::distribute_sparsity_pattern(
    dsp,
    this->dof_handler.locally_owned_dofs(),
    this->mpi_communicator,
    this->locally_relevant_dofs);
  system_matrix.reinit(this->locally_owned_dofs,
                       this->locally_owned_dofs,
                       dsp,
                       this->mpi_communicator);

#endif

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
GLSNitscheNavierStokesSolver<dim>::assembleGLS()
{
  if (assemble_matrix)
    system_matrix = 0;
  this->system_rhs = 0;

  double         viscosity_ = this->nsparam.physical_properties.viscosity;
  Function<dim> *l_forcing_function = this->forcing_function;

  QGauss<dim>                      quadrature_formula(this->degreeQuadrature_);
  const MappingQ<dim>              mapping(this->degreeVelocity_,
                              this->nsparam.fem_parameters.qmapping_all);
  FEValues<dim>                    fe_values(mapping,
                          this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients |
                            update_hessians);
  const unsigned int               dofs_per_cell = this->fe.dofs_per_cell;
  const unsigned int               n_q_points    = quadrature_formula.size();
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
    this->simulationControl->get_time_steps_vector();

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

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf1)
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
          fe_values.reinit(cell);

          if (dim == 2)
            h = std::sqrt(4. * cell->measure() / M_PI) / this->degreeVelocity_;
          else if (dim == 3)
            h =
              pow(6 * cell->measure() / M_PI, 1. / 3.) / this->degreeVelocity_;

          local_matrix = 0;
          local_rhs    = 0;

          // Gather velocity (values, gradient and laplacian)
          fe_values[velocities].get_function_values(this->evaluation_point,
                                                    present_velocity_values);
          fe_values[velocities].get_function_gradients(
            this->evaluation_point, present_velocity_gradients);
          fe_values[velocities].get_function_laplacians(
            this->evaluation_point, present_velocity_laplacians);

          // Gather pressure (values, gradient)
          fe_values[pressure].get_function_values(this->evaluation_point,
                                                  present_pressure_values);
          fe_values[pressure].get_function_gradients(
            this->evaluation_point, present_pressure_gradients);

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
              // Calculation of the magnitude of the velocity for the
              // stabilization parameter
              const double u_mag = std::max(present_velocity_values[q].norm(),
                                            1e-12 * GLS_u_scale);

              // Calculation of the GLS stabilization parameter. The
              // stabilization parameter used is different if the simulation is
              // steady or unsteady. In the unsteady case it includes the value
              // of the time-step
              double tau;
              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::steady)
                tau = 1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                                     9 * std::pow(4 * viscosity_ / (h * h), 2));
              else
                tau = 1. /
                      std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                                9 * std::pow(4 * viscosity_ / (h * h), 2));

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

              // Calculate the divergence of the velocity
              double present_velocity_divergence =
                trace(present_velocity_gradients[q]);

              // Calculate the strong residual for GLS stabilization
              auto strong_residual =
                present_velocity_gradients[q] * present_velocity_values[q] +
                present_pressure_gradients[q] -
                viscosity_ * present_velocity_laplacians[q] - force;

              if (velocity_source ==
                  Parameters::VelocitySource::VelocitySourceType::srf)
                {
                  if (dim == 2)
                    {
                      strong_residual +=
                        2 * omega_z * (-1.) *
                        cross_product_2d(present_velocity_values[q]);
                      auto centrifugal =
                        omega_z * (-1.) *
                        cross_product_2d(
                          omega_z * (-1.) *
                          cross_product_2d(quadrature_points[q]));
                      strong_residual += centrifugal;
                    }
                  else if (dim == 3)
                    {
                      strong_residual +=
                        2 * cross_product_3d(omega_vector,
                                             present_velocity_values[q]);
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
                  Parameters::SimulationControl::TimeSteppingMethod::bdf1)
                strong_residual += bdf_coefs[0] * present_velocity_values[q] +
                                   bdf_coefs[1] * p1_velocity_values[q];

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                strong_residual += bdf_coefs[0] * present_velocity_values[q] +
                                   bdf_coefs[1] * p1_velocity_values[q] +
                                   bdf_coefs[2] * p2_velocity_values[q];

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                strong_residual += bdf_coefs[0] * present_velocity_values[q] +
                                   bdf_coefs[1] * p1_velocity_values[q] +
                                   bdf_coefs[2] * p2_velocity_values[q] +
                                   bdf_coefs[3] * p3_velocity_values[q];


              if (is_sdirk_step1(scheme))
                strong_residual +=
                  sdirk_coefs[0][0] * present_velocity_values[q] +
                  sdirk_coefs[0][1] * p1_velocity_values[q];

              if (is_sdirk_step2(scheme))
                {
                  strong_residual +=
                    sdirk_coefs[1][0] * present_velocity_values[q] +
                    sdirk_coefs[1][1] * p1_velocity_values[q] +
                    sdirk_coefs[1][2] * p2_velocity_values[q];
                }

              if (is_sdirk_step3(scheme))
                {
                  strong_residual +=
                    sdirk_coefs[2][0] * present_velocity_values[q] +
                    sdirk_coefs[2][1] * p1_velocity_values[q] +
                    sdirk_coefs[2][2] * p2_velocity_values[q] +
                    sdirk_coefs[2][3] * p3_velocity_values[q];
                }

              // Matrix assembly
              if (assemble_matrix)
                {
                  // We loop over the column first to prevent recalculation of
                  // the strong jacobian in the inner loop
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      auto strong_jac =
                        (present_velocity_gradients[q] * phi_u[j] +
                         grad_phi_u[j] * present_velocity_values[q] +
                         grad_phi_p[j] - viscosity_ * laplacian_phi_u[j]);

                      if (is_bdf(scheme))
                        strong_jac += phi_u[j] * bdf_coefs[0];
                      if (is_sdirk(scheme))
                        strong_jac += phi_u[j] * sdirk_coefs[0][0];

                      if (velocity_source ==
                          Parameters::VelocitySource::VelocitySourceType::srf)
                        {
                          if (dim == 2)
                            strong_jac +=
                              2 * omega_z * (-1.) * cross_product_2d(phi_u[j]);
                          else if (dim == 3)
                            strong_jac +=
                              2 * cross_product_3d(omega_vector, phi_u[j]);
                        }

                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                          local_matrix(i, j) +=
                            (
                              // Momentum terms
                              viscosity_ *
                                scalar_product(grad_phi_u[j], grad_phi_u[i]) +
                              present_velocity_gradients[q] * phi_u[j] *
                                phi_u[i] +
                              grad_phi_u[j] * present_velocity_values[q] *
                                phi_u[i] -
                              div_phi_u[i] * phi_p[j] +
                              // Continuity
                              phi_p[i] * div_phi_u[j]) *
                            fe_values.JxW(q);

                          // Mass matrix
                          if (is_bdf(scheme))
                            local_matrix(i, j) += phi_u[j] * phi_u[i] *
                                                  bdf_coefs[0] *
                                                  fe_values.JxW(q);

                          if (is_sdirk(scheme))
                            local_matrix(i, j) += phi_u[j] * phi_u[i] *
                                                  sdirk_coefs[0][0] *
                                                  fe_values.JxW(q);

                          // PSPG GLS term
                          local_matrix(i, j) +=
                            tau * strong_jac * grad_phi_p[i] * fe_values.JxW(q);

                          if (velocity_source == Parameters::VelocitySource::
                                                   VelocitySourceType::srf)
                            {
                              if (dim == 2)
                                local_matrix(i, j) +=
                                  2 * omega_z * (-1.) *
                                  cross_product_2d(phi_u[j]) * phi_u[i] *
                                  fe_values.JxW(q);

                              else if (dim == 3)
                                local_matrix(i, j) +=
                                  2 * cross_product_3d(omega_vector, phi_u[j]) *
                                  phi_u[i] * fe_values.JxW(q);
                            }


                          // PSPG TAU term is currently disabled because it does
                          // not alter the matrix sufficiently
                          // local_matrix(i, j) +=
                          //  -tau * tau * tau * 4 / h / h *
                          //  (present_velocity_values[q] * phi_u[j]) *
                          //  strong_residual * grad_phi_p[i] *
                          //  fe_values.JxW(q);

                          // Jacobian is currently incomplete
                          if (SUPG)
                            {
                              local_matrix(i, j) +=
                                tau *
                                (strong_jac * (grad_phi_u[i] *
                                               present_velocity_values[q]) +
                                 strong_residual * (grad_phi_u[i] * phi_u[j])) *
                                fe_values.JxW(q);

                              // SUPG TAU term is currently disabled because it
                              // does not alter the matrix sufficiently
                              // local_matrix(i, j)
                              // +=
                              //   -strong_residual
                              //   * (grad_phi_u[i]
                              //   *
                              //   present_velocity_values[q])
                              //   * tau * tau *
                              //   tau * 4 / h / h
                              //   *
                              //   (present_velocity_values[q]
                              //   * phi_u[j]) *
                              //   fe_values.JxW(q);
                            }
                        }
                    }
                }

              // Assembly of the right-hand side
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Navier-Stokes Residual
                  local_rhs(i) +=
                    (
                      // Momentum
                      -viscosity_ *
                        scalar_product(present_velocity_gradients[q],
                                       grad_phi_u[i]) -
                      present_velocity_gradients[q] *
                        present_velocity_values[q] * phi_u[i] +
                      present_pressure_values[q] * div_phi_u[i] +
                      force * phi_u[i] -
                      // Continuity
                      present_velocity_divergence * phi_p[i]) *
                    fe_values.JxW(q);

                  // Residual associated with BDF schemes
                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf1)
                    local_rhs(i) -=
                      bdf_coefs[0] *
                      (present_velocity_values[q] - p1_velocity_values[q]) *
                      phi_u[i] * fe_values.JxW(q);

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                    local_rhs(i) -=
                      (bdf_coefs[0] * (present_velocity_values[q] * phi_u[i]) +
                       bdf_coefs[1] * (p1_velocity_values[q] * phi_u[i]) +
                       bdf_coefs[2] * (p2_velocity_values[q] * phi_u[i])) *
                      fe_values.JxW(q);

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                    local_rhs(i) -=
                      (bdf_coefs[0] * (present_velocity_values[q] * phi_u[i]) +
                       bdf_coefs[1] * (p1_velocity_values[q] * phi_u[i]) +
                       bdf_coefs[2] * (p2_velocity_values[q] * phi_u[i]) +
                       bdf_coefs[3] * (p3_velocity_values[q] * phi_u[i])) *
                      fe_values.JxW(q);

                  // Residuals associated with SDIRK schemes
                  if (is_sdirk_step1(scheme))
                    local_rhs(i) -=
                      (sdirk_coefs[0][0] *
                         (present_velocity_values[q] * phi_u[i]) +
                       sdirk_coefs[0][1] * (p1_velocity_values[q] * phi_u[i])) *
                      fe_values.JxW(q);

                  if (is_sdirk_step2(scheme))
                    {
                      local_rhs(i) -=
                        (sdirk_coefs[1][0] *
                           (present_velocity_values[q] * phi_u[i]) +
                         sdirk_coefs[1][1] *
                           (p1_velocity_values[q] * phi_u[i]) +
                         sdirk_coefs[1][2] *
                           (p2_velocity_values[q] * phi_u[i])) *
                        fe_values.JxW(q);
                    }

                  if (is_sdirk_step3(scheme))
                    {
                      local_rhs(i) -=
                        (sdirk_coefs[2][0] *
                           (present_velocity_values[q] * phi_u[i]) +
                         sdirk_coefs[2][1] *
                           (p1_velocity_values[q] * phi_u[i]) +
                         sdirk_coefs[2][2] *
                           (p2_velocity_values[q] * phi_u[i]) +
                         sdirk_coefs[2][3] *
                           (p3_velocity_values[q] * phi_u[i])) *
                        fe_values.JxW(q);
                    }

                  if (velocity_source ==
                      Parameters::VelocitySource::VelocitySourceType::srf)
                    {
                      if (dim == 2)
                        {
                          local_rhs(i) +=
                            -2 * omega_z * (-1.) *
                            cross_product_2d(present_velocity_values[q]) *
                            phi_u[i] * fe_values.JxW(q);
                          auto centrifugal =
                            omega_z * (-1.) *
                            cross_product_2d(
                              omega_z * (-1.) *
                              cross_product_2d(quadrature_points[q]));
                          local_rhs(i) +=
                            -centrifugal * phi_u[i] * fe_values.JxW(q);
                        }
                      else if (dim == 3)
                        {
                          local_rhs(i) +=
                            -2 *
                            cross_product_3d(omega_vector,
                                             present_velocity_values[q]) *
                            phi_u[i] * fe_values.JxW(q);
                          local_rhs(i) +=
                            -cross_product_3d(
                              omega_vector,
                              cross_product_3d(omega_vector,
                                               quadrature_points[q])) *
                            phi_u[i] * fe_values.JxW(q);
                        }
                    }

                  // PSPG GLS term
                  local_rhs(i) +=
                    -tau * (strong_residual * grad_phi_p[i]) * fe_values.JxW(q);

                  // SUPG GLS term
                  if (SUPG)
                    {
                      local_rhs(i) +=
                        -tau *
                        (strong_residual *
                         (grad_phi_u[i] * present_velocity_values[q])) *
                        fe_values.JxW(q);
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
                                                          this->system_rhs);
            }
          else
            {
              constraints_used.distribute_local_to_global(local_rhs,
                                                          local_dof_indices,
                                                          this->system_rhs);
            }
        }
    }
  if (assemble_matrix)
    system_matrix.compress(VectorOperation::add);
  this->system_rhs.compress(VectorOperation::add);
}

/**
 * Set the initial condition using a L2 or a viscous solver
 **/
template <int dim>
void
GLSNitscheNavierStokesSolver<dim>::set_initial_condition(
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
      this->present_solution = this->newton_update;
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
GLSNitscheNavierStokesSolver<dim>::assemble_L2_projection()
{
  system_matrix    = 0;
  this->system_rhs = 0;
  QGauss<dim>                 quadrature_formula(this->degreeQuadrature_);
  const MappingQ<dim>         mapping(this->degreeVelocity_,
                              this->nsparam.fem_parameters.qmapping_all);
  FEValues<dim>               fe_values(mapping,
                          this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  const unsigned int          dofs_per_cell = this->fe.dofs_per_cell;
  const unsigned int          n_q_points    = quadrature_formula.size();
  FullMatrix<double>          local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>              local_rhs(dofs_per_cell);
  std::vector<Vector<double>> initial_velocity(n_q_points,
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
          const AffineConstraints<double> &constraints_used =
            this->nonzero_constraints;
          constraints_used.distribute_local_to_global(local_matrix,
                                                      local_rhs,
                                                      local_dof_indices,
                                                      system_matrix,
                                                      this->system_rhs);
        }
    }
  system_matrix.compress(VectorOperation::add);
  this->system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
GLSNitscheNavierStokesSolver<dim>::assemble_matrix_and_rhs(
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
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::none>();
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
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<true,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
    }
}
template <int dim>
void
GLSNitscheNavierStokesSolver<dim>::assemble_rhs(
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
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
                    Parameters::VelocitySource::VelocitySourceType::none>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::none>();
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
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
      else if (time_stepping_method ==
               Parameters::SimulationControl::TimeSteppingMethod::steady)
        assembleGLS<false,
                    Parameters::SimulationControl::TimeSteppingMethod::steady,
                    Parameters::VelocitySource::VelocitySourceType::srf>();
    }
}

template <int dim>
void
GLSNitscheNavierStokesSolver<dim>::solve_linear_system(const bool initial_step,
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
  else
    throw(std::runtime_error("This solver is not allowed"));
}

template <int dim>
void
GLSNitscheNavierStokesSolver<dim>::setup_ILU()
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
GLSNitscheNavierStokesSolver<dim>::setup_AMG()
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
  if (this->degreeVelocity_ > 1)
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
GLSNitscheNavierStokesSolver<dim>::solve_system_GMRES(const bool   initial_step,
                                               const double absolute_residual,
                                               const double relative_residual,
                                               const bool   renewed_matrix)
{
  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);

  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << std::setprecision(
                       this->nsparam.linear_solver.residual_precision)
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);
  TrilinosWrappers::SolverGMRES solver(solver_control);

  if (renewed_matrix || !ilu_preconditioner)
    setup_ILU();

  {
    TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

    solver.solve(system_matrix,
                 completely_distributed_solution,
                 this->system_rhs,
                 *ilu_preconditioner);

    if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }
  }
  constraints_used.distribute(completely_distributed_solution);
  this->newton_update = completely_distributed_solution;
}

template <int dim>
void
GLSNitscheNavierStokesSolver<dim>::solve_system_BiCGStab(
  const bool   initial_step,
  const double absolute_residual,
  const double relative_residual,
  const bool   renewed_matrix)
{
  TimerOutput::Scope t(this->computing_timer, "solve");

  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);
  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << std::setprecision(
                       this->nsparam.linear_solver.residual_precision)
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
                 this->system_rhs,
                 *ilu_preconditioner);

    if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }
    constraints_used.distribute(completely_distributed_solution);
    this->newton_update = completely_distributed_solution;
  }
}

template <int dim>
void
GLSNitscheNavierStokesSolver<dim>::solve_system_AMG(const bool   initial_step,
                                             const double absolute_residual,
                                             const double relative_residual,
                                             const bool   renewed_matrix)
{
  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;

  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);
  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << std::setprecision(
                       this->nsparam.linear_solver.residual_precision)
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);
  TrilinosWrappers::SolverGMRES solver(solver_control);

  if (renewed_matrix || !amg_preconditioner)
    setup_AMG();

  {
    TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

    solver.solve(system_matrix,
                 completely_distributed_solution,
                 this->system_rhs,
                 *amg_preconditioner);

    if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }

    constraints_used.distribute(completely_distributed_solution);

    this->newton_update = completely_distributed_solution;
  }
}

template <int dim>
void
GLSNitscheNavierStokesSolver<dim>::solve()
{
  read_mesh_and_manifolds(this->triangulation,
                          this->nsparam.mesh,
                          this->nsparam.manifolds_parameters,
                          this->nsparam.boundary_conditions);

  this->setup_dofs();
  this->set_initial_condition(this->nsparam.initial_condition->type,
                              this->nsparam.restart_parameters.restart);

  while (this->simulationControl->integrate())
    {
      this->simulationControl->print_progression(this->pcout);
      if (this->simulationControl->is_at_start())
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
template class GLSNitscheNavierStokesSolver<2>;
template class GLSNitscheNavierStokesSolver<3>;
