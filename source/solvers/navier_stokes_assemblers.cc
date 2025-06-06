// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/simulation_control.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/navier_stokes_assemblers.h>
#include <solvers/stabilization.h>

template <int dim>
void
PSPGSUPGNavierStokesAssemblerCore<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &viscosity_vector =
    scratch_data.kinematic_viscosity;

  const std::vector<double> &viscosity_for_stabilization_vector =
    scratch_data.kinematic_viscosity_for_stabilization;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix        = copy_data.local_matrix;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Pressure scaling factor
  const double pressure_scaling_factor = scratch_data.pressure_scaling_factor;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const double          kinematic_viscosity = viscosity_vector[q];
      const Tensor<1, dim> &velocity = scratch_data.velocity_values[q];
      const Tensor<2, dim> &velocity_gradient =
        scratch_data.velocity_gradients[q];
      const Tensor<1, dim> &velocity_laplacian =
        scratch_data.velocity_laplacians[q];

      const Tensor<1, dim> &pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Forcing term
      const Tensor<1, dim> &force       = scratch_data.force[q];
      double                mass_source = scratch_data.mass_source[q];

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(
            u_mag, viscosity_for_stabilization_vector[q], h) :
          calculate_navier_stokes_gls_tau_transient(
            u_mag, viscosity_for_stabilization_vector[q], h, sdt);

      // Calculate the strong residual for GLS stabilization
      auto strong_residual = velocity_gradient * velocity + pressure_gradient -
                             kinematic_viscosity * velocity_laplacian - force +
                             mass_source * velocity + strong_residual_vec[q];

      std::vector<Tensor<1, dim>> grad_phi_u_j_x_velocity(n_dofs);
      std::vector<Tensor<1, dim>> velocity_gradient_x_phi_u_j(n_dofs);


      // We loop over the column first to prevent recalculation
      // of the strong jacobian in the inner loop
      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto &phi_u_j           = scratch_data.phi_u[q][j];
          const auto &grad_phi_u_j      = scratch_data.grad_phi_u[q][j];
          const auto &laplacian_phi_u_j = scratch_data.laplacian_phi_u[q][j];

          const auto &grad_phi_p_j =
            pressure_scaling_factor * scratch_data.grad_phi_p[q][j];

          strong_jacobian_vec[q][j] +=
            (velocity_gradient * phi_u_j + grad_phi_u_j * velocity +
             grad_phi_p_j - kinematic_viscosity * laplacian_phi_u_j +
             mass_source * phi_u_j);

          // Store these temporary products in auxiliary variables for speed
          grad_phi_u_j_x_velocity[j]     = grad_phi_u_j * velocity;
          velocity_gradient_x_phi_u_j[j] = velocity_gradient * phi_u_j;
        }



      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const unsigned int component_i = scratch_data.components[i];

          const auto &phi_u_i      = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &div_phi_u_i  = scratch_data.div_phi_u[q][i];
          const auto &phi_p_i      = scratch_data.phi_p[q][i];
          const auto &grad_phi_p_i = scratch_data.grad_phi_p[q][i];


          // Store these temporary products in auxiliary variables for speed
          const auto grad_phi_u_i_x_velocity = grad_phi_u_i * velocity;
          const auto strong_residual_x_grad_phi_u_i =
            strong_residual * grad_phi_u_i;

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const unsigned int component_j = scratch_data.components[j];

              const auto &phi_u_j = scratch_data.phi_u[q][j];

              const auto &phi_p_j =
                pressure_scaling_factor * scratch_data.phi_p[q][j];

              const auto &strong_jac = strong_jacobian_vec[q][j];

              double local_matrix_ij =
                component_j == dim ? -div_phi_u_i * phi_p_j : 0;
              if (component_i == dim)
                {
                  const auto &div_phi_u_j = scratch_data.div_phi_u[q][j];

                  local_matrix_ij += phi_p_i * div_phi_u_j;

                  // PSPG GLS term
                  local_matrix_ij += tau * (strong_jac * grad_phi_p_i);
                }

              if (component_i < dim && component_j < dim)
                {
                  const auto &grad_phi_u_j = scratch_data.grad_phi_u[q][j];

                  local_matrix_ij += velocity_gradient_x_phi_u_j[j] * phi_u_i +
                                     grad_phi_u_j_x_velocity[j] * phi_u_i;

                  if (component_i == component_j)
                    {
                      local_matrix_ij +=
                        kinematic_viscosity * (grad_phi_u_j[component_j] *
                                               grad_phi_u_i[component_i]) +
                        mass_source * phi_u_j * phi_u_i;
                    }
                }
              if (component_i < dim)
                {
                  // The jacobian matrix for the SUPG formulation
                  // currently does not include the jacobian of the
                  // stabilization parameter tau. Our experience has shown that
                  // does not alter the number of newton iteration for
                  // convergence, but greatly simplifies assembly.

                  local_matrix_ij +=
                    tau * (strong_jac * grad_phi_u_i_x_velocity +
                           strong_residual_x_grad_phi_u_i * phi_u_j);
                }

              local_matrix_ij *= JxW;
              local_matrix(i, j) += local_matrix_ij;
            }
        }
    }
}

template <int dim>
void
PSPGSUPGNavierStokesAssemblerCore<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &viscosity_vector =
    scratch_data.kinematic_viscosity;

  const std::vector<double> &viscosity_for_stabilization_vector =
    scratch_data.kinematic_viscosity_for_stabilization;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs           = copy_data.local_rhs;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Physical properties
      const double kinematic_viscosity = viscosity_vector[q];

      // Velocity
      const Tensor<1, dim> &velocity   = scratch_data.velocity_values[q];
      const double velocity_divergence = scratch_data.velocity_divergences[q];
      const Tensor<2, dim> &velocity_gradient =
        scratch_data.velocity_gradients[q];
      const Tensor<1, dim> &velocity_laplacian =
        scratch_data.velocity_laplacians[q];

      // Pressure
      const double          pressure = scratch_data.pressure_values[q];
      const Tensor<1, dim> &pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Forcing term
      const Tensor<1, dim> &force       = scratch_data.force[q];
      double                mass_source = scratch_data.mass_source[q];
      // Calculation of the magnitude of the
      // velocity for the stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(
            u_mag, viscosity_for_stabilization_vector[q], h) :
          calculate_navier_stokes_gls_tau_transient(
            u_mag, viscosity_for_stabilization_vector[q], h, sdt);



      // Calculate the strong residual for GLS stabilization
      auto strong_residual = velocity_gradient * velocity + pressure_gradient -
                             kinematic_viscosity * velocity_laplacian - force +
                             mass_source * velocity + strong_residual_vec[q];

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_u_i      = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &phi_p_i      = scratch_data.phi_p[q][i];
          const auto &grad_phi_p_i = scratch_data.grad_phi_p[q][i];
          const auto &div_phi_u_i  = scratch_data.div_phi_u[q][i];

          double local_rhs_i = 0;

          // Navier-Stokes Residual
          local_rhs_i +=
            (
              // Momentum
              -kinematic_viscosity *
                scalar_product(velocity_gradient, grad_phi_u_i) -
              velocity_gradient * velocity * phi_u_i + pressure * div_phi_u_i +
              force * phi_u_i - mass_source * velocity * phi_u_i -
              // Continuity
              velocity_divergence * phi_p_i + mass_source * phi_p_i) *
            JxW;

          // PSPG GLS term
          local_rhs_i += -tau * (strong_residual * grad_phi_p_i) * JxW;

          // SUPG GLS term
          local_rhs_i +=
            -tau * (strong_residual * (grad_phi_u_i * velocity)) * JxW;

          local_rhs(i) += local_rhs_i;
        }
    }
}



template class PSPGSUPGNavierStokesAssemblerCore<2>;
template class PSPGSUPGNavierStokesAssemblerCore<3>;

template <int dim>
void
GLSNavierStokesAssemblerCore<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &viscosity_vector =
    scratch_data.kinematic_viscosity;

  const std::vector<double> &viscosity_for_stabilization_vector =
    scratch_data.kinematic_viscosity_for_stabilization;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix        = copy_data.local_matrix;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Pressure scaling factor
  const double pressure_scaling_factor = scratch_data.pressure_scaling_factor;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const double          kinematic_viscosity = viscosity_vector[q];
      const Tensor<1, dim> &velocity = scratch_data.velocity_values[q];
      const Tensor<2, dim> &velocity_gradient =
        scratch_data.velocity_gradients[q];
      const Tensor<1, dim> &velocity_laplacian =
        scratch_data.velocity_laplacians[q];

      const Tensor<1, dim> &pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Forcing term
      const Tensor<1, dim> &force       = scratch_data.force[q];
      double                mass_source = scratch_data.mass_source[q];

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(
            u_mag, viscosity_for_stabilization_vector[q], h) :
          calculate_navier_stokes_gls_tau_transient(
            u_mag, viscosity_for_stabilization_vector[q], h, sdt);

      // LSIC stabilization term
      const double tau_lsic = 0.5 * u_mag * h;


      // Calculate the strong residual for GLS stabilization
      auto strong_residual = velocity_gradient * velocity + pressure_gradient -
                             kinematic_viscosity * velocity_laplacian - force +
                             mass_source * velocity + strong_residual_vec[q];

      std::vector<Tensor<1, dim>> grad_phi_u_j_x_velocity(n_dofs);
      std::vector<Tensor<1, dim>> velocity_gradient_x_phi_u_j(n_dofs);


      // We loop over the column first to prevent recalculation
      // of the strong jacobian in the inner loop
      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto &phi_u_j           = scratch_data.phi_u[q][j];
          const auto &grad_phi_u_j      = scratch_data.grad_phi_u[q][j];
          const auto &laplacian_phi_u_j = scratch_data.laplacian_phi_u[q][j];

          const auto &grad_phi_p_j =
            pressure_scaling_factor * scratch_data.grad_phi_p[q][j];

          strong_jacobian_vec[q][j] +=
            (velocity_gradient * phi_u_j + grad_phi_u_j * velocity +
             grad_phi_p_j - kinematic_viscosity * laplacian_phi_u_j +
             mass_source * phi_u_j);

          // Store these temporary products in auxiliary variables for speed
          grad_phi_u_j_x_velocity[j]     = grad_phi_u_j * velocity;
          velocity_gradient_x_phi_u_j[j] = velocity_gradient * phi_u_j;
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const unsigned int component_i = scratch_data.components[i];

          const auto &phi_u_i           = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i      = scratch_data.grad_phi_u[q][i];
          const auto &div_phi_u_i       = scratch_data.div_phi_u[q][i];
          const auto &phi_p_i           = scratch_data.phi_p[q][i];
          const auto &grad_phi_p_i      = scratch_data.grad_phi_p[q][i];
          const auto &laplacian_phi_u_i = scratch_data.laplacian_phi_u[q][i];



          // Store these temporary products in auxiliary variables for speed
          const auto grad_phi_u_i_x_velocity = grad_phi_u_i * velocity;
          const auto strong_residual_x_grad_phi_u_i =
            strong_residual * grad_phi_u_i;

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const unsigned int component_j = scratch_data.components[j];

              const auto &phi_u_j = scratch_data.phi_u[q][j];

              const auto &phi_p_j =
                pressure_scaling_factor * scratch_data.phi_p[q][j];

              const auto &strong_jac = strong_jacobian_vec[q][j];

              double local_matrix_ij =
                component_j == dim ? -div_phi_u_i * phi_p_j : 0;
              if (component_i == dim)
                {
                  const auto &div_phi_u_j = scratch_data.div_phi_u[q][j];
                  local_matrix_ij += phi_p_i * div_phi_u_j; // continuity

                  // PSPG GLS term
                  local_matrix_ij += tau * (strong_jac * grad_phi_p_i);
                }

              if (component_i < dim && component_j < dim)
                {
                  const auto &grad_phi_u_j = scratch_data.grad_phi_u[q][j];

                  local_matrix_ij += velocity_gradient_x_phi_u_j[j] * phi_u_i +
                                     grad_phi_u_j_x_velocity[j] * phi_u_i;

                  // LSIC GLS term
                  const auto &div_phi_u_j = scratch_data.div_phi_u[q][j];
                  local_matrix_ij += tau_lsic * (div_phi_u_i * div_phi_u_j);

                  if (component_i == component_j)
                    {
                      local_matrix_ij +=
                        kinematic_viscosity * (grad_phi_u_j[component_j] *
                                               grad_phi_u_i[component_i]) +
                        mass_source * phi_u_j * phi_u_i;
                    }
                }
              if (component_i < dim)
                {
                  // The jacobian matrix for the SUPG formulation
                  // currently does not include the jacobian of the
                  // stabilization parameter tau. Our experience has shown that
                  // does not alter the number of newton iteration for
                  // convergence, but greatly simplifies assembly.
                  local_matrix_ij +=
                    tau *
                    (strong_jac * (grad_phi_u_i_x_velocity -
                                   kinematic_viscosity * laplacian_phi_u_i) +
                     strong_residual_x_grad_phi_u_i * phi_u_j);
                }

              local_matrix_ij *= JxW;
              local_matrix(i, j) += local_matrix_ij;
            }
        }
    }
}

template <int dim>
void
GLSNavierStokesAssemblerCore<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &viscosity_vector =
    scratch_data.kinematic_viscosity;

  const std::vector<double> &viscosity_for_stabilization_vector =
    scratch_data.kinematic_viscosity_for_stabilization;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs           = copy_data.local_rhs;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Physical properties
      const double kinematic_viscosity = viscosity_vector[q];

      // Velocity
      const Tensor<1, dim> &velocity   = scratch_data.velocity_values[q];
      const double velocity_divergence = scratch_data.velocity_divergences[q];
      const Tensor<2, dim> &velocity_gradient =
        scratch_data.velocity_gradients[q];
      const Tensor<1, dim> &velocity_laplacian =
        scratch_data.velocity_laplacians[q];

      // Pressure
      const double          pressure = scratch_data.pressure_values[q];
      const Tensor<1, dim> &pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Forcing term
      const Tensor<1, dim> &force       = scratch_data.force[q];
      double                mass_source = scratch_data.mass_source[q];
      // Calculation of the magnitude of the
      // velocity for the stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(
            u_mag, viscosity_for_stabilization_vector[q], h) :
          calculate_navier_stokes_gls_tau_transient(
            u_mag, viscosity_for_stabilization_vector[q], h, sdt);


      // LSIC stabilization term
      const double tau_lsic = u_mag * h / 2;


      // Calculate the strong residual for GLS stabilization
      auto strong_residual = velocity_gradient * velocity + pressure_gradient -
                             kinematic_viscosity * velocity_laplacian - force +
                             mass_source * velocity + strong_residual_vec[q];

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_u_i           = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i      = scratch_data.grad_phi_u[q][i];
          const auto &phi_p_i           = scratch_data.phi_p[q][i];
          const auto &grad_phi_p_i      = scratch_data.grad_phi_p[q][i];
          const auto &div_phi_u_i       = scratch_data.div_phi_u[q][i];
          const auto &laplacian_phi_u_i = scratch_data.laplacian_phi_u[q][i];


          double local_rhs_i = 0;

          // Navier-Stokes Residual
          local_rhs_i +=
            (
              // Momentum
              -kinematic_viscosity *
                scalar_product(velocity_gradient, grad_phi_u_i) -
              velocity_gradient * velocity * phi_u_i + pressure * div_phi_u_i +
              force * phi_u_i - mass_source * velocity * phi_u_i -
              // Continuity
              velocity_divergence * phi_p_i + mass_source * phi_p_i) *
            JxW;

          // PSPG GLS term
          local_rhs_i += -tau * (strong_residual * grad_phi_p_i) * JxW;

          // LSIC GLS term
          local_rhs_i += -tau_lsic * (div_phi_u_i * velocity_divergence) * JxW;

          // SUPG GLS term
          local_rhs_i +=
            -tau *
            (strong_residual * (grad_phi_u_i * velocity -
                                kinematic_viscosity * laplacian_phi_u_i)) *
            JxW;

          local_rhs(i) += local_rhs_i;
        }
    }
}



template class GLSNavierStokesAssemblerCore<2>;
template class GLSNavierStokesAssemblerCore<3>;

template <int dim>
void
GLSNavierStokesAssemblerNonNewtonianCore<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &viscosity_vector =
    scratch_data.kinematic_viscosity;

  const std::vector<double> &viscosity_for_stabilization_vector =
    scratch_data.kinematic_viscosity_for_stabilization;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix        = copy_data.local_matrix;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Pressure scaling factor
  const double pressure_scaling_factor = scratch_data.pressure_scaling_factor;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const double kinematic_viscosity = viscosity_vector[q];

      const Tensor<1, dim> &velocity = scratch_data.velocity_values[q];
      const Tensor<2, dim> &velocity_gradient =
        scratch_data.velocity_gradients[q];
      const Tensor<1, dim> &velocity_laplacian =
        scratch_data.velocity_laplacians[q];
      const Tensor<3, dim> &velocity_hessian =
        scratch_data.velocity_hessians[q];
      const Tensor<1, dim> &pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Calculate shear rate (at each q)
      const Tensor<2, dim> shear_rate =
        velocity_gradient + transpose(velocity_gradient);

      // Calculate the shear rate magnitude
      double shear_rate_magnitude = calculate_shear_rate_magnitude(shear_rate);
      // Set the shear rate magnitude to 1e-12 if it is too close to zero,
      // since the viscosity gradient is undefined for shear_rate_magnitude = 0
      shear_rate_magnitude =
        shear_rate_magnitude > 1e-12 ? shear_rate_magnitude : 1e-12;

      // Calculate viscosity gradient
      const Tensor<1, dim> viscosity_gradient =
        this->get_kinematic_viscosity_gradient(
          velocity_gradient,
          velocity_hessian,
          shear_rate_magnitude,
          scratch_data.grad_kinematic_viscosity_shear_rate[q]);

      // Forcing term
      const Tensor<1, dim> force       = scratch_data.force[q];
      double               mass_source = scratch_data.mass_source[q];

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(
            u_mag, viscosity_for_stabilization_vector[q], h) :
          calculate_navier_stokes_gls_tau_transient(
            u_mag, viscosity_for_stabilization_vector[q], h, sdt);

      // Calculate the strong residual for GLS stabilization
      auto strong_residual = velocity_gradient * velocity + pressure_gradient -
                             shear_rate * viscosity_gradient -
                             kinematic_viscosity * velocity_laplacian - force +
                             mass_source * velocity + strong_residual_vec[q];

      std::vector<Tensor<1, dim>> grad_phi_u_j_x_velocity(n_dofs);
      std::vector<Tensor<1, dim>> velocity_gradient_x_phi_u_j(n_dofs);


      // We loop over the column first to prevent recalculation
      // of the strong jacobian in the inner loop
      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto &phi_u_j           = scratch_data.phi_u[q][j];
          const auto &grad_phi_u_j      = scratch_data.grad_phi_u[q][j];
          const auto &laplacian_phi_u_j = scratch_data.laplacian_phi_u[q][j];

          const auto &grad_phi_p_j =
            pressure_scaling_factor * scratch_data.grad_phi_p[q][j];

          const auto &grad_phi_u_j_non_newtonian =
            grad_phi_u_j + transpose(grad_phi_u_j);

          strong_jacobian_vec[q][j] +=
            (velocity_gradient * phi_u_j + grad_phi_u_j * velocity +
             grad_phi_p_j - kinematic_viscosity * laplacian_phi_u_j -
             grad_phi_u_j_non_newtonian * viscosity_gradient +
             mass_source * phi_u_j);

          // Store these temporary products in auxiliary variables for speed
          grad_phi_u_j_x_velocity[j]     = grad_phi_u_j * velocity;
          velocity_gradient_x_phi_u_j[j] = velocity_gradient * phi_u_j;
        }

      shear_rate_magnitude =
        shear_rate_magnitude > 1e-3 ? shear_rate_magnitude : 1e-3;

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_u_i      = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &div_phi_u_i  = scratch_data.div_phi_u[q][i];
          const auto &phi_p_i      = scratch_data.phi_p[q][i];
          const auto &grad_phi_p_i = scratch_data.grad_phi_p[q][i];

          // Store these temporary products in auxiliary variables for speed
          const auto grad_phi_u_i_x_velocity = grad_phi_u_i * velocity;
          const auto strong_residual_x_grad_phi_u_i =
            strong_residual * grad_phi_u_i;

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto &phi_u_j      = scratch_data.phi_u[q][j];
              const auto &grad_phi_u_j = scratch_data.grad_phi_u[q][j];
              const auto &div_phi_u_j  = scratch_data.div_phi_u[q][j];

              const auto &grad_phi_u_j_non_newtonian =
                grad_phi_u_j + transpose(grad_phi_u_j);

              const auto &phi_p_j =
                pressure_scaling_factor * scratch_data.phi_p[q][j];

              const auto &strong_jac = strong_jacobian_vec[q][j];

              double local_matrix_ij =
                kinematic_viscosity *
                  scalar_product(grad_phi_u_j_non_newtonian, grad_phi_u_i) +
                0.5 * scratch_data.grad_kinematic_viscosity_shear_rate[q] /
                  shear_rate_magnitude *
                  scalar_product(grad_phi_u_j_non_newtonian, shear_rate) *
                  scalar_product(shear_rate, grad_phi_u_i) +
                velocity_gradient_x_phi_u_j[j] * phi_u_i +
                grad_phi_u_j_x_velocity[j] * phi_u_i - div_phi_u_i * phi_p_j +
                mass_source * phi_u_j * phi_u_i +
                // Continuity
                phi_p_i * div_phi_u_j;

              // PSPG GLS term
              local_matrix_ij += tau * (strong_jac * grad_phi_p_i);

              // The jacobian matrix for the SUPG formulation
              // currently does not include the jacobian of the stabilization
              // parameter tau. Our experience has shown that does not alter the
              // number of newton iteration for convergence, but greatly
              // simplifies assembly.
              if (SUPG)
                {
                  local_matrix_ij +=
                    tau * (strong_jac * grad_phi_u_i_x_velocity +
                           strong_residual_x_grad_phi_u_i * phi_u_j);
                }
              local_matrix_ij *= JxW;
              local_matrix(i, j) += local_matrix_ij;
            }
        }
    }
}

template <int dim>
void
GLSNavierStokesAssemblerNonNewtonianCore<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &viscosity_vector =
    scratch_data.kinematic_viscosity;

  const std::vector<double> &viscosity_for_stabilization_vector =
    scratch_data.kinematic_viscosity_for_stabilization;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs           = copy_data.local_rhs;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Physical properties
      const double kinematic_viscosity = viscosity_vector[q];

      // Velocity
      const Tensor<1, dim> &velocity   = scratch_data.velocity_values[q];
      const double velocity_divergence = scratch_data.velocity_divergences[q];
      const Tensor<2, dim> &velocity_gradient =
        scratch_data.velocity_gradients[q];
      const Tensor<1, dim> &velocity_laplacian =
        scratch_data.velocity_laplacians[q];
      const Tensor<3, dim> &velocity_hessian =
        scratch_data.velocity_hessians[q];

      // Calculate shear rate (at each q)
      const Tensor<2, dim> shear_rate =
        velocity_gradient + transpose(velocity_gradient);

      // Calculate the shear rate magnitude
      double shear_rate_magnitude = calculate_shear_rate_magnitude(shear_rate);

      shear_rate_magnitude =
        shear_rate_magnitude > 1e-12 ? shear_rate_magnitude : 1e-12;

      // Calculate kinematic viscosity gradient
      const Tensor<1, dim> kinematic_viscosity_gradient =
        this->get_kinematic_viscosity_gradient(
          velocity_gradient,
          velocity_hessian,
          shear_rate_magnitude,
          scratch_data.grad_kinematic_viscosity_shear_rate[q]);

      // Pressure
      const double         pressure = scratch_data.pressure_values[q];
      const Tensor<1, dim> pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Forcing term
      const Tensor<1, dim> force       = scratch_data.force[q];
      double               mass_source = scratch_data.mass_source[q];
      // Calculation of the magnitude of the
      // velocity for the stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(
            u_mag, viscosity_for_stabilization_vector[q], h) :
          calculate_navier_stokes_gls_tau_transient(
            u_mag, viscosity_for_stabilization_vector[q], h, sdt);


      // Calculate the strong residual for GLS stabilization
      auto strong_residual = velocity_gradient * velocity + pressure_gradient -
                             shear_rate * kinematic_viscosity_gradient -
                             kinematic_viscosity * velocity_laplacian - force +
                             mass_source * velocity + strong_residual_vec[q];

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_u_i      = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &phi_p_i      = scratch_data.phi_p[q][i];
          const auto &grad_phi_p_i = scratch_data.grad_phi_p[q][i];
          const auto &div_phi_u_i  = scratch_data.div_phi_u[q][i];

          double local_rhs_i = 0;

          // Navier-Stokes Residual
          local_rhs_i +=
            (
              // Momentum
              -kinematic_viscosity * scalar_product(shear_rate, grad_phi_u_i) -
              velocity_gradient * velocity * phi_u_i + pressure * div_phi_u_i +
              force * phi_u_i - mass_source * velocity * phi_u_i -
              // Continuity
              velocity_divergence * phi_p_i + mass_source * phi_p_i) *
            JxW;

          // PSPG GLS term
          local_rhs_i += -tau * (strong_residual * grad_phi_p_i) * JxW;

          // SUPG GLS term
          if (SUPG)
            {
              local_rhs_i +=
                -tau * (strong_residual * (grad_phi_u_i * velocity)) * JxW;
            }
          local_rhs(i) += local_rhs_i;
        }
    }
}

template class GLSNavierStokesAssemblerNonNewtonianCore<2>;
template class GLSNavierStokesAssemblerNonNewtonianCore<3>;


template <int dim>
void
GLSNavierStokesAssemblerSRF<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &strong_jacobian = copy_data.strong_jacobian;
  auto &local_matrix    = copy_data.local_matrix;

  // SRF Source term
  //----------------------------------
  // Angular velocity of the rotating frame. This is always a 3D vector even
  // in 2D.
  Tensor<1, dim> omega_vector;

  double omega_z  = velocity_sources.omega_z;
  omega_vector[0] = velocity_sources.omega_x;
  omega_vector[1] = velocity_sources.omega_y;
  if (dim == 3)
    omega_vector[2] = velocity_sources.omega_z;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Velocity
      const std::vector<Tensor<1, dim>> velocity = {
        scratch_data.velocity_values[q]};

      if (dim == 2)
        {
          strong_residual[q] +=
            2 * omega_z * (-1.) * cross_product_2d(velocity[0]);
          auto centrifugal =
            omega_z * (-1.) *
            cross_product_2d(
              omega_z * (-1.) *
              cross_product_2d(scratch_data.quadrature_points[q]));
          strong_residual[q] += centrifugal;
        }
      else // dim == 3
        {
          strong_residual[q] += 2 * cross_product_3d(omega_vector, velocity[0]);
          strong_residual[q] += cross_product_3d(
            omega_vector,
            cross_product_3d(omega_vector, scratch_data.quadrature_points[q]));
        }

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto &phi_u_j = scratch_data.phi_u[q][j];
          if (dim == 2)
            strong_jacobian[q][j] +=
              2 * omega_z * (-1.) * cross_product_2d(phi_u_j);
          else if (dim == 3)
            strong_jacobian[q][j] +=
              2 * cross_product_3d(omega_vector, phi_u_j);
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i = scratch_data.phi_u[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto &phi_u_j = scratch_data.phi_u[q][j];

              if (dim == 2)
                local_matrix(i, j) += 2 * omega_z * (-1.) *
                                      cross_product_2d(phi_u_j) * phi_u_i *
                                      JxW[q];

              else if (dim == 3)
                local_matrix(i, j) += 2 *
                                      cross_product_3d(omega_vector, phi_u_j) *
                                      phi_u_i * JxW[q];
            }
        }
    }
}

template <int dim>
void
GLSNavierStokesAssemblerSRF<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW               = scratch_data.JxW;
  const auto        &quadrature_points = scratch_data.quadrature_points;
  const unsigned int n_q_points        = scratch_data.n_q_points;
  const unsigned int n_dofs            = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &local_rhs       = copy_data.local_rhs;

  // SRF Source term
  //----------------------------------
  // Angular velocity of the rotating frame. This is always a 3D vector even
  // in 2D.
  Tensor<1, dim> omega_vector;

  double omega_z  = velocity_sources.omega_z;
  omega_vector[0] = velocity_sources.omega_x;
  omega_vector[1] = velocity_sources.omega_y;
  if (dim == 3)
    omega_vector[2] = velocity_sources.omega_z;


  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Velocity
      const std::vector<Tensor<1, dim>> velocity = {
        scratch_data.velocity_values[q]};

      if (dim == 2)
        {
          strong_residual[q] +=
            2 * omega_z * (-1.) * cross_product_2d(velocity[0]);
          auto centrifugal =
            omega_z * (-1.) *
            cross_product_2d(
              omega_z * (-1.) *
              cross_product_2d(scratch_data.quadrature_points[q]));
          strong_residual[q] += centrifugal;
        }
      else // dim == 3
        {
          strong_residual[q] += 2 * cross_product_3d(omega_vector, velocity[0]);
          strong_residual[q] += cross_product_3d(
            omega_vector,
            cross_product_3d(omega_vector, scratch_data.quadrature_points[q]));
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i = scratch_data.phi_u[q][i];

          if (dim == 2)
            {
              local_rhs(i) += -2 * omega_z * (-1.) *
                              cross_product_2d(velocity[0]) * phi_u_i * JxW[q];
              auto centrifugal =
                omega_z * (-1.) *
                cross_product_2d(omega_z * (-1.) *
                                 cross_product_2d(quadrature_points[q]));
              local_rhs(i) += -centrifugal * phi_u_i * JxW[q];
            }
          else if (dim == 3)
            {
              local_rhs(i) += -2 * cross_product_3d(omega_vector, velocity[0]) *
                              phi_u_i * JxW[q];
              local_rhs(i) +=
                -cross_product_3d(omega_vector,
                                  cross_product_3d(omega_vector,
                                                   quadrature_points[q])) *
                phi_u_i * JxW[q];
            }
        }
    }
}



template class GLSNavierStokesAssemblerSRF<2>;
template class GLSNavierStokesAssemblerSRF<3>;

template <int dim>
void
GLSNavierStokesAssemblerBDF<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &strong_jacobian = copy_data.strong_jacobian;
  auto &local_matrix    = copy_data.local_matrix;

  // Time stepping information
  const auto method = this->simulation_control->get_assembly_method();
  // Vector for the BDF coefficients
  const Vector<double> &bdf_coefs =
    this->simulation_control->get_bdf_coefficients();

  std::vector<Tensor<1, dim>> velocity(1 +
                                       number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      velocity[0] = scratch_data.velocity_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        velocity[p + 1] = scratch_data.previous_velocity_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += bdf_coefs[p] * velocity[p];
        }

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian[q][j] += bdf_coefs[0] * scratch_data.phi_u[q][j];
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const Tensor<1, dim> &phi_u_i = scratch_data.phi_u[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const Tensor<1, dim> &phi_u_j = scratch_data.phi_u[q][j];

              local_matrix(i, j) += phi_u_j * phi_u_i * bdf_coefs[0] * JxW[q];
            }
        }
    }
}

template <int dim>
void
GLSNavierStokesAssemblerBDF<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &local_rhs       = copy_data.local_rhs;

  // Time stepping information
  const auto method = this->simulation_control->get_assembly_method();
  // Vector for the BDF coefficients
  const Vector<double> &bdf_coefs =
    this->simulation_control->get_bdf_coefficients();

  std::vector<Tensor<1, dim>> velocity(1 +
                                       number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      velocity[0] = scratch_data.velocity_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        velocity[p + 1] = scratch_data.previous_velocity_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += (bdf_coefs[p] * velocity[p]);
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i     = scratch_data.phi_u[q][i];
          double     local_rhs_i = 0;
          for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
               ++p)
            {
              local_rhs_i -= bdf_coefs[p] * (velocity[p] * phi_u_i);
            }
          local_rhs(i) += local_rhs_i * JxW[q];
        }
    }
}

template class GLSNavierStokesAssemblerBDF<2>;
template class GLSNavierStokesAssemblerBDF<3>;

template <int dim>
void
BlockNavierStokesAssemblerNonNewtonianCore<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &kinematic_viscosity =
    scratch_data.kinematic_viscosity;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_matrix = copy_data.local_matrix;

  // Local variables to reuse multiplications
  std::vector<Tensor<1, dim>> grad_phi_u_j_x_velocity(n_dofs);
  std::vector<Tensor<1, dim>> velocity_gradient_x_phi_u_j(n_dofs);

  // Pressure scaling factor
  const double pressure_scaling_factor = scratch_data.pressure_scaling_factor;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const Tensor<1, dim> &velocity = scratch_data.velocity_values[q];
      const Tensor<2, dim> &velocity_gradient =
        scratch_data.velocity_gradients[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];



      // We loop over the column first to prevent recalculation
      // of the strong jacobian in the inner loop
      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto &phi_u_j      = scratch_data.phi_u[q][j];
          const auto &grad_phi_u_j = scratch_data.grad_phi_u[q][j];

          // Store these temporary products in auxiliary variables for speed
          grad_phi_u_j_x_velocity[j]     = grad_phi_u_j * velocity;
          velocity_gradient_x_phi_u_j[j] = velocity_gradient * phi_u_j;
        }



      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_u_i      = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &div_phi_u_i  = scratch_data.div_phi_u[q][i];
          const auto &phi_p_i      = scratch_data.phi_p[q][i];


          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto &grad_phi_u_j = scratch_data.grad_phi_u[q][j];
              const auto &div_phi_u_j  = scratch_data.div_phi_u[q][j];

              const auto &phi_p_j =
                pressure_scaling_factor * scratch_data.phi_p[q][j];

              const auto &grad_phi_u_j_non_newtonian =
                grad_phi_u_j + transpose(grad_phi_u_j);

              double local_matrix_ij =
                kinematic_viscosity[q] *
                  scalar_product(grad_phi_u_j_non_newtonian, grad_phi_u_i) +
                velocity_gradient_x_phi_u_j[j] * phi_u_i +
                grad_phi_u_j_x_velocity[j] * phi_u_i - div_phi_u_i * phi_p_j -
                // Continuity
                phi_p_i * div_phi_u_j + gamma * div_phi_u_j * div_phi_u_i +
                // Mass matrix for S block in schur complement
                phi_p_i * phi_p_j;

              local_matrix_ij *= JxW;
              local_matrix(i, j) += local_matrix_ij;
            }
        }
    }
}

template <int dim>
void
BlockNavierStokesAssemblerNonNewtonianCore<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &kinematic_viscosity_vector =
    scratch_data.kinematic_viscosity;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_rhs = copy_data.local_rhs;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Physical properties
      const double kinematic_viscosity = kinematic_viscosity_vector[q];

      // Velocity
      const Tensor<1, dim> &velocity   = scratch_data.velocity_values[q];
      const double velocity_divergence = scratch_data.velocity_divergences[q];
      const Tensor<2, dim> &velocity_gradient =
        scratch_data.velocity_gradients[q];

      // Calculate shear rate (at each q)
      const Tensor<2, dim> shear_rate =
        velocity_gradient + transpose(velocity_gradient);

      // Pressure
      const double pressure = scratch_data.pressure_values[q];

      // Forcing term
      const Tensor<1, dim> &force = scratch_data.force[q];

      // Store JxW in local variable for faster access
      const double JxW = JxW_vec[q];

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i      = scratch_data.phi_u[q][i];
          const auto grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto phi_p_i      = scratch_data.phi_p[q][i];
          const auto div_phi_u_i  = scratch_data.div_phi_u[q][i];

          double local_rhs_i = 0;

          // Navier-Stokes Residual
          local_rhs_i +=
            (
              // Momentum
              -kinematic_viscosity * scalar_product(shear_rate, grad_phi_u_i) -
              velocity_gradient * velocity * phi_u_i + pressure * div_phi_u_i +
              force * phi_u_i +
              // Continuity
              velocity_divergence * phi_p_i -
              gamma * velocity_divergence * div_phi_u_i) *
            JxW;

          local_rhs(i) += local_rhs_i;
        }
    }
}

template class BlockNavierStokesAssemblerNonNewtonianCore<2>;
template class BlockNavierStokesAssemblerNonNewtonianCore<3>;


template <int dim>
void
BlockNavierStokesAssemblerCore<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &kinematic_viscosity_vector =
    scratch_data.kinematic_viscosity;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_matrix = copy_data.local_matrix;

  // Local variables to reuse multiplications
  std::vector<Tensor<1, dim>> grad_phi_u_j_x_velocity(n_dofs);
  std::vector<Tensor<1, dim>> velocity_gradient_x_phi_u_j(n_dofs);

  // Pressure scaling factor
  const double pressure_scaling_factor = scratch_data.pressure_scaling_factor;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Physical properties
      const double kinematic_viscosity = kinematic_viscosity_vector[q];


      // Gather into local variables the relevant fields
      const Tensor<1, dim> &velocity = scratch_data.velocity_values[q];
      const Tensor<2, dim> &velocity_gradient =
        scratch_data.velocity_gradients[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // We loop over the column first to prevent recalculation
      // of the strong jacobian in the inner loop
      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto &phi_u_j      = scratch_data.phi_u[q][j];
          const auto &grad_phi_u_j = scratch_data.grad_phi_u[q][j];

          // Store these temporary products in auxiliary variables for speed
          grad_phi_u_j_x_velocity[j]     = grad_phi_u_j * velocity;
          velocity_gradient_x_phi_u_j[j] = velocity_gradient * phi_u_j;
        }



      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_u_i      = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &div_phi_u_i  = scratch_data.div_phi_u[q][i];
          const auto &phi_p_i      = scratch_data.phi_p[q][i];


          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto &grad_phi_u_j = scratch_data.grad_phi_u[q][j];
              const auto &div_phi_u_j  = scratch_data.div_phi_u[q][j];

              const auto &phi_p_j =
                pressure_scaling_factor * scratch_data.phi_p[q][j];

              double local_matrix_ij =
                kinematic_viscosity *
                  scalar_product(grad_phi_u_j, grad_phi_u_i) +
                velocity_gradient_x_phi_u_j[j] * phi_u_i +
                grad_phi_u_j_x_velocity[j] * phi_u_i - div_phi_u_i * phi_p_j -
                // Continuity
                phi_p_i * div_phi_u_j + gamma * div_phi_u_j * div_phi_u_i +
                // Mass matrix for S block in schur complement
                phi_p_i * phi_p_j;

              local_matrix_ij *= JxW;
              local_matrix(i, j) += local_matrix_ij;
            }
        }
    }
}

template <int dim>
void
BlockNavierStokesAssemblerCore<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &kinematic_viscosity_vector =
    scratch_data.kinematic_viscosity;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_rhs = copy_data.local_rhs;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Physical properties
      const double kinematic_viscosity = kinematic_viscosity_vector[q];


      // Velocity
      const Tensor<1, dim> &velocity   = scratch_data.velocity_values[q];
      const double velocity_divergence = scratch_data.velocity_divergences[q];
      const Tensor<2, dim> &velocity_gradient =
        scratch_data.velocity_gradients[q];

      // Pressure
      const double pressure = scratch_data.pressure_values[q];

      // Forcing term
      const Tensor<1, dim> &force = scratch_data.force[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_u_i      = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &phi_p_i      = scratch_data.phi_p[q][i];
          const auto &div_phi_u_i  = scratch_data.div_phi_u[q][i];

          double local_rhs_i = 0;

          // Navier-Stokes Residual
          local_rhs_i += (
                           // Momentum
                           -kinematic_viscosity *
                             scalar_product(velocity_gradient, grad_phi_u_i) -
                           velocity_gradient * velocity * phi_u_i +
                           pressure * div_phi_u_i + force * phi_u_i +
                           // Continuity
                           velocity_divergence * phi_p_i -
                           gamma * velocity_divergence * div_phi_u_i) *
                         JxW;

          local_rhs(i) += local_rhs_i;
        }
    }
}



template class BlockNavierStokesAssemblerCore<2>;
template class BlockNavierStokesAssemblerCore<3>;

template <int dim>
void
LaplaceAssembly<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW_vec       = scratch_data.JxW;
  const unsigned int n_q_points    = scratch_data.n_q_points;
  const unsigned int n_dofs        = scratch_data.n_dofs;
  const double kinematic_viscosity = scratch_data.kinematic_viscosity_scale;
  const double h                   = scratch_data.cell_size;
  // Copy data elements
  auto &local_matrix = copy_data.local_matrix;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &grad_phi_p_i = scratch_data.grad_phi_p[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto &grad_phi_u_j = scratch_data.grad_phi_u[q][j];
              const auto &grad_phi_p_j = scratch_data.grad_phi_p[q][j];

              // Units are not consistent for the following line, but it is
              // better that way for the problem scaling and it doesn't affect
              // the end result for our purposes
              // Laplacian on the velocity terms
              double local_matrix_ij =
                kinematic_viscosity *
                scalar_product(grad_phi_u_j, grad_phi_u_i);

              // Laplacian on the pressure terms
              local_matrix_ij += 1 / kinematic_viscosity * h *
                                 scalar_product(grad_phi_p_j, grad_phi_p_i);

              // The jacobian matrix for the SUPG formulation
              // currently does not include the jacobian of the stabilization
              // parameter tau. Our experience has shown that does not alter the
              // number of newton iteration for convergence, but greatly
              // simplifies assembly.

              local_matrix_ij *= JxW;
              local_matrix(i, j) += local_matrix_ij;
            }
        }
    }
}

template <int dim>
void
LaplaceAssembly<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW_vec       = scratch_data.JxW;
  const unsigned int n_q_points    = scratch_data.n_q_points;
  const unsigned int n_dofs        = scratch_data.n_dofs;
  auto              &local_rhs     = copy_data.local_rhs;
  const double kinematic_viscosity = scratch_data.kinematic_viscosity_scale;
  const double h                   = scratch_data.cell_size;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Velocity
      const Tensor<2, dim> &velocity_gradient =
        scratch_data.velocity_gradients[q];

      // Pressure
      const Tensor<1, dim> &pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];


      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &grad_phi_p_i = scratch_data.grad_phi_p[q][i];


          double local_rhs_i = 0;

          // Laplacian on the velocity terms
          local_rhs_i += -kinematic_viscosity *
                         scalar_product(velocity_gradient, grad_phi_u_i) * JxW;


          // Laplacian on the pressure terms
          local_rhs_i += -1 / kinematic_viscosity * h *
                         scalar_product(pressure_gradient, grad_phi_p_i) * JxW;

          if (scratch_data.components[i] == dim)
            local_rhs_i /= scratch_data.pressure_scaling_factor;
          local_rhs(i) += local_rhs_i;
        }
    }
}

template class LaplaceAssembly<2>;
template class LaplaceAssembly<3>;


template <int dim>
void
BuoyancyAssembly<dim>::assemble_matrix(
  const NavierStokesScratchData<dim> & /*scratch_data*/,
  StabilizedMethodsTensorCopyData<dim> & /*copy_data*/)
{}



template <int dim>
void
BuoyancyAssembly<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  auto &local_rhs       = copy_data.local_rhs;
  auto &strong_residual = copy_data.strong_residual;


  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Forcing term (gravity)
      const Tensor<1, dim> &force = scratch_data.force[q];

      const double thermal_expansion = scratch_data.thermal_expansion[q];


      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Current temperature values
      double current_temperature = scratch_data.temperature_values[q];

      strong_residual[q] += force * thermal_expansion *
                            (current_temperature - reference_temperature);

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i = scratch_data.phi_u[q][i];

          // Laplacian on the velocity terms
          local_rhs(i) -= force * thermal_expansion *
                          (current_temperature - reference_temperature) *
                          phi_u_i * JxW;
        }
    }
}

template class BuoyancyAssembly<2>;
template class BuoyancyAssembly<3>;

template <int dim>
void
PhaseChangeDarcyAssembly<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW_vec      = scratch_data.JxW;
  const unsigned int n_q_points   = scratch_data.n_q_points;
  const unsigned int n_dofs       = scratch_data.n_dofs;
  const auto        &velocities   = scratch_data.velocity_values;
  const auto        &temperatures = scratch_data.temperature_values;

  auto &local_matrix    = copy_data.local_matrix;
  auto &strong_residual = copy_data.strong_residual;
  auto &strong_jacobian = copy_data.strong_jacobian;


  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];
      // Current temperature values
      double current_temperature = temperatures[q];
      // Depending on the value of the temperature, calculate the Darcy
      // permeability. This min(max) approach is used to calculate the liquid
      // fraction without requiring if conditions.
      const double liquid_fraction =
        std::min(1.,
                 std::max((current_temperature -
                           phase_change_parameters.T_solidus) /
                            (phase_change_parameters.T_liquidus -
                             phase_change_parameters.T_solidus),
                          0.));
      double darcy_penalty =
        phase_change_parameters.penalty_l * liquid_fraction +
        (1. - liquid_fraction) * phase_change_parameters.penalty_s;

      strong_residual[q] += darcy_penalty * velocities[q];

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian[q][j] += darcy_penalty * scratch_data.phi_u[q][j];
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_u_i = scratch_data.phi_u[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto &phi_u_j = scratch_data.phi_u[q][j];
              local_matrix(i, j) += darcy_penalty * phi_u_i * phi_u_j * JxW;
            }
        }
    }
}

template <int dim>
void
PhaseChangeDarcyAssembly<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW_vec      = scratch_data.JxW;
  const unsigned int n_q_points   = scratch_data.n_q_points;
  const unsigned int n_dofs       = scratch_data.n_dofs;
  const auto        &velocities   = scratch_data.velocity_values;
  const auto        &temperatures = scratch_data.temperature_values;


  auto &local_rhs       = copy_data.local_rhs;
  auto &strong_residual = copy_data.strong_residual;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];
      // Current temperature values
      double current_temperature = temperatures[q];
      // Depending on the value of the temperature, calculate the Darcy
      // permeability. This min(max) approach is used to calculate the liquid
      // fraction without requiring if conditions.
      const double liquid_fraction =
        std::min(1.,
                 std::max((current_temperature -
                           phase_change_parameters.T_solidus) /
                            (phase_change_parameters.T_liquidus -
                             phase_change_parameters.T_solidus),
                          0.));
      double darcy_penalty =
        phase_change_parameters.penalty_l * liquid_fraction +
        (1. - liquid_fraction) * phase_change_parameters.penalty_s;

      strong_residual[q] += darcy_penalty * velocities[q];

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i = scratch_data.phi_u[q][i];

          // Laplacian on the velocity terms
          local_rhs(i) -= darcy_penalty * velocities[q] * phi_u_i * JxW;
        }
    }
}

template class PhaseChangeDarcyAssembly<2>;
template class PhaseChangeDarcyAssembly<3>;


template <int dim>
void
PressureBoundaryCondition<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  // Scheme and physical properties
  // To generalize for dependent viscosity
  const double kinematic_viscosity = scratch_data.kinematic_viscosity_scale;

  // Loop and quadrature information
  Tensor<2, dim> identity;
  for (unsigned int d = 0; d < dim; ++d)
    {
      identity[d][d] = 1;
    }

  std::vector<std::vector<std::vector<Tensor<1, dim>>>> gn_j;
  gn_j = std::vector<std::vector<std::vector<Tensor<1, dim>>>>(
    scratch_data.n_faces,
    std::vector<std::vector<Tensor<1, dim>>>(scratch_data.n_faces_q_points,
                                             std::vector<Tensor<1, dim>>(
                                               scratch_data.n_dofs)));

  auto &local_matrix = copy_data.local_matrix;

  // Loop over the faces of the cell.
  for (unsigned int f = 0; f < scratch_data.n_faces; ++f)
    {
      // Check if the face is on a boundary
      if (scratch_data.is_boundary_face[f])
        {
          types::boundary_id boundary_id = scratch_data.boundary_face_id[f];
          // Check if the face is part of the boundary that as a
          // pressure BC.
          if (this->pressure_boundary_conditions.type.at(boundary_id) ==
              BoundaryConditions::BoundaryType::pressure)
            {
              // Assemble the matrix of the BC
              for (unsigned int q = 0; q < scratch_data.n_faces_q_points; ++q)
                {
                  const double JxW = scratch_data.face_JxW[f][q];
                  for (const unsigned int j :
                       scratch_data.fe_face_values.dof_indices())
                    {
                      gn_j[f][q][j] =
                        (-kinematic_viscosity *
                         (scratch_data.face_grad_phi_u[f][q][j])) *
                        scratch_data.face_normal[f][q];
                    }
                  for (const unsigned int i :
                       scratch_data.fe_face_values.dof_indices())
                    {
                      for (const unsigned int j :
                           scratch_data.fe_face_values.dof_indices())
                        {
                          local_matrix[i][j] +=
                            -scratch_data.face_phi_u[f][q][i] * gn_j[f][q][j] *
                            JxW;
                        }
                    }
                }
            }
        }
    }
}

template <int dim>
void
PressureBoundaryCondition<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  // Scheme and physical properties
  const double kinematic_viscosity = scratch_data.kinematic_viscosity_scale;

  // Loop and quadrature information
  Tensor<2, dim> identity;
  for (unsigned int d = 0; d < dim; ++d)
    {
      identity[d][d] = 1;
    }
  std::vector<std::vector<double>> prescribed_pressure_values;
  prescribed_pressure_values = std::vector<std::vector<double>>(
    scratch_data.n_faces, std::vector<double>(scratch_data.n_faces_q_points));

  std::vector<std::vector<Tensor<1, dim>>> gn_bc =
    std::vector<std::vector<Tensor<1, dim>>>(scratch_data.n_faces,
                                             std::vector<Tensor<1, dim>>(
                                               scratch_data.n_faces_q_points));


  auto &local_rhs = copy_data.local_rhs;

  // Pressure boundary condition, loop on faces
  // implementation similar to deal.ii step-22
  // Loop over the BCs
  // Loop over the faces of the cell.
  for (unsigned int f = 0; f < scratch_data.n_faces; ++f)
    {
      // Check if the face is on a boundary
      if (scratch_data.is_boundary_face[f])
        {
          types::boundary_id boundary_id = scratch_data.boundary_face_id[f];

          NavierStokesPressureFunctionDefined<dim> function_p(
            &pressure_boundary_conditions.navier_stokes_functions
               .at(boundary_id)
               ->p);
          // Check if the face is part of the boundary that as a
          // pressure BC.
          if (this->pressure_boundary_conditions.type.at(boundary_id) ==
              BoundaryConditions::BoundaryType::pressure)
            {
              // Assemble the rhs of the BC
              for (unsigned int q = 0; q < scratch_data.n_faces_q_points; ++q)
                {
                  const double JxW = scratch_data.face_JxW[f][q];
                  prescribed_pressure_values[f][q] =
                    function_p.value(scratch_data.face_quadrature_points[f][q],
                                     dim);
                  gn_bc[f][q] =
                    (-prescribed_pressure_values[f][q] * identity -
                     kinematic_viscosity *
                       (scratch_data.face_velocity_gradients[f][q])) *
                    scratch_data.face_normal[f][q];
                  for (const unsigned int i :
                       scratch_data.fe_face_values.dof_indices())
                    {
                      local_rhs(i) -=
                        -scratch_data.face_phi_u[f][q][i] * (gn_bc[f][q]) * JxW;
                    }
                }
            }
        }
    }
}

template class PressureBoundaryCondition<2>;
template class PressureBoundaryCondition<3>;

template <int dim>
void
WeakDirichletBoundaryCondition<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  // Scheme and physical properties
  const double kinematic_viscosity = scratch_data.kinematic_viscosity_scale;
  // Loop and quadrature information
  Tensor<2, dim> identity;
  for (unsigned int d = 0; d < dim; ++d)
    {
      identity[d][d] = 1;
    }
  std::vector<std::vector<Tensor<1, dim>>> prescribed_velocity_values;
  prescribed_velocity_values =
    std::vector<std::vector<Tensor<1, dim>>>(scratch_data.n_faces,
                                             std::vector<Tensor<1, dim>>(
                                               scratch_data.n_faces_q_points));
  const FiniteElement<dim> &fe = scratch_data.fe_face_values.get_fe();

  const double penalty_parameter =
    1. / std::pow(scratch_data.cell_size, fe.degree + 1);
  auto &local_matrix = copy_data.local_matrix;

  // Loop over the faces of the cell.
  for (unsigned int f = 0; f < scratch_data.n_faces; ++f)
    {
      // Check if the face is on a boundary
      if (scratch_data.is_boundary_face[f])
        {
          types::boundary_id boundary_id = scratch_data.boundary_face_id[f];
          // Check if the face is part of the boundary that as a
          // pressure BC.
          if (this->boundary_conditions.type.at(boundary_id) ==
              BoundaryConditions::BoundaryType::function_weak)
            {
              const double beta = boundary_conditions.beta.at(boundary_id);

              // Assemble the matrix of the BC
              for (unsigned int q = 0; q < scratch_data.n_faces_q_points; ++q)
                {
                  const double JxW = scratch_data.face_JxW[f][q];
                  for (const unsigned int i :
                       scratch_data.fe_face_values.dof_indices())
                    {
                      const auto comp_i = fe.system_to_component_index(i).first;
                      if (comp_i < dim)
                        {
                          for (const unsigned int j :
                               scratch_data.fe_face_values.dof_indices())
                            {
                              const auto comp_j =
                                fe.system_to_component_index(j).first;
                              if (comp_i == comp_j)
                                {
                                  double beta_terms =
                                    penalty_parameter * beta *
                                    (scratch_data.face_phi_u[f][q][j][comp_i]) *
                                    scratch_data.face_phi_u[f][q][i][comp_i] *
                                    JxW;
                                  double grad_phi_terms =
                                    ((kinematic_viscosity *
                                      scratch_data.face_phi_u[f][q][j]) *
                                     (scratch_data.face_grad_phi_u[f][q][i] *
                                      scratch_data.face_normal[f][q])) *
                                    JxW;
                                  double surface_stress_term =
                                    (kinematic_viscosity *
                                     scratch_data.face_grad_phi_u[f][q][j] *
                                     scratch_data.face_normal[f][q]) *
                                    scratch_data.face_phi_u[f][q][i] * JxW;

                                  local_matrix(i, j) += +beta_terms -
                                                        grad_phi_terms -
                                                        surface_stress_term;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template <int dim>
void
WeakDirichletBoundaryCondition<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  // Scheme and physical properties
  const double kinematic_viscosity = scratch_data.kinematic_viscosity_scale;

  // Loop and quadrature information
  Tensor<2, dim> identity;
  for (unsigned int d = 0; d < dim; ++d)
    {
      identity[d][d] = 1;
    }
  std::vector<std::vector<Tensor<1, dim>>> prescribed_velocity_values =
    std::vector<std::vector<Tensor<1, dim>>>(scratch_data.n_faces,
                                             std::vector<Tensor<1, dim>>(
                                               scratch_data.n_faces_q_points));

  const FiniteElement<dim> &fe = scratch_data.fe_face_values.get_fe();

  const double penalty_parameter =
    1. / std::pow(scratch_data.cell_size, fe.degree + 1);
  auto &local_rhs = copy_data.local_rhs;
  // Loop over the faces of the cell.
  for (unsigned int f = 0; f < scratch_data.n_faces; ++f)
    {
      // Check if the face is on a boundary
      if (scratch_data.is_boundary_face[f])
        {
          types::boundary_id boundary_id = scratch_data.boundary_face_id[f];
          // Check if the face is part of the boundary that has a
          // weakly imposed Dirichlet BC.
          if (this->boundary_conditions.type.at(boundary_id) ==
              BoundaryConditions::BoundaryType::function_weak)
            {
              const double beta = boundary_conditions.beta.at(boundary_id);

              NavierStokesFunctionDefined<dim> function_v(
                &boundary_conditions.navier_stokes_functions.at(boundary_id)->u,
                &boundary_conditions.navier_stokes_functions.at(boundary_id)->v,
                &boundary_conditions.navier_stokes_functions.at(boundary_id)
                   ->w);
              for (unsigned int q = 0; q < scratch_data.n_faces_q_points; ++q)
                {
                  const double JxW = scratch_data.face_JxW[f][q];
                  for (unsigned int d = 0; d < dim; ++d)
                    {
                      prescribed_velocity_values[f][q][d] = function_v.value(
                        scratch_data.face_quadrature_points[f][q], d);
                    }
                  for (const unsigned int i :
                       scratch_data.fe_face_values.dof_indices())
                    {
                      const auto comp_i = fe.system_to_component_index(i).first;
                      if (comp_i < dim)
                        {
                          double beta_terms =
                            penalty_parameter * beta *
                            (scratch_data.face_velocity_values[f][q][comp_i] -
                             prescribed_velocity_values[f][q][comp_i]) *
                            scratch_data.face_phi_u[f][q][i][comp_i] * JxW;
                          double grad_phi_terms =
                            ((kinematic_viscosity *
                              (scratch_data.face_velocity_values[f][q] -
                               prescribed_velocity_values[f][q])) *
                             (scratch_data.face_grad_phi_u[f][q][i] *
                              scratch_data.face_normal[f][q])) *
                            JxW;
                          double surface_stress_term =
                            (kinematic_viscosity *
                             scratch_data.face_velocity_gradients[f][q] *
                             scratch_data.face_normal[f][q] *
                             scratch_data.face_phi_u[f][q][i]) *
                            JxW;

                          local_rhs(i) +=
                            -beta_terms + grad_phi_terms + surface_stress_term;
                        }
                    }
                }
            }
        }
    }
}

template class WeakDirichletBoundaryCondition<2>;
template class WeakDirichletBoundaryCondition<3>;

template <int dim>
void
PartialSlipDirichletBoundaryCondition<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  // Scheme and physical properties
  const double kinematic_viscosity = scratch_data.kinematic_viscosity[0];

  // Loop and quadrature information
  Tensor<2, dim> identity;
  for (unsigned int d = 0; d < dim; ++d)
    {
      identity[d][d] = 1;
    }
  std::vector<std::vector<Tensor<1, dim>>> prescribed_velocity_values;
  prescribed_velocity_values =
    std::vector<std::vector<Tensor<1, dim>>>(scratch_data.n_faces,
                                             std::vector<Tensor<1, dim>>(
                                               scratch_data.n_faces_q_points));
  const FiniteElement<dim> &fe = scratch_data.fe_face_values.get_fe();

  const double penalty_parameter =
    1. / std::pow(scratch_data.cell_size, fe.degree + 1);
  auto &local_matrix = copy_data.local_matrix;

  // Loop over the faces of the cell.
  for (unsigned int f = 0; f < scratch_data.n_faces; ++f)
    {
      // Check if the face is on a boundary
      if (scratch_data.is_boundary_face[f])
        {
          types::boundary_id boundary_id = scratch_data.boundary_face_id[f];
          // Check if the face is part of the boundary that as a
          // pressure BC.
          if (this->boundary_conditions.type.at(boundary_id) ==
              BoundaryConditions::BoundaryType::partial_slip)
            {
              const double beta = boundary_conditions.beta.at(boundary_id);
              const double boundary_layer_thickness =
                boundary_conditions.boundary_layer_thickness.at(boundary_id);
              const double beta_tangent =
                kinematic_viscosity / boundary_layer_thickness;

              // Assemble the matrix of the BC
              for (unsigned int q = 0; q < scratch_data.n_faces_q_points; ++q)
                {
                  const double JxW = scratch_data.face_JxW[f][q];
                  for (const unsigned int i :
                       scratch_data.fe_face_values.dof_indices())
                    {
                      const auto comp_i = fe.system_to_component_index(i).first;
                      if (comp_i < dim)
                        {
                          for (const unsigned int j :
                               scratch_data.fe_face_values.dof_indices())
                            {
                              const auto comp_j =
                                fe.system_to_component_index(j).first;
                              if (comp_i == comp_j)
                                {
                                  double beta_terms_normal =
                                    penalty_parameter * beta *
                                    scratch_data.face_normal[f][q][comp_i] *
                                    (scratch_data.face_normal[f][q][comp_i] *
                                     scratch_data.face_phi_u[f][q][j][comp_i]) *
                                    scratch_data.face_phi_u[f][q][i][comp_i] *
                                    JxW;
                                  double beta_terms_tangent =
                                    beta_tangent *
                                    (scratch_data.face_phi_u[f][q][j][comp_i] -
                                     scratch_data.face_normal[f][q][comp_i] *
                                       (scratch_data.face_normal[f][q][comp_i] *
                                        scratch_data
                                          .face_phi_u[f][q][j][comp_i])) *
                                    scratch_data.face_phi_u[f][q][i][comp_i] *
                                    JxW;

                                  /* We don't use these terms has its
                                   * seems they don't affect the
                                   * solution for the normal component
                                   */
                                  /*
                                                                            double
                                     grad_phi_terms =
                                                                              scratch_data.face_normal[f][q][comp_i]*((kinematic_viscosity
                                     * scratch_data .face_phi_u[f][q][j]) *
                                                                               (scratch_data
                                                                                  .face_grad_phi_u[f][q][i] *
                                                                                scratch_data.face_normal[f][q])) *
                                                                              JxW;
                                                                            double
                                     surface_stress_term =
                                                                              scratch_data.face_normal[f][q][comp_i]*
                                     (kinematic_viscosity * scratch_data
                                                                                 .face_grad_phi_u[f][q][j] *
                                                                               scratch_data.face_normal[f][q]) *
                                                                              scratch_data.face_phi_u[f][q][i]
                                     * JxW;
                                  */

                                  local_matrix(i, j) +=
                                    +beta_terms_normal +
                                    beta_terms_tangent; //-
                                                        // grad_phi_terms
                                                        //-
                                  // surface_stress_term;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template <int dim>
void
PartialSlipDirichletBoundaryCondition<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  // Scheme and physical properties
  const double kinematic_viscosity = scratch_data.kinematic_viscosity[0];

  // Loop and quadrature information
  Tensor<2, dim> identity;
  for (unsigned int d = 0; d < dim; ++d)
    {
      identity[d][d] = 1;
    }
  std::vector<std::vector<Tensor<1, dim>>> prescribed_velocity_values =
    std::vector<std::vector<Tensor<1, dim>>>(scratch_data.n_faces,
                                             std::vector<Tensor<1, dim>>(
                                               scratch_data.n_faces_q_points));

  const FiniteElement<dim> &fe = scratch_data.fe_face_values.get_fe();

  const double penalty_parameter =
    1. / std::pow(scratch_data.cell_size, fe.degree + 1);
  auto &local_rhs = copy_data.local_rhs;
  // Loop over the faces of the cell.
  for (unsigned int f = 0; f < scratch_data.n_faces; ++f)
    {
      // Check if the face is on a boundary
      if (scratch_data.is_boundary_face[f])
        {
          types::boundary_id boundary_id = scratch_data.boundary_face_id[f];

          // Check if the face is part of the boundary that has a
          // weakly imposed Dirichlet BC.
          if (this->boundary_conditions.type.at(boundary_id) ==
              BoundaryConditions::BoundaryType::partial_slip)
            {
              const double beta = boundary_conditions.beta.at(boundary_id);
              const double boundary_layer_thickness =
                boundary_conditions.boundary_layer_thickness.at(boundary_id);
              const double beta_tangent =
                kinematic_viscosity / boundary_layer_thickness;
              NavierStokesFunctionDefined<dim> function_v(
                &boundary_conditions.navier_stokes_functions.at(boundary_id)->u,
                &boundary_conditions.navier_stokes_functions.at(boundary_id)->v,
                &boundary_conditions.navier_stokes_functions.at(boundary_id)
                   ->w);
              for (unsigned int q = 0; q < scratch_data.n_faces_q_points; ++q)
                {
                  const double JxW = scratch_data.face_JxW[f][q];
                  for (unsigned int d = 0; d < dim; ++d)
                    {
                      prescribed_velocity_values[f][q][d] = function_v.value(
                        scratch_data.face_quadrature_points[f][q], d);
                    }
                  for (const unsigned int i :
                       scratch_data.fe_face_values.dof_indices())
                    {
                      const auto comp_i = fe.system_to_component_index(i).first;
                      if (comp_i < dim)
                        {
                          double beta_terms_normal =
                            scratch_data.face_normal[f][q][comp_i] *
                            penalty_parameter * beta *
                            (scratch_data.face_normal[f][q][comp_i] *
                               scratch_data.face_velocity_values[f][q][comp_i] -
                             scratch_data.face_normal[f][q][comp_i] *
                               prescribed_velocity_values[f][q][comp_i]) *
                            scratch_data.face_phi_u[f][q][i][comp_i] * JxW;
                          double beta_terms_tangent =
                            beta_tangent *
                            (scratch_data.face_velocity_values[f][q][comp_i] -
                             scratch_data.face_normal[f][q][comp_i] *
                               (scratch_data.face_normal[f][q][comp_i] *
                                scratch_data
                                  .face_velocity_values[f][q][comp_i]) -
                             (prescribed_velocity_values[f][q][comp_i] -
                              scratch_data.face_normal[f][q][comp_i] *
                                (scratch_data.face_normal[f][q][comp_i] *
                                 prescribed_velocity_values[f][q][comp_i]))) *
                            scratch_data.face_phi_u[f][q][i][comp_i] * JxW;
                          /* We don't use these terms has its seems they
                           * don't affect the solution for the normal
                           * component */
                          /*
                          double grad_phi_terms =
                            scratch_data.face_normal[f][q][comp_i]*((kinematic_viscosity
                          * (scratch_data.face_velocity_values[f][q] -
                               prescribed_velocity_values[f][q])) *
                             (scratch_data.face_grad_phi_u[f][q][i] *
                              scratch_data.face_normal[f][q])) *
                            JxW;
                          double surface_stress_term =
                            scratch_data.face_normal[f][q][comp_i]*(kinematic_viscosity
                          * scratch_data .face_velocity_gradients[f][q]
                          * scratch_data.face_normal[f][q] *
                             scratch_data.face_phi_u[f][q][i]) *
                            JxW;*/



                          local_rhs(i) +=
                            -beta_terms_normal -
                            beta_terms_tangent; //+ grad_phi_terms +
                                                // surface_stress_term;
                        }
                    }
                }
            }
        }
    }
}

template class PartialSlipDirichletBoundaryCondition<2>;
template class PartialSlipDirichletBoundaryCondition<3>;

template <int dim>
void
OutletBoundaryCondition<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  const FiniteElement<dim> &fe = scratch_data.fe_face_values.get_fe();

  const double penalty_parameter =
    1. / std::pow(scratch_data.cell_size, fe.degree + 1);
  auto &local_matrix = copy_data.local_matrix;

  // Loop over the faces of the cell.
  for (unsigned int f = 0; f < scratch_data.n_faces; ++f)
    {
      // Check if the face is on a boundary
      if (scratch_data.is_boundary_face[f])
        {
          types::boundary_id boundary_id = scratch_data.boundary_face_id[f];
          // Check if the face is part of the boundary that as a
          // pressure BC.
          if (this->boundary_conditions.type.at(boundary_id) ==
              BoundaryConditions::BoundaryType::outlet)
            {
              const double beta = boundary_conditions.beta.at(boundary_id);

              // Assemble the matrix of the BC
              for (unsigned int q = 0; q < scratch_data.n_faces_q_points; ++q)
                {
                  const double JxW = scratch_data.face_JxW[f][q];
                  for (const unsigned int i :
                       scratch_data.fe_face_values.dof_indices())
                    {
                      double normal_outflux =
                        std::min(0.,
                                 scratch_data.face_velocity_values[f][q] *
                                   scratch_data.face_normal[f][q]);

                      const auto comp_i = fe.system_to_component_index(i).first;
                      if (comp_i < dim)
                        {
                          for (const unsigned int j :
                               scratch_data.fe_face_values.dof_indices())
                            {
                              const auto comp_j =
                                fe.system_to_component_index(j).first;
                              if (comp_i == comp_j)
                                {
                                  double beta_terms =
                                    penalty_parameter * beta * normal_outflux *
                                    (scratch_data.face_phi_u[f][q][j][comp_i] *
                                     scratch_data.face_phi_u[f][q][i][comp_i]) *
                                    JxW;

                                  local_matrix(i, j) += -beta_terms;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template <int dim>
void
OutletBoundaryCondition<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  const FiniteElement<dim> &fe = scratch_data.fe_face_values.get_fe();

  const double penalty_parameter =
    1. / std::pow(scratch_data.cell_size, fe.degree + 1);
  auto &local_rhs = copy_data.local_rhs;

  // Loop over the faces of the cell.
  for (unsigned int f = 0; f < scratch_data.n_faces; ++f)
    {
      // Check if the face is on a boundary
      if (scratch_data.is_boundary_face[f])
        {
          types::boundary_id boundary_id = scratch_data.boundary_face_id[f];

          // Check if the face is part of the boundary that has a
          // weakly imposed Dirichlet BC.
          if (this->boundary_conditions.type.at(boundary_id) ==
              BoundaryConditions::BoundaryType::outlet)
            {
              const double beta = boundary_conditions.beta.at(boundary_id);

              for (unsigned int q = 0; q < scratch_data.n_faces_q_points; ++q)
                {
                  const double JxW = scratch_data.face_JxW[f][q];
                  for (const unsigned int i :
                       scratch_data.fe_face_values.dof_indices())
                    {
                      // Calculate beta term depending on the
                      // value of  u*n. If it is positive (outgoing
                      // flow) then
                      double normal_outflux =
                        std::min(0.,
                                 (scratch_data.face_velocity_values[f][q] *
                                  scratch_data.face_normal[f][q]));

                      const auto comp_i = fe.system_to_component_index(i).first;
                      if (comp_i < dim)
                        {
                          double beta_terms =
                            penalty_parameter * beta * normal_outflux *
                            (scratch_data.face_velocity_values[f][q][comp_i] *
                             scratch_data.face_phi_u[f][q][i][comp_i]) *
                            JxW;

                          local_rhs(i) += +beta_terms;
                        }
                    }
                }
            }
        }
    }
}

template class OutletBoundaryCondition<2>;
template class OutletBoundaryCondition<3>;


template <int dim>
void
NavierStokesAssemblerALE<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  /// Loop and quadrature informations
  const auto                   &JxW_vec    = scratch_data.JxW;
  const unsigned int            n_q_points = scratch_data.n_q_points;
  const std::vector<Point<dim>> quadrature_points =
    scratch_data.quadrature_points;
  const unsigned int n_dofs = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix        = copy_data.local_matrix;

  // ALE components
  Tensor<1, dim>                                  velocity_ale;
  std::shared_ptr<Functions::ParsedFunction<dim>> velocity_ale_function =
    ale.velocity;
  Vector<double> velocity_ale_vector(dim);

  // assembling local matrix and right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      velocity_ale_function->vector_value(quadrature_points[q],
                                          velocity_ale_vector);
      for (unsigned int d = 0; d < dim; ++d)
        velocity_ale[d] = velocity_ale_vector[d];

      // Store JxW in local variable for faster access
      const double JxW = JxW_vec[q];

      // Calculate strong residual vector
      strong_residual_vec[q] +=
        -scratch_data.velocity_gradients[q] * velocity_ale;

      // Strong residual jacobian calculation
      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian_vec[q][j] +=
            -scratch_data.grad_phi_u[q][j] * velocity_ale;
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i = scratch_data.phi_u[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const Tensor<2, dim> grad_phi_u_j = scratch_data.grad_phi_u[q][j];

              // Weak form for : -u_ALE * gradu
              local_matrix(i, j) +=
                -phi_u_i * (grad_phi_u_j * velocity_ale) * JxW;
            }
        }

    } // end loop on quadrature points
}

template <int dim>
void
NavierStokesAssemblerALE<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature informations
  const auto                   &JxW_vec    = scratch_data.JxW;
  const unsigned int            n_q_points = scratch_data.n_q_points;
  const std::vector<Point<dim>> quadrature_points =
    scratch_data.quadrature_points;
  const unsigned int n_dofs = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs           = copy_data.local_rhs;

  // ALE components
  Tensor<1, dim>                                  velocity_ale;
  std::shared_ptr<Functions::ParsedFunction<dim>> velocity_ale_function =
    ale.velocity;
  Vector<double> velocity_ale_vector(dim);

  // assembling local matrix and right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      velocity_ale_function->vector_value(quadrature_points[q],
                                          velocity_ale_vector);
      for (unsigned int d = 0; d < dim; ++d)
        velocity_ale[d] = velocity_ale_vector[d];

      // Store JxW in local variable for faster access
      const double JxW = JxW_vec[q];

      // Calculate strong residual vector
      strong_residual_vec[q] +=
        -scratch_data.velocity_gradients[q] * velocity_ale;

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          local_rhs[i] +=
            (scratch_data.phi_u[q][i] *
             (scratch_data.velocity_gradients[q] * velocity_ale)) *
            JxW;
        }
    } // end loop on quadrature points
}

template class NavierStokesAssemblerALE<3>;
template class NavierStokesAssemblerALE<2>;
