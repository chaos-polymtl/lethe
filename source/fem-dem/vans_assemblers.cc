// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/simulation_control.h>
#include <core/utilities.h>

#include <solvers/stabilization.h>

#include <fem-dem/vans_assemblers.h>

#include <deal.II/base/tensor.h>

template <int dim>
void
VANSAssemblerCoreModelB<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Kinematic viscosity at Gauss points
  const std::vector<double> &viscosity_vector =
    scratch_data.kinematic_viscosity;

  // Loop and quadrature informations
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

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const double         kinematic_viscosity = viscosity_vector[q];
      const Tensor<1, dim> velocity = scratch_data.velocity_values[q];
      const Tensor<1, dim> previous_velocity =
        scratch_data.previous_velocity_values[0][q];
      const Tensor<2, dim> velocity_gradient =
        scratch_data.velocity_gradients[q];
      const Tensor<1, dim> velocity_laplacian =
        scratch_data.velocity_laplacians[q];

      const Tensor<1, dim> pressure_gradient =
        scratch_data.pressure_gradients[q];

      const double         void_fraction = scratch_data.void_fraction_values[q];
      const Tensor<1, dim> void_fraction_gradients =
        scratch_data.void_fraction_gradient_values[q];

      // Forcing term
      const Tensor<1, dim> force       = scratch_data.force[q];
      double               mass_source = scratch_data.mass_source[q];

      double u_mag = 0;
      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      if (cfd_dem.implicit_stabilization)
        u_mag = std::max(velocity.norm(), 1e-12);
      else
        {
          if (this->simulation_control->get_current_time() ==
              this->simulation_control->get_time_step())
            u_mag = std::max(velocity.norm(), 1e-12);
          else
            u_mag = std::max(previous_velocity.norm(), 1e-12);
        }

      // Grad-div weight factor
      const double gamma =
        calculate_gamma(u_mag, kinematic_viscosity, h, cfd_dem.cstar);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(u_mag,
                                                 kinematic_viscosity,
                                                 h) :
          calculate_navier_stokes_gls_tau_transient(u_mag,
                                                    kinematic_viscosity,
                                                    h,
                                                    sdt);

      // Calculate the strong residual for GLS stabilization
      auto strong_residual = velocity_gradient * velocity * void_fraction +
                             // Mass Source
                             mass_source * velocity
                             // Pressure
                             + pressure_gradient -
                             // Kinematic viscosity and Force
                             kinematic_viscosity * velocity_laplacian -
                             force * void_fraction + strong_residual_vec[q];

      // Pressure scaling factor
      const double pressure_scaling_factor =
        scratch_data.pressure_scaling_factor;

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
            (velocity_gradient * phi_u_j * void_fraction +
             grad_phi_u_j * velocity * void_fraction +
             // Mass Source
             mass_source * phi_u_j +
             // Pressure
             grad_phi_p_j
             // Kinematic viscosity
             - kinematic_viscosity * laplacian_phi_u_j);
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const int   component_i  = scratch_data.components[i];
          const auto &phi_u_i      = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &div_phi_u_i  = scratch_data.div_phi_u[q][i];
          const auto &phi_p_i      = scratch_data.phi_p[q][i];
          const auto &grad_phi_p_i = scratch_data.grad_phi_p[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const int   component_j  = scratch_data.components[j];
              const auto &phi_u_j      = scratch_data.phi_u[q][j];
              const auto &grad_phi_u_j = scratch_data.grad_phi_u[q][j];
              const auto &div_phi_u_j  = scratch_data.div_phi_u[q][j];

              const auto &phi_p_j =
                pressure_scaling_factor * scratch_data.phi_p[q][j];

              const auto &strong_jac = strong_jacobian_vec[q][j];

              // Pressure
              double local_matrix_ij =
                component_j == dim ? -div_phi_u_i * phi_p_j : 0;

              if (component_i == dim)
                {
                  // Continuity
                  local_matrix_ij +=
                    phi_p_i * ((void_fraction * div_phi_u_j) +
                               (phi_u_j * void_fraction_gradients));

                  // PSPG GLS term
                  local_matrix_ij += tau * (strong_jac * grad_phi_p_i);
                }

              if (component_i < dim && component_j < dim)
                {
                  // Convection
                  local_matrix_ij +=
                    ((phi_u_j * void_fraction * velocity_gradient * phi_u_i) +
                     (grad_phi_u_j * void_fraction * velocity * phi_u_i));

                  if (component_i == component_j)
                    {
                      local_matrix_ij +=
                        // Deviatoric stress tensor
                        kinematic_viscosity * (grad_phi_u_j[component_j] *
                                               grad_phi_u_i[component_i])
                        // Mass source
                        + mass_source * phi_u_j * phi_u_i;
                    }
                }

              // The jacobian matrix for the SUPG formulation
              // currently does not include the jacobian of the stabilization
              // parameter tau. Our experience has shown that does not alter the
              // number of newton iteration for convergence, but greatly
              // simplifies assembly.
              if (SUPG && component_i < dim)
                {
                  local_matrix_ij +=
                    tau * (strong_jac * grad_phi_u_i * velocity +
                           strong_residual * grad_phi_u_i * phi_u_j);
                }

              // Grad-div stabilization
              if (cfd_dem.grad_div == true)
                {
                  local_matrix_ij += gamma *
                                     (div_phi_u_j * void_fraction +
                                      phi_u_j * void_fraction_gradients) *
                                     div_phi_u_i;
                }

              local_matrix_ij *= JxW;
              local_matrix(i, j) += local_matrix_ij;
            }
        }
    }
}

template <int dim>
void
VANSAssemblerCoreModelB<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &viscosity_vector =
    scratch_data.kinematic_viscosity;

  // Loop and quadrature informations
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
      // Physical property
      const double kinematic_viscosity = viscosity_vector[q];

      // Velocity
      const Tensor<1, dim> velocity = scratch_data.velocity_values[q];
      const Tensor<1, dim> previous_velocity =
        scratch_data.previous_velocity_values[0][q];
      const double velocity_divergence = scratch_data.velocity_divergences[q];
      const Tensor<2, dim> velocity_gradient =
        scratch_data.velocity_gradients[q];
      const Tensor<1, dim> velocity_laplacian =
        scratch_data.velocity_laplacians[q];

      // Pressure
      const double         pressure = scratch_data.pressure_values[q];
      const Tensor<1, dim> pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Void Fraction
      const double         void_fraction = scratch_data.void_fraction_values[q];
      const Tensor<1, dim> void_fraction_gradients =
        scratch_data.void_fraction_gradient_values[q];

      // Forcing term
      const Tensor<1, dim> force       = scratch_data.force[q];
      double               mass_source = scratch_data.mass_source[q];

      double u_mag = 0;
      // Calculation of the magnitude of the
      // velocity for the stabilization parameter
      if (cfd_dem.implicit_stabilization)
        u_mag = std::max(velocity.norm(), 1e-12);
      else
        {
          if (this->simulation_control->get_current_time() ==
              this->simulation_control->get_time_step())
            u_mag = std::max(velocity.norm(), 1e-12);
          else
            u_mag = std::max(previous_velocity.norm(), 1e-12);
        }

      // Grad-div weight factor
      const double gamma =
        calculate_gamma(u_mag, kinematic_viscosity, h, cfd_dem.cstar);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(u_mag,
                                                 kinematic_viscosity,
                                                 h) :
          calculate_navier_stokes_gls_tau_transient(u_mag,
                                                    kinematic_viscosity,
                                                    h,
                                                    sdt);

      // Calculate the strong residual for GLS stabilization
      auto strong_residual = velocity_gradient * velocity * void_fraction +
                             // Mass Source
                             mass_source * velocity
                             // Pressure
                             + pressure_gradient -
                             // Kinematic viscosity and Force
                             kinematic_viscosity * velocity_laplacian -
                             force * void_fraction + strong_residual_vec[q];

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const int  component_i  = scratch_data.components[i];
          const auto phi_u_i      = scratch_data.phi_u[q][i];
          const auto grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto phi_p_i      = scratch_data.phi_p[q][i];
          const auto grad_phi_p_i = scratch_data.grad_phi_p[q][i];
          const auto div_phi_u_i  = scratch_data.div_phi_u[q][i];

          double local_rhs_i = 0;

          // Navier-Stokes Residual
          if (component_i < dim)
            local_rhs_i +=
              // Momentum
              -kinematic_viscosity *
                scalar_product(velocity_gradient, grad_phi_u_i) -
              velocity_gradient * velocity * void_fraction * phi_u_i
              // Mass Source
              - mass_source * velocity * phi_u_i
              // Pressure and Force
              + pressure * div_phi_u_i + force * void_fraction * phi_u_i;

          if (component_i == dim)
            // Continuity
            local_rhs_i += -(velocity_divergence * void_fraction +
                             velocity * void_fraction_gradients - mass_source) *
                           phi_p_i;

          // PSPG GLS term
          local_rhs_i += -tau * (strong_residual * grad_phi_p_i);

          // SUPG GLS term
          if (SUPG && component_i < dim)
            {
              local_rhs_i +=
                -tau * (strong_residual * (grad_phi_u_i * velocity));
            }

          // Grad-div stabilization
          if (cfd_dem.grad_div == true)
            {
              local_rhs_i -= gamma *
                             (void_fraction * velocity_divergence +
                              velocity * void_fraction_gradients) *
                             div_phi_u_i;
            }

          local_rhs(i) += local_rhs_i * JxW;
        }
    }
}

template class VANSAssemblerCoreModelB<2>;
template class VANSAssemblerCoreModelB<3>;

template <int dim>
void
VANSAssemblerCoreModelA<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &viscosity_vector =
    scratch_data.kinematic_viscosity;

  // Loop and quadrature informations
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

  // Grad-div weight factor

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Physical properties
      const double kinematic_viscosity = viscosity_vector[q];

      // Gather into local variables the relevant fields
      const Tensor<1, dim> velocity = scratch_data.velocity_values[q];
      const Tensor<1, dim> previous_velocity =
        scratch_data.previous_velocity_values[0][q];
      const Tensor<2, dim> velocity_gradient =
        scratch_data.velocity_gradients[q];
      const Tensor<1, dim> velocity_laplacian =
        scratch_data.velocity_laplacians[q];

      const Tensor<1, dim> pressure_gradient =
        scratch_data.pressure_gradients[q];

      const double         void_fraction = scratch_data.void_fraction_values[q];
      const Tensor<1, dim> void_fraction_gradients =
        scratch_data.void_fraction_gradient_values[q];

      // Forcing term
      const Tensor<1, dim> force       = scratch_data.force[q];
      double               mass_source = scratch_data.mass_source[q];

      double u_mag = 0;
      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      if (cfd_dem.implicit_stabilization)
        u_mag = std::max(velocity.norm(), 1e-12);
      else
        {
          if (this->simulation_control->get_current_time() ==
              this->simulation_control->get_time_step())
            u_mag = std::max(velocity.norm(), 1e-12);
          else
            u_mag = std::max(previous_velocity.norm(), 1e-12);
        }

      // Grad-div weight factor
      const double gamma =
        calculate_gamma(u_mag, kinematic_viscosity, h, cfd_dem.cstar);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(u_mag,
                                                 kinematic_viscosity,
                                                 h) :
          calculate_navier_stokes_gls_tau_transient(u_mag,
                                                    kinematic_viscosity,
                                                    h,
                                                    sdt);

      // Calculate the strong residual for GLS stabilization
      auto strong_residual =
        velocity_gradient * velocity * void_fraction +
        // Mass Source
        mass_source * velocity
        // Pressure
        + void_fraction * pressure_gradient -
        // Kinematic viscosity and Force
        void_fraction * kinematic_viscosity * velocity_laplacian -
        force * void_fraction + strong_residual_vec[q];

      // Pressure scaling factor
      const double pressure_scaling_factor =
        scratch_data.pressure_scaling_factor;

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
            (velocity_gradient * phi_u_j * void_fraction +
             grad_phi_u_j * velocity * void_fraction +
             // Mass Source
             mass_source * phi_u_j +
             // Pressure
             void_fraction * grad_phi_p_j
             // Kinematic viscosity
             - void_fraction * kinematic_viscosity * laplacian_phi_u_j);
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const int   component_i  = scratch_data.components[i];
          const auto &phi_u_i      = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &div_phi_u_i  = scratch_data.div_phi_u[q][i];
          const auto &phi_p_i      = scratch_data.phi_p[q][i];
          const auto &grad_phi_p_i = scratch_data.grad_phi_p[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const int   component_j  = scratch_data.components[j];
              const auto &phi_u_j      = scratch_data.phi_u[q][j];
              const auto &grad_phi_u_j = scratch_data.grad_phi_u[q][j];
              const auto &div_phi_u_j  = scratch_data.div_phi_u[q][j];
              const auto &phi_p_j =
                pressure_scaling_factor * scratch_data.phi_p[q][j];

              const auto &strong_jac = strong_jacobian_vec[q][j];

              // Pressure
              double local_matrix_ij =
                -(div_phi_u_i * phi_p_j * void_fraction +
                  phi_p_j * void_fraction_gradients * phi_u_i);

              if (component_i == dim)
                {
                  local_matrix_ij +=
                    phi_p_i * ((void_fraction * div_phi_u_j) +
                               (phi_u_j * void_fraction_gradients));

                  // PSPG GLS term
                  local_matrix_ij += tau * (strong_jac * grad_phi_p_i);
                }

              if (component_i < dim && component_j < dim)
                {
                  // Convection
                  local_matrix_ij +=
                    ((phi_u_j * void_fraction * velocity_gradient * phi_u_i) +
                     (grad_phi_u_j * void_fraction * velocity * phi_u_i));

                  if (component_i == component_j)
                    {
                      local_matrix_ij += void_fraction * kinematic_viscosity *
                                           (grad_phi_u_j[component_j] *
                                            grad_phi_u_i[component_i]) +
                                         kinematic_viscosity * grad_phi_u_j *
                                           void_fraction_gradients * phi_u_i +
                                         mass_source * phi_u_j * phi_u_i;
                    }
                }

              // The jacobian matrix for the SUPG formulation
              // currently does not include the jacobian of the stabilization
              // parameter tau. Our experience has shown that does not alter the
              // number of newton iteration for convergence, but greatly
              // simplifies assembly.
              if (SUPG && component_i < dim)
                {
                  local_matrix_ij +=
                    tau * (strong_jac * grad_phi_u_i * velocity +
                           strong_residual * grad_phi_u_i * phi_u_j);
                }

              // Grad-div stabilization
              if (cfd_dem.grad_div == true)
                {
                  local_matrix_ij += gamma *
                                     (div_phi_u_j * void_fraction +
                                      phi_u_j * void_fraction_gradients) *
                                     div_phi_u_i;
                }

              local_matrix_ij *= JxW;
              local_matrix(i, j) += local_matrix_ij;
            }
        }
    }
}

template <int dim>
void
VANSAssemblerCoreModelA<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &viscosity_vector =
    scratch_data.kinematic_viscosity;

  // Loop and quadrature informations
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
      const Tensor<1, dim> velocity = scratch_data.velocity_values[q];
      const Tensor<1, dim> previous_velocity =
        scratch_data.previous_velocity_values[0][q];
      const double velocity_divergence = scratch_data.velocity_divergences[q];
      const Tensor<2, dim> velocity_gradient =
        scratch_data.velocity_gradients[q];
      const Tensor<1, dim> velocity_laplacian =
        scratch_data.velocity_laplacians[q];

      // Pressure
      const double         pressure = scratch_data.pressure_values[q];
      const Tensor<1, dim> pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Void Fraction
      const double         void_fraction = scratch_data.void_fraction_values[q];
      const Tensor<1, dim> void_fraction_gradients =
        scratch_data.void_fraction_gradient_values[q];

      // Forcing term
      const Tensor<1, dim> force       = scratch_data.force[q];
      double               mass_source = scratch_data.mass_source[q];

      double u_mag = 0;
      // Calculation of the magnitude of the
      // velocity for the stabilization parameter
      if (cfd_dem.implicit_stabilization)
        u_mag = std::max(velocity.norm(), 1e-12);
      else
        {
          if (this->simulation_control->get_current_time() ==
              this->simulation_control->get_time_step())
            u_mag = std::max(velocity.norm(), 1e-12);
          else
            u_mag = std::max(previous_velocity.norm(), 1e-12);
        }

      // Grad-div weight factor
      const double gamma =
        calculate_gamma(u_mag, kinematic_viscosity, h, cfd_dem.cstar);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(u_mag,
                                                 kinematic_viscosity,
                                                 h) :
          calculate_navier_stokes_gls_tau_transient(u_mag,
                                                    kinematic_viscosity,
                                                    h,
                                                    sdt);

      // Calculate the strong residual for GLS stabilization
      auto strong_residual =
        velocity_gradient * velocity * void_fraction +
        // Mass Source
        mass_source * velocity
        // Pressure
        + void_fraction * pressure_gradient -
        // Kinematic viscosity and Force
        void_fraction * kinematic_viscosity * velocity_laplacian -
        force * void_fraction + strong_residual_vec[q];

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const int component_i = scratch_data.components[i];

          const auto phi_u_i      = scratch_data.phi_u[q][i];
          const auto grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto phi_p_i      = scratch_data.phi_p[q][i];
          const auto grad_phi_p_i = scratch_data.grad_phi_p[q][i];
          const auto div_phi_u_i  = scratch_data.div_phi_u[q][i];

          // Navier-Stokes Residual
          double local_rhs_i = 0;

          if (component_i < dim)
            local_rhs_i += (
              // Momentum
              -(void_fraction * kinematic_viscosity *
                  scalar_product(velocity_gradient, grad_phi_u_i) +
                kinematic_viscosity * velocity_gradient *
                  void_fraction_gradients * phi_u_i) -
              velocity_gradient * velocity * void_fraction * phi_u_i
              // Mass Source
              - mass_source * velocity * phi_u_i
              // Pressure and Force
              + (void_fraction * pressure * div_phi_u_i +
                 pressure * void_fraction_gradients * phi_u_i) +
              force * void_fraction * phi_u_i);

          if (component_i == dim)
            // Continuity
            local_rhs_i -= (velocity_divergence * void_fraction +
                            velocity * void_fraction_gradients - mass_source) *
                           phi_p_i;

          // PSPG GLS term
          local_rhs_i += -tau * (strong_residual * grad_phi_p_i);

          // SUPG GLS term
          if (SUPG)
            {
              local_rhs_i +=
                -tau * (strong_residual * (grad_phi_u_i * velocity));
            }

          // Grad-div stabilization
          if (cfd_dem.grad_div == true)
            {
              local_rhs_i -= gamma *
                             (void_fraction * velocity_divergence +
                              velocity * void_fraction_gradients) *
                             div_phi_u_i;
            }

          local_rhs_i *= JxW;
          local_rhs(i) += local_rhs_i;
        }
    }
}

template class VANSAssemblerCoreModelA<2>;
template class VANSAssemblerCoreModelA<3>;

template <int dim>
void
VANSAssemblerBDF<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature informations
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

      double void_fraction = scratch_data.void_fraction_values[q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += bdf_coefs[p] * velocity[p] * void_fraction;
        }

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian[q][j] +=
            void_fraction * bdf_coefs[0] * scratch_data.phi_u[q][j];
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const Tensor<1, dim> &phi_u_i = scratch_data.phi_u[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const Tensor<1, dim> &phi_u_j = scratch_data.phi_u[q][j];

              local_matrix(i, j) +=
                void_fraction * phi_u_j * phi_u_i * bdf_coefs[0] * JxW[q];
            }
        }
    }
}

template <int dim>
void
VANSAssemblerBDF<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Physical properties
  const std::vector<double> &kinematic_viscosity =
    scratch_data.kinematic_viscosity;

  // Loop and quadrature informations
  const double       h          = scratch_data.cell_size;
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
  std::vector<double> void_fraction(1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      velocity[0] = scratch_data.velocity_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        velocity[p + 1] = scratch_data.previous_velocity_values[p][q];

      void_fraction[0] = scratch_data.void_fraction_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        void_fraction[p + 1] = scratch_data.previous_void_fraction_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += bdf_coefs[p] * velocity[p] * void_fraction[0];
        }
      double u_mag = 0;

      if (cfd_dem.implicit_stabilization)
        u_mag = std::max(velocity[0].norm(), 1e-12);
      else
        {
          // Grad-div weight factor used to be constant 0.1
          if (this->simulation_control->get_current_time() ==
              this->simulation_control->get_time_step())
            u_mag = std::max(velocity[0].norm(), 1e-12);
          else
            u_mag = std::max(velocity[1].norm(), 1e-12);
        }

      // Grad-div weight factor
      const double gamma =
        calculate_gamma(u_mag, kinematic_viscosity[q], h, cfd_dem.cstar);

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i     = scratch_data.phi_u[q][i];
          const auto div_phi_u_i = scratch_data.div_phi_u[q][i];
          double     local_rhs_i = 0;

          for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
               ++p)
            {
              local_rhs_i -=
                (bdf_coefs[p] * void_fraction[0] * (velocity[p] * phi_u_i));

              if (cfd_dem.void_fraction_time_derivative == true)
                {
                  local_rhs_i -= bdf_coefs[p] *
                                 (void_fraction[p] * scratch_data.phi_p[q][i]);

                  if (cfd_dem.grad_div == true)
                    {
                      local_rhs_i -=
                        gamma * (bdf_coefs[p] * void_fraction[p]) * div_phi_u_i;
                    }
                }
            }
          local_rhs(i) += local_rhs_i * JxW[q];
        }
    }
}

template class VANSAssemblerBDF<2>;
template class VANSAssemblerBDF<3>;

template <int dim>
void
VANSAssemblerDiFelice<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "VANSAssemblerDiFelice<dim>::calculate_particle_fluid_interactions"));

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "VANSAssemblerDiFelice<dim>::calculate_particle_fluid_interactions"));

  double      cell_void_fraction = 0;
  double      C_d                = 0;
  const auto &relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p      = scratch_data.Re_particle;
  const auto &density   = scratch_data.density_at_particle_location;
  auto       &beta_drag = scratch_data.beta_drag;

  // The explicit_particle_volumetric_acceleration_on_fluid variable is used to
  // store the drag if the coupling is explicit.
  auto &explicit_particle_volumetric_acceleration_on_fluid =
    scratch_data.explicit_particle_volumetric_acceleration_on_fluid;

  Tensor<1, dim> drag_force;

  const auto pic = scratch_data.pic;
  beta_drag      = 0;
  // i_particle is an index that runs from 0 to n_particles_in_cell.
  // It is used to access vectors of particle-specific variables and dependent
  // quantities needed to compute the fluid’s effect on each particle in the
  // cell.

  unsigned int i_particle = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      cell_void_fraction =
        std::min(scratch_data.cell_void_fraction[i_particle], 1.0);

      // Di Felice Drag Model CD Calculation
      C_d = Utilities::fixed_power<2>((0.63 + 4.8 / sqrt(Re_p[i_particle]))) *
            pow(cell_void_fraction,
                2 - (3.7 - 0.65 * exp(-Utilities::fixed_power<2>(
                                        1.5 - log10(Re_p[i_particle])) *
                                      0.5)));

      double momentum_transfer_coefficient =
        (0.5 * C_d * M_PI *
         Utilities::fixed_power<2>(
           particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp]) /
         4) *
        relative_velocity[i_particle].norm();

      if (cfd_dem.drag_coupling == Parameters::DragCoupling::fully_implicit ||
          cfd_dem.drag_coupling == Parameters::DragCoupling::semi_implicit)
        {
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                momentum_transfer_coefficient] =
            momentum_transfer_coefficient;
          beta_drag += momentum_transfer_coefficient;
        }
      else // Explicit drag model.
        {
          explicit_particle_volumetric_acceleration_on_fluid -=
            momentum_transfer_coefficient * relative_velocity[i_particle] /
            scratch_data.cell_volume;
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                momentum_transfer_coefficient] = 0;
        }

      drag_force = density[i_particle] * momentum_transfer_coefficient *
                   relative_velocity[i_particle];

      for (int d = 0; d < dim; ++d)
        {
          particle_properties
            [DEM::CFDDEMProperties::PropertiesIndex::fem_drag_x + d] +=
            drag_force[d];
        }

      i_particle += 1;
    }

  beta_drag = beta_drag / scratch_data.cell_volume;
}

template class VANSAssemblerDiFelice<2>;
template class VANSAssemblerDiFelice<3>;

template <int dim>
void
VANSAssemblerRong<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // Physical Properties
  Assert(!scratch_data.properties_manager.is_non_newtonian(),
         RequiresConstantViscosity(
           "VANSAssemblerRong<dim>::calculate_particle_fluid_interactions"));

  Assert(scratch_data.properties_manager.density_is_constant(),
         RequiresConstantDensity(
           "VANSAssemblerRong<dim>::calculate_particle_fluid_interactions"));

  double      cell_void_fraction = 0;
  double      C_d                = 0;
  const auto &relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p      = scratch_data.Re_particle;
  const auto &density   = scratch_data.density_at_particle_location;
  auto       &beta_drag = scratch_data.beta_drag;

  // The explicit_particle_volumetric_acceleration_on_fluid variable is used to
  // store the drag if the coupling is explicit.
  auto &explicit_particle_volumetric_acceleration_on_fluid =
    scratch_data.explicit_particle_volumetric_acceleration_on_fluid;

  Tensor<1, dim> drag_force;

  const auto pic          = scratch_data.pic;
  beta_drag               = 0;
  unsigned int i_particle = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      cell_void_fraction =
        std::min(scratch_data.cell_void_fraction[i_particle], 1.0);

      // Rong Drag Model CD Calculation
      C_d =
        Utilities::fixed_power<2>((0.63 + 4.8 / sqrt(Re_p[i_particle]))) *
        pow(cell_void_fraction,
            2 -
              (2.65 * (cell_void_fraction + 1) -
               (5.3 - (3.5 * cell_void_fraction)) *
                 Utilities::fixed_power<2>(cell_void_fraction) *
                 exp(-Utilities::fixed_power<2>(1.5 - log10(Re_p[i_particle])) *
                     0.5)));

      double momentum_transfer_coefficient =
        (0.5 * C_d * M_PI *
         Utilities::fixed_power<2>(
           particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp]) *
         0.25) *
        relative_velocity[i_particle].norm();

      if (cfd_dem.drag_coupling == Parameters::DragCoupling::fully_implicit ||
          cfd_dem.drag_coupling == Parameters::DragCoupling::semi_implicit)
        {
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                momentum_transfer_coefficient] =
            momentum_transfer_coefficient;
          beta_drag += momentum_transfer_coefficient;
        }
      else // Explicit drag model.
        {
          explicit_particle_volumetric_acceleration_on_fluid -=
            momentum_transfer_coefficient * relative_velocity[i_particle] /
            scratch_data.cell_volume;
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                momentum_transfer_coefficient] = 0;
        }
      drag_force = density[i_particle] * momentum_transfer_coefficient *
                   relative_velocity[i_particle];

      for (int d = 0; d < dim; ++d)
        {
          particle_properties
            [DEM::CFDDEMProperties::PropertiesIndex::fem_drag_x + d] +=
            drag_force[d];
        }

      i_particle += 1;
    }

  beta_drag = beta_drag / scratch_data.cell_volume;
}

template class VANSAssemblerRong<2>;
template class VANSAssemblerRong<3>;

template <int dim>
void
VANSAssemblerDallavalle<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "VANSAssemblerDallavalle<dim>::calculate_particle_fluid_interactions"));

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "VANSAssemblerDallavalle<dim>::calculate_particle_fluid_interactions"));

  double      C_d = 0;
  const auto &relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p      = scratch_data.Re_particle;
  const auto &density   = scratch_data.density_at_particle_location;
  auto       &beta_drag = scratch_data.beta_drag;

  // The explicit_particle_volumetric_acceleration_on_fluid variable is used to
  // store the drag if the coupling is explicit.
  auto &explicit_particle_volumetric_acceleration_on_fluid =
    scratch_data.explicit_particle_volumetric_acceleration_on_fluid;

  Tensor<1, dim> drag_force;

  const auto pic          = scratch_data.pic;
  beta_drag               = 0;
  unsigned int i_particle = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      // Dallavalle Drag Model CD Calculation
      C_d = Utilities::fixed_power<2>(0.63 + 4.8 / sqrt(Re_p[i_particle]));

      double momentum_transfer_coefficient =
        (0.5 * C_d * M_PI *
         Utilities::fixed_power<2>(
           particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp]) *
         0.25) *
        relative_velocity[i_particle].norm();

      if (cfd_dem.drag_coupling == Parameters::DragCoupling::fully_implicit ||
          cfd_dem.drag_coupling == Parameters::DragCoupling::semi_implicit)
        {
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                momentum_transfer_coefficient] =
            momentum_transfer_coefficient;
          beta_drag += momentum_transfer_coefficient;
        }
      else // Explicit drag model.
        {
          explicit_particle_volumetric_acceleration_on_fluid -=
            momentum_transfer_coefficient * relative_velocity[i_particle] /
            scratch_data.cell_volume;
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                momentum_transfer_coefficient] = 0;
        }

      drag_force = density[i_particle] * momentum_transfer_coefficient *
                   relative_velocity[i_particle];

      for (int d = 0; d < dim; ++d)
        {
          particle_properties
            [DEM::CFDDEMProperties::PropertiesIndex::fem_drag_x + d] +=
            drag_force[d];
        }

      i_particle += 1;
    }

  beta_drag = beta_drag / scratch_data.cell_volume;
}

template class VANSAssemblerDallavalle<2>;
template class VANSAssemblerDallavalle<3>;

template <int dim>
void
VANSAssemblerKochHill<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "VANSAssemblerKochHill<dim>::calculate_particle_fluid_interactions"));

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "VANSAssemblerKochHill<dim>::calculate_particle_fluid_interactions"));

  double      cell_void_fraction = 0;
  const auto &relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p    = scratch_data.Re_particle;
  const auto &density = scratch_data.density_at_particle_location;
  const auto &kinematic_viscosity =
    scratch_data.kinematic_viscosity_at_particle_location;
  auto &beta_drag = scratch_data.beta_drag;

  // The explicit_particle_volumetric_acceleration_on_fluid variable is used to
  // store the drag if the coupling is explicit.
  auto &explicit_particle_volumetric_acceleration_on_fluid =
    scratch_data.explicit_particle_volumetric_acceleration_on_fluid;

  Tensor<1, dim> drag_force;

  const auto pic          = scratch_data.pic;
  beta_drag               = 0;
  unsigned int i_particle = 0;

  double f0 = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      cell_void_fraction =
        std::min(scratch_data.cell_void_fraction[i_particle], 1.0);

      // Koch and Hill Drag Model Calculation
      if ((1 - cell_void_fraction) < 0.4)
        {
          f0 = (1 + 3 * sqrt((1 - cell_void_fraction + DBL_MIN) / 2) +
                (135.0 / 64) * (1 - cell_void_fraction) *
                  log(1 - cell_void_fraction + DBL_MIN) +
                16.14 * (1 - cell_void_fraction)) /
               (1 + 0.681 * (1 - cell_void_fraction) -
                8.48 * Utilities::fixed_power<2>(1 - cell_void_fraction) +
                8.16 * Utilities::fixed_power<3>(1 - cell_void_fraction));
        }
      else if ((1 - cell_void_fraction) >= 0.4)
        {
          f0 = 10 * (1 - cell_void_fraction) /
               Utilities::fixed_power<3>(cell_void_fraction);
        }

      double f3 = 0.0673 + 0.212 * (1 - cell_void_fraction) +
                  0.0232 / Utilities::fixed_power<5>(cell_void_fraction);

      double momentum_transfer_coefficient =
        ((18 * kinematic_viscosity[i_particle] *
          Utilities::fixed_power<2>(cell_void_fraction) *
          (1 - cell_void_fraction)) /
         Utilities::fixed_power<2>(
           particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp])) *
        (f0 + 0.5 * f3 * Re_p[i_particle]) *
        (M_PI *
         pow(particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp],
             dim) /
         (2 * dim)) /
        (1 - cell_void_fraction + DBL_MIN);

      if (cfd_dem.drag_coupling == Parameters::DragCoupling::fully_implicit ||
          cfd_dem.drag_coupling == Parameters::DragCoupling::semi_implicit)
        {
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                momentum_transfer_coefficient] =
            momentum_transfer_coefficient;
          beta_drag += momentum_transfer_coefficient;
        }
      else // Explicit drag model.
        {
          explicit_particle_volumetric_acceleration_on_fluid -=
            momentum_transfer_coefficient * relative_velocity[i_particle] /
            scratch_data.cell_volume;
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                momentum_transfer_coefficient] = 0;
        }

      drag_force = density[i_particle] * momentum_transfer_coefficient *
                   relative_velocity[i_particle];

      for (int d = 0; d < dim; ++d)
        {
          particle_properties
            [DEM::CFDDEMProperties::PropertiesIndex::fem_drag_x + d] +=
            drag_force[d];
        }

      i_particle += 1;
    }

  beta_drag = beta_drag / scratch_data.cell_volume;
}

template class VANSAssemblerKochHill<2>;
template class VANSAssemblerKochHill<3>;

template <int dim>
void
VANSAssemblerBeetstra<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "VANSAssemblerBeetstra<dim>::calculate_particle_fluid_interactions"));

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "VANSAssemblerBeetstra<dim>::calculate_particle_fluid_interactions"));

  double      cell_void_fraction = 0;
  double      F0                 = 0;
  const auto &relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p    = scratch_data.Re_particle;
  const auto &density = scratch_data.density_at_particle_location;
  const auto &kinematic_viscosity =
    scratch_data.kinematic_viscosity_at_particle_location;
  auto &beta_drag = scratch_data.beta_drag;

  // The explicit_particle_volumetric_acceleration_on_fluid variable is used to
  // store the drag if the coupling is explicit.
  auto &explicit_particle_volumetric_acceleration_on_fluid =
    scratch_data.explicit_particle_volumetric_acceleration_on_fluid;

  Tensor<1, dim> drag_force;

  const auto pic          = scratch_data.pic;
  beta_drag               = 0;
  unsigned int i_particle = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      cell_void_fraction =
        std::min(scratch_data.cell_void_fraction[i_particle], 1.0);

      // Beetstra drag coefficient
      F0 = 10 * (1 - cell_void_fraction) /
             Utilities::fixed_power<2>(cell_void_fraction) +
           Utilities::fixed_power<2>(cell_void_fraction) *
             (1 + 1.5 * sqrt((1 - cell_void_fraction + DBL_MIN))) +
           0.413 * (Re_p[i_particle]) /
             (24 * Utilities::fixed_power<2>(cell_void_fraction)) *
             ((1 / cell_void_fraction) +
              3 * (1 - cell_void_fraction) * cell_void_fraction +
              8.4 * pow((Re_p[i_particle]), -0.343)) /
             (1 + pow(10, 3 * (1 - cell_void_fraction)) *
                    pow((Re_p[i_particle]),
                        -(1 + 4 * (1 - cell_void_fraction)) * 0.5));

      double momentum_transfer_coefficient =
        F0 * 3 * M_PI * kinematic_viscosity[i_particle] * cell_void_fraction *
        particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp];

      if (cfd_dem.drag_coupling == Parameters::DragCoupling::fully_implicit ||
          cfd_dem.drag_coupling == Parameters::DragCoupling::semi_implicit)
        {
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                momentum_transfer_coefficient] =
            momentum_transfer_coefficient;
          beta_drag += momentum_transfer_coefficient;
        }
      else // Explicit drag model.
        {
          explicit_particle_volumetric_acceleration_on_fluid -=
            momentum_transfer_coefficient * relative_velocity[i_particle] /
            scratch_data.cell_volume;
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                momentum_transfer_coefficient] = 0;
        }

      drag_force = density[i_particle] * momentum_transfer_coefficient *
                   relative_velocity[i_particle];

      for (int d = 0; d < dim; ++d)
        {
          particle_properties
            [DEM::CFDDEMProperties::PropertiesIndex::fem_drag_x + d] +=
            drag_force[d];
        }

      i_particle += 1;
    }

  beta_drag = beta_drag / scratch_data.cell_volume;
}

template class VANSAssemblerBeetstra<2>;
template class VANSAssemblerBeetstra<3>;

template <int dim>
void
VANSAssemblerGidaspow<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "VANSAssemblerGidaspow<dim>::calculate_particle_fluid_interactions"));

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "VANSAssemblerGidaspow<dim>::calculate_particle_fluid_interactions"));

  double      cell_void_fraction = 0;
  const auto &relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p    = scratch_data.Re_particle;
  const auto &density = scratch_data.density_at_particle_location;
  const auto &kinematic_viscosity =
    scratch_data.kinematic_viscosity_at_particle_location;
  auto &beta_drag = scratch_data.beta_drag;

  // The explicit_particle_volumetric_acceleration_on_fluid variable is used to
  // store the drag if the coupling is explicit.
  auto &explicit_particle_volumetric_acceleration_on_fluid =
    scratch_data.explicit_particle_volumetric_acceleration_on_fluid;


  Tensor<1, dim> drag_force;

  const auto pic                           = scratch_data.pic;
  double     momentum_transfer_coefficient = 0;
  beta_drag                                = 0;
  unsigned int i_particle                  = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      cell_void_fraction =
        std::min(scratch_data.cell_void_fraction[i_particle], 1.0);

      // Gidaspow Drag Model
      double particle_density =
        particle_properties[DEM::CFDDEMProperties::PropertiesIndex::mass] /
        (M_PI *
         Utilities::fixed_power<dim>(
           particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp]) /
         (2.0 * dim));

      if (cell_void_fraction > 0.8)
        {
          momentum_transfer_coefficient =
            (18 * pow(cell_void_fraction, -3.65) *
             (1 + 0.15 * pow(Re_p[i_particle], 0.687))) *
            (particle_properties[DEM::CFDDEMProperties::PropertiesIndex::mass] *
             kinematic_viscosity[i_particle] /
             (Utilities::fixed_power<2, double>(
                particle_properties
                  [DEM::CFDDEMProperties::PropertiesIndex::dp]) *
              particle_density));
        }
      else
        {
          // Assuming the sphericity of particles = 1
          momentum_transfer_coefficient =
            (150 * (1 - cell_void_fraction) /
               Utilities::fixed_power<2, double>(cell_void_fraction) +
             1.75 * Re_p[i_particle] /
               Utilities::fixed_power<2, double>(cell_void_fraction)) *
            (particle_properties[DEM::CFDDEMProperties::PropertiesIndex::mass] *
             kinematic_viscosity[i_particle] /
             (Utilities::fixed_power<2, double>(
                particle_properties
                  [DEM::CFDDEMProperties::PropertiesIndex::dp]) *
              particle_density));
        }

      if (cfd_dem.drag_coupling == Parameters::DragCoupling::fully_implicit ||
          cfd_dem.drag_coupling == Parameters::DragCoupling::semi_implicit)
        {
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                momentum_transfer_coefficient] =
            momentum_transfer_coefficient;
          beta_drag += momentum_transfer_coefficient;
        }
      else // Explicit drag model.
        {
          explicit_particle_volumetric_acceleration_on_fluid -=
            momentum_transfer_coefficient * relative_velocity[i_particle] /
            scratch_data.cell_volume;
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                momentum_transfer_coefficient] = 0;
        }

      particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                            momentum_transfer_coefficient] =
        momentum_transfer_coefficient;

      drag_force = density[i_particle] * momentum_transfer_coefficient *
                   relative_velocity[i_particle];


      for (int d = 0; d < dim; ++d)
        {
          particle_properties
            [DEM::CFDDEMProperties::PropertiesIndex::fem_drag_x + d] +=
            drag_force[d];
        }

      i_particle += 1;
    }

  beta_drag = beta_drag / scratch_data.cell_volume;
}

template class VANSAssemblerGidaspow<2>;
template class VANSAssemblerGidaspow<3>;

template <int dim>
void
VANSAssemblerSaffmanMei<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "VANSAssemblerSaffmanMei<dim>::calculate_particle_fluid_interactions"));

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "VANSAssemblerSaffmanMei<dim>::calculate_particle_fluid_interactions"));

  // This implementation follows the formulation in the book "Multiphase Flows
  // with Droplets and Particles" by Crowe et al. (2011) and the brief
  // communication article "An approximate expression for the shear lift force
  // on a spherical particle at finite reynolds number" by Mei (1992)
  double C_s   = 0;
  double alpha = 0;

  const auto &relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p    = scratch_data.Re_particle;
  const auto &density = scratch_data.density_at_particle_location;
  const auto &kinematic_viscosity =
    scratch_data.kinematic_viscosity_at_particle_location;

  auto &vorticity_2d =
    scratch_data.fluid_velocity_curls_at_particle_location_2d;
  auto &vorticity_3d =
    scratch_data.fluid_velocity_curls_at_particle_location_3d;
  auto &explicit_particle_volumetric_acceleration_on_fluid =
    scratch_data.explicit_particle_volumetric_acceleration_on_fluid;

  Tensor<1, dim> lift_force;

  const auto   pic        = scratch_data.pic;
  unsigned int i_particle = 0;

  if constexpr (dim == 2)
    {
      for (auto &particle : pic)
        {
          auto particle_properties = particle.get_properties();

          // Saffman-Mei coefficient
          alpha =
            0.5 *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] /
            relative_velocity[i_particle].norm() *
            abs(vorticity_2d[i_particle][0] + 1e-12);

          if (Re_p[i_particle] <= 40)
            C_s = (1 - 0.3314 * sqrt(alpha)) * exp(-0.1 * Re_p[i_particle]) +
                  0.3314 * sqrt(alpha);
          else if (Re_p[i_particle] > 40)
            C_s = 0.0524 * sqrt(alpha * Re_p[i_particle]);

          // Vorticity vector
          Tensor<1, 2> vorticity;
          vorticity[0] = vorticity_2d[i_particle][0];
          vorticity[1] = vorticity_2d[i_particle][1];

          // Saffman Lift force
          lift_force[0] =
            C_s * 1.61 *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] *
            density[i_particle] *
            sqrt(kinematic_viscosity[i_particle] + DBL_MIN) /
            sqrt(vorticity_2d[i_particle].norm()) *
            (relative_velocity[i_particle][0] * vorticity[1]);

          lift_force[1] =
            C_s * 1.61 *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] *
            density[i_particle] *
            sqrt(kinematic_viscosity[i_particle] + DBL_MIN) /
            sqrt(vorticity.norm() + 1e-12) *
            (relative_velocity[i_particle][1] * vorticity[0]);

          for (int d = 0; d < dim; ++d)
            {
              // Apply lift force on the particle
              particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                    fem_force_two_way_coupling_x +
                                  d] += lift_force[d];

              // Apply lift force on the fluid
              explicit_particle_volumetric_acceleration_on_fluid[d] -=
                lift_force[d] /
                (density[i_particle] * scratch_data.cell_volume);
            }

          i_particle += 1;
        }
    }
  else if constexpr (dim == 3)
    {
      for (auto &particle : pic)
        {
          auto particle_properties = particle.get_properties();

          // Saffman-Mei coefficient C_s
          alpha =
            0.5 *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] /
            relative_velocity[i_particle].norm() *
            (vorticity_3d[i_particle].norm() + 1e-12);

          if (Re_p[i_particle] <= 40)
            C_s = (1 - 0.3314 * sqrt(alpha)) * exp(-0.1 * Re_p[i_particle]) +
                  0.3314 * sqrt(alpha);
          else if (Re_p[i_particle] > 40)
            C_s = 0.0524 * sqrt(alpha * Re_p[i_particle]);

          // Vorticity tensor
          Tensor<1, 3> vorticity = vorticity_3d[i_particle];

          lift_force =
            C_s * 1.61 *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] *
            density[i_particle] *
            sqrt(kinematic_viscosity[i_particle] + DBL_MIN) /
            sqrt(vorticity.norm() + 1e-12) *
            (cross_product_3d(relative_velocity[i_particle], vorticity));

          for (int d = 0; d < dim; ++d)
            {
              // Apply lift force on the particle
              particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                    fem_force_two_way_coupling_x +
                                  d] += lift_force[d];

              // Apply lift force on the fluid
              explicit_particle_volumetric_acceleration_on_fluid[d] -=
                lift_force[d] /
                (density[i_particle] * scratch_data.cell_volume);
            }
          i_particle += 1;
        }
    }
}

template class VANSAssemblerSaffmanMei<2>;
template class VANSAssemblerSaffmanMei<3>;

template <int dim>
void
VANSAssemblerMagnus<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // Physical Properties
  Assert(scratch_data.properties_manager.density_is_constant(),
         RequiresConstantDensity(
           "VANSAssemblerMagnus<dim>::calculate_particle_fluid_interactions"));

  // This implementation follows the formulation in the book "Multiphase Flows
  // with Droplets and Particles" by Crowe et al. (2011).
  double C_m = 0.;

  const auto &relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p    = scratch_data.Re_particle;
  const auto &density = scratch_data.density_at_particle_location;

  auto &explicit_particle_volumetric_acceleration_on_fluid =
    scratch_data.explicit_particle_volumetric_acceleration_on_fluid;

  Tensor<1, dim> lift_force;

  const auto   pic        = scratch_data.pic;
  unsigned int i_particle = 0;

  if constexpr (dim == 2)
    {
      for (auto &particle : pic)
        {
          auto particle_properties = particle.get_properties();

          // Add factor to rotational velocity to avoid divisions by zero
          double omega_z = particle_properties
                             [DEM::CFDDEMProperties::PropertiesIndex::omega_z] +
                           1e-15;

          double omega_norm = abs(omega_z);

          // Spin parameter
          double spin_parameter =
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] *
            omega_norm / (2.0 * relative_velocity[i_particle].norm());

          // Magnus lift coefficient
          if (spin_parameter > 1.0 && spin_parameter < 6.0 &&
              Re_p[i_particle] > 10.0 && Re_p[i_particle] < 140.0)
            {
              // Oesterlé and Dinh (1998)
              C_m = 0.45 + (2 * spin_parameter - 0.45) *
                             exp(-0.075 * pow(spin_parameter, 0.4) *
                                 pow(Re_p[i_particle], 0.7));
            }
          else
            C_m = 2.0 * spin_parameter;

          // Magnus Lift force
          lift_force[0] =
            0.5 * C_m *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] *
            density[i_particle] * relative_velocity[i_particle].norm() *
            (omega_z / omega_norm * relative_velocity[i_particle][1]);

          lift_force[1] =
            0.5 * C_m *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] *
            density[i_particle] * relative_velocity[i_particle].norm() *
            (omega_z / omega_norm * relative_velocity[i_particle][0]);

          for (int d = 0; d < dim; ++d)
            {
              // Apply lift force on the particle
              particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                    fem_force_two_way_coupling_x +
                                  d] += lift_force[d];

              // Apply lift force on the fluid
              explicit_particle_volumetric_acceleration_on_fluid[d] -=
                lift_force[d] /
                (density[i_particle] * scratch_data.cell_volume);
            }
          i_particle += 1;
        }
    }

  else if constexpr (dim == 3)
    {
      for (auto &particle : pic)
        {
          auto particle_properties = particle.get_properties();

          Tensor<1, dim> omega;

          for (int d = 0; d < dim; ++d)
            {
              omega[d] =
                particle_properties
                  [DEM::CFDDEMProperties::PropertiesIndex::omega_x + d] +
                1e-15;
            }

          // Spin parameter
          double spin_parameter =
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] *
            omega.norm() / (2.0 * relative_velocity[i_particle].norm());

          // Magnus lift coefficient
          if (spin_parameter > 1.0 && spin_parameter < 6.0 &&
              Re_p[i_particle] > 10.0 && Re_p[i_particle] < 140.0)
            {
              // Oesterlé and Dinh (1998)
              C_m = 0.45 + (2 * spin_parameter - 0.45) *
                             exp(-0.075 * pow(spin_parameter, 0.4) *
                                 pow(Re_p[i_particle], 0.7));
            }
          else
            C_m = 2.0 * spin_parameter;

          Tensor<1, dim> rotational_vector = omega / omega.norm();

          // Magnus Lift force
          lift_force =
            0.125 * M_PI *
            Utilities::fixed_power<2>(
              particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp]) *
            C_m * density[i_particle] * relative_velocity[i_particle].norm() *
            (cross_product_3d(rotational_vector,
                              relative_velocity[i_particle]));

          for (int d = 0; d < dim; ++d)
            {
              // Apply lift force on the particle
              particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                    fem_force_two_way_coupling_x +
                                  d] += lift_force[d];

              // Apply lift force on the fluid
              explicit_particle_volumetric_acceleration_on_fluid[d] -=
                lift_force[d] /
                (density[i_particle] * scratch_data.cell_volume);
            }
          i_particle += 1;
        }
    }
}

template class VANSAssemblerMagnus<2>;
template class VANSAssemblerMagnus<3>;

template <int dim>
void
VANSAssemblerViscousTorque<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "VANSAssemblerViscousTorque<dim>::calculate_particle_fluid_interactions"));

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "VANSAssemblerViscousTorque<dim>::calculate_particle_fluid_interactions"));

  const auto &density = scratch_data.density_at_particle_location;
  const auto &kinematic_viscosity =
    scratch_data.kinematic_viscosity_at_particle_location;
  const auto pic = scratch_data.pic;

  unsigned int i_particle = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      // Extract constant factor to avoid remultiplying it during the loop
      const double factor =
        M_PI *
        Utilities::fixed_power<3, double>(
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp]) *
        kinematic_viscosity[i_particle] * density[i_particle] * 0.5;

      for (int d = 0; d < dim; d++)
        {
          particle_properties
            [DEM::CFDDEMProperties::PropertiesIndex::fem_torque_x + d] -=
            factor * particle_properties
                       [DEM::CFDDEMProperties::PropertiesIndex::omega_x + d];
        }
      i_particle += 1;
    }
}

template class VANSAssemblerViscousTorque<2>;
template class VANSAssemblerViscousTorque<3>;

template <int dim>
void
VANSAssemblerVorticalTorque<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "VANSAssemblerVorticalTorque<dim>::calculate_particle_fluid_interactions"));

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "VANSAssemblerVorticalTorque<dim>::calculate_particle_fluid_interactions"));

  const auto &density = scratch_data.density_at_particle_location;
  const auto &kinematic_viscosity =
    scratch_data.kinematic_viscosity_at_particle_location;

  auto &vorticity_3d =
    scratch_data.fluid_velocity_curls_at_particle_location_3d;

  const auto pic = scratch_data.pic;

  // Local index used to access the local fields calculated at the particle
  // location
  unsigned int i_particle = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      // Extract constant factor to avoid remultiplying it during the loop
      const double factor =
        M_PI *
        Utilities::fixed_power<3, double>(
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp]) *
        kinematic_viscosity[i_particle] * density[i_particle] * 0.5;

      for (int d = 0; d < dim; d++)
        {
          // Calculate and apply viscous torque
          particle_properties
            [DEM::CFDDEMProperties::PropertiesIndex::fem_torque_x + d] +=
            factor * vorticity_3d[i_particle][d];
        }
      i_particle += 1;
    }
}

template class VANSAssemblerVorticalTorque<2>;
template class VANSAssemblerVorticalTorque<3>;

template <int dim>
void
VANSAssemblerBuoyancy<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // Physical Properties
  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "VANSAssemblerBuoyancy<dim>::calculate_particle_fluid_interactions"));

  const auto   pic = scratch_data.pic;
  Tensor<1, 3> buoyancy_force;

  const auto &density = scratch_data.density_at_particle_location;

  unsigned int i_particle = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      // Buoyancy Force
      buoyancy_force =
        -gravity * (4.0 / 3) * M_PI *
        Utilities::fixed_power<3>(
          (particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] /
           2.0));

      for (int d = 0; d < dim; ++d)
        {
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                fem_force_one_way_coupling_x +
                              d] += buoyancy_force[d] * density[i_particle];
        }
      i_particle += 1;
    }
}

template class VANSAssemblerBuoyancy<2>;
template class VANSAssemblerBuoyancy<3>;

template <int dim>
void
VANSAssemblerPressureForce<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // Physical Properties
  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "VANSAssemblerPressureForce<dim>::calculate_particle_fluid_interactions"));

  const auto pic = scratch_data.pic;
  auto      &explicit_particle_volumetric_acceleration_on_fluid =
    scratch_data.explicit_particle_volumetric_acceleration_on_fluid;
  auto pressure_gradients =
    scratch_data.fluid_pressure_gradients_at_particle_location;
  Tensor<1, dim> pressure_force;

  unsigned int i_particle = 0;

  const auto &density = scratch_data.density_at_particle_location;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      // Pressure Force
      pressure_force =
        -(M_PI *
          Utilities::fixed_power<dim>(
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp]) /
          (2 * dim)) *
        pressure_gradients[i_particle];

      for (int d = 0; d < dim; ++d)
        {
          // Apply pressure force to the particles only, when we are solving
          // model A of the VANS. When we are solving Model B, apply the
          // pressure force back on the fluid by lumping it in the
          // explicit_particle_volumetric_acceleration_on_fluid.
          if (cfd_dem.vans_model == Parameters::VANSModel::modelA)
            {
              particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                    fem_force_one_way_coupling_x +
                                  d] += pressure_force[d] * density[i_particle];
            }
          if (cfd_dem.vans_model == Parameters::VANSModel::modelB)
            {
              particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                    fem_force_two_way_coupling_x +
                                  d] += pressure_force[d] * density[i_particle];
              explicit_particle_volumetric_acceleration_on_fluid[d] -=
                pressure_force[d] / scratch_data.cell_volume;
            }
        }

      i_particle += 1;
    }
}

template class VANSAssemblerPressureForce<2>;
template class VANSAssemblerPressureForce<3>;

template <int dim>
void
VANSAssemblerShearForce<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // Kinematic viscosity and density are currently assumed constant within the
  // same fluid phase. Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "VANSAssemblerShearForce<dim>::calculate_particle_fluid_interactions"));

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "VANSAssemblerShearForce<dim>::calculate_particle_fluid_interactions"));

  const auto pic = scratch_data.pic;
  auto      &explicit_particle_volumetric_acceleration_on_fluid =
    scratch_data.explicit_particle_volumetric_acceleration_on_fluid;
  auto &velocity_laplacians =
    scratch_data.fluid_velocity_laplacian_at_particle_location;
  Tensor<1, dim> shear_force;

  unsigned int i_particle = 0;

  const auto &density = scratch_data.density_at_particle_location;
  const std::vector<double> &kinematic_viscosity =
    scratch_data.kinematic_viscosity_at_particle_location;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      // Shear Force
      shear_force =
        -(M_PI *
          Utilities::fixed_power<dim>(
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp]) /
          (2 * dim)) *
        kinematic_viscosity[i_particle] * velocity_laplacians[i_particle];

      for (int d = 0; d < dim; ++d)
        {
          // Apply shear force to the particles only, when we are solving
          // model A of the VANS. When we are solving Model B, apply the shear
          // force back on the fluid by lumping it in the
          // explicit_particle_volumetric_acceleration_on_fluid.

          if (cfd_dem.vans_model == Parameters::VANSModel::modelA)
            {
              particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                    fem_force_one_way_coupling_x +
                                  d] += shear_force[d] * density[i_particle];
            }
          if (cfd_dem.vans_model == Parameters::VANSModel::modelB)
            {
              particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                    fem_force_two_way_coupling_x +
                                  d] += shear_force[d] * density[i_particle];
              explicit_particle_volumetric_acceleration_on_fluid[d] -=
                shear_force[d] / scratch_data.cell_volume;
            }
        }

      i_particle += 1;
    }
}

template class VANSAssemblerShearForce<2>;
template class VANSAssemblerShearForce<3>;

template <int dim>
void
VANSAssemblerFPI<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &strong_jacobian = copy_data.strong_jacobian;
  auto &local_matrix    = copy_data.local_matrix;
  auto  beta_drag       = scratch_data.beta_drag;
  auto &explicit_particle_volumetric_acceleration_on_fluid =
    scratch_data.explicit_particle_volumetric_acceleration_on_fluid;
  const Tensor<1, dim> average_particles_velocity =
    scratch_data.average_particle_velocity;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const Tensor<1, dim> velocity = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Subtraction of forces applied on fluid from the residual for GLS
      // stabilization
      strong_residual[q] -=
        // drag applied on fluid
        (-beta_drag * (velocity - average_particles_velocity)
         // other two-way coupling forces applied on the fluid
         + explicit_particle_volumetric_acceleration_on_fluid);

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto &phi_u_j = scratch_data.phi_u[q][j];
          strong_jacobian[q][j] +=
            // Drag Force
            beta_drag * phi_u_j;
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_u_i = scratch_data.phi_u[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto &phi_u_j = scratch_data.phi_u[q][j];

              local_matrix(i, j) += // Drag Force
                beta_drag * phi_u_j * phi_u_i * JxW;
            }
        }
    }
}

template <int dim>
void
VANSAssemblerFPI<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &local_rhs       = copy_data.local_rhs;
  auto  beta_drag       = scratch_data.beta_drag;
  auto &explicit_particle_volumetric_acceleration_on_fluid =
    scratch_data.explicit_particle_volumetric_acceleration_on_fluid;
  const Tensor<1, dim> average_particles_velocity =
    scratch_data.average_particle_velocity;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Velocity
      const Tensor<1, dim> velocity = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculate the strong residual for GLS stabilization
      // No need to have separate codes for models A and B here, as the
      // correct forces are included in
      // explicit_particle_volumetric_acceleration_on_fluid in the
      // particle-fluid interaction assemblers, depending on the VANS model.

      // Subtraction of forces applied on fluid from the residual for GLS
      // stabilization
      strong_residual[q] -=
        // drag applied on fluid
        (-beta_drag * (velocity - average_particles_velocity)
         // other two-way coupling forces applied on the fluid
         + explicit_particle_volumetric_acceleration_on_fluid);

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i = scratch_data.phi_u[q][i];

          local_rhs(i) +=
            // + drag applied on fluid
            (-beta_drag * (velocity - average_particles_velocity) +
             // + other two-way coupling forces applied on the fluid
             explicit_particle_volumetric_acceleration_on_fluid) *
            phi_u_i * JxW;
        }
    }
}
template class VANSAssemblerFPI<2>;
template class VANSAssemblerFPI<3>;

template <int dim>
void
VANSAssemblerFPIProjection<dim>::assemble_matrix(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  // The CFD-DEM solver only supports constant density for now
  const double density = scratch_data.density_scale;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &strong_jacobian = copy_data.strong_jacobian;
  auto &local_matrix    = copy_data.local_matrix;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      const Tensor<1, dim> &two_way_coupling_force =
        -scratch_data.particle_two_way_coupling_force_values[q] / density;
      const Tensor<1, dim> &fluid_drag =
        -scratch_data.particle_drag_values[q] / density;
      const double &beta_drag =
        scratch_data.particle_momentum_transfer_coefficient_values[q];
      const Tensor<1, dim> particles_velocity =
        scratch_data.particle_velocity_values[q];
      const Tensor<1, dim> velocity = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Subtraction of forces applied on fluid from the residual for GLS
      // stabilization
      // If the coupling is explicit, beta_drag is zero. Otherwise, fluid_drag
      // is zero.
      strong_residual[q] -=
        (fluid_drag - beta_drag * (velocity - particles_velocity) +
         two_way_coupling_force);

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto &phi_u_j = scratch_data.phi_u[q][j];
          strong_jacobian[q][j] +=
            // Drag Force
            beta_drag * phi_u_j;
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_u_i = scratch_data.phi_u[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto &phi_u_j = scratch_data.phi_u[q][j];

              local_matrix(i, j) += // Drag Force
                beta_drag * phi_u_j * phi_u_i * JxW;
            }
        }
    }
}

template <int dim>
void
VANSAssemblerFPIProjection<dim>::assemble_rhs(
  const NavierStokesScratchData<dim>   &scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  // The CFD-DEM solver only supports constant density for now
  const double density = scratch_data.density_scale;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &local_rhs       = copy_data.local_rhs;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      const Tensor<1, dim> &two_way_coupling_force =
        -scratch_data.particle_two_way_coupling_force_values[q] / density;
      const Tensor<1, dim> &fluid_drag =
        -scratch_data.particle_drag_values[q] / density;
      const double &beta_drag =
        scratch_data.particle_momentum_transfer_coefficient_values[q];
      const Tensor<1, dim> particles_velocity =
        scratch_data.particle_velocity_values[q];
      const Tensor<1, dim> velocity = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculate the strong residual for GLS stabilization
      // Drag Force and other two-way coupling forces
      // If the coupling is explicit, beta_drag is zero. Otherwise, fluid_drag
      // is zero.
      strong_residual[q] -=
        (fluid_drag - beta_drag * (velocity - particles_velocity) +
         two_way_coupling_force);

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i = scratch_data.phi_u[q][i];
          // Drag Force
          // The distinction between Model A and B of the VANS equations is
          // made in the shear and pressure forces assemblers.
          // If the coupling is explicit, beta_drag is zero. Otherwise,
          // fluid_drag is zero.
          local_rhs(i) +=
            (fluid_drag - beta_drag * (velocity - particles_velocity) +
             two_way_coupling_force) *
            phi_u_i * JxW;
        }
    }
}
template class VANSAssemblerFPIProjection<2>;
template class VANSAssemblerFPIProjection<3>;
