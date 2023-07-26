#include <core/bdf.h>
#include <core/sdirk.h>
#include <core/simulation_control.h>

#include <solvers/isothermal_compressible_navier_stokes_assembler.h>
#include <solvers/stabilization.h>


template <int dim>
void
GLSIsothermalCompressibleNavierStokesAssemblerCore<dim>::assemble_matrix(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Physical properties
  const std::vector<double> &viscosity_vector = scratch_data.viscosity;
  const std::vector<double> &density_vector   = scratch_data.density;
  const double               density_psi      = scratch_data.density_psi;
  const double               density_ref      = scratch_data.density_ref;

  // Loop and quadrature information
  const auto &       JxW_vec    = scratch_data.JxW;
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
      const double viscosity         = viscosity_vector[q];
      const double density           = density_vector[q];
      const double dynamic_viscosity = density_ref * viscosity;

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
          calculate_navier_stokes_gls_tau_steady(u_mag, viscosity, h) :
          calculate_navier_stokes_gls_tau_transient(u_mag, viscosity, h, sdt);

      // Calculate the strong residual for GLS stabilization
      auto strong_residual = density * velocity_gradient * velocity +
                             pressure_gradient -
                             dynamic_viscosity * velocity_laplacian - force +
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
            (density * velocity_gradient * phi_u_j + grad_phi_u_j * velocity +
             grad_phi_p_j - dynamic_viscosity * laplacian_phi_u_j +
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
              const auto &grad_phi_p_j =
                pressure_scaling_factor * scratch_data.grad_phi_p[q][j];

              const auto &strong_jac = strong_jacobian_vec[q][j];

              // $$\nabla p$$ from the momentum equation
              double local_matrix_ij =
                component_j == dim ? -div_phi_u_i * phi_p_j : 0;

              // continuity equation
              if (component_i == dim)
                {
                  const auto &div_phi_u_j = scratch_data.div_phi_u[q][j];

                  local_matrix_ij +=
                    // $$\psi \mathbf{u} \nabla p$$
                    phi_p_i * density_psi * phi_u_j * pressure_gradient +
                    phi_p_i * density_psi * velocity * grad_phi_p_j +

                    // $$\rho \nabla \cdot \mathbf{u}$$
                    phi_p_i * density * div_phi_u_j;

                  // PSPG GLS term
                  local_matrix_ij += tau * (strong_jac * grad_phi_p_i);
                }

              // momentum equation
              if (component_i < dim && component_j < dim)
                {
                  local_matrix_ij +=
                    // $$\rho (\mathbf{u} \cdot \nabla) \mathbf{u}$$
                    phi_u_i * density * velocity_gradient_x_phi_u_j[j] +
                    phi_u_i * density * grad_phi_u_j_x_velocity[j];

                  if (component_i == component_j)
                    {
                      const auto &grad_phi_u_j = scratch_data.grad_phi_u[q][j];

                      local_matrix_ij +=
                        // weak form of $$\mu \nabla^{2} \mathbf{u}$$
                        dynamic_viscosity * (grad_phi_u_j[component_j] *
                                             grad_phi_u_i[component_i]) +
                        // volume source
                        phi_u_i * mass_source * phi_u_j;
                    }
                }

              // The jacobian matrix for the SUPG formulation
              // currently does not include the jacobian of the
              // stabilization parameter tau. Our experience has shown that
              // does not alter the number of newton iteration for
              // convergence, but greatly simplifies assembly.
              local_matrix_ij +=
                tau * (strong_jac * (grad_phi_u_i_x_velocity) +
                       strong_residual_x_grad_phi_u_i * phi_u_j);

              local_matrix_ij *= JxW;
              local_matrix(i, j) += local_matrix_ij;
            }
        }
    }
}

template <int dim>
void
GLSIsothermalCompressibleNavierStokesAssemblerCore<dim>::assemble_rhs(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Physical properties
  const std::vector<double> &viscosity_vector = scratch_data.viscosity;
  const std::vector<double> &density_vector   = scratch_data.density;
  const double               density_psi      = scratch_data.density_psi;
  const double               density_ref      = scratch_data.density_ref;

  // Loop and quadrature information
  const auto &       JxW_vec    = scratch_data.JxW;
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
      const double viscosity         = viscosity_vector[q];
      const double density           = density_vector[q];
      const double dynamic_viscosity = density_ref * viscosity;

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
          calculate_navier_stokes_gls_tau_steady(u_mag, viscosity, h) :
          calculate_navier_stokes_gls_tau_transient(u_mag, viscosity, h, sdt);

      // Calculate the strong residual for GLS stabilization
      auto strong_residual = density * velocity_gradient * velocity +
                             pressure_gradient -
                             dynamic_viscosity * velocity_laplacian - force +
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
              // Continuity
              -phi_p_i * density_psi * velocity * pressure_gradient -
              phi_p_i * density * velocity_divergence - phi_p_i * mass_source -
              // Momentum
              phi_u_i * density * velocity_gradient * velocity +
              div_phi_u_i * pressure -
              dynamic_viscosity *
                scalar_product(velocity_gradient, grad_phi_u_i) +
              phi_u_i * force - mass_source * velocity * phi_u_i) *
            JxW;

          // PSPG GLS term
          local_rhs_i +=
            -tau / density * (strong_residual * grad_phi_p_i) * JxW;

          // SUPG GLS term
          local_rhs_i +=
            -tau * (strong_residual * (grad_phi_u_i * velocity)) * JxW;

          local_rhs(i) += local_rhs_i;
        }
    }
}

template class GLSIsothermalCompressibleNavierStokesAssemblerCore<2>;
template class GLSIsothermalCompressibleNavierStokesAssemblerCore<3>;

template <int dim>
void
GLSIsothermalCompressibleNavierStokesAssemblerBDF<dim>::assemble_matrix(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Physical properties
  const std::vector<double> &density_vector = scratch_data.density;
  const double               density_psi    = scratch_data.density_psi;

  // Loop and quadrature information
  const auto &       JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &strong_jacobian = copy_data.strong_jacobian;
  auto &local_matrix    = copy_data.local_matrix;

  // Time stepping information
  const auto          method = this->simulation_control->get_assembly_method();
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Vector for the BDF coefficients
  Vector<double>      bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<double> pressure(1 + number_of_previous_solutions(method));
  std::vector<Tensor<1, dim>> velocity(1 +
                                       number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Physical properties
      const double density = density_vector[q];

      velocity[0] = scratch_data.velocity_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        velocity[p + 1] = scratch_data.previous_velocity_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        strong_residual[q] += bdf_coefs[p] * density * velocity[p];

      for (unsigned int j = 0; j < n_dofs; ++j)
        strong_jacobian[q][j] +=
          bdf_coefs[0] * density * scratch_data.phi_u[q][j];

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double          phi_p_i = scratch_data.phi_p[q][i];
          const Tensor<1, dim> &phi_u_i = scratch_data.phi_u[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const double &        phi_p_j = scratch_data.phi_p[q][j];
              const Tensor<1, dim> &phi_u_j = scratch_data.phi_u[q][j];

              local_matrix(i, j) +=
                (phi_p_i * density_psi * phi_p_j * bdf_coefs[0] +
                 phi_u_i * density * phi_u_j * bdf_coefs[0]) *
                JxW[q];
            }
        }
    }
}

template <int dim>
void
GLSIsothermalCompressibleNavierStokesAssemblerBDF<dim>::assemble_rhs(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Physical properties
  const std::vector<double> &density_vector = scratch_data.density;
  const double               density_psi    = scratch_data.density_psi;

  // Loop and quadrature information
  const auto &       JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &local_rhs       = copy_data.local_rhs;

  // Time stepping information
  const auto          method = this->simulation_control->get_assembly_method();
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Vector for the BDF coefficients
  Vector<double>      bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<double> pressure(1 + number_of_previous_solutions(method));
  std::vector<Tensor<1, dim>> velocity(1 +
                                       number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Physical properties
      const double density = density_vector[q];

      pressure[0] = scratch_data.pressure_values[q];
      velocity[0] = scratch_data.velocity_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        {
          pressure[p + 1] = scratch_data.previous_pressure_values[p][q];
          velocity[p + 1] = scratch_data.previous_velocity_values[p][q];
        }

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        strong_residual[q] += (bdf_coefs[p] * density * velocity[p]);


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double          phi_p_i     = scratch_data.phi_p[q][i];
          const Tensor<1, dim> &phi_u_i     = scratch_data.phi_u[q][i];
          double                local_rhs_i = 0;
          for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
               ++p)
            {
              local_rhs_i -=
                bdf_coefs[p] * (phi_p_i * density_psi * pressure[p] +
                                phi_u_i * density * velocity[p]);
            }
          local_rhs(i) += local_rhs_i * JxW[q];
        }
    }
}

template class GLSIsothermalCompressibleNavierStokesAssemblerBDF<2>;
template class GLSIsothermalCompressibleNavierStokesAssemblerBDF<3>;
