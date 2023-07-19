#include <core/bdf.h>
#include <core/sdirk.h>
#include <core/simulation_control.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/navier_stokes_vof_assemblers.h>
#include <solvers/stabilization.h>

template <int dim>
void
GLSNavierStokesVOFAssemblerCore<dim>::assemble_matrix(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
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

  std::vector<double> &phase_values = scratch_data.phase_values;

  Assert(scratch_data.properties_manager.density_is_constant(),
         RequiresConstantDensity(
           "GLSNavierStokesVOFAssemblerCore<dim>::assemble_matrix"));

  // Phase cutoff to limit continuity application on non-conservative fluid
  const double phase_cutoff = 1e-6;

  // Determine whether continuity condition is solved in this cell.
  // Removing the conservation condition on the lowest density fluid
  // can improve the wetting mechanism in the framework of incompressible
  // fluids. See documentation for more details.
  auto max_phase_cell =
    std::max_element(std::begin(phase_values), std::end(phase_values));
  bool solve_continuity(true);

  if (vof_parameters.conservation.conservative_fluid ==
      Parameters::FluidIndicator::fluid0)
    {
      if (*max_phase_cell > 1. - phase_cutoff)
        solve_continuity = false;
    }
  else if (vof_parameters.conservation.conservative_fluid ==
           Parameters::FluidIndicator::fluid1)
    {
      if (*max_phase_cell < phase_cutoff)
        solve_continuity = false;
    }

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const Tensor<1, dim> velocity = scratch_data.velocity_values[q];
      const Tensor<2, dim> velocity_gradient =
        scratch_data.velocity_gradients[q];
      const Tensor<1, dim> velocity_laplacian =
        scratch_data.velocity_laplacians[q];
      const Tensor<3, dim> &velocity_hessian =
        scratch_data.velocity_hessians[q];

      // From hessian, calculate grad (div (u)) term needed for VOF problems
      Tensor<1, dim> grad_div_velocity;
      for (int d = 0; d < dim; ++d)
        {
          for (int k = 0; k < dim; ++k)
            {
              // hessian[c][i][j] is the (i,j)th component of the matrix of
              // second derivatives of the cth vector component of the vector
              // field at quadrature point q of the current cell
              grad_div_velocity[d] += velocity_hessian[k][d][k];
            }
        }

      const Tensor<1, dim> pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Forcing term
      Tensor<1, dim> force = scratch_data.force[q];

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the equivalent properties at the quadrature point
      double       density_eq           = scratch_data.density[q];
      double       viscosity_eq         = scratch_data.viscosity[q];
      const double dynamic_viscosity_eq = density_eq * viscosity_eq;

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(u_mag, viscosity_eq, h, 1) :
          calculate_navier_stokes_gls_tau_transient(
            u_mag, viscosity_eq, h, sdt, 1);

      // Calculate the strong residual for GLS stabilization
      auto strong_residual = density_eq * velocity_gradient * velocity +
                             pressure_gradient -
                             dynamic_viscosity_eq * velocity_laplacian -
                             dynamic_viscosity_eq * grad_div_velocity -
                             density_eq * force + strong_residual_vec[q];

      std::vector<Tensor<1, dim>> grad_phi_u_j_x_velocity(n_dofs);
      std::vector<Tensor<1, dim>> velocity_gradient_x_phi_u_j(n_dofs);

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
            scratch_data.grad_phi_p[q][j] * pressure_scaling_factor;

          strong_jacobian_vec[q][j] +=
            (density_eq * velocity_gradient * phi_u_j +
             density_eq * grad_phi_u_j * velocity + grad_phi_p_j -
             dynamic_viscosity_eq * laplacian_phi_u_j);

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
              const auto &grad_shear_rate_j =
                grad_phi_u_j + transpose(grad_phi_u_j);

              const auto &phi_p_j =
                scratch_data.phi_p[q][j] * pressure_scaling_factor;

              const auto &strong_jac = strong_jacobian_vec[q][j];

              double local_matrix_ij =
                dynamic_viscosity_eq *
                  scalar_product(grad_shear_rate_j, grad_phi_u_i) +
                density_eq * velocity_gradient_x_phi_u_j[j] * phi_u_i +
                density_eq * grad_phi_u_j_x_velocity[j] * phi_u_i -
                div_phi_u_i * phi_p_j;

              if (solve_continuity)
                {
                  // Continuity
                  local_matrix_ij += phi_p_i * div_phi_u_j;

                  // PSPG GLS term
                  local_matrix_ij +=
                    tau / density_eq * (strong_jac * grad_phi_p_i);
                }
              else
                {
                  // assemble Jacobian corresponding to p = 0
                  local_matrix_ij += phi_p_i * phi_p_j;
                }

              // Jacobian is currently incomplete
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
GLSNavierStokesVOFAssemblerCore<dim>::assemble_rhs(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
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

  std::vector<double> &phase_values = scratch_data.phase_values;

  // Phase cutoff to limit continuity application on non-conservative fluid
  const double phase_cutoff = 1e-6;

  // Determine whether continuity condition is solved in this cell.
  // Removing the conservation condition on the lowest density fluid
  // can improve the wetting mechanism in the framework of incompressible
  // fluids. See documentation for more details.
  auto max_phase_cell =
    std::max_element(std::begin(phase_values), std::end(phase_values));
  bool solve_continuity(true);

  if (vof_parameters.conservation.conservative_fluid ==
      Parameters::FluidIndicator::fluid0)
    {
      if (*max_phase_cell > 1. - phase_cutoff)
        solve_continuity = false;
    }
  else if (vof_parameters.conservation.conservative_fluid ==
           Parameters::FluidIndicator::fluid1)
    {
      if (*max_phase_cell < phase_cutoff)
        solve_continuity = false;
    }

  Assert(scratch_data.properties_manager.density_is_constant(),
         RequiresConstantDensity(
           "GLSNavierStokesVOFAssemblerCore<dim>::assemble_matrix"));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Velocity
      const Tensor<1, dim> velocity    = scratch_data.velocity_values[q];
      const double velocity_divergence = scratch_data.velocity_divergences[q];
      const Tensor<2, dim> velocity_gradient =
        scratch_data.velocity_gradients[q];
      const Tensor<1, dim> velocity_laplacian =
        scratch_data.velocity_laplacians[q];


      const Tensor<3, dim> &velocity_hessian =
        scratch_data.velocity_hessians[q];
      // From hessian, calculate grad (div (u)) term needed for VOF problems
      Tensor<1, dim> grad_div_velocity;
      for (int d = 0; d < dim; ++d)
        {
          for (int k = 0; k < dim; ++k)
            {
              // hessian[c][i][j] is the (i,j)th component of the matrix of
              // second derivatives of the cth vector component of the vector
              // field at quadrature point q of the current cell
              grad_div_velocity[d] += velocity_hessian[k][d][k];
            }
        }


      // Calculate shear rate
      const Tensor<2, dim> shear_rate =
        velocity_gradient + transpose(velocity_gradient);

      // Pressure
      const double         pressure = scratch_data.pressure_values[q];
      const Tensor<1, dim> pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Forcing term
      Tensor<1, dim> force = scratch_data.force[q];

      // Calculation of the magnitude of the
      // velocity for the stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the equivalent properties at the quadrature point
      double       density_eq           = scratch_data.density[q];
      double       viscosity_eq         = scratch_data.viscosity[q];
      const double dynamic_viscosity_eq = density_eq * viscosity_eq;

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(u_mag, viscosity_eq, h, 1.) :
          calculate_navier_stokes_gls_tau_transient(
            u_mag, viscosity_eq, h, sdt, 1.);


      // Calculate the strong residual for GLS stabilization
      auto strong_residual = density_eq * velocity_gradient * velocity +
                             pressure_gradient -
                             dynamic_viscosity_eq * velocity_laplacian -
                             dynamic_viscosity_eq * grad_div_velocity -
                             density_eq * force + strong_residual_vec[q];

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i      = scratch_data.phi_u[q][i];
          const auto grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto phi_p_i      = scratch_data.phi_p[q][i];
          const auto grad_phi_p_i = scratch_data.grad_phi_p[q][i];
          const auto div_phi_u_i  = scratch_data.div_phi_u[q][i];


          // Navier-Stokes Residual
          // Momentum
          local_rhs(i) +=
            (-dynamic_viscosity_eq * scalar_product(shear_rate, grad_phi_u_i) -
             density_eq * velocity_gradient * velocity * phi_u_i +
             pressure * div_phi_u_i + density_eq * force * phi_u_i) *
            JxW;

          if (solve_continuity)
            {
              // Continuity
              local_rhs(i) += -(velocity_divergence * phi_p_i) * JxW;

              // PSPG GLS term
              local_rhs(i) +=
                -tau / density_eq * (strong_residual * grad_phi_p_i) * JxW;
            }
          else
            {
              // assemble RHS for p = 0
              local_rhs(i) += -phi_p_i * pressure * JxW;
            }

          // SUPG GLS term
          if (SUPG)
            {
              local_rhs(i) +=
                -tau * (strong_residual * (grad_phi_u_i * velocity)) * JxW;
            }
        }
    }
}



template class GLSNavierStokesVOFAssemblerCore<2>;
template class GLSNavierStokesVOFAssemblerCore<3>;

template <int dim>
void
GLSNavierStokesVOFAssemblerBDF<dim>::assemble_matrix(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
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
  Vector<double> bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<Tensor<1, dim>> velocity(1 +
                                       number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      velocity[0] = scratch_data.velocity_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        velocity[p + 1] = scratch_data.previous_velocity_values[p][q];

      const double density = scratch_data.density[q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += density * bdf_coefs[p] * velocity[p];
        }

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian[q][j] +=
            density * bdf_coefs[0] * scratch_data.phi_u[q][j];
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const Tensor<1, dim> &phi_u_i = scratch_data.phi_u[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const Tensor<1, dim> &phi_u_j = scratch_data.phi_u[q][j];

              local_matrix(i, j) +=
                phi_u_j * phi_u_i * density * bdf_coefs[0] * JxW[q];
            }
        }
    }
}

template <int dim>
void
GLSNavierStokesVOFAssemblerBDF<dim>::assemble_rhs(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
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
  Vector<double> bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<Tensor<1, dim>> velocity(1 +
                                       number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      velocity[0] = scratch_data.velocity_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        velocity[p + 1] = scratch_data.previous_velocity_values[p][q];

      const double density = scratch_data.density[q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += (density * bdf_coefs[p] * velocity[p]);
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i     = scratch_data.phi_u[q][i];
          double     local_rhs_i = 0;
          for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
               ++p)
            {
              local_rhs_i -= density * bdf_coefs[p] * (velocity[p] * phi_u_i);
            }
          local_rhs(i) += local_rhs_i * JxW[q];
        }
    }
}

template class GLSNavierStokesVOFAssemblerBDF<2>;
template class GLSNavierStokesVOFAssemblerBDF<3>;


template <int dim>
void
GLSNavierStokesVOFAssemblerSTF<dim>::assemble_matrix(
  NavierStokesScratchData<dim> & /*scratch_data*/,
  StabilizedMethodsTensorCopyData<dim> & /*copy_data*/)
{}

template <int dim>
void
GLSNavierStokesVOFAssemblerSTF<dim>::assemble_rhs(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Densities of phases
  Assert(scratch_data.properties_manager.density_is_constant(),
         RequiresConstantDensity(
           "GLSNavierStokesVOFAssemblerCore<dim>::assemble_matrix"));

  // Loop and quadrature information
  const auto &       JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &local_rhs       = copy_data.local_rhs;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Surface tension coefficient
      const double surface_tension_coef = scratch_data.surface_tension[q];

      // Gather pfg and curvature values
      const double &        curvature_value = scratch_data.curvature_values[q];
      const Tensor<1, dim> &phase_gradient_value =
        scratch_data.filtered_phase_gradient_values[q];
      const double JxW_value = JxW[q];

      const Tensor<1, dim> surface_tension_force =
        -surface_tension_coef * curvature_value * phase_gradient_value;

      strong_residual[q] += surface_tension_force;

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i     = scratch_data.phi_u[q][i];
          double     local_rhs_i = 0;

          local_rhs_i -= surface_tension_force * phi_u_i;
          local_rhs(i) += local_rhs_i * JxW_value;
        }
    }
}

template class GLSNavierStokesVOFAssemblerSTF<2>;
template class GLSNavierStokesVOFAssemblerSTF<3>;


template <int dim>
void
GLSNavierStokesVOFAssemblerMarangoni<dim>::assemble_matrix(
  NavierStokesScratchData<dim> & /*scratch_data*/,
  StabilizedMethodsTensorCopyData<dim> & /*copy_data*/)
{}

template <int dim>
void
GLSNavierStokesVOFAssemblerMarangoni<dim>::assemble_rhs(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Surface tension gradient
  const double surface_tension_gradient =
    STF_properties.surface_tension_gradient;

  // Densities of phases
  Assert(scratch_data.properties_manager.density_is_constant(),
         RequiresConstantDensity(
           "GLSNavierStokesVOFAssemblerCore<dim>::assemble_matrix"));

  // Loop and quadrature information
  const auto &       JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &local_rhs       = copy_data.local_rhs;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Surface tension coefficient
      const double surface_tension_coef = scratch_data.surface_tension[q];

      const double &curvature_value = scratch_data.curvature_values[q];

      // Gather phase fraction gradient
      const Tensor<1, dim> &phase_gradient_value =
        scratch_data.filtered_phase_gradient_values[q];

      const double phase_gradient_norm = phase_gradient_value.norm();

      const Tensor<1, dim> normalized_phase_fraction_gradient =
        phase_gradient_value / (phase_gradient_norm + DBL_MIN);

      // Gather temperature
      const double temperature = scratch_data.temperature_values[q];

      // Gather temperature gradient
      const Tensor<1, dim> temperature_gradient =
        scratch_data.temperature_gradients[q];

      const double JxW_value = JxW[q];


      const Tensor<1, dim> surface_tension_force =
        -(surface_tension_coef + surface_tension_gradient * temperature) *
        curvature_value * phase_gradient_value;

      const Tensor<1, dim> marangoni_effect =
        -surface_tension_gradient *
        (temperature_gradient -
         normalized_phase_fraction_gradient *
           (normalized_phase_fraction_gradient * temperature_gradient)) *
        phase_gradient_norm;

      strong_residual[q] += marangoni_effect + surface_tension_force;

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i     = scratch_data.phi_u[q][i];
          double     local_rhs_i = 0;

          local_rhs_i -= (marangoni_effect + surface_tension_force) * phi_u_i;
          local_rhs(i) += local_rhs_i * JxW_value;
        }
    }
}

template class GLSNavierStokesVOFAssemblerMarangoni<2>;
template class GLSNavierStokesVOFAssemblerMarangoni<3>;


template <int dim>
void
GLSNavierStokesVOFAssemblerNonNewtonianCore<dim>::assemble_matrix(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
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

  // Phase values and limiters
  std::vector<double> &phase_values = scratch_data.phase_values;
  // std::vector<double> &phase_values_m1 =
  // scratch_data.previous_phase_values[0];
  // std::vector<Tensor<1, dim>> &phase_gradient_values =
  // scratch_data.phase_gradient_values;

  // Phase cutoff to limit continuity application on non-conservative fluid
  const double phase_cutoff = 1e-6;

  // Determine whether continuity condition is solved in this cell.
  // Removing the conservation condition on the lowest density fluid
  // can improve the wetting mechanism in the framework of incompressible
  // fluids. See documentation for more details.
  auto max_phase_cell =
    std::max_element(std::begin(phase_values), std::end(phase_values));
  bool solve_continuity(true);

  if (vof_parameters.conservation.conservative_fluid ==
      Parameters::FluidIndicator::fluid0)
    {
      if (*max_phase_cell > 1. - phase_cutoff)
        solve_continuity = false;
    }
  else if (vof_parameters.conservation.conservative_fluid ==
           Parameters::FluidIndicator::fluid1)
    {
      if (*max_phase_cell < phase_cutoff)
        solve_continuity = false;
    }

  Assert(scratch_data.properties_manager.density_is_constant(),
         RequiresConstantDensity(
           "GLSNavierStokesVOFAssemblerCore<dim>::assemble_matrix"));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
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
        this->get_viscosity_gradient(velocity_gradient,
                                     velocity_hessian,
                                     shear_rate_magnitude,
                                     scratch_data.grad_viscosity_shear_rate[q]);

      // Forcing term
      Tensor<1, dim> force = scratch_data.force[q];

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the equivalent properties at the quadrature point
      double               density_eq           = scratch_data.density[q];
      double               viscosity_eq         = scratch_data.viscosity[q];
      const double         dynamic_viscosity_eq = density_eq * viscosity_eq;
      const Tensor<1, dim> dynamic_viscosity_gradient =
        density_eq * viscosity_gradient;

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(u_mag,
                                                 viscosity_eq,
                                                 h,
                                                 density_eq) :
          calculate_navier_stokes_gls_tau_transient(
            u_mag, viscosity_eq, h, sdt, density_eq);

      // Calculate the strong residual for GLS stabilization
      auto strong_residual = density_eq * velocity_gradient * velocity +
                             pressure_gradient -
                             shear_rate * dynamic_viscosity_gradient -
                             dynamic_viscosity_eq * velocity_laplacian -
                             density_eq * force + strong_residual_vec[q];

      std::vector<Tensor<1, dim>> grad_phi_u_j_x_velocity(n_dofs);
      std::vector<Tensor<1, dim>> velocity_gradient_x_phi_u_j(n_dofs);

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
            scratch_data.grad_phi_p[q][j] * pressure_scaling_factor;

          const auto &grad_phi_u_j_non_newtonian =
            grad_phi_u_j + transpose(grad_phi_u_j);

          strong_jacobian_vec[q][j] +=
            density_eq * velocity_gradient * phi_u_j +
            density_eq * grad_phi_u_j * velocity + grad_phi_p_j -
            dynamic_viscosity_eq * laplacian_phi_u_j -
            grad_phi_u_j_non_newtonian * dynamic_viscosity_gradient;

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
                scratch_data.phi_p[q][j] * pressure_scaling_factor;

              const auto &strong_jac = strong_jacobian_vec[q][j];

              double local_matrix_ij =
                dynamic_viscosity_eq *
                  scalar_product(grad_phi_u_j_non_newtonian, grad_phi_u_i) +
                0.5 * scratch_data.grad_viscosity_shear_rate[q] /
                  shear_rate_magnitude *
                  scalar_product(grad_phi_u_j_non_newtonian, shear_rate) *
                  scalar_product(shear_rate, grad_phi_u_i) +
                density_eq * velocity_gradient_x_phi_u_j[j] * 0.5 * phi_u_i +
                density_eq * grad_phi_u_j_x_velocity[j] * phi_u_i -
                div_phi_u_i * phi_p_j;

              if (solve_continuity)
                {
                  // Continuity
                  local_matrix_ij += phi_p_i * div_phi_u_j;

                  // PSPG GLS term
                  local_matrix_ij += tau * (strong_jac * grad_phi_p_i);
                }
              else
                {
                  // assemble Jacobian corresponding to p = 0
                  local_matrix_ij += phi_p_i * phi_p_j;
                }

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
GLSNavierStokesVOFAssemblerNonNewtonianCore<dim>::assemble_rhs(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
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

  // Phase values and limiters
  std::vector<double> &phase_values = scratch_data.phase_values;
  // std::vector<double> &phase_values_m1 =
  // scratch_data.previous_phase_values[0];
  // std::vector<Tensor<1, dim>> &phase_gradient_values =
  // scratch_data.phase_gradient_values;

  // Phase cutoff to limit continuity application on non-conservative fluid
  const double phase_cutoff = 1e-6;

  // Determine whether continuity condition is solved in this cell.
  // Removing the conservation condition on the lowest density fluid
  // can improve the wetting mechanism in the framework of incompressible
  // fluids. See documentation for more details.
  auto max_phase_cell =
    std::max_element(std::begin(phase_values), std::end(phase_values));
  bool solve_continuity(true);

  if (vof_parameters.conservation.conservative_fluid ==
      Parameters::FluidIndicator::fluid0)
    {
      if (*max_phase_cell > 1. - phase_cutoff)
        solve_continuity = false;
    }
  else if (vof_parameters.conservation.conservative_fluid ==
           Parameters::FluidIndicator::fluid1)
    {
      if (*max_phase_cell < phase_cutoff)
        solve_continuity = false;
    }

  Assert(scratch_data.properties_manager.density_is_constant(),
         RequiresConstantDensity(
           "GLSNavierStokesVOFAssemblerCore<dim>::assemble_matrix"));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
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

      // Calculate viscosity gradient
      const Tensor<1, dim> viscosity_gradient =
        this->get_viscosity_gradient(velocity_gradient,
                                     velocity_hessian,
                                     shear_rate_magnitude,
                                     scratch_data.grad_viscosity_shear_rate[q]);

      // Pressure
      const double         pressure = scratch_data.pressure_values[q];
      const Tensor<1, dim> pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Forcing term
      Tensor<1, dim> force = scratch_data.force[q];

      // Calculation of the magnitude of the
      // velocity for the stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the equivalent properties at the quadrature point
      double               density_eq           = scratch_data.density[q];
      double               viscosity_eq         = scratch_data.viscosity[q];
      const double         dynamic_viscosity_eq = density_eq * viscosity_eq;
      const Tensor<1, dim> dynamic_viscosity_gradient =
        density_eq * viscosity_gradient;

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          calculate_navier_stokes_gls_tau_steady(u_mag,
                                                 viscosity_eq,
                                                 h,
                                                 density_eq) :
          calculate_navier_stokes_gls_tau_transient(
            u_mag, viscosity_eq, h, sdt, density_eq);


      // Calculate the strong residual for GLS stabilization
      auto strong_residual = density_eq * velocity_gradient * velocity +
                             pressure_gradient -
                             shear_rate * dynamic_viscosity_gradient -
                             dynamic_viscosity_eq * velocity_laplacian -
                             density_eq * force + strong_residual_vec[q];

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_u_i      = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &phi_p_i      = scratch_data.phi_p[q][i];
          const auto &grad_phi_p_i = scratch_data.grad_phi_p[q][i];
          const auto &div_phi_u_i  = scratch_data.div_phi_u[q][i];

          // Navier-Stokes Residual
          local_rhs(i) +=
            (
              // Momentum
              -dynamic_viscosity_eq * scalar_product(shear_rate, grad_phi_u_i) -
              density_eq * velocity_gradient * velocity * phi_u_i +
              pressure * div_phi_u_i + density_eq * force * phi_u_i) *
            JxW;

          if (solve_continuity)
            {
              // Continuity
              local_rhs(i) += -(velocity_divergence * phi_p_i) * JxW;

              // PSPG GLS term
              local_rhs(i) += -tau * (strong_residual * grad_phi_p_i) * JxW;
            }
          else
            {
              // assemble RHS for p = 0
              local_rhs(i) += -phi_p_i * pressure * JxW;
            }

          // SUPG GLS term
          if (SUPG)
            {
              local_rhs(i) +=
                -tau * (strong_residual * (grad_phi_u_i * velocity)) * JxW;
            }
        }
    }
}

template class GLSNavierStokesVOFAssemblerNonNewtonianCore<2>;
template class GLSNavierStokesVOFAssemblerNonNewtonianCore<3>;
