#include <core/bdf.h>
#include <core/sdirk.h>
#include <core/simulation_control.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/navier_stokes_vof_assemblers.h>

template <int dim>
void
GLSNavierStokesVOFAssemblerCore<dim>::assemble_matrix(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature informations
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

  // Limit force application : not applied if cell density is below
  // density_ratio of the maximum density (e.g. when one of the fluids is air)
  const double density_ratio      = 2;
  double       phase_force_cutoff = 0;

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "GLSVansAssemblerDiFelice<dim>::calculate_particle_fluid_interactions"));

  if (scratch_data.density_0[0] < scratch_data.density_1[0] &&
      scratch_data.density_0[0] * density_ratio < scratch_data.density_1[0])
    {
      // gravity not will be applied for phase < phase_force_cutoff
      phase_force_cutoff = 1e-6;
    }
  if (scratch_data.density_1[0] < scratch_data.density_0[0] &&
      scratch_data.density_1[0] * density_ratio < scratch_data.density_0[0])
    {
      // gravity not will be applied for phase > phase_force_cutoff
      phase_force_cutoff = 1 - 1e-6;
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

      const Tensor<1, dim> pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Forcing term
      Tensor<1, dim> force = scratch_data.force[q];

      if (phase_force_cutoff < 0.5 && phase_values[q] < phase_force_cutoff)
        force = 0;
      // Gravity not applied on phase 1
      if (phase_force_cutoff > 0.5 && phase_values[q] > phase_force_cutoff)
        force = 0;

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the equivalent density at the quadrature point
      double density_eq = scratch_data.density[q];

      double viscosity_eq = scratch_data.viscosity[q];

      const double dynamic_viscosity_eq = density_eq * viscosity_eq;

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          1. / std::sqrt(std::pow(2. * density_eq * u_mag / h, 2) +
                         9 * std::pow(4 * dynamic_viscosity_eq / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) +
                         std::pow(2. * density_eq * u_mag / h, 2) +
                         9 * std::pow(4 * dynamic_viscosity_eq / (h * h), 2));

      // Calculate the strong residual for GLS stabilization
      auto strong_residual = density_eq * velocity_gradient * velocity +
                             pressure_gradient -
                             dynamic_viscosity_eq * velocity_laplacian -
                             density_eq * force + strong_residual_vec[q];

      std::vector<Tensor<1, dim>> grad_phi_u_j_x_velocity(n_dofs);
      std::vector<Tensor<1, dim>> velocity_gradient_x_phi_u_j(n_dofs);


      // We loop over the column first to prevent recalculation
      // of the strong jacobian in the inner loop
      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto &phi_u_j           = scratch_data.phi_u[q][j];
          const auto &grad_phi_u_j      = scratch_data.grad_phi_u[q][j];
          const auto &laplacian_phi_u_j = scratch_data.laplacian_phi_u[q][j];

          const auto &grad_phi_p_j = scratch_data.grad_phi_p[q][j];

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

              const auto &phi_p_j = scratch_data.phi_p[q][j];

              const auto &strong_jac = strong_jacobian_vec[q][j];

              double local_matrix_ij =
                dynamic_viscosity_eq *
                  scalar_product(grad_phi_u_j, grad_phi_u_i) +
                density_eq * velocity_gradient_x_phi_u_j[j] * phi_u_i +
                density_eq * grad_phi_u_j_x_velocity[j] * phi_u_i -
                div_phi_u_i * phi_p_j +
                // Continuity
                phi_p_i * div_phi_u_j;

              // PSPG GLS term
              local_matrix_ij += tau * (strong_jac * grad_phi_p_i);

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
  // Loop and quadrature informations
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

  // Limit force application : not applied if cell density is below
  // density_ratio of the maximum density (e.g. when one of the fluids is air)
  const double density_ratio      = 2;
  double       phase_force_cutoff = 0;

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "GLSVansAssemblerDiFelice<dim>::calculate_particle_fluid_interactions"));

  if (scratch_data.density_0[0] < scratch_data.density_1[0] &&
      scratch_data.density_0[0] * density_ratio < scratch_data.density_1[0])
    {
      // gravity not will be applied for phase < phase_force_cutoff
      phase_force_cutoff = 1e-6;
    }
  if (scratch_data.density_1[0] < scratch_data.density_0[0] &&
      scratch_data.density_1[0] * density_ratio < scratch_data.density_0[0])
    {
      // gravity not will be applied for phase > phase_force_cutoff
      phase_force_cutoff = 1 - 1e-6;
    }

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

      // Pressure
      const double         pressure = scratch_data.pressure_values[q];
      const Tensor<1, dim> pressure_gradient =
        scratch_data.pressure_gradients[q];

      // Forcing term
      Tensor<1, dim> force = scratch_data.force[q];

      if (phase_force_cutoff < 0.5 && phase_values[q] < phase_force_cutoff)
        force = 0;
      // Gravity not applied on phase 1
      if (phase_force_cutoff > 0.5 && phase_values[q] > phase_force_cutoff)
        force = 0;

      // Calculation of the magnitude of the
      // velocity for the stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the equivalent density at the quadrature point
      double density_eq = scratch_data.density[q];

      double       viscosity_eq         = scratch_data.viscosity[q];
      const double dynamic_viscosity_eq = density_eq * viscosity_eq;

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          1. / std::sqrt(std::pow(2. * density_eq * u_mag / h, 2) +
                         9 * std::pow(4 * dynamic_viscosity_eq / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) +
                         std::pow(2. * density_eq * u_mag / h, 2) +
                         9 * std::pow(4 * dynamic_viscosity_eq / (h * h), 2));


      // Calculate the strong residual for GLS stabilization
      auto strong_residual = density_eq * velocity_gradient * velocity +
                             pressure_gradient -
                             dynamic_viscosity_eq * velocity_laplacian -
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
          local_rhs(i) +=
            (
              // Momentum
              -dynamic_viscosity_eq *
                scalar_product(velocity_gradient, grad_phi_u_i) -
              density_eq * velocity_gradient * velocity * phi_u_i +
              pressure * div_phi_u_i + density_eq * force * phi_u_i -
              // Continuity
              velocity_divergence * phi_p_i) *
            JxW;

          // PSPG GLS term
          local_rhs(i) += -tau * (strong_residual * grad_phi_p_i) * JxW;

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
  // Loop and quadrature informations
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

      std::vector<double> densities(number_of_previous_solutions(method) + 1);

      densities[0] = scratch_data.density[q];

      for (unsigned int p = 1; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          densities[p] = scratch_data.previous_density[p - 1][q];
        }

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += densities[p] * bdf_coefs[p] * velocity[p];
        }

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian[q][j] +=
            densities[0] * bdf_coefs[0] * scratch_data.phi_u[q][j];
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const Tensor<1, dim> &phi_u_i = scratch_data.phi_u[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const Tensor<1, dim> &phi_u_j = scratch_data.phi_u[q][j];

              local_matrix(i, j) +=
                phi_u_j * phi_u_i * densities[0] * bdf_coefs[0] * JxW[q];
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
  // Loop and quadrature informations
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


      std::vector<double> densities(number_of_previous_solutions(method) + 1);

      densities[0] = scratch_data.density[q];

      for (unsigned int p = 1; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          densities[p] = scratch_data.previous_density[p - 1][q];
        }

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += (densities[p] * bdf_coefs[p] * velocity[p]);
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i     = scratch_data.phi_u[q][i];
          double     local_rhs_i = 0;
          for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
               ++p)
            {
              local_rhs_i -=
                densities[p] * bdf_coefs[p] * (velocity[p] * phi_u_i);
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
  // Surface tension coefficient
  const double surface_tension_coef = STF_properties.surface_tension_coef;

  // Loop and quadrature informations
  const auto &       JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &local_rhs       = copy_data.local_rhs;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather pfg and curvature values
      const double &        curvature_value = scratch_data.curvature_values[q];
      const Tensor<1, dim> &pfg_value       = scratch_data.pfg_values[q];
      const double          JxW_value       = JxW[q];
      const Tensor<1, dim>  temp_STF =
        -2.0 * surface_tension_coef * curvature_value * pfg_value;

      strong_residual[q] += temp_STF;

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i     = scratch_data.phi_u[q][i];
          double     local_rhs_i = 0;

          local_rhs_i -= temp_STF * phi_u_i;
          local_rhs(i) += local_rhs_i * JxW_value;
        }
    }
}

template class GLSNavierStokesVOFAssemblerSTF<2>;
template class GLSNavierStokesVOFAssemblerSTF<3>;
