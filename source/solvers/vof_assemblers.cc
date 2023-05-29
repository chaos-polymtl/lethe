#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/copy_data.h>
#include <solvers/vof_assemblers.h>



template <int dim>
void
VOFAssemblerCore<dim>::assemble_matrix(VOFScratchData<dim>       &scratch_data,
                                       StabilizedMethodsCopyData &copy_data)
{
  // Scheme and physical properties
  const auto method = this->simulation_control->get_assembly_method();

  // Loop and quadrature informations
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Add a small diffusivity, used in the context of the wetting mechanism
  const double diffusivity = this->vof_parameters.diffusivity;

  // Copy data elements
  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix        = copy_data.local_matrix;

  // assembling local matrix and right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const Tensor<1, dim> velocity = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      const Tensor<1, dim> phase_gradient = scratch_data.phase_gradients[q];
      const double         phase_gradient_norm = phase_gradient.norm();

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter and the compression term for the phase
      // indicator
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Implementation of a DCDD shock capturing scheme.
      // For more information see
      // Tezduyar, T. E., & Park, Y. J. (1986). Discontinuity-capturing
      // finite element formulations for nonlinear
      // convection-diffusion-reaction equations. Computer methods in
      // applied mechanics and engineering, 59(3), 307-325.

      // Gather the order of the VOF interpolation
      const double order = this->fem_parameters.VOF_order;

      // Calculate the artificial viscosity of the shock capture
      const double vdcdd = (0.5 * h) * (velocity.norm());

      const double tolerance = 1e-12;

      Tensor<1, dim> gradient_unit_vector =
        phase_gradient / (phase_gradient_norm + tolerance);

      // We neglect to remove the diffusion aligned with the velocity
      // as is done in the original article. This term generates poorer
      // results than just using the vdcdd approach.
      // Tensor<1, dim> s = velocity / (velocity.norm() + 1e-12);
      // const Tensor<2, dim> k_corr      = (r * s) * outer_product(s,
      // s);
      // const Tensor<2, dim> k_corr      = (r * s) * outer_product(s, s);
      const Tensor<2, dim> gradient_unit_tensor =
        outer_product(gradient_unit_vector, gradient_unit_vector);
      const Tensor<2, dim> dcdd_factor = gradient_unit_tensor; // - k_corr;

      // We neglect the gradient of the shock capturing viscosity in the
      // jacobian matrix since it tends to destabilize the matrix Gradient of
      // the shock capturing viscosity for the assembly of the jacobian matrix
      // const double d_vdcdd = order * (0.5 * h * h) *
      //                       (velocity.norm() * velocity.norm()) *
      //                       pow(phase_gradient_norm * h, order - 1);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step. Hypothesis : advection dominated problem
      // (Pe>3) [Bochev et al., Stability of the SUPG finite element
      // method for transient advection-diffusion problems, CMAME 2004]
      const double tau =
        is_steady(method) ?
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * diffusivity / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * diffusivity / (h * h), 2));

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto grad_phi_phase_j      = scratch_data.grad_phi[q][j];
          const auto laplacian_phi_phase_j = scratch_data.laplacian_phi[q][j];

          // Strong Jacobian associated with the GLS
          // stabilization
          strong_jacobian_vec[q][j] +=
            velocity * grad_phi_phase_j - diffusivity * laplacian_phi_phase_j;

          // DCDD shock capturing
          if (DCDD)
            strong_jacobian_vec[q][j] += -vdcdd * laplacian_phi_phase_j;
        }



      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_phase_i      = scratch_data.phi[q][i];
          const auto grad_phi_phase_i = scratch_data.grad_phi[q][i];


          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto grad_phi_phase_j = scratch_data.grad_phi[q][j];

              // Weak form for advection-diffusion:
              // u * grad(phase) - diffusivity * laplacian(phase) = 0
              local_matrix(i, j) +=
                (phi_phase_i * velocity * grad_phi_phase_j +
                 diffusivity * grad_phi_phase_i * grad_phi_phase_j) *
                JxW;

              // Addition to the cell matrix for GLS stabilization
              local_matrix(i, j) += tau * strong_jacobian_vec[q][j] *
                                    (grad_phi_phase_i * velocity) * JxW;

              // DCDD shock capturing
              if (DCDD)
                {
                  local_matrix(i, j) +=
                    (vdcdd * scalar_product(grad_phi_phase_j,
                                            grad_phi_phase_i) //+
                     // d_vdcdd * grad_phi_phase_j.norm() *
                     //   scalar_product(phase_gradient,
                     //                  dcdd_factor * grad_phi_phase_i)
                     ) *
                    JxW;
                }
            }
        }
    } // end loop on quadrature points
}



template <int dim>
void
VOFAssemblerCore<dim>::assemble_rhs(VOFScratchData<dim>       &scratch_data,
                                    StabilizedMethodsCopyData &copy_data)
{
  // Scheme and physical properties
  const auto method = this->simulation_control->get_assembly_method();

  // Add a small diffusivity, used in the context of the wetting mechanism
  const double diffusivity = this->vof_parameters.diffusivity;

  // Loop and quadrature informations
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs           = copy_data.local_rhs;

  // assembling local matrix and right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const Tensor<1, dim> phase_gradient   = scratch_data.phase_gradients[q];
      const double         phase_laplacians = scratch_data.phase_laplacians[q];
      const Tensor<1, dim> velocity         = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW                 = JxW_vec[q];
      const double phase_gradient_norm = phase_gradient.norm();

      // Implementation of a DCDD shock capturing scheme.
      // For more information see
      // Tezduyar, T. E., & Park, Y. J. (1986). Discontinuity-capturing
      // finite element formulations for nonlinear
      // convection-diffusion-reaction equations. Computer methods in
      // applied mechanics and engineering, 59(3), 307-325.

      // Gather the order of the VOF interpolation
      const double order = this->fem_parameters.VOF_order;

      // Calculate the artificial viscosity of the shock capture
      const double vdcdd = (0.5 * h) * (velocity.norm());

      const double tolerance = 1e-12;

      Tensor<1, dim> gradient_unit_vector =
        phase_gradient / (phase_gradient_norm + tolerance);

      // We neglect to remove the diffusion aligned with the velocity
      // as is done in the original article. This term generates poorer
      // results than just using the vdcdd approach.
      // Tensor<1, dim> s = velocity / (velocity.norm() + 1e-12);
      // const Tensor<2, dim> k_corr      = (r * s) * outer_product(s,
      // s);
      // const Tensor<2, dim> k_corr      = (r * s) * outer_product(s, s);
      const Tensor<2, dim> gradient_unit_tensor =
        outer_product(gradient_unit_vector, gradient_unit_vector);
      const Tensor<2, dim> dcdd_factor = gradient_unit_tensor; // - k_corr;

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), tolerance);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step. Hypothesis : advection dominated problem
      // (Pe>3) [Bochev et al., Stability of the SUPG finite element
      // method for transient advection-diffusion problems, CMAME 2004]
      const double tau =
        is_steady(method) ?
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * diffusivity / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * diffusivity / (h * h), 2));

      // Calculate the strong residual for GLS stabilization
      strong_residual_vec[q] +=
        velocity * phase_gradient - diffusivity * phase_laplacians;

      // DCDD shock capturing
      if (DCDD)
        strong_residual_vec[q] += -vdcdd * phase_laplacians;

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_phase_i      = scratch_data.phi[q][i];
          const auto grad_phi_phase_i = scratch_data.grad_phi[q][i];


          // rhs for: u * grad(phase) - diffusivity * laplacian(phase) = 0
          local_rhs(i) -= (phi_phase_i * velocity * phase_gradient +
                           diffusivity * grad_phi_phase_i * phase_gradient) *
                          JxW;

          local_rhs(i) -=
            tau * (strong_residual_vec[q] * (grad_phi_phase_i * velocity)) *
            JxW;

          // DCDD shock capturing
          if (DCDD)
            {
              local_rhs(i) +=
                -vdcdd * scalar_product(phase_gradient, grad_phi_phase_i) * JxW;
            }
        }
    }
}

template class VOFAssemblerCore<2>;
template class VOFAssemblerCore<3>;

template <int dim>
void
VOFAssemblerBDF<dim>::assemble_matrix(VOFScratchData<dim>       &scratch_data,
                                      StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature informations
  const auto        &JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_jacobian = copy_data.strong_jacobian;
  auto &strong_residual = copy_data.strong_residual;
  auto &local_matrix    = copy_data.local_matrix;

  // Time stepping information
  const auto          method = this->simulation_control->get_assembly_method();
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Vector for the BDF coefficients
  Vector<double>      bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<double> phase_value(1 + number_of_previous_solutions(method));

  // Extrapolate velocity to t+dt using the BDF scheme
  std::vector<double> time_vector = simulation_control->get_simulation_times();
  bdf_extrapolate(time_vector,
                  scratch_data.previous_velocity_values,
                  number_of_previous_solutions(method),
                  scratch_data.velocity_values);

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      phase_value[0] = scratch_data.present_phase_values[q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        phase_value[p + 1] = scratch_data.previous_phase_values[p][q];


      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += (bdf_coefs[p] * phase_value[p]);
        }

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian[q][j] += bdf_coefs[0] * scratch_data.phi[q][j];
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_phase_i = scratch_data.phi[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto phi_phase_j = scratch_data.phi[q][j];

              local_matrix(i, j) +=
                phi_phase_j * phi_phase_i * bdf_coefs[0] * JxW[q];
            }
        }
    }
}

template <int dim>
void
VOFAssemblerBDF<dim>::assemble_rhs(VOFScratchData<dim>       &scratch_data,
                                   StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature informations
  const auto        &JxW        = scratch_data.JxW;
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
  std::vector<double> phase_value(1 + number_of_previous_solutions(method));

  // Extrapolate velocity to t+dt using the BDF scheme
  std::vector<double> time_vector = simulation_control->get_simulation_times();
  bdf_extrapolate(time_vector,
                  scratch_data.previous_velocity_values,
                  number_of_previous_solutions(method),
                  scratch_data.velocity_values);

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      phase_value[0] = scratch_data.present_phase_values[q];


      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        phase_value[p + 1] = scratch_data.previous_phase_values[p][q];



      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += (bdf_coefs[p] * phase_value[p]);
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double phi_phase_i = scratch_data.phi[q][i];
          double       local_rhs_i = 0;
          for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
               ++p)
            {
              local_rhs_i -= bdf_coefs[p] * (phase_value[p] * phi_phase_i);
            }
          local_rhs(i) += local_rhs_i * JxW[q];
        }
    }
}

template class VOFAssemblerBDF<2>;
template class VOFAssemblerBDF<3>;
