#include <core/bdf.h>
#include <core/sdirk.h>
#include <core/simulation_control.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <fem-dem/vans_assemblers.h>

#include <deal.II/base/tensor.h>

template <int dim>
void
GLSVansAssemblerCoreModelB<dim>::assemble_matrix(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Viscosity at Gauss points
  const std::vector<double> &viscosity_vector = scratch_data.viscosity;

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

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const double         viscosity = viscosity_vector[q];
      const Tensor<1, dim> velocity  = scratch_data.velocity_values[q];
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
          if (simulation_control->get_current_time() ==
              simulation_control->get_time_step())
            u_mag = std::max(velocity.norm(), 1e-12);
          else
            u_mag = std::max(previous_velocity.norm(), 1e-12);
        }

      // Grad-div weight factor
      const double gamma = calculate_gamma(u_mag, viscosity, h, cfd_dem.cstar);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * viscosity / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * viscosity / (h * h), 2));

      // Calculate the strong residual for GLS stabilization
      auto strong_residual = velocity_gradient * velocity * void_fraction +
                             // Mass Source
                             mass_source * velocity
                             // Pressure
                             + pressure_gradient -
                             // Viscosity and Force
                             viscosity * velocity_laplacian -
                             force * void_fraction + strong_residual_vec[q];

      // We loop over the column first to prevent recalculation
      // of the strong jacobian in the inner loop
      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto &phi_u_j           = scratch_data.phi_u[q][j];
          const auto &grad_phi_u_j      = scratch_data.grad_phi_u[q][j];
          const auto &laplacian_phi_u_j = scratch_data.laplacian_phi_u[q][j];

          const auto &grad_phi_p_j = scratch_data.grad_phi_p[q][j];

          strong_jacobian_vec[q][j] +=
            (velocity_gradient * phi_u_j * void_fraction +
             grad_phi_u_j * velocity * void_fraction +
             // Mass Source
             mass_source * phi_u_j +
             // Pressure
             grad_phi_p_j
             // Viscosity
             - viscosity * laplacian_phi_u_j);
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_u_i      = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &div_phi_u_i  = scratch_data.div_phi_u[q][i];
          const auto &phi_p_i      = scratch_data.phi_p[q][i];
          const auto &grad_phi_p_i = scratch_data.grad_phi_p[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto &phi_u_j      = scratch_data.phi_u[q][j];
              const auto &grad_phi_u_j = scratch_data.grad_phi_u[q][j];
              const auto &div_phi_u_j  = scratch_data.div_phi_u[q][j];

              const auto &phi_p_j = scratch_data.phi_p[q][j];

              const auto &strong_jac = strong_jacobian_vec[q][j];

              double local_matrix_ij =
                viscosity * scalar_product(grad_phi_u_j, grad_phi_u_i) +
                // Convection
                ((phi_u_j * void_fraction * velocity_gradient * phi_u_i) +
                 (grad_phi_u_j * void_fraction * velocity * phi_u_i))
                // Mass source term
                + mass_source * phi_u_j * phi_u_i
                // Pressure
                - (div_phi_u_i * phi_p_j) +
                // Continuity
                phi_p_i * ((void_fraction * div_phi_u_j) +
                           (phi_u_j * void_fraction_gradients));

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
GLSVansAssemblerCoreModelB<dim>::assemble_rhs(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &viscosity_vector = scratch_data.viscosity;

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

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Physical property
      const double viscosity = viscosity_vector[q];

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
          if (simulation_control->get_current_time() ==
              simulation_control->get_time_step())
            u_mag = std::max(velocity.norm(), 1e-12);
          else
            u_mag = std::max(previous_velocity.norm(), 1e-12);
        }

      // Grad-div weight factor
      const double gamma = calculate_gamma(u_mag, viscosity, h, cfd_dem.cstar);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * viscosity / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * viscosity / (h * h), 2));

      // Calculate the strong residual for GLS stabilization
      auto strong_residual = velocity_gradient * velocity * void_fraction +
                             // Mass Source
                             mass_source * velocity
                             // Pressure
                             + pressure_gradient -
                             // Viscosity and Force
                             viscosity * velocity_laplacian -
                             force * void_fraction + strong_residual_vec[q];

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i      = scratch_data.phi_u[q][i];
          const auto grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto phi_p_i      = scratch_data.phi_p[q][i];
          const auto grad_phi_p_i = scratch_data.grad_phi_p[q][i];
          const auto div_phi_u_i  = scratch_data.div_phi_u[q][i];

          double local_rhs_i = 0;

          // Navier-Stokes Residual
          local_rhs_i +=
            (
              // Momentum
              -viscosity * scalar_product(velocity_gradient, grad_phi_u_i) -
              velocity_gradient * velocity * void_fraction * phi_u_i
              // Mass Source
              - mass_source * velocity * phi_u_i
              // Pressure and Force
              + pressure * div_phi_u_i +
              force * void_fraction * phi_u_i
              // Continuity
              - (velocity_divergence * void_fraction +
                 velocity * void_fraction_gradients - mass_source) *
                  phi_p_i) *
            JxW;

          // PSPG GLS term
          local_rhs_i += -tau * (strong_residual * grad_phi_p_i) * JxW;

          // SUPG GLS term
          if (SUPG)
            {
              local_rhs_i +=
                -tau * (strong_residual * (grad_phi_u_i * velocity)) * JxW;
            }

          // Grad-div stabilization
          if (cfd_dem.grad_div == true)
            {
              local_rhs_i -= gamma *
                             (void_fraction * velocity_divergence +
                              velocity * void_fraction_gradients) *
                             div_phi_u_i * JxW;
            }

          local_rhs(i) += local_rhs_i;
        }
    }
}


template class GLSVansAssemblerCoreModelB<2>;

template class GLSVansAssemblerCoreModelB<3>;

template <int dim>
void
GLSVansAssemblerCoreModelA<dim>::assemble_matrix(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &viscosity_vector = scratch_data.viscosity;

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

  // Grad-div weight factor

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Physical properties
      const double viscosity = viscosity_vector[q];


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
          if (simulation_control->get_current_time() ==
              simulation_control->get_time_step())
            u_mag = std::max(velocity.norm(), 1e-12);
          else
            u_mag = std::max(previous_velocity.norm(), 1e-12);
        }

      // Grad-div weight factor
      const double gamma = calculate_gamma(u_mag, viscosity, h, cfd_dem.cstar);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * viscosity / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * viscosity / (h * h), 2));

      // Calculate the strong residual for GLS stabilization
      auto strong_residual = velocity_gradient * velocity * void_fraction +
                             // Mass Source
                             mass_source * velocity
                             // Pressure
                             + void_fraction * pressure_gradient -
                             // Viscosity and Force
                             void_fraction * viscosity * velocity_laplacian -
                             force * void_fraction + strong_residual_vec[q];

      // We loop over the column first to prevent recalculation
      // of the strong jacobian in the inner loop
      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto &phi_u_j           = scratch_data.phi_u[q][j];
          const auto &grad_phi_u_j      = scratch_data.grad_phi_u[q][j];
          const auto &laplacian_phi_u_j = scratch_data.laplacian_phi_u[q][j];

          const auto &grad_phi_p_j = scratch_data.grad_phi_p[q][j];

          strong_jacobian_vec[q][j] +=
            (velocity_gradient * phi_u_j * void_fraction +
             grad_phi_u_j * velocity * void_fraction +
             // Mass Source
             mass_source * phi_u_j +
             // Pressure
             void_fraction * grad_phi_p_j
             // Viscosity
             - void_fraction * viscosity * laplacian_phi_u_j);
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_u_i      = scratch_data.phi_u[q][i];
          const auto &grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto &div_phi_u_i  = scratch_data.div_phi_u[q][i];
          const auto &phi_p_i      = scratch_data.phi_p[q][i];
          const auto &grad_phi_p_i = scratch_data.grad_phi_p[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto &phi_u_j      = scratch_data.phi_u[q][j];
              const auto &grad_phi_u_j = scratch_data.grad_phi_u[q][j];
              const auto &div_phi_u_j  = scratch_data.div_phi_u[q][j];

              const auto &phi_p_j = scratch_data.phi_p[q][j];

              const auto &strong_jac = strong_jacobian_vec[q][j];

              double local_matrix_ij =
                (void_fraction * viscosity *
                   scalar_product(grad_phi_u_j, grad_phi_u_i) +
                 viscosity * grad_phi_u_j * void_fraction_gradients * phi_u_i) +
                // Convection
                ((phi_u_j * void_fraction * velocity_gradient * phi_u_i) +
                 (grad_phi_u_j * void_fraction * velocity * phi_u_i))
                // Mass source term
                + mass_source * phi_u_j * phi_u_i
                // Pressure
                - (div_phi_u_i * phi_p_j * void_fraction +
                   phi_p_j * void_fraction_gradients * phi_u_i) +
                // Continuity
                phi_p_i * ((void_fraction * div_phi_u_j) +
                           (phi_u_j * void_fraction_gradients));

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
GLSVansAssemblerCoreModelA<dim>::assemble_rhs(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &viscosity_vector = scratch_data.viscosity;

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

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Physical properties
      const double viscosity = viscosity_vector[q];

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
          if (simulation_control->get_current_time() ==
              simulation_control->get_time_step())
            u_mag = std::max(velocity.norm(), 1e-12);
          else
            u_mag = std::max(previous_velocity.norm(), 1e-12);
        }

      // Grad-div weight factor
      const double gamma = calculate_gamma(u_mag, viscosity, h, cfd_dem.cstar);

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation
      // is steady or unsteady. In the unsteady case it includes the
      // value of the time-step
      const double tau =
        this->simulation_control->get_assembly_method() ==
            Parameters::SimulationControl::TimeSteppingMethod::steady ?
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * viscosity / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * viscosity / (h * h), 2));

      // Calculate the strong residual for GLS stabilization
      auto strong_residual = velocity_gradient * velocity * void_fraction +
                             // Mass Source
                             mass_source * velocity
                             // Pressure
                             + void_fraction * pressure_gradient -
                             // Viscosity and Force
                             void_fraction * viscosity * velocity_laplacian -
                             force * void_fraction + strong_residual_vec[q];

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i      = scratch_data.phi_u[q][i];
          const auto grad_phi_u_i = scratch_data.grad_phi_u[q][i];
          const auto phi_p_i      = scratch_data.phi_p[q][i];
          const auto grad_phi_p_i = scratch_data.grad_phi_p[q][i];
          const auto div_phi_u_i  = scratch_data.div_phi_u[q][i];

          double local_rhs_i = 0;

          // Navier-Stokes Residual
          local_rhs_i +=
            (
              // Momentum
              -(void_fraction * viscosity *
                  scalar_product(velocity_gradient, grad_phi_u_i) +
                viscosity * velocity_gradient * void_fraction_gradients *
                  phi_u_i) -
              velocity_gradient * velocity * void_fraction * phi_u_i
              // Mass Source
              - mass_source * velocity * phi_u_i
              // Pressure and Force
              + (void_fraction * pressure * div_phi_u_i +
                 pressure * void_fraction_gradients * phi_u_i) +
              force * void_fraction * phi_u_i
              // Continuity
              - (velocity_divergence * void_fraction +
                 velocity * void_fraction_gradients - mass_source) *
                  phi_p_i) *
            JxW;

          // PSPG GLS term
          local_rhs_i += -tau * (strong_residual * grad_phi_p_i) * JxW;

          // SUPG GLS term
          if (SUPG)
            {
              local_rhs_i +=
                -tau * (strong_residual * (grad_phi_u_i * velocity)) * JxW;
            }

          // Grad-div stabilization
          if (cfd_dem.grad_div == true)
            {
              local_rhs_i -= gamma *
                             (void_fraction * velocity_divergence +
                              velocity * void_fraction_gradients) *
                             div_phi_u_i * JxW;
            }

          local_rhs(i) += local_rhs_i;
        }
    }
}

template class GLSVansAssemblerCoreModelA<2>;

template class GLSVansAssemblerCoreModelA<3>;

template <int dim>
void
GLSVansAssemblerBDF<dim>::assemble_matrix(
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
GLSVansAssemblerBDF<dim>::assemble_rhs(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Physical properties
  const std::vector<double> &viscosity = scratch_data.viscosity;

  // Loop and quadrature informations
  const double       h          = scratch_data.cell_size;
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
          if (simulation_control->get_current_time() ==
              simulation_control->get_time_step())
            u_mag = std::max(velocity[0].norm(), 1e-12);
          else
            u_mag = std::max(velocity[1].norm(), 1e-12);
        }

      // Grad-div weight factor
      const double gamma =
        calculate_gamma(u_mag, viscosity[q], h, cfd_dem.cstar);


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

template class GLSVansAssemblerBDF<2>;

template class GLSVansAssemblerBDF<3>;

template <int dim>
void
GLSVansAssemblerDiFelice<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // particle_number is an increment that goes from 0 to n_particles_in_cell.
  // It is incremented at the end of the loop over particles and is used to
  // point to the element of the vectors relative_velocity and
  // fluid_velocity_at_particle_location corresponding to the particle being
  // looped over.
  unsigned int particle_number;
  double       cell_void_fraction = 0;
  double       C_d                = 0;
  const auto & relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p      = scratch_data.Re_particle;
  auto &      beta_drag = scratch_data.beta_drag;

  Tensor<1, dim> drag_force;


  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "GLSVansAssemblerDiFelice<dim>::calculate_particle_fluid_interactions"));

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "GLSVansAssemblerDiFelice<dim>::calculate_particle_fluid_interactions"));
  const double density = scratch_data.properties_manager.density_scale;

  const auto pic  = scratch_data.pic;
  beta_drag       = 0;
  particle_number = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      cell_void_fraction = scratch_data.cell_void_fraction[particle_number];

      // Di Felice Drag Model CD Calculation
      C_d = pow((0.63 + 4.8 / sqrt(Re_p[particle_number])), 2) *
            pow(cell_void_fraction,
                2 - (3.7 -
                     0.65 *
                       exp(-pow((1.5 - log10(Re_p[particle_number])), 2) / 2)));

      double momentum_transfer_coefficient =
        (0.5 * C_d * M_PI *
         pow(particle_properties[DEM::PropertiesIndex::dp], 2) / 4) *
        relative_velocity[particle_number].norm();

      beta_drag += momentum_transfer_coefficient;

      drag_force = density * momentum_transfer_coefficient *
                   relative_velocity[particle_number];

      for (int d = 0; d < dim; ++d)
        {
          particle_properties[DEM::PropertiesIndex::fem_force_x + d] +=
            drag_force[d];
        }

      particle_number += 1;
    }

  beta_drag = beta_drag / scratch_data.cell_volume;
}

template class GLSVansAssemblerDiFelice<2>;

template class GLSVansAssemblerDiFelice<3>;

template <int dim>
void
GLSVansAssemblerRong<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  unsigned int particle_number;
  double       cell_void_fraction = 0;
  double       C_d                = 0;
  const auto & relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p      = scratch_data.Re_particle;
  auto &      beta_drag = scratch_data.beta_drag;

  Tensor<1, dim> drag_force;

  // Physical Properties
  Assert(!scratch_data.properties_manager.is_non_newtonian(),
         RequiresConstantViscosity(
           "GLSVansAssemblerRong<dim>::calculate_particle_fluid_interactions"));

  Assert(scratch_data.properties_manager.density_is_constant(),
         RequiresConstantDensity(
           "GLSVansAssemblerRong<dim>::calculate_particle_fluid_interactions"));
  const double density = scratch_data.properties_manager.density_scale;

  const auto pic  = scratch_data.pic;
  beta_drag       = 0;
  particle_number = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      cell_void_fraction = scratch_data.cell_void_fraction[particle_number];

      // Rong Drag Model CD Calculation
      C_d = pow((0.63 + 4.8 / sqrt(Re_p[particle_number])), 2) *
            pow(cell_void_fraction,
                2 - (2.65 * (cell_void_fraction + 1) -
                     (5.3 - (3.5 * cell_void_fraction)) *
                       pow(cell_void_fraction, 2) *
                       exp(-pow(1.5 - log10(Re_p[particle_number]), 2) / 2)));

      double momentum_transfer_coefficient =
        (0.5 * C_d * M_PI *
         pow(particle_properties[DEM::PropertiesIndex::dp], 2) / 4) *
        relative_velocity[particle_number].norm();

      beta_drag += momentum_transfer_coefficient;

      drag_force = density * momentum_transfer_coefficient *
                   relative_velocity[particle_number];

      for (int d = 0; d < dim; ++d)
        {
          particle_properties[DEM::PropertiesIndex::fem_force_x + d] +=
            drag_force[d];
        }

      particle_number += 1;
    }

  beta_drag = beta_drag / scratch_data.cell_volume;
}

template class GLSVansAssemblerRong<2>;

template class GLSVansAssemblerRong<3>;

template <int dim>
void
GLSVansAssemblerDallavalle<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  unsigned int particle_number;
  double       C_d = 0;
  const auto & relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p      = scratch_data.Re_particle;
  auto &      beta_drag = scratch_data.beta_drag;

  Tensor<1, dim> drag_force;

  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "GLSVansAssemblerDallavalle<dim>::calculate_particle_fluid_interactions"));

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "GLSVansAssemblerDallavalle<dim>::calculate_particle_fluid_interactions"));
  const double density = scratch_data.properties_manager.density_scale;

  const auto pic  = scratch_data.pic;
  beta_drag       = 0;
  particle_number = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      // Dallavalle Drag Model CD Calculation
      C_d = pow((0.63 + 4.8 / sqrt(Re_p[particle_number])), 2);

      double momentum_transfer_coefficient =
        (0.5 * C_d * M_PI *
         pow(particle_properties[DEM::PropertiesIndex::dp], 2) / 4) *
        relative_velocity[particle_number].norm();

      beta_drag += momentum_transfer_coefficient;

      drag_force = density * momentum_transfer_coefficient *
                   relative_velocity[particle_number];

      for (int d = 0; d < dim; ++d)
        {
          particle_properties[DEM::PropertiesIndex::fem_force_x + d] +=
            drag_force[d];
        }

      particle_number += 1;
    }

  beta_drag = beta_drag / scratch_data.cell_volume;
}

template class GLSVansAssemblerDallavalle<2>;

template class GLSVansAssemblerDallavalle<3>;

template <int dim>
void
GLSVansAssemblerKochHill<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  unsigned int particle_number;
  double       cell_void_fraction = 0;
  const auto & relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p      = scratch_data.Re_particle;
  auto &      beta_drag = scratch_data.beta_drag;

  Tensor<1, dim> drag_force;

  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "GLSVansAssemblerKochHill<dim>::calculate_particle_fluid_interactions"));
  const double viscosity = scratch_data.properties_manager.viscosity_scale;

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "GLSVansAssemblerKochHill<dim>::calculate_particle_fluid_interactions"));
  const double density = scratch_data.properties_manager.density_scale;

  const auto pic  = scratch_data.pic;
  beta_drag       = 0;
  particle_number = 0;

  double f0 = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      cell_void_fraction = scratch_data.cell_void_fraction[particle_number];

      // Koch and Hill Drag Model Calculation
      if ((1 - cell_void_fraction) < 0.4)
        {
          f0 = (1 + 3 * sqrt((1 - cell_void_fraction) / 2) +
                (135.0 / 64) * (1 - cell_void_fraction) *
                  log(1 - cell_void_fraction) +
                16.14 * (1 - cell_void_fraction)) /
               (1 + 0.681 * (1 - cell_void_fraction) -
                8.48 * pow(1 - cell_void_fraction, 2) +
                8.16 * pow(1 - cell_void_fraction, 3));
        }
      else if ((1 - cell_void_fraction) >= 0.4)
        {
          f0 = 10 * (1 - cell_void_fraction) / pow(cell_void_fraction, 3);
        }

      double f3 = 0.0673 + 0.212 * (1 - cell_void_fraction) +
                  0.0232 / pow(cell_void_fraction, 5);

      double momentum_transfer_coefficient =
        ((18 * viscosity * pow(cell_void_fraction, 2) *
          (1 - cell_void_fraction)) /
         pow(particle_properties[DEM::PropertiesIndex::dp], 2)) *
        (f0 + 0.5 * f3 * Re_p[particle_number]) *
        (M_PI * pow(particle_properties[DEM::PropertiesIndex::dp], dim) /
         (2 * dim)) /
        (1 - cell_void_fraction);

      beta_drag += momentum_transfer_coefficient;

      drag_force = density * momentum_transfer_coefficient *
                   relative_velocity[particle_number];

      for (int d = 0; d < dim; ++d)
        {
          particle_properties[DEM::PropertiesIndex::fem_force_x + d] +=
            drag_force[d];
        }

      particle_number += 1;
    }

  beta_drag = beta_drag / scratch_data.cell_volume;
}

template class GLSVansAssemblerKochHill<2>;

template class GLSVansAssemblerKochHill<3>;

template <int dim>
void
GLSVansAssemblerBeetstra<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  unsigned int particle_number;
  double       cell_void_fraction = 0;
  double       F0                 = 0;
  const auto & relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p      = scratch_data.Re_particle;
  auto &      beta_drag = scratch_data.beta_drag;


  Tensor<1, dim> drag_force;


  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "GLSVansAssemblerBeetstra<dim>::calculate_particle_fluid_interactions"));
  const double viscosity = scratch_data.properties_manager.viscosity_scale;

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "GLSVansAssemblerBeetstra<dim>::calculate_particle_fluid_interactions"));
  const double density = scratch_data.properties_manager.density_scale;

  const auto pic  = scratch_data.pic;
  beta_drag       = 0;
  particle_number = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      cell_void_fraction = scratch_data.cell_void_fraction[particle_number];

      // Beetstra drag coefficient
      F0 = 10 * (1 - cell_void_fraction) / pow(cell_void_fraction, 2) +
           pow(cell_void_fraction, 2) *
             (1 + 1.5 * sqrt((1 - cell_void_fraction))) +
           0.413 * (Re_p[particle_number]) / (24 * pow(cell_void_fraction, 2)) *
             ((1 / cell_void_fraction) +
              3 * (1 - cell_void_fraction) * cell_void_fraction +
              8.4 * pow((Re_p[particle_number]), -0.343)) /
             (1 + pow(10, 3 * (1 - cell_void_fraction)) *
                    pow((Re_p[particle_number]),
                        -(1 + 4 * (1 - cell_void_fraction)) * 0.5));

      double momentum_transfer_coefficient =
        F0 * 3 * M_PI * viscosity * cell_void_fraction *
        particle_properties[DEM::PropertiesIndex::dp];

      beta_drag += momentum_transfer_coefficient;

      drag_force = density * momentum_transfer_coefficient *
                   relative_velocity[particle_number];

      for (int d = 0; d < dim; ++d)
        {
          particle_properties[DEM::PropertiesIndex::fem_force_x + d] +=
            drag_force[d];
        }

      particle_number += 1;
    }

  beta_drag = beta_drag / scratch_data.cell_volume;
}

template class GLSVansAssemblerBeetstra<2>;

template class GLSVansAssemblerBeetstra<3>;

template <int dim>
void
GLSVansAssemblerGidaspow<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // particle_number is an increment that goes from 0 to n_particles_in_cell.
  // It is incremented at the end of the loop over particles and is used to
  // point to the element of the vectors relative_velocity and
  // fluid_velocity_at_particle_location corresponding to the particle being
  // looped over.
  unsigned int particle_number;
  double       cell_void_fraction = 0;
  double       C_d                = 0;
  const auto & relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &average_relative_velocity =
    scratch_data.average_fluid_particle_relative_velocity;
  const auto &Re_p      = scratch_data.Re_particle;
  auto &      beta_drag = scratch_data.beta_drag;

  Tensor<1, dim> drag_force;


  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "GLSVansAssemblerGidaspow<dim>::calculate_particle_fluid_interactions"));
  const double viscosity = scratch_data.properties_manager.viscosity_scale;

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "GLSVansAssemblerGidaspow<dim>::calculate_particle_fluid_interactions"));
  const double density = scratch_data.properties_manager.density_scale;


  const auto pic                           = scratch_data.pic;
  double     momentum_transfer_coefficient = 0;
  beta_drag                                = 0;
  particle_number                          = 0;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      cell_void_fraction = scratch_data.cell_void_fraction[particle_number];

      // Gidaspow Drag Model CD Calculation
      if (Re_p[particle_number] < 1000)
        {
          C_d = 24 / Re_p[particle_number] *
                (1 + 0.15 * pow(Re_p[particle_number], 0.687));
        }
      else
        {
          C_d = 0.44;
        }

      if (cell_void_fraction >= 0.8)
        {
          momentum_transfer_coefficient =
            0.75 * C_d * cell_void_fraction * average_relative_velocity.norm() *
            (1 - cell_void_fraction) * pow(cell_void_fraction, -2.65) /
            particle_properties[DEM::PropertiesIndex::dp];
        }
      else
        {
          // Assuming the sphericity of particles = 1
          momentum_transfer_coefficient =
            150 * pow((1 - cell_void_fraction), 2) * viscosity /
              (cell_void_fraction *
               pow(particle_properties[DEM::PropertiesIndex::dp], 2)) +
            1.75 * (1 - cell_void_fraction) * average_relative_velocity.norm() /
              particle_properties[DEM::PropertiesIndex::dp];
        }

      beta_drag += particle_properties[DEM::PropertiesIndex::mass] *
                   momentum_transfer_coefficient / density;

      drag_force = particle_properties[DEM::PropertiesIndex::mass] *
                   momentum_transfer_coefficient *
                   relative_velocity[particle_number];

      for (int d = 0; d < dim; ++d)
        {
          particle_properties[DEM::PropertiesIndex::fem_force_x + d] +=
            drag_force[d];
        }

      particle_number += 1;
    }

  beta_drag = beta_drag / scratch_data.cell_volume;
}

template class GLSVansAssemblerGidaspow<2>;

template class GLSVansAssemblerGidaspow<3>;

template <int dim>
void
GLSVansAssemblerSaffmanMei<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // particle_number is an increment that goes from 0 to n_particles_in_cell.
  // It is incremented at the end of the loop over particles and is used to
  // point to the element of the vectors relative_velocity and
  // fluid_velocity_at_particle_location corresponding to the particle being
  // looped over.

  // This implementation follows the formulation in the book "Multiphase Flows
  // with Droplets and Particles" by Crowe et al. (2011) and the brief
  // communication article "An approximate expression for the shear lift force
  // on a spherical particle at finite reynolds number" by Mei (1992)
  unsigned int particle_number;
  double       C_s   = 0;
  double       alpha = 0;

  const auto &relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p = scratch_data.Re_particle;

  auto &velocity_curls_2d =
    scratch_data.fluid_velocity_curls_at_particle_location_2d;
  auto &velocity_curls_3d =
    scratch_data.fluid_velocity_curls_at_particle_location_3d;
  auto &undisturbed_flow_force = scratch_data.undisturbed_flow_force;

  Tensor<1, dim> lift_force;

  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "GLSVansAssemblerSaffmanMei<dim>::calculate_particle_fluid_interactions"));
  const double viscosity = scratch_data.properties_manager.viscosity_scale;

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "GLSVansAssemblerSaffmanMei<dim>::calculate_particle_fluid_interactions"));
  const double density = scratch_data.properties_manager.density_scale;

  const auto pic  = scratch_data.pic;
  particle_number = 0;

  if constexpr (dim == 2)
    {
      for (auto &particle : pic)
        {
          auto particle_properties = particle.get_properties();

          // Saffman-Mei coefficient
          alpha = 0.5 * particle_properties[DEM::PropertiesIndex::dp] /
                  relative_velocity[particle_number].norm() *
                  abs(velocity_curls_2d[particle_number][0] + 1e-9);

          if (Re_p[particle_number] <= 40)
            C_s =
              (1 - 0.3314 * sqrt(alpha)) * exp(-0.1 * Re_p[particle_number]) +
              0.3314 * sqrt(alpha);
          else if (Re_p[particle_number] > 40)
            C_s = 0.0524 * sqrt(alpha * Re_p[particle_number]);

          // Vorticity vector
          Tensor<1, 2> vorticity;
          vorticity[0] = velocity_curls_2d[particle_number][0];
          vorticity[1] = -velocity_curls_2d[particle_number][1];

          // Saffman Lift force
          lift_force[0] =
            C_s * 1.61 * particle_properties[DEM::PropertiesIndex::dp] *
            particle_properties[DEM::PropertiesIndex::dp] * density *
            sqrt(viscosity + DBL_MIN) /
            sqrt(velocity_curls_2d[particle_number].norm()) *
            (relative_velocity[particle_number][0] * vorticity[1]);

          lift_force[1] =
            C_s * 1.61 * particle_properties[DEM::PropertiesIndex::dp] *
            particle_properties[DEM::PropertiesIndex::dp] * density *
            sqrt(viscosity + DBL_MIN) / sqrt(vorticity.norm() + 1e-9) *
            (relative_velocity[particle_number][1] * vorticity[0]);


          for (int d = 0; d < dim; ++d)
            {
              // Apply lift force on the particle
              particle_properties[DEM::PropertiesIndex::fem_force_x + d] +=
                lift_force[d];

              // Apply lift force on the fluid
              undisturbed_flow_force[d] +=
                lift_force[d] / scratch_data.cell_volume;
            }

          particle_number += 1;
        }
    }
  else if constexpr (dim == 3)
    {
      for (auto &particle : pic)
        {
          auto particle_properties = particle.get_properties();

          // Saffman-Mei coefficient C_s
          alpha = 0.5 * particle_properties[DEM::PropertiesIndex::dp] /
                  relative_velocity[particle_number].norm() *
                  (velocity_curls_3d[particle_number].norm() + 1e-9);

          if (Re_p[particle_number] <= 40)
            C_s =
              (1 - 0.3314 * sqrt(alpha)) * exp(-0.1 * Re_p[particle_number]) +
              0.3314 * sqrt(alpha);
          else if (Re_p[particle_number] > 40)
            C_s = 0.0524 * sqrt(alpha * Re_p[particle_number]);

          // Vorticity tensor
          Tensor<1, 3> vorticity = velocity_curls_3d[particle_number];

          lift_force =
            C_s * 1.61 * particle_properties[DEM::PropertiesIndex::dp] *
            particle_properties[DEM::PropertiesIndex::dp] * density *
            sqrt(viscosity + DBL_MIN) / sqrt(vorticity.norm() + 1e-9) *
            (cross_product_3d(relative_velocity[particle_number], vorticity));

          for (int d = 0; d < dim; ++d)
            {
              // Apply lift force on the particle
              particle_properties[DEM::PropertiesIndex::fem_force_x + d] +=
                lift_force[d];

              // Apply lift force on the fluid
              undisturbed_flow_force[d] +=
                lift_force[d] / scratch_data.cell_volume;
            }
          particle_number += 1;
        }
    }
}

template class GLSVansAssemblerSaffmanMei<2>;

template class GLSVansAssemblerSaffmanMei<3>;

template <int dim>
void
GLSVansAssemblerMagnus<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  // particle_number is an increment that goes from 0 to n_particles_in_cell.
  // It is incremented at the end of the loop over particles and is used to
  // point to the element of the vectors relative_velocity and
  // fluid_velocity_at_particle_location corresponding to the particle being
  // looped over.

  // This implementation follows the formulation in the book "Multiphase Flows
  // with Droplets and Particles" by Crowe et al. (2011).
  unsigned int particle_number;
  double       C_m = 0;

  const auto &relative_velocity =
    scratch_data.fluid_particle_relative_velocity_at_particle_location;
  const auto &Re_p = scratch_data.Re_particle;

  auto &undisturbed_flow_force = scratch_data.undisturbed_flow_force;

  Tensor<1, dim> lift_force;

  // Physical Properties
  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "GLSVansAssemblerSaffmanMei<dim>::calculate_particle_fluid_interactions"));
  const double density = scratch_data.properties_manager.density_scale;

  const auto pic  = scratch_data.pic;
  particle_number = 0;

  if constexpr (dim == 2)
    {
      for (auto &particle : pic)
        {
          auto particle_properties = particle.get_properties();

          // Spin parameter
          double spin_parameter =
            particle_properties[DEM::PropertiesIndex::dp] *
            abs(particle_properties[DEM::PropertiesIndex::omega_z]) /
            (2.0 * relative_velocity[particle_number].norm());

          // Magnus lift coefficient
          if (spin_parameter > 1.0 && spin_parameter < 6.0 &&
              Re_p[particle_number] > 10.0 && Re_p[particle_number] < 140.0)
            {
              // Oesterlé and Dinh (1998)
              C_m = 0.45 + (2 * spin_parameter - 0.45) *
                             exp(-0.075 * pow(spin_parameter, 0.4) *
                                 pow(Re_p[particle_number], 0.7));
            }
          else
            C_m = 2.0 * spin_parameter;

          // Magnus Lift force
          lift_force[0] =
            0.125 * M_PI *
            pow(particle_properties[DEM::PropertiesIndex::dp], 2.0) * density *
            relative_velocity[particle_number].norm() *
            (particle_properties[DEM::PropertiesIndex::omega_z] /
             abs(particle_properties[DEM::PropertiesIndex::omega_z]) *
             relative_velocity[particle_number][1]);

          lift_force[1] =
            0.125 * M_PI *
            pow(particle_properties[DEM::PropertiesIndex::dp], 2.0) * density *
            relative_velocity[particle_number].norm() *
            (particle_properties[DEM::PropertiesIndex::omega_z] /
             abs(particle_properties[DEM::PropertiesIndex::omega_z]) *
             relative_velocity[particle_number][0]);


          for (int d = 0; d < dim; ++d)
            {
              // Apply lift force on the particle
              particle_properties[DEM::PropertiesIndex::fem_force_x + d] +=
                lift_force[d];

              // Apply lift force on the fluid
              undisturbed_flow_force[d] +=
                lift_force[d] / scratch_data.cell_volume;
            }
          particle_number += 1;
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
              omega[d] = particle_properties[DEM::PropertiesIndex::omega_x + d];
            }

          // Spin parameter
          double spin_parameter =
            particle_properties[DEM::PropertiesIndex::dp] * omega.norm() /
            (2.0 * relative_velocity[particle_number].norm());

          // Magnus lift coefficient
          if (spin_parameter > 1.0 && spin_parameter < 6.0 &&
              Re_p[particle_number] > 10.0 && Re_p[particle_number] < 140.0)
            {
              // Oesterlé and Dinh (1998)
              C_m = 0.45 + (2 * spin_parameter - 0.45) *
                             exp(-0.075 * pow(spin_parameter, 0.4) *
                                 pow(Re_p[particle_number], 0.7));
            }
          else
            C_m = 2.0 * spin_parameter;

          Tensor<1, dim> rotational_vector = omega / omega.norm();

          // Magnus Lift force
          lift_force = 0.125 * M_PI *
                       pow(particle_properties[DEM::PropertiesIndex::dp], 2.0) *
                       density * relative_velocity[particle_number].norm() *
                       (cross_product_3d(rotational_vector,
                                         relative_velocity[particle_number]));


          for (int d = 0; d < dim; ++d)
            {
              // Apply lift force on the particle
              particle_properties[DEM::PropertiesIndex::fem_force_x + d] +=
                lift_force[d];

              // Apply lift force on the fluid
              undisturbed_flow_force[d] +=
                lift_force[d] / scratch_data.cell_volume;
            }
          particle_number += 1;
        }
    }
}

template class GLSVansAssemblerMagnus<2>;
template class GLSVansAssemblerMagnus<3>;


template <int dim>
void
GLSVansAssemblerBuoyancy<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  const auto   pic = scratch_data.pic;
  Tensor<1, 3> buoyancy_force;

  // Physical Properties
  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "GLSVansAssemblerBuoyancy<dim>::calculate_particle_fluid_interactions"));

  const double density = scratch_data.properties_manager.density_scale;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      // Buoyancy Force
      buoyancy_force =
        -lagrangian_physical_properties.g * (4.0 / 3) * M_PI *
        pow((particle_properties[DEM::PropertiesIndex::dp] / 2.0), 3);

      for (int d = 0; d < dim; ++d)
        {
          particle_properties[DEM::PropertiesIndex::fem_force_x + d] +=
            buoyancy_force[d] * density;
        }
    }
}

template class GLSVansAssemblerBuoyancy<2>;

template class GLSVansAssemblerBuoyancy<3>;

template <int dim>
void
GLSVansAssemblerPressureForce<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  unsigned int particle_number;
  const auto   pic                    = scratch_data.pic;
  auto &       undisturbed_flow_force = scratch_data.undisturbed_flow_force;
  auto         pressure_gradients =
    scratch_data.fluid_pressure_gradients_at_particle_location;
  Tensor<1, dim> pressure_force;

  particle_number = 0;

  // Physical Properties
  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "GLSVansAssemblerPressureForce<dim>::calculate_particle_fluid_interactions"));


  const double density = scratch_data.properties_manager.density_scale;
  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      // Pressure Force
      pressure_force =
        -(M_PI * pow(particle_properties[DEM::PropertiesIndex::dp], dim) /
          (2 * dim)) *
        pressure_gradients[particle_number];

      for (int d = 0; d < dim; ++d)
        {
          particle_properties[DEM::PropertiesIndex::fem_force_x + d] +=
            pressure_force[d] * density;

          // Apply pressure force to the particles only, when we are solving
          // model A of the VANS. When we are solving Model B, apply the
          // pressure force back on the fluid by lumping it in the
          // undisturbed_flow_force.
          if (cfd_dem.vans_model == Parameters::VANSModel::modelB)
            {
              undisturbed_flow_force[d] +=
                pressure_force[d] / scratch_data.cell_volume;
            }
        }

      particle_number += 1;
    }
}

template class GLSVansAssemblerPressureForce<2>;

template class GLSVansAssemblerPressureForce<3>;

template <int dim>
void
GLSVansAssemblerShearForce<dim>::calculate_particle_fluid_interactions(
  NavierStokesScratchData<dim> &scratch_data)
{
  unsigned int particle_number;
  const auto   pic                    = scratch_data.pic;
  auto &       undisturbed_flow_force = scratch_data.undisturbed_flow_force;
  auto &       velocity_laplacians =
    scratch_data.fluid_velocity_laplacian_at_particle_location;
  Tensor<1, dim> shear_force;

  particle_number = 0;

  // Viscosity and density are currently assumed constant from the particle
  // point of view.
  // Physical Properties
  Assert(
    !scratch_data.properties_manager.is_non_newtonian(),
    RequiresConstantViscosity(
      "GLSVansAssemblerDallavalle<dim>::calculate_particle_fluid_interactions"));
  const double viscosity = scratch_data.properties_manager.viscosity_scale;

  Assert(
    scratch_data.properties_manager.density_is_constant(),
    RequiresConstantDensity(
      "GLSVansAssemblerDallavalle<dim>::calculate_particle_fluid_interactions"));
  const double density = scratch_data.properties_manager.density_scale;

  // Loop over particles in cell
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      // Shear Force
      shear_force =
        -(M_PI * pow(particle_properties[DEM::PropertiesIndex::dp], dim) /
          (2 * dim)) *
        viscosity * velocity_laplacians[particle_number];

      for (int d = 0; d < dim; ++d)
        {
          particle_properties[DEM::PropertiesIndex::fem_force_x + d] +=
            shear_force[d] * density;

          // Apply shear force to the particles only, when we are solving
          // model A of the VANS. When we are solving Model B, apply the shear
          // force back on the fluid by lumping it in the
          // undisturbed_flow_force.
          if (cfd_dem.vans_model == Parameters::VANSModel::modelB)
            {
              undisturbed_flow_force[d] +=
                shear_force[d] / scratch_data.cell_volume;
            }
        }

      particle_number += 1;
    }
}

template class GLSVansAssemblerShearForce<2>;

template class GLSVansAssemblerShearForce<3>;


template <int dim>
void
GLSVansAssemblerFPI<dim>::assemble_matrix(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual        = copy_data.strong_residual;
  auto &strong_jacobian        = copy_data.strong_jacobian;
  auto &local_matrix           = copy_data.local_matrix;
  auto  beta_drag              = scratch_data.beta_drag;
  auto &undisturbed_flow_force = scratch_data.undisturbed_flow_force;
  const Tensor<1, dim> average_particles_velocity =
    scratch_data.average_particle_velocity;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const Tensor<1, dim> velocity = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Calculate the strong residual for GLS stabilization
      if (cfd_dem.vans_model == Parameters::VANSModel::modelB)
        {
          strong_residual[q] += // Drag Force
            (beta_drag * (velocity - average_particles_velocity) +
             undisturbed_flow_force);
        }
      else if (cfd_dem.vans_model == Parameters::VANSModel::modelA)
        {
          strong_residual[q] += // Drag Force
            beta_drag * (velocity - average_particles_velocity);
        }

      // We loop over the column first to prevent recalculation
      // of the strong jacobian in the inner loop
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
GLSVansAssemblerFPI<dim>::assemble_rhs(
  NavierStokesScratchData<dim> &        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual        = copy_data.strong_residual;
  auto &local_rhs              = copy_data.local_rhs;
  auto  beta_drag              = scratch_data.beta_drag;
  auto &undisturbed_flow_force = scratch_data.undisturbed_flow_force;
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
      if (cfd_dem.vans_model == Parameters::VANSModel::modelB)
        {
          strong_residual[q] += // Drag Force
            (beta_drag * (velocity - average_particles_velocity) +
             undisturbed_flow_force);
        }
      else if (cfd_dem.vans_model == Parameters::VANSModel::modelA)
        {
          strong_residual[q] += // Drag Force
            beta_drag * (velocity - average_particles_velocity);
        }

      // Assembly of the right-hand side
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_u_i = scratch_data.phi_u[q][i];
          // Drag Force
          //  Model B of the VANS
          if (cfd_dem.vans_model == Parameters::VANSModel::modelB)
            {
              local_rhs(i) -=
                (beta_drag * (velocity - average_particles_velocity) +
                 undisturbed_flow_force) *
                phi_u_i * JxW;
            }
          //  Model A of the VANS
          if (cfd_dem.vans_model == Parameters::VANSModel::modelA)
            {
              local_rhs(i) -=
                (beta_drag * (velocity - average_particles_velocity)) *
                phi_u_i * JxW;
            }
        }
    }
}

template class GLSVansAssemblerFPI<2>;

template class GLSVansAssemblerFPI<3>;