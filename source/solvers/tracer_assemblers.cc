#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/copy_data.h>
#include <solvers/tracer_assemblers.h>


template <int dim>
void
TracerAssemblerCore<dim>::assemble_matrix(TracerScratchData<dim> &scratch_data,
                                          StabilizedMethodsCopyData &copy_data)
{
  // Scheme and physical properties
  const double diffusivity = physical_properties.tracer_diffusivity;
  const auto   method      = this->simulation_control->get_assembly_method();

  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;


  // Copy data elements
  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix        = copy_data.local_matrix;

  // assembling local matrix and right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const Tensor<1, dim> tracer_gradient = scratch_data.tracer_gradients[q];
      const Tensor<1, dim> velocity        = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Shock capturing viscosity term
      const double order = scratch_data.fe_values_tracer.get_fe().degree;

      const double vdcdd = (0.5 * h) * (velocity.norm() * velocity.norm()) *
                           pow(tracer_gradient.norm() * h, order);

      const double   tolerance = 1e-12;
      Tensor<1, dim> s         = velocity / (velocity.norm() + tolerance);
      Tensor<1, dim> r = tracer_gradient / (tracer_gradient.norm() + tolerance);

      const Tensor<2, dim> k_corr      = (r * s) * outer_product(s, s);
      const Tensor<2, dim> rr          = outer_product(r, r);
      const Tensor<2, dim> dcdd_factor = rr - k_corr;


      const double d_vdcdd = order * (0.5 * h * h) *
                             (velocity.norm() * velocity.norm()) *
                             pow(tracer_gradient.norm() * h, order - 1);



      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), tolerance);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step
      const double tau =
        is_steady(method) ?
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * diffusivity / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * diffusivity / (h * h), 2));

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const Tensor<1, dim> grad_phi_T_j = scratch_data.grad_phi[q][j];
          const double laplacian_phi_T_j    = scratch_data.laplacian_phi[q][j];
          strong_jacobian_vec[q][j] +=
            velocity * grad_phi_T_j - diffusivity * laplacian_phi_T_j;

          if (DCDD)
            strong_jacobian_vec[q][j] += -vdcdd * laplacian_phi_T_j;
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_T_i      = scratch_data.phi[q][i];
          const auto grad_phi_T_i = scratch_data.grad_phi[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const Tensor<1, dim> grad_phi_T_j = scratch_data.grad_phi[q][j];

              // Weak form : - D * laplacian T +  u * gradT - f=0
              local_matrix(i, j) += (diffusivity * grad_phi_T_i * grad_phi_T_j +
                                     phi_T_i * velocity * grad_phi_T_j) *
                                    JxW;



              local_matrix(i, j) += tau * strong_jacobian_vec[q][j] *
                                    (grad_phi_T_i * velocity) * JxW;

              if (DCDD)
                {
                  local_matrix(i, j) +=
                    (vdcdd * scalar_product(grad_phi_T_j,
                                            dcdd_factor * grad_phi_T_i) +
                     d_vdcdd * grad_phi_T_j.norm() *
                       scalar_product(tracer_gradient,
                                      dcdd_factor * grad_phi_T_i)) *
                    JxW;
                }
            }
        }
    } // end loop on quadrature points
}


template <int dim>
void
TracerAssemblerCore<dim>::assemble_rhs(TracerScratchData<dim> &   scratch_data,
                                       StabilizedMethodsCopyData &copy_data)
{
  // Scheme and physical properties
  const double diffusivity = physical_properties.tracer_diffusivity;
  const auto   method      = this->simulation_control->get_assembly_method();

  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
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
      const Tensor<1, dim> tracer_gradient  = scratch_data.tracer_gradients[q];
      const double         tracer_laplacian = scratch_data.tracer_laplacians[q];
      const Tensor<1, dim> velocity         = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Shock capturing viscosity term
      const double order = scratch_data.fe_values_tracer.get_fe().degree;


      const double vdcdd = (0.5 * h) * (velocity.norm() * velocity.norm()) *
                           pow(tracer_gradient.norm() * h, order);

      const double   tolerance = 1e-12;
      Tensor<1, dim> s         = velocity / (velocity.norm() + tolerance);
      Tensor<1, dim> r = tracer_gradient / (tracer_gradient.norm() + tolerance);

      const Tensor<2, dim> k_corr      = (r * s) * outer_product(s, s);
      const Tensor<2, dim> rr          = outer_product(r, r);
      const Tensor<2, dim> dcdd_factor = rr - k_corr;

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), tolerance);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step
      const double tau =
        is_steady(method) ?
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * diffusivity / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * diffusivity / (h * h), 2));

      // Calculate the strong residual for GLS stabilization
      strong_residual_vec[q] +=
        velocity * tracer_gradient - diffusivity * tracer_laplacian;

      if (DCDD)
        strong_residual_vec[q] += -vdcdd * tracer_laplacian;

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_T_i      = scratch_data.phi[q][i];
          const auto grad_phi_T_i = scratch_data.grad_phi[q][i];

          // rhs for : - D * laplacian T +  u * grad T - f=0
          local_rhs(i) -= (diffusivity * grad_phi_T_i * tracer_gradient +
                           phi_T_i * velocity * tracer_gradient -
                           scratch_data.source[q] * phi_T_i) *
                          JxW;

          local_rhs(i) -=
            tau * (strong_residual_vec[q] * (grad_phi_T_i * velocity)) * JxW;

          if (DCDD)
            {
              local_rhs(i) +=
                -vdcdd *
                scalar_product(tracer_gradient, dcdd_factor * grad_phi_T_i) *
                JxW;
            }
        }
    } // end loop on quadrature points
}

template class TracerAssemblerCore<2>;
template class TracerAssemblerCore<3>;

template <int dim>
void
TracerAssemblerBDF<dim>::assemble_matrix(TracerScratchData<dim> &scratch_data,
                                         StabilizedMethodsCopyData &copy_data)
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
  Vector<double>      bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<double> tracer(1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      tracer[0] = scratch_data.tracer_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        tracer[p + 1] = scratch_data.previous_tracer_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += bdf_coefs[p] * tracer[p];
        }

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian[q][j] += bdf_coefs[0] * scratch_data.phi[q][j];
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double phi_u_i = scratch_data.phi[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const double phi_u_j = scratch_data.phi[q][j];

              local_matrix(i, j) += phi_u_j * phi_u_i * bdf_coefs[0] * JxW[q];
            }
        }
    }
}

template <int dim>
void
TracerAssemblerBDF<dim>::assemble_rhs(TracerScratchData<dim> &   scratch_data,
                                      StabilizedMethodsCopyData &copy_data)
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
  Vector<double>      bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<double> tracer(1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      tracer[0] = scratch_data.tracer_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        tracer[p + 1] = scratch_data.previous_tracer_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += (bdf_coefs[p] * tracer[p]);
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double phi_u_i     = scratch_data.phi[q][i];
          double       local_rhs_i = 0;
          for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
               ++p)
            {
              local_rhs_i -= bdf_coefs[p] * (tracer[p] * phi_u_i);
            }
          local_rhs(i) += local_rhs_i * JxW[q];
        }
    }
}

template class TracerAssemblerBDF<2>;
template class TracerAssemblerBDF<3>;
