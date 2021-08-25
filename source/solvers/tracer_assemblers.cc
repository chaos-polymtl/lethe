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
      const double order =
        this->simulation_parameters.fem_parameters.tracer_order;

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

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_T_i      = scratch_data.phi[i];
          const auto grad_phi_T_i = scratch_data.grad_phi[i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto phi_T_j           = scratch_data.phi[j];
              const auto grad_phi_T_j      = scratch_data.grad_phi[j];
              const auto laplacian_phi_T_j = scratch_data.laplacian_phi[j];


              // Weak form : - D * laplacian T +  u * gradT - f=0
              local_matrix(i, j) += (diffusivity * grad_phi_T_i * grad_phi_T_j +
                                     phi_T_i * velocity * grad_phi_T_j) *
                                    JxW;

              strong_jacobian_vec[q][j] +=
                velocity * grad_phi_T_j - diffusivity * laplacian_phi_T_j;

              if (DCDD)
                strong_jacobian_vec[q][j] += -vdcdd * laplacian_phi_T_j;


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
      const double order =
        this->simulation_parameters.fem_parameters.tracer_order;

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

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_T_i      = scratch_data.phi[i];
          const auto grad_phi_T_i = scratch_data.grad_phi[i];

          // rhs for : - D * laplacian T +  u * grad T - f=0
          local_rhs(i) -= (diffusivity * grad_phi_T_i * tracer_gradient +
                           phi_T_i * velocity * tracer_gradient -
                           scratch_data.source[q] * phi_T_i) *
                          JxW;

          // Calculate the strong residual for GLS stabilization
          strong_residual_vec[q] =
            velocity * tracer_gradient - diffusivity * tracer_laplacian;

          if (DCDD)
            strong_residual_vec[q] += -vdcdd * tracer_laplacian;

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
