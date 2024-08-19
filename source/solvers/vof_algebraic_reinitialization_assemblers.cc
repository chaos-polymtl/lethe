#include <solvers/vof_algebraic_reinitialization_assemblers.h>

template <int dim>
void
VOFAlgebraicReinitializationAssemblerCore<dim>::assemble_matrix(
  const VOFAlgebraicReinitializationScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData                          &copy_data)
{
  // Time stepping scheme
  //  const auto method = this->simulation_control->get_assembly_method();

  // Loop and quadrature information
  const std::vector<double> &JxW_vec    = scratch_data.JxW;
  const unsigned int         n_q_points = scratch_data.n_q_points;
  const unsigned int         n_dofs     = scratch_data.n_dofs;
  const double               h          = scratch_data.cell_size;

  // TODO AMISHGA check if stabilization is necessary

  //  // Time steps and inverse time steps which is used for stabilization
  //  constant std::vector<double> time_steps_vector =
  //    this->simulation_control->get_time_steps_vector();
  //  const double dt  = time_steps_vector[0];
  //  const double sdt = 1. / dt;

  // Diffusivity coefficient
  const double diffusivity_coefficient =
    2 * std::pow(h, 0.9); // TODO AMISHGA make it with parameters

  // Copy data elements
  //  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix = copy_data.local_matrix;

  // Assemble local matrix
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const double vof_phase_value = scratch_data.present_vof_phase_values[q];
      const Tensor<1, dim> vof_phase_gradient =
        scratch_data.present_vof_phase_gradients[q];

      // Store JxW in a local variable for faster access
      const double JxW = JxW_vec[q];

      // Compute normal vector with the VOF solution
      const Tensor<1, dim> interface_normal =
        vof_phase_gradient / vof_phase_gradient.norm();

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto grad_ksi_phase_i = scratch_data.grad_ksi[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto ksi_phase_j      = scratch_data.ksi[q][j];
              const auto grad_ksi_phase_j = scratch_data.grad_ksi[q][j];

              // Weak form (compressible and diffusive terms):
              // \f$ \int_\Omega \nabla v \cdot \left[\epsilon (\nabla
              // \phi_\mathrm{reinit} \cdot \mathbf{n}) \mathbf{n}
              // - (\phi_\mathrm{reinit} (1-\phi_\mathrm{reinit) \mathbf{n})
              // \right] \mathrm{d} \Omega \f$
              local_matrix(i, j) +=
                (grad_ksi_phase_i *
                 ((ksi_phase_j - 2 * vof_phase_value * ksi_phase_j) *
                    interface_normal -
                  diffusivity_coefficient *
                    (grad_ksi_phase_j * interface_normal) * interface_normal)) *
                JxW;
            } // end j loop on dofs
        }     // end i loop on dofs
    }         // end loop on quadrature points
}

template <int dim>
void
VOFAlgebraicReinitializationAssemblerCore<dim>::assemble_rhs(
  const VOFAlgebraicReinitializationScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData                          &copy_data)
{
  // Scheme and physical properties
  //  const auto method = this->simulation_control->get_assembly_method();

  // Loop and quadrature information
  const std::vector<double> &JxW_vec    = scratch_data.JxW;
  const unsigned int         n_q_points = scratch_data.n_q_points;
  const unsigned int         n_dofs     = scratch_data.n_dofs;
  const double               h          = scratch_data.cell_size;

  // TODO AMISHGA check if stabilization is necessary

  //  // Time steps and inverse time steps which is used for stabilization
  //  constant std::vector<double> time_steps_vector =
  //    this->simulation_control->get_time_steps_vector();
  //  const double dt  = time_steps_vector[0];
  //  const double sdt = 1. / dt;

  // Diffusivity coefficient
  const double diffusivity_coefficient =
    2 * std::pow(h, 0.9); // TODO AMISHGA make it with parameters

  // Copy data elements
  //  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs = copy_data.local_rhs;

  // Assemble local right-hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const double vof_phase_value = scratch_data.present_vof_phase_values[q];
      const Tensor<1, dim> vof_phase_gradient =
        scratch_data.present_vof_phase_gradients[q];

      // Store JxW in a local variable for faster access
      const double JxW = JxW_vec[q];

      // Compute normal vector with the VOF solution
      const Tensor<1, dim> interface_normal =
        vof_phase_gradient / vof_phase_gradient.norm();

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto grad_phi_phase_i = scratch_data.grad_ksi[q][i];

          local_rhs(i) -=
            (grad_phi_phase_i *
             ((vof_phase_value - Utilities::fixed_power<2>(vof_phase_value)) *
                interface_normal -
              diffusivity_coefficient *
                (vof_phase_gradient * interface_normal) * interface_normal)) *
            JxW;
        } // end i loop on dofs
    }     // end loop on quadrature points
}

template class VOFAlgebraicReinitializationAssemblerCore<2>;
template class VOFAlgebraicReinitializationAssemblerCore<3>;


template <int dim>
void
VOFAlgebraicReinitializationAssemblerBDF<dim>::assemble_matrix(
  const VOFAlgebraicReinitializationScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData                          &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  //  auto &strong_jacobian = copy_data.strong_jacobian;
  //  auto &strong_residual = copy_data.strong_residual;
  auto &local_matrix = copy_data.local_matrix;

  // Time stepping information
  const auto method = this->simulation_control->get_assembly_method();
  // Vector for the BDF coefficients
  const Vector<double> &bdf_coefs =
    this->simulation_control->get_bdf_coefficients();
  std::vector<double> phase_value(1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      phase_value[0] = scratch_data.present_phase_algebraic_reinit_values[q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        phase_value[p + 1] =
          scratch_data.previous_phase_algebraic_reinit_values[p][q];

      //      for (unsigned int p = 0; p < number_of_previous_solutions(method)
      //      + 1;
      //           ++p)
      //        {
      //          strong_residual[q] += (bdf_coefs[p] * phase_value[p]);
      //        }
      //
      //      for (unsigned int j = 0; j < n_dofs; ++j)
      //        {
      //          strong_jacobian[q][j] += bdf_coefs[0] *
      //          scratch_data.phi[q][j];
      //        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_phase_i = scratch_data.ksi[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto phi_phase_j = scratch_data.ksi[q][j];

              local_matrix(i, j) +=
                phi_phase_j * phi_phase_i * bdf_coefs[0] * JxW[q];
            }
        }
    }
}

template <int dim>
void
VOFAlgebraicReinitializationAssemblerBDF<dim>::assemble_rhs(
  const VOFAlgebraicReinitializationScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData                          &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  //  auto &strong_residual = copy_data.strong_residual;
  auto &local_rhs = copy_data.local_rhs;

  // Time stepping information
  const auto method = this->simulation_control->get_assembly_method();
  // Vector for the BDF coefficients
  const Vector<double> &bdf_coefs =
    this->simulation_control->get_bdf_coefficients();
  std::vector<double> phase_value(1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      phase_value[0] = scratch_data.present_phase_algebraic_reinit_values[q];


      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        phase_value[p + 1] =
          scratch_data.previous_phase_algebraic_reinit_values[p][q];



      //      for (unsigned int p = 0; p < number_of_previous_solutions(method)
      //      + 1;
      //           ++p)
      //        {
      //          strong_residual[q] += (bdf_coefs[p] * phase_value[p]);
      //        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double phi_phase_i = scratch_data.ksi[q][i];
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

template class VOFAlgebraicReinitializationAssemblerBDF<2>;
template class VOFAlgebraicReinitializationAssemblerBDF<3>;
