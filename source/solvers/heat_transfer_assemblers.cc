#include <core/bdf.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/copy_data.h>
#include <solvers/heat_transfer.h>
#include <solvers/heat_transfer_assemblers.h>

template <int dim>
void
HeatTransferAssemblerCore<dim>::assemble_matrix(
  HeatTransferScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{
  // Gather physical properties in case of mono fluids simulations (to be
  // modified by cell in case of multiple fluids simulations)
  const std::vector<double> &density       = scratch_data.density;
  const std::vector<double> &specific_heat = scratch_data.specific_heat;
  const std::vector<double> &thermal_conductivity =
    scratch_data.thermal_conductivity;


  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const double       h          = scratch_data.cell_size;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Vector for the BDF coefficients
  // The coefficients are stored in the following fashion :
  // 0 - n+1
  // 1 - n
  // 2 - n-1
  // 3 - n-2
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Time steps and inverse time steps which is used for numerous calculations
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Copy data elements
  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix        = copy_data.local_matrix;

  // assembling local matrix and right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      const double rho_cp = density[q] * specific_heat[q];
      const double alpha  = thermal_conductivity[q] / rho_cp;

      const auto method = this->simulation_control->get_assembly_method();

      // Store JxW in local variable for faster access
      const double JxW = JxW_vec[q];

      const auto velocity = scratch_data.velocity_values[q];

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step
      const double tau =
        is_steady(method) ?
          1. / std::sqrt(std::pow(2. * rho_cp * u_mag / h, 2) +
                         9 * std::pow(4 * alpha / (h * h), 2)) :
          1. /
            std::sqrt(std::pow(sdt, 2) + std::pow(2. * rho_cp * u_mag / h, 2) +
                      9 * std::pow(4 * alpha / (h * h), 2));

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_T_i      = scratch_data.phi_T[q][i];
          const auto grad_phi_T_i = scratch_data.grad_phi_T[q][i];


          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const Tensor<1, dim> grad_phi_T_j = scratch_data.grad_phi_T[q][j];


              // Weak form for : - k * laplacian T + rho * cp *
              //                  u * gradT - f -
              //                  tau:grad(u) =0
              // Hypothesis : incompressible newtonian fluid
              // so tau:grad(u) =
              // mu*(grad(u)+transpose(grad(u)).transpose(grad(u))
              local_matrix(i, j) +=
                (thermal_conductivity[q] * grad_phi_T_i * grad_phi_T_j +
                 rho_cp * phi_T_i * velocity * grad_phi_T_j) *
                JxW;

              local_matrix(i, j) += tau * strong_jacobian_vec[q][j] *
                                    (grad_phi_T_i * velocity) * JxW;
            }
        }

    } // end loop on quadrature points
}

template <int dim>
void
HeatTransferAssemblerCore<dim>::assemble_rhs(
  HeatTransferScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{
  const auto         method = this->simulation_control->get_assembly_method();
  const unsigned int n_q_points = scratch_data.n_q_points;
  const double       h          = scratch_data.cell_size;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  const std::vector<double> &density       = scratch_data.density;
  const std::vector<double> &specific_heat = scratch_data.specific_heat;
  const std::vector<double> &thermal_conductivity =
    scratch_data.thermal_conductivity;


  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs           = copy_data.local_rhs;

  // assembling right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      const Tensor<1, dim> temperature_gradient =
        scratch_data.temperature_gradients[q];
      const double temperature_laplacian =
        scratch_data.present_temperature_laplacians[q];

      // Gather physical properties in case of mono fluids simulations (to be
      // modified by cell in case of multiple fluids simulations)
      double rho_cp = density[q] * specific_heat[q];
      double alpha  = thermal_conductivity[q] / rho_cp;

      // Store JxW in local variable for faster access
      const double JxW = scratch_data.fe_values_T.JxW(q);

      const auto velocity = scratch_data.velocity_values[q];

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step
      const double tau =
        is_steady(method) ?
          1. / std::sqrt(std::pow(2. * rho_cp * u_mag / h, 2) +
                         9 * std::pow(4 * alpha / (h * h), 2)) :
          1. /
            std::sqrt(std::pow(sdt, 2) + std::pow(2. * rho_cp * u_mag / h, 2) +
                      9 * std::pow(4 * alpha / (h * h), 2));

      // Calculate the strong residual for GLS stabilization
      strong_residual_vec[q] += rho_cp * velocity * temperature_gradient -
                                thermal_conductivity[q] * temperature_laplacian;


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_T_i      = scratch_data.phi_T[q][i];
          const auto grad_phi_T_i = scratch_data.grad_phi_T[q][i];

          // rhs for : - k * laplacian T + rho * cp * u * grad T - f
          // -grad(u)*grad(u) = 0
          local_rhs(i) -=
            (thermal_conductivity[q] * grad_phi_T_i * temperature_gradient +
             rho_cp * phi_T_i * velocity * temperature_gradient -
             scratch_data.source[q] * phi_T_i) *
            JxW;

          local_rhs(i) -=
            tau * (strong_residual_vec[q] * (grad_phi_T_i * velocity)) * JxW;
        }

    } // end loop on quadrature points
}

template class HeatTransferAssemblerCore<2>;
template class HeatTransferAssemblerCore<3>;

template <int dim>
void
HeatTransferAssemblerBDF<dim>::assemble_matrix(
  HeatTransferScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{
  // Loop and quadrature informations
  const auto &       JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  const std::vector<double> &density       = scratch_data.density;
  const std::vector<double> &specific_heat = scratch_data.specific_heat;
  const std::vector<double> &grad_specific_heat_temperature =
    scratch_data.grad_specific_heat_temperature;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &strong_jacobian = copy_data.strong_jacobian;
  auto &local_matrix    = copy_data.local_matrix;

  // Time stepping information
  const auto          method = this->simulation_control->get_assembly_method();
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Gather physical properties in case of mono fluids simulations (to be
  // modified by cell in case of multiple fluids simulations)

  const double h = scratch_data.cell_size;


  // Vector for the BDF coefficients
  Vector<double>      bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<double> temperature(1 + number_of_previous_solutions(method));


  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      const double rho    = density[q];
      const double rho_cp = rho * specific_heat[q];
      double       dT_dt  = 0;

      temperature[0] = scratch_data.present_temperature_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        temperature[p + 1] = scratch_data.previous_temperature_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += rho_cp * bdf_coefs[p] * temperature[p];
          dT_dt += bdf_coefs[p] * temperature[p];
        }

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian[q][j] +=
            rho_cp * bdf_coefs[0] * scratch_data.phi_T[q][j];
        }

      const double tau_ggls =
        std::pow(h, scratch_data.fe_values_T.get_fe().degree + 1) / 6. / rho_cp;

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double phi_T_i      = scratch_data.phi_T[q][i];
          const auto   grad_phi_T_i = scratch_data.grad_phi_T[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const double phi_T_j      = scratch_data.phi_T[q][j];
              const auto   grad_phi_T_j = scratch_data.grad_phi_T[q][j];


              local_matrix(i, j) +=
                rho_cp * phi_T_j * phi_T_i * bdf_coefs[0] * JxW[q];

              local_matrix(i, j) += rho * grad_specific_heat_temperature[q] *
                                    phi_T_j * phi_T_i * dT_dt * JxW[q];


              if (GGLS)
                {
                  local_matrix(i, j) += rho_cp * rho_cp * tau_ggls *
                                        (grad_phi_T_i * grad_phi_T_j) *
                                        bdf_coefs[0] * JxW[q];
                }
            }
        }
    } // end loop on quadrature points
}

template <int dim>
void
HeatTransferAssemblerBDF<dim>::assemble_rhs(
  HeatTransferScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{
  // Gather physical properties in case of mono fluids simulations (to be
  // modified by cell in case of multiple fluids simulations)
  const std::vector<double> &density       = scratch_data.density;
  const std::vector<double> &specific_heat = scratch_data.specific_heat;


  // Loop and quadrature informations
  const auto &       JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &local_rhs       = copy_data.local_rhs;

  // Time stepping information
  const auto          method = this->simulation_control->get_assembly_method();
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Vector for the BDF coefficients
  Vector<double>      bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<double> temperature(1 + number_of_previous_solutions(method));
  std::vector<Tensor<1, dim>> temperature_gradient(
    1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      const double rho_cp = density[q] * specific_heat[q];

      const double tau_ggls =
        std::pow(h, scratch_data.fe_values_T.get_fe().degree + 1) / 6. / rho_cp;


      temperature[0]          = scratch_data.present_temperature_values[q];
      temperature_gradient[0] = scratch_data.temperature_gradients[q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        {
          temperature[p + 1] = scratch_data.previous_temperature_values[p][q];
          temperature_gradient[p + 1] =
            scratch_data.previous_temperature_gradients[p][q];
        }

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += rho_cp * bdf_coefs[p] * temperature[p];
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double phi_T_i      = scratch_data.phi_T[q][i];
          const auto   grad_phi_T_i = scratch_data.grad_phi_T[q][i];
          double       local_rhs_i  = 0;
          for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
               ++p)
            {
              local_rhs_i -= rho_cp * bdf_coefs[p] * (temperature[p] * phi_T_i);

              if (GGLS)
                {
                  local_rhs_i -= rho_cp * rho_cp * tau_ggls * grad_phi_T_i *
                                 bdf_coefs[p] * temperature_gradient[p];
                }
            }
          local_rhs(i) += local_rhs_i * JxW[q];
        }
    } // end loop on quadrature points
}

template class HeatTransferAssemblerBDF<2>;
template class HeatTransferAssemblerBDF<3>;


template <int dim>
void
HeatTransferAssemblerRobinBC<dim>::assemble_matrix(
  HeatTransferScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;
  auto &local_matrix = copy_data.local_matrix;

  // Robin boundary condition, loop on faces (Newton's cooling law)
  // implementation similar to deal.ii step-7
  for (unsigned int i_bc = 0; i_bc < this->boundary_conditions_ht.size; ++i_bc)
    {
      if (this->boundary_conditions_ht.type[i_bc] ==
          BoundaryConditions::BoundaryType::convection)
        {
          const double h = this->boundary_conditions_ht.h[i_bc];
          for (unsigned int f = 0; f < scratch_data.n_faces; ++f)
            {
              if (scratch_data.boundary_face_id[f] ==
                  this->boundary_conditions_ht.id[i_bc])
                {
                  for (unsigned int q = 0; q < scratch_data.n_faces_q_points;
                       ++q)
                    {
                      const double JxW = scratch_data.face_JxW[f][q];
                      for (unsigned int i = 0; i < scratch_data.n_dofs; ++i)
                        {
                          const double phi_face_T_i =
                            scratch_data.phi_face_T[f][q][i];

                          for (unsigned int j = 0; j < scratch_data.n_dofs; ++j)
                            {
                              const double phi_face_T_j =
                                scratch_data.phi_face_T[f][q][j];
                              local_matrix(i, j) +=
                                phi_face_T_i * phi_face_T_j * h * JxW;
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
HeatTransferAssemblerRobinBC<dim>::assemble_rhs(
  HeatTransferScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  auto &local_rhs = copy_data.local_rhs;


  // Robin boundary condition, loop on faces (Newton's cooling law)
  // implementation similar to deal.ii step-7
  for (unsigned int i_bc = 0; i_bc < this->boundary_conditions_ht.size; ++i_bc)
    {
      if (this->boundary_conditions_ht.type[i_bc] ==
          BoundaryConditions::BoundaryType::convection)
        {
          const double h     = this->boundary_conditions_ht.h[i_bc];
          const double T_inf = this->boundary_conditions_ht.Tinf[i_bc];

          for (unsigned int f = 0; f < scratch_data.n_faces; ++f)
            {
              if (scratch_data.boundary_face_id[f] ==
                  this->boundary_conditions_ht.id[i_bc])
                {
                  for (unsigned int q = 0; q < scratch_data.n_faces_q_points;
                       ++q)
                    {
                      const double T_face =
                        scratch_data.temperature_face_value[f][q];
                      const double JxW = scratch_data.face_JxW[f][q];
                      for (unsigned int i = 0; i < scratch_data.n_dofs; ++i)
                        {
                          const double phi_face_T_i =
                            scratch_data.phi_face_T[f][q][i];
                          local_rhs(i) -=
                            phi_face_T_i * h * (T_face - T_inf) * JxW;
                        }
                    }
                }
            }
        }
    }
}

template class HeatTransferAssemblerRobinBC<2>;
template class HeatTransferAssemblerRobinBC<3>;


template <int dim>
void
HeatTransferAssemblerViscousDissipation<dim>::assemble_matrix(
  HeatTransferScratchData<dim> & /*scratch_data*/,
  StabilizedMethodsCopyData &
  /*copy_data*/)
{}

template <int dim>
void
HeatTransferAssemblerViscousDissipation<dim>::assemble_rhs(
  HeatTransferScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  const std::vector<double> &density   = scratch_data.density;
  const std::vector<double> &viscosity = scratch_data.viscosity;


  // Time steps and inverse time steps which is used for stabilization
  // constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Copy data elements
  auto &local_rhs = copy_data.local_rhs;

  // assembling right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather physical properties in case of mono fluids simulations (to be
      // modified by cell in case of multiple fluids simulations)
      const double dynamic_viscosity = viscosity[q] * density[q];

      // Store JxW in local variable for faster access
      const double JxW = scratch_data.fe_values_T.JxW(q);

      const auto velocity_gradient = scratch_data.velocity_gradient_values[q];


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_T_i = scratch_data.phi_T[q][i];

          local_rhs(i) -=
            (-dynamic_viscosity * phi_T_i *
             scalar_product(velocity_gradient + transpose(velocity_gradient),
                            transpose(velocity_gradient))) *
            JxW;
        }
    }
}

template class HeatTransferAssemblerViscousDissipation<2>;
template class HeatTransferAssemblerViscousDissipation<3>;
