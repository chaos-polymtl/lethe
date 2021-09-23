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
  double density              = physical_properties.density;
  double specific_heat        = physical_properties.specific_heat;
  double thermal_conductivity = physical_properties.thermal_conductivity;
  double rho_cp               = density * specific_heat;
  double alpha                = thermal_conductivity / rho_cp;

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
    simulation_control->get_time_steps_vector();

  // Time steps and inverse time steps which is used for numerous calculations
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Copy data elements
  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix        = copy_data.local_matrix;

  // assembling local matrix and right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      if (simulation_parameters.multiphysics.free_surface)
        {
          // Calculation of the equivalent physical properties at the
          // quadrature point
          density =
            calculate_point_property(scratch_data.phase_values[q],
                                     physical_properties.fluids[0].density,
                                     physical_properties.fluids[1].density);

          specific_heat = calculate_point_property(
            scratch_data.phase_values[q],
            physical_properties.fluids[0].specific_heat,
            physical_properties.fluids[1].specific_heat);

          thermal_conductivity = calculate_point_property(
            scratch_data.phase_values[q],
            physical_properties.fluids[0].thermal_conductivity,
            physical_properties.fluids[1].thermal_conductivity);

          // Useful definitions
          rho_cp = density * specific_heat;
          alpha  = thermal_conductivity / rho_cp;
        }

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


      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const Tensor<1, dim> grad_phi_T_j = scratch_data.grad_phi_T[q][j];
          const double laplacian_phi_T_j = scratch_data.laplacian_phi_T[q][j];
          strong_jacobian_vec[q][j] += rho_cp * velocity * grad_phi_T_j -
                                       thermal_conductivity * laplacian_phi_T_j;
        }


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
                (thermal_conductivity * grad_phi_T_i * grad_phi_T_j +
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
  std::cout << "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh " << std::endl;

  const auto         method = this->simulation_control->get_assembly_method();
  const unsigned int n_q_points = scratch_data.n_q_points;
  const double       h          = scratch_data.cell_size;
  const unsigned int n_dofs     = scratch_data.n_dofs;


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
      double density              = physical_properties.density;
      double specific_heat        = physical_properties.specific_heat;
      double thermal_conductivity = physical_properties.thermal_conductivity;
      double viscosity            = physical_properties.viscosity;

      double dynamic_viscosity = viscosity * density;
      double rho_cp            = density * specific_heat;
      double alpha             = thermal_conductivity / rho_cp;

      if (this->simulation_parameters.multiphysics.free_surface)
        {
          // Calculation of the equivalent physical properties at the
          // quadrature point
          density =
            calculate_point_property(scratch_data.phase_values[q],
                                     physical_properties.fluids[0].density,
                                     physical_properties.fluids[1].density);

          viscosity =
            calculate_point_property(scratch_data.phase_values[q],
                                     physical_properties.fluids[0].viscosity,
                                     physical_properties.fluids[1].viscosity);

          specific_heat = calculate_point_property(
            scratch_data.phase_values[q],
            physical_properties.fluids[0].specific_heat,
            physical_properties.fluids[1].specific_heat);

          thermal_conductivity = calculate_point_property(
            scratch_data.phase_values[q],
            physical_properties.fluids[0].thermal_conductivity,
            physical_properties.fluids[1].thermal_conductivity);

          // Useful definitions
          dynamic_viscosity = viscosity * density;
          rho_cp            = density * specific_heat;
          alpha             = thermal_conductivity / rho_cp;
        }

      // Store JxW in local variable for faster access
      const double JxW = scratch_data.fe_values_T.JxW(q);

      const auto velocity          = scratch_data.velocity_values[q];
      const auto velocity_gradient = scratch_data.velocity_gradient_values[q];

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
                                thermal_conductivity * temperature_laplacian;


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_T_i      = scratch_data.phi_T[q][i];
          const auto grad_phi_T_i = scratch_data.grad_phi_T[q][i];

          // rhs for : - k * laplacian T + rho * cp * u * grad T - f
          // -grad(u)*grad(u) = 0
          local_rhs(i) -=
            (thermal_conductivity * grad_phi_T_i * temperature_gradient +
             rho_cp * phi_T_i * velocity * temperature_gradient -
             scratch_data.source[q] * phi_T_i) *
            JxW;

          if (this->simulation_parameters.multiphysics.viscous_dissipation)
            {
              local_rhs(i) -= (-dynamic_viscosity * phi_T_i *
                               scalar_product(velocity_gradient +
                                                transpose(velocity_gradient),
                                              transpose(velocity_gradient))) *
                              JxW;
            }

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

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &strong_jacobian = copy_data.strong_jacobian;
  auto &local_matrix    = copy_data.local_matrix;

  // Time stepping information
  const auto          method = this->simulation_control->get_assembly_method();
  std::vector<double> time_steps_vector =
    simulation_control->get_time_steps_vector();

  // Gather physical properties in case of mono fluids simulations (to be
  // modified by cell in case of multiple fluids simulations)
  double density       = physical_properties.density;
  double specific_heat = physical_properties.specific_heat;
  double rho_cp        = density * specific_heat;

  const double h = scratch_data.cell_size;


  // Vector for the BDF coefficients
  Vector<double>      bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<double> temperature(1 + number_of_previous_solutions(method));


  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      if (this->simulation_parameters.multiphysics.free_surface)
        {
          // Calculation of the equivalent physical properties at the
          // quadrature point
          density =
            calculate_point_property(scratch_data.phase_values[q],
                                     physical_properties.fluids[0].density,
                                     physical_properties.fluids[1].density);

          specific_heat = calculate_point_property(
            scratch_data.phase_values[q],
            physical_properties.fluids[0].specific_heat,
            physical_properties.fluids[1].specific_heat);

          // Useful definitions
          rho_cp = density * specific_heat;
        }


      temperature[0] = scratch_data.present_temperature_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        temperature[p + 1] = scratch_data.previous_temperature_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += rho_cp * bdf_coefs[p] * temperature[p];
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
  double density       = physical_properties.density;
  double specific_heat = physical_properties.specific_heat;
  double rho_cp        = density * specific_heat;


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
    simulation_control->get_time_steps_vector();

  // Vector for the BDF coefficients
  Vector<double>      bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<double> temperature(1 + number_of_previous_solutions(method));
  std::vector<Tensor<1, dim>> temperature_gradient(
    1 + number_of_previous_solutions(method));


  const double tau_ggls =
    std::pow(h, scratch_data.fe_values_T.get_fe().degree + 1) / 6. / rho_cp;

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      if (this->simulation_parameters.multiphysics.free_surface)
        {
          // Calculation of the equivalent physical properties at the
          // quadrature point
          density =
            calculate_point_property(scratch_data.phase_values[q],
                                     physical_properties.fluids[0].density,
                                     physical_properties.fluids[1].density);

          specific_heat = calculate_point_property(
            scratch_data.phase_values[q],
            physical_properties.fluids[0].specific_heat,
            physical_properties.fluids[1].specific_heat);

          // Useful definitions
          rho_cp = density * specific_heat;
        }

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
HeatTransferAssemblerRBC<dim>::assemble_matrix(
  HeatTransferScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{
  auto &local_matrix = copy_data.local_matrix;

  std::vector<double> face_temperature_values =
    scratch_data.present_face_temperature_values;

  // Robin boundary condition, loop on faces (Newton's cooling law)
  // implementation similar to deal.ii step-7
  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions_ht.size;
       ++i_bc)
    {
      if (this->simulation_parameters.boundary_conditions_ht.type[i_bc] ==
          BoundaryConditions::BoundaryType::convection)
        {
          const double h =
            this->simulation_parameters.boundary_conditions_ht.h[i_bc];


          if (scratch_data.cell->is_locally_owned())
            {
              for (unsigned int face = 0;
                   face < GeometryInfo<dim>::faces_per_cell;
                   face++)
                {
                  if (scratch_data.cell->face(face)->at_boundary() &&
                      (scratch_data.cell->face(face)->boundary_id() ==
                       this->simulation_parameters.boundary_conditions_ht
                         .id[i_bc]))
                    {
                      scratch_data.fe_face_values_ht.reinit(scratch_data.cell,
                                                            face);
                      scratch_data.fe_face_values_ht.get_function_values(
                        scratch_data.evaluation_point, face_temperature_values);
                      {
                        for (const unsigned int q :
                             scratch_data.fe_face_values_ht
                               .quadrature_point_indices())
                          {
                            const double JxW =
                              scratch_data.fe_face_values_ht.JxW(q);
                            for (unsigned int k :
                                 scratch_data.fe_values_T.dof_indices())
                              scratch_data.phi_face_T[k] =
                                scratch_data.fe_face_values_ht.shape_value(k,
                                                                           q);

                            for (const unsigned int i :
                                 scratch_data.fe_values_T.dof_indices())
                              {
                                for (const unsigned int j :
                                     scratch_data.fe_values_T.dof_indices())
                                  {
                                    // Weak form modification
                                    local_matrix(i, j) +=
                                      scratch_data.phi_face_T[i] *
                                      scratch_data.phi_face_T[j] * h * JxW;
                                  }
                              }
                          }
                      }
                    }
                }
            }
        }
    } // end loop for Robin condition
}

template <int dim>
void
HeatTransferAssemblerRBC<dim>::assemble_rhs(
  HeatTransferScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{
  std::vector<double> face_temperature_values =
    scratch_data.present_face_temperature_values;

  auto &local_rhs = copy_data.local_rhs;

  // Robin boundary condition, loop on faces (Newton's cooling law)
  // implementation similar to deal.ii step-7
  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions_ht.size;
       ++i_bc)
    {
      if (this->simulation_parameters.boundary_conditions_ht.type[i_bc] ==
          BoundaryConditions::BoundaryType::convection)
        {
          const double h =
            this->simulation_parameters.boundary_conditions_ht.h[i_bc];
          const double T_inf =
            this->simulation_parameters.boundary_conditions_ht.Tinf[i_bc];


          if (scratch_data.cell->is_locally_owned())
            {
              for (unsigned int face = 0;
                   face < GeometryInfo<dim>::faces_per_cell;
                   face++)
                {
                  if (scratch_data.cell->face(face)->at_boundary() &&
                      (scratch_data.cell->face(face)->boundary_id() ==
                       this->simulation_parameters.boundary_conditions_ht
                         .id[i_bc]))
                    {
                      scratch_data.fe_face_values_ht.reinit(scratch_data.cell,
                                                            face);
                      scratch_data.fe_face_values_ht.get_function_values(
                        scratch_data.evaluation_point,
                        scratch_data.present_face_temperature_values);
                      {
                        for (const unsigned int q :
                             scratch_data.fe_face_values_ht
                               .quadrature_point_indices())
                          {
                            const double JxW =
                              scratch_data.fe_face_values_ht.JxW(q);
                            for (unsigned int k :
                                 scratch_data.fe_values_T.dof_indices())
                              scratch_data.phi_face_T[k] =
                                scratch_data.fe_face_values_ht.shape_value(k,
                                                                           q);

                            for (const unsigned int i :
                                 scratch_data.fe_values_T.dof_indices())
                              {
                                // Residual
                                local_rhs(i) -=
                                  scratch_data.phi_face_T[i] * h *
                                  (scratch_data
                                     .present_face_temperature_values[q] -
                                   T_inf) *
                                  JxW;
                              }
                          }
                      }
                    }
                }
            }
        }
    } // end loop for Robin condition
}

template class HeatTransferAssemblerRBC<2>;
template class HeatTransferAssemblerRBC<3>;
