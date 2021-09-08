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
  double viscosity            = physical_properties.viscosity;

  double dynamic_viscosity = viscosity * density;
  double rho_cp            = density * specific_heat;
  double alpha             = thermal_conductivity / rho_cp;

  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

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

  auto &local_matrix = copy_data.local_matrix;

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

      const auto method = this->simulation_control->get_assembly_method();

      // Store JxW in local variable for faster access
      const double JxW = JxW_vec[q];

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
      const double tau_ggls =
        std::pow(h, scratch_data.fe_values_T.get_fe().degree + 1) / 6. / rho_cp;

      // Gather the shape functions and their gradient
      for (unsigned int k : scratch_data.fe_values_T.dof_indices())
        {
          scratch_data.phi_T[k] = scratch_data.fe_values_T.shape_value(k, q);
          scratch_data.grad_phi_T[k] =
            scratch_data.fe_values_T.shape_grad(k, q);
          scratch_data.hess_phi_T[k] =
            scratch_data.fe_values_T.shape_hessian(k, q);

          scratch_data.laplacian_phi_T[k] = trace(scratch_data.hess_phi_T[k]);
        }



      for (const unsigned int i : scratch_data.fe_values_T.dof_indices())
        {
          const auto phi_T_i      = scratch_data.phi_T[i];
          const auto grad_phi_T_i = scratch_data.grad_phi_T[i];


          for (const unsigned int j : scratch_data.fe_values_T.dof_indices())
            {
              const auto phi_T_j           = scratch_data.phi_T[j];
              const auto grad_phi_T_j      = scratch_data.grad_phi_T[j];
              const auto laplacian_phi_T_j = scratch_data.laplacian_phi_T[j];


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

              auto strong_jacobian = rho_cp * velocity * grad_phi_T_j -
                                     thermal_conductivity * laplacian_phi_T_j;

              local_matrix(i, j) +=
                tau * strong_jacobian * (grad_phi_T_i * velocity) * JxW;
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
  const auto   method = this->simulation_control->get_assembly_method();
  const double h      = scratch_data.cell_size;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs           = copy_data.local_rhs;

  // assembling local matrix and right hand side
  for (const unsigned int q :
       scratch_data.fe_values_T.quadrature_point_indices())
    {
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
      const double tau_ggls =
        std::pow(h, scratch_data.fe_values_T.get_fe().degree + 1) / 6. / rho_cp;

      // Gather the shape functions and their gradient
      for (unsigned int k : scratch_data.fe_values_T.dof_indices())
        {
          scratch_data.phi_T[k] = scratch_data.fe_values_T.shape_value(k, q);
          scratch_data.grad_phi_T[k] =
            scratch_data.fe_values_T.shape_grad(k, q);
          scratch_data.hess_phi_T[k] =
            scratch_data.fe_values_T.shape_hessian(k, q);

          scratch_data.laplacian_phi_T[k] = trace(scratch_data.hess_phi_T[k]);
        }



      for (const unsigned int i : scratch_data.fe_values_T.dof_indices())
        {
          const auto phi_T_i      = scratch_data.phi_T[i];
          const auto grad_phi_T_i = scratch_data.grad_phi_T[i];

          // rhs for : - k * laplacian T + rho * cp * u * grad T - f
          // -grad(u)*grad(u) = 0
          local_rhs(i) -= (thermal_conductivity * grad_phi_T_i *
                             scratch_data.temperature_gradients[q] +
                           rho_cp * phi_T_i * velocity *
                             scratch_data.temperature_gradients[q] -
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
  const auto method = this->simulation_control->get_assembly_method();

  // Gather physical properties in case of mono fluids simulations (to be
  // modified by cell in case of multiple fluids simulations)
  double density              = physical_properties.density;
  double specific_heat        = physical_properties.specific_heat;
  double thermal_conductivity = physical_properties.thermal_conductivity;
  double viscosity            = physical_properties.viscosity;

  double dynamic_viscosity = viscosity * density;
  double rho_cp            = density * specific_heat;
  double alpha             = thermal_conductivity / rho_cp;

  const double h = scratch_data.cell_size;



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

  // Vector for the BDF coefficients
  Vector<double> bdf_coefs = bdf_coefficients(method, time_steps_vector);

  // Copy data elements
  auto &strong_jacobian = copy_data.strong_jacobian;
  auto &local_matrix    = copy_data.local_matrix;

  // assembling local matrix and right hand side
  for (const unsigned int q :
       scratch_data.fe_values_T.quadrature_point_indices())
    {
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

      const double tau_ggls =
        std::pow(h, scratch_data.fe_values_T.get_fe().degree + 1) / 6. / rho_cp;

      // Gather the shape functions and their gradient
      for (unsigned int k : scratch_data.fe_values_T.dof_indices())
        {
          scratch_data.phi_T[k] = scratch_data.fe_values_T.shape_value(k, q);
          scratch_data.grad_phi_T[k] =
            scratch_data.fe_values_T.shape_grad(k, q);
          scratch_data.hess_phi_T[k] =
            scratch_data.fe_values_T.shape_hessian(k, q);

          scratch_data.laplacian_phi_T[k] = trace(scratch_data.hess_phi_T[k]);
        }



      for (const unsigned int i : scratch_data.fe_values_T.dof_indices())
        {
          const auto phi_T_i      = scratch_data.phi_T[i];
          const auto grad_phi_T_i = scratch_data.grad_phi_T[i];


          for (const unsigned int j : scratch_data.fe_values_T.dof_indices())
            {
              const auto phi_T_j      = scratch_data.phi_T[j];
              const auto grad_phi_T_j = scratch_data.grad_phi_T[j];

              local_matrix(i, j) +=
                rho_cp * phi_T_j * phi_T_i * bdf_coefs[0] * JxW;

              strong_jacobian[q][j] += rho_cp * phi_T_j * bdf_coefs[0];

              if (GGLS)
                {
                  local_matrix(i, j) += rho_cp * rho_cp * tau_ggls *
                                        (grad_phi_T_i * grad_phi_T_j) *
                                        bdf_coefs[0] * JxW;
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
  const auto   method = this->simulation_control->get_assembly_method();
  const double h      = scratch_data.cell_size;



  std::vector<double> time_steps_vector =
    simulation_control->get_time_steps_vector();
  // Time steps and inverse time steps which is used for numerous calculations
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;


  // Vector for the BDF coefficients
  Vector<double> bdf_coefs = bdf_coefficients(method, time_steps_vector);

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs           = copy_data.local_rhs;


  // assembling local matrix and right hand side
  for (const unsigned int q :
       scratch_data.fe_values_T.quadrature_point_indices())
    {
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
      const double tau_ggls =
        std::pow(h, scratch_data.fe_values_T.get_fe().degree + 1) / 6. / rho_cp;

      // Gather the shape functions and their gradient
      for (unsigned int k : scratch_data.fe_values_T.dof_indices())
        {
          scratch_data.phi_T[k] = scratch_data.fe_values_T.shape_value(k, q);
          scratch_data.grad_phi_T[k] =
            scratch_data.fe_values_T.shape_grad(k, q);
          scratch_data.hess_phi_T[k] =
            scratch_data.fe_values_T.shape_hessian(k, q);

          scratch_data.laplacian_phi_T[k] = trace(scratch_data.hess_phi_T[k]);
        }

      for (const unsigned int i : scratch_data.fe_values_T.dof_indices())
        {
          const auto phi_T_i      = scratch_data.phi_T[i];
          const auto grad_phi_T_i = scratch_data.grad_phi_T[i];


          // Calculate the strong residual for GLS stabilization
          auto strong_residual =
            rho_cp * scratch_data.velocity_values[q] *
              scratch_data.temperature_gradients[q] -
            thermal_conductivity *
              scratch_data.present_temperature_laplacians[q];


          // Residual associated with BDF schemes
          for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
               ++p)
            {
              local_rhs(i) -= rho_cp * bdf_coefs[p] *
                              scratch_data.present_temperature_values[p] *
                              phi_T_i * JxW;

              strong_residual += rho_cp * bdf_coefs[p] *
                                 scratch_data.present_temperature_values[p];

              if (GGLS)
                {
                  local_rhs(i) -= rho_cp * rho_cp * tau_ggls * grad_phi_T_i *
                                  bdf_coefs[p] *
                                  scratch_data.temperature_gradients[p] * JxW;
                }
            }
        }
    } // end loop on quadrature points
}

template class HeatTransferAssemblerBDF<2>;
template class HeatTransferAssemblerBDF<3>;
