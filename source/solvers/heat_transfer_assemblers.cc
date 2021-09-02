#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/copy_data.h>
#include <solvers/heat_transfer_assemblers.h>

template <int dim>
void
HeatTransferAssemblerCore<dim>::assemble_matrix(
  HeatTransferScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{
  // assembling local matrix and right hand side
  for (const unsigned int q : fe_values_ht.quadrature_point_indices())
    {
      if (this->simulation_parameters.multiphysics.free_surface)
        {
          // Calculation of the equivalent physical properties at the
          // quadrature point
          density =
            calculate_point_property(phase_values[q],
                                     physical_properties.fluids[0].density,
                                     physical_properties.fluids[1].density);

          viscosity =
            calculate_point_property(phase_values[q],
                                     physical_properties.fluids[0].viscosity,
                                     physical_properties.fluids[1].viscosity);

          specific_heat = calculate_point_property(
            phase_values[q],
            physical_properties.fluids[0].specific_heat,
            physical_properties.fluids[1].specific_heat);

          thermal_conductivity = calculate_point_property(
            phase_values[q],
            physical_properties.fluids[0].thermal_conductivity,
            physical_properties.fluids[1].thermal_conductivity);

          // Useful definitions
          dynamic_viscosity = viscosity * density;
          rho_cp            = density * specific_heat;
          alpha             = thermal_conductivity / rho_cp;
        }

      // Store JxW in local variable for faster access
      const double JxW = fe_values_ht.JxW(q);

      const auto velocity          = velocity_values[q];
      const auto velocity_gradient = velocity_gradient_values[q];


      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step
      const double tau =
        is_steady(time_stepping_method) ?
          1. / std::sqrt(std::pow(2. * rho_cp * u_mag / h, 2) +
                         9 * std::pow(4 * alpha / (h * h), 2)) :
          1. /
            std::sqrt(std::pow(sdt, 2) + std::pow(2. * rho_cp * u_mag / h, 2) +
                      9 * std::pow(4 * alpha / (h * h), 2));
      const double tau_ggls = std::pow(h, fe->degree + 1) / 6. / rho_cp;

      // Gather the shape functions and their gradient
      for (unsigned int k : fe_values_ht.dof_indices())
        {
          phi_T[k]      = fe_values_ht.shape_value(k, q);
          grad_phi_T[k] = fe_values_ht.shape_grad(k, q);
          hess_phi_T[k] = fe_values_ht.shape_hessian(k, q);

          laplacian_phi_T[k] = trace(hess_phi_T[k]);
        }



      for (const unsigned int i : fe_values_ht.dof_indices())
        {
          const auto phi_T_i      = phi_T[i];
          const auto grad_phi_T_i = grad_phi_T[i];


          if (assemble_matrix)
            {
              for (const unsigned int j : fe_values_ht.dof_indices())
                {
                  const auto phi_T_j           = phi_T[j];
                  const auto grad_phi_T_j      = grad_phi_T[j];
                  const auto laplacian_phi_T_j = laplacian_phi_T[j];


                  // Weak form for : - k * laplacian T + rho * cp *
                  //                  u * gradT - f -
                  //                  tau:grad(u) =0
                  // Hypothesis : incompressible newtonian fluid
                  // so tau:grad(u) =
                  // mu*(grad(u)+transpose(grad(u)).transpose(grad(u))
                  cell_matrix(i, j) +=
                    (thermal_conductivity * grad_phi_T_i * grad_phi_T_j +
                     rho_cp * phi_T_i * velocity * grad_phi_T_j) *
                    JxW;

                  auto strong_jacobian =
                    rho_cp * velocity * grad_phi_T_j -
                    thermal_conductivity * laplacian_phi_T_j;

                  // Mass matrix for transient simulation
                  if (is_bdf(time_stepping_method))
                    {
                      cell_matrix(i, j) +=
                        rho_cp * phi_T_j * phi_T_i * bdf_coefs[0] * JxW;

                      strong_jacobian += rho_cp * phi_T_j * bdf_coefs[0];

                      if (GGLS)
                        {
                          cell_matrix(i, j) += rho_cp * rho_cp * tau_ggls *
                                               (grad_phi_T_i * grad_phi_T_j) *
                                               bdf_coefs[0] * JxW;
                        }
                    }

                  cell_matrix(i, j) +=
                    tau * strong_jacobian * (grad_phi_T_i * velocity) * JxW;
                }
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
  // assembling local matrix and right hand side
  for (const unsigned int q : fe_values_ht.quadrature_point_indices())
    {
      if (this->simulation_parameters.multiphysics.free_surface)
        {
          // Calculation of the equivalent physical properties at the
          // quadrature point
          density =
            calculate_point_property(phase_values[q],
                                     physical_properties.fluids[0].density,
                                     physical_properties.fluids[1].density);

          viscosity =
            calculate_point_property(phase_values[q],
                                     physical_properties.fluids[0].viscosity,
                                     physical_properties.fluids[1].viscosity);

          specific_heat = calculate_point_property(
            phase_values[q],
            physical_properties.fluids[0].specific_heat,
            physical_properties.fluids[1].specific_heat);

          thermal_conductivity = calculate_point_property(
            phase_values[q],
            physical_properties.fluids[0].thermal_conductivity,
            physical_properties.fluids[1].thermal_conductivity);

          // Useful definitions
          dynamic_viscosity = viscosity * density;
          rho_cp            = density * specific_heat;
          alpha             = thermal_conductivity / rho_cp;
        }

      // Store JxW in local variable for faster access
      const double JxW = fe_values_ht.JxW(q);

      const auto velocity          = velocity_values[q];
      const auto velocity_gradient = velocity_gradient_values[q];


      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step
      const double tau =
        is_steady(time_stepping_method) ?
          1. / std::sqrt(std::pow(2. * rho_cp * u_mag / h, 2) +
                         9 * std::pow(4 * alpha / (h * h), 2)) :
          1. /
            std::sqrt(std::pow(sdt, 2) + std::pow(2. * rho_cp * u_mag / h, 2) +
                      9 * std::pow(4 * alpha / (h * h), 2));
      const double tau_ggls = std::pow(h, fe->degree + 1) / 6. / rho_cp;

      // Gather the shape functions and their gradient
      for (unsigned int k : fe_values_ht.dof_indices())
        {
          phi_T[k]      = fe_values_ht.shape_value(k, q);
          grad_phi_T[k] = fe_values_ht.shape_grad(k, q);
          hess_phi_T[k] = fe_values_ht.shape_hessian(k, q);

          laplacian_phi_T[k] = trace(hess_phi_T[k]);
        }



      for (const unsigned int i : fe_values_ht.dof_indices())
        {
          const auto phi_T_i      = phi_T[i];
          const auto grad_phi_T_i = grad_phi_T[i];

          // rhs for : - k * laplacian T + rho * cp * u * grad T - f
          // -grad(u)*grad(u) = 0
          cell_rhs(i) -=
            (thermal_conductivity * grad_phi_T_i * temperature_gradients[q] +
             rho_cp * phi_T_i * velocity * temperature_gradients[q] -
             source_term_values[q] * phi_T_i) *
            JxW;

          if (this->simulation_parameters.multiphysics.viscous_dissipation)
            {
              cell_rhs(i) -= (-dynamic_viscosity * phi_T_i *
                              scalar_product(velocity_gradient +
                                               transpose(velocity_gradient),
                                             transpose(velocity_gradient))) *
                             JxW;
            }

          // Calculate the strong residual for GLS stabilization
          auto strong_residual =
            rho_cp * velocity_values[q] * temperature_gradients[q] -
            thermal_conductivity * present_temperature_laplacians[q];



          // Residual associated with BDF schemes
          if (time_stepping_method ==
                Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
              time_stepping_method ==
                Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
            {
              cell_rhs(i) -= rho_cp *
                             (bdf_coefs[0] * present_temperature_values[q] +
                              bdf_coefs[1] * p1_temperature_values[q]) *
                             phi_T_i * JxW;

              strong_residual +=
                rho_cp * (bdf_coefs[0] * present_temperature_values[q] +
                          bdf_coefs[1] * p1_temperature_values[q]);

              if (GGLS)
                {
                  cell_rhs(i) -= rho_cp * rho_cp * tau_ggls * grad_phi_T_i *
                                 (bdf_coefs[0] * temperature_gradients[q] +
                                  bdf_coefs[1] * p1_temperature_gradients[q]) *
                                 JxW;
                }
            }

          if (time_stepping_method ==
              Parameters::SimulationControl::TimeSteppingMethod::bdf2)
            {
              cell_rhs(i) -= rho_cp *
                             (bdf_coefs[0] * present_temperature_values[q] +
                              bdf_coefs[1] * p1_temperature_values[q] +
                              bdf_coefs[2] * p2_temperature_values[q]) *
                             phi_T_i * JxW;

              strong_residual +=
                rho_cp * (bdf_coefs[0] * present_temperature_values[q] +
                          bdf_coefs[1] * p1_temperature_values[q] +
                          bdf_coefs[2] * p2_temperature_values[q]);

              if (GGLS)
                {
                  cell_rhs(i) -= rho_cp * rho_cp * tau_ggls * grad_phi_T_i *
                                 (bdf_coefs[0] * temperature_gradients[q] +
                                  bdf_coefs[1] * p1_temperature_gradients[q] +
                                  bdf_coefs[2] * p2_temperature_gradients[q]) *
                                 JxW;
                }
            }

          if (time_stepping_method ==
              Parameters::SimulationControl::TimeSteppingMethod::bdf3)
            {
              cell_rhs(i) -= rho_cp *
                             (bdf_coefs[0] * present_temperature_values[q] +
                              bdf_coefs[1] * p1_temperature_values[q] +
                              bdf_coefs[2] * p2_temperature_values[q] +
                              bdf_coefs[3] * p3_temperature_values[q]) *
                             phi_T_i * JxW;

              strong_residual +=
                rho_cp * (bdf_coefs[0] * present_temperature_values[q] +
                          bdf_coefs[1] * p1_temperature_values[q] +
                          bdf_coefs[2] * p2_temperature_values[q] +
                          bdf_coefs[3] * p3_temperature_values[q]);

              if (GGLS)
                {
                  cell_rhs(i) -= rho_cp * rho_cp * tau_ggls * grad_phi_T_i *
                                 (bdf_coefs[0] * temperature_gradients[q] +
                                  bdf_coefs[1] * p1_temperature_gradients[q] +
                                  bdf_coefs[2] * p2_temperature_gradients[q] +
                                  bdf_coefs[3] * p3_temperature_gradients[q]) *
                                 JxW;
                }
            }


          cell_rhs(i) -=
            tau * (strong_residual * (grad_phi_T_i * velocity_values[q])) * JxW;
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
{}

template <int dim>
void
HeatTransferAssemblerBDF<dim>::assemble_rhs(
  HeatTransferScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{}

template class HeatTransferAssemblerBDF<2>;
template class HeatTransferAssemblerBDF<3>;

template <int dim>
void
HeatTransferAssemblerRBC<dim>::assemble_matrix(
  HeatTransferScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{
  // Robin boundary condition, loop on faces (Newton's cooling law)
  // implementation similar to deal.ii step-7
  for (unsigned int i_bc = 0;
       i_bc < simulation_parameters.boundary_conditions_ht.size;
       ++i_bc)
    {
      if (this->simulation_parameters.boundary_conditions_ht.type[i_bc] ==
          BoundaryConditions::BoundaryType::convection)
        {
          const double h = simulation_parameters.boundary_conditions_ht.h[i_bc];
          const double T_inf =
            simulation_parameters.boundary_conditions_ht.Tinf[i_bc];
          std::vector<double> phi_face_T(dofs_per_cell);

          if (cell->is_locally_owned())
            {
              for (unsigned int face = 0;
                   face < GeometryInfo<dim>::faces_per_cell;
                   face++)
                {
                  if (cell->face(face)->at_boundary() &&
                      (cell->face(face)->boundary_id() ==
                       simulation_parameters.boundary_conditions_ht.id[i_bc]))
                    {
                      fe_face_values_ht.reinit(cell, face);
                      fe_face_values_ht.get_function_values(
                        evaluation_point, present_face_temperature_values);
                      {
                        for (const unsigned int q :
                             fe_face_values_ht.quadrature_point_indices())
                          {
                            const double JxW = fe_face_values_ht.JxW(q);
                            for (unsigned int k : fe_values_ht.dof_indices())
                              phi_face_T[k] =
                                fe_face_values_ht.shape_value(k, q);

                            for (const unsigned int i :
                                 fe_values_ht.dof_indices())
                              {
                                if (assemble_matrix)
                                  {
                                    for (const unsigned int j :
                                         fe_values_ht.dof_indices())
                                      {
                                        // Weak form modification
                                        cell_matrix(i, j) += phi_face_T[i] *
                                                             phi_face_T[j] * h *
                                                             JxW;
                                      }
                                  }
                                // Residual
                                cell_rhs(i) -=
                                  phi_face_T[i] * h *
                                  (present_face_temperature_values[q] - T_inf) *
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

template <int dim>
void
HeatTransferAssemblerRBC<dim>::assemble_rhs(
  HeatTransferScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{}

template class HeatTransferAssemblerRBC<2>;
template class HeatTransferAssemblerRBC<3>;
