// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/utilities.h>

#include <solvers/rans_turbulence_scratch_data.h>


template <int dim>
void
RANSTurbulenceScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_rans.get_quadrature().size();
  this->n_dofs     = fe_values_rans.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW = std::vector<double>(n_q_points);

  // Forcing term array
  this->source = std::vector<double>(n_q_points);

  // Initialize arrays related to velocity and pressure
  // Velocity
  this->velocity_values          = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_gradient_values = std::vector<Tensor<2, dim>>(n_q_points);
  this->shear_rate_values        = std::vector<double>(n_q_points);

  // Initial arrays for physical properties
  this->turbulent_viscosity              = std::vector<double>(n_q_points);

  // Velocity for BDF schemes
  this->previous_turbulent_viscosity =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));

  // Initialize arrays related to shape functions
  // Velocity shape functions
  this->phi_T =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));
  this->grad_phi_T = std::vector<std::vector<Tensor<1, dim>>>(
    n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
  this->hess_phi_T = std::vector<std::vector<Tensor<2, dim>>>(
    n_q_points, std::vector<Tensor<2, dim>>(n_dofs));
  this->laplacian_phi_T =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));

  // Physical properties
  fields.insert(
    std::pair<field, std::vector<double>>(field::shear_rate, n_q_points));
  fields.insert(
    std::pair<field, std::vector<double>>(field::turbulent_viscosity, n_q_points));
}

template <int dim>
void
RANSTurbulenceScratchData<dim>::calculate_physical_properties()
{
  set_field_vector(field::temperature,
                   this->present_temperature_values,
                   this->fields);

  if (properties_manager.field_is_required(field::temperature_p1))
    set_field_vector(field::temperature_p1,
                     this->previous_temperature_values[0],
                     this->fields);

  if (properties_manager.field_is_required(field::temperature_p2))
    set_field_vector(field::temperature_p2,
                     this->previous_temperature_values[1],
                     this->fields);

  if (properties_manager.field_is_required(field::shear_rate))
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // Calculate shear rate (at each q)
          const Tensor<2, dim> shear_rate_tensor =
            velocity_gradient_values[q] +
            transpose(velocity_gradient_values[q]);

          // Calculate the shear rate magnitude
          shear_rate_values[q] =
            calculate_shear_rate_magnitude(shear_rate_tensor);
        }

      set_field_vector(field::shear_rate, shear_rate_values, this->fields);
    }

  if (properties_manager.field_is_required(field::pressure))
    {
      set_field_vector(field::pressure, this->pressure_values, this->fields);
    }

  if (material_id < 1 || properties_manager.get_number_of_solids() < 1)
    {
      switch (properties_manager.get_number_of_fluids())
        {
          case 1:
            {
              const auto density_model = properties_manager.get_density();
              const auto specific_heat_model =
                properties_manager.get_specific_heat();
              const auto thermal_conductivity_model =
                properties_manager.get_thermal_conductivity();
              const auto rheology_model = properties_manager.get_rheology();

              density_model->vector_value(fields, density);
              specific_heat_model->vector_value(fields, specific_heat);
              specific_heat_model->vector_jacobian(
                fields, field::temperature, grad_specific_heat_temperature);
              thermal_conductivity_model->vector_value(fields,
                                                       thermal_conductivity);
              rheology_model->get_dynamic_viscosity_vector(
                density_model->get_density_ref(), fields, dynamic_viscosity);

              break;
            }
            break;
          default:
            throw std::runtime_error("Unsupported number of fluids (>2)");
        }
    }
  else
    {
      const auto density_model = properties_manager.get_density(0, material_id);
      const auto specific_heat_model =
        properties_manager.get_specific_heat(0, material_id);
      const auto thermal_conductivity_model =
        properties_manager.get_thermal_conductivity(0, material_id);
      const auto rheology_model =
        properties_manager.get_rheology(0, material_id);

      density_model->vector_value(fields, density);
      specific_heat_model->vector_value(fields, specific_heat);
      specific_heat_model->vector_jacobian(fields,
                                           field::temperature,
                                           grad_specific_heat_temperature);
      thermal_conductivity_model->vector_value(fields, thermal_conductivity);
      rheology_model->get_dynamic_viscosity_vector(
        density_model->get_density_ref(), fields, dynamic_viscosity);
    }
}


template class RANSTurbulenceScratchData<2>;
template class RANSTurbulenceScratchData<3>;
