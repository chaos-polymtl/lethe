// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/utilities.h>

#include <solvers/heat_transfer_scratch_data.h>

template <int dim>
void
HeatTransferScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_T.get_quadrature().size();
  this->n_dofs     = fe_values_T.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW.reinit(n_q_points);

  // Forcing term array
  this->source = std::vector<double>(n_q_points);

  // Initialize arrays related to velocity and pressure
  // Velocity
  this->velocity_values          = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_gradient_values = std::vector<Tensor<2, dim>>(n_q_points);
  this->shear_rate_values        = std::vector<double>(n_q_points);
  // Pressure
  this->pressure_values = std::vector<double>(n_q_points);

  // Initialize arrays related to temperature
  this->present_temperature_values = std::vector<double>(n_q_points);
  this->temperature_gradients      = std::vector<Tensor<1, dim>>(n_q_points);
  this->present_temperature_laplacians = std::vector<double>(n_q_points);

  // Initial arrays for physical properties
  this->specific_heat                  = std::vector<double>(n_q_points);
  this->grad_specific_heat_temperature = std::vector<double>(n_q_points);
  this->thermal_conductivity           = std::vector<double>(n_q_points);
  this->density                        = std::vector<double>(n_q_points);
  this->dynamic_viscosity              = std::vector<double>(n_q_points);

  // Velocity for BDF schemes
  this->previous_temperature_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));
  this->previous_temperature_gradients =
    std::vector<std::vector<Tensor<1, dim>>>(
      maximum_number_of_previous_solutions(),
      std::vector<Tensor<1, dim>>(n_q_points));

  // Initialize arrays related to shape functions
  // Velocity shape functions
  this->phi_T.reinit(TableIndices<2>(n_q_points, n_dofs));
  this->grad_phi_T.reinit(TableIndices<2>(n_q_points, n_dofs));
  this->hess_phi_T.reinit(TableIndices<2>(n_q_points, n_dofs));
  this->laplacian_phi_T.reinit(TableIndices<2>(n_q_points, n_dofs));

  // Physical properties
  fields.insert(
    std::pair<field, std::vector<double>>(field::temperature, n_q_points));
  fields.insert(
    std::pair<field, std::vector<double>>(field::temperature_p1, n_q_points));
  fields.insert(
    std::pair<field, std::vector<double>>(field::temperature_p2, n_q_points));
  fields.insert(
    std::pair<field, std::vector<double>>(field::shear_rate, n_q_points));
  fields.insert(
    std::pair<field, std::vector<double>>(field::pressure, n_q_points));
}

template <int dim>
void
HeatTransferScratchData<dim>::enable_cls(
  const FiniteElement<dim>          &fe,
  const Quadrature<dim>             &quadrature,
  const Mapping<dim>                &mapping,
  const Parameters::CLS_PhaseFilter &phase_filter_parameters)
{
  gather_cls    = true;
  fe_values_cls = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Allocate CLS values
  filtered_phase_values = std::vector<double>(this->n_q_points);
  filtered_phase_gradient_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // Allocate physical properties
  specific_heat_0                  = std::vector<double>(n_q_points);
  density_0                        = std::vector<double>(n_q_points);
  thermal_conductivity_0           = std::vector<double>(n_q_points);
  dynamic_viscosity_0              = std::vector<double>(n_q_points);
  grad_specific_heat_temperature_0 = std::vector<double>(n_q_points);

  specific_heat_1                  = std::vector<double>(n_q_points);
  density_1                        = std::vector<double>(n_q_points);
  thermal_conductivity_1           = std::vector<double>(n_q_points);
  dynamic_viscosity_1              = std::vector<double>(n_q_points);
  grad_specific_heat_temperature_1 = std::vector<double>(n_q_points);

  // Create filter
  filter = ConservativeLevelSetFilterBase::model_cast(phase_filter_parameters);
}

template <int dim>
void
HeatTransferScratchData<dim>::enable_cls(
  const FiniteElement<dim>                              &fe,
  const Quadrature<dim>                                 &quadrature,
  const Mapping<dim>                                    &mapping,
  const std::shared_ptr<ConservativeLevelSetFilterBase> &filter)
{
  gather_cls    = true;
  fe_values_cls = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Allocate CLS values
  filtered_phase_values = std::vector<double>(this->n_q_points);
  filtered_phase_gradient_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // Allocate physical properties
  specific_heat_0                  = std::vector<double>(n_q_points);
  density_0                        = std::vector<double>(n_q_points);
  thermal_conductivity_0           = std::vector<double>(n_q_points);
  dynamic_viscosity_0              = std::vector<double>(n_q_points);
  grad_specific_heat_temperature_0 = std::vector<double>(n_q_points);

  specific_heat_1                  = std::vector<double>(n_q_points);
  density_1                        = std::vector<double>(n_q_points);
  thermal_conductivity_1           = std::vector<double>(n_q_points);
  dynamic_viscosity_1              = std::vector<double>(n_q_points);
  grad_specific_heat_temperature_1 = std::vector<double>(n_q_points);

  // Create filter
  this->filter = filter;
}

template <int dim>
void
HeatTransferScratchData<dim>::enable_time_harmonic_maxwell(
  const FiniteElement<dim> &fe,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_time_harmonic_maxwell = true;
  fe_values_time_harmonic_maxwell =
    std::make_shared<FEValues<dim>>(mapping, fe, quadrature, update_values);

  // Allocate time-harmonic Maxwell values
  this->electric_field_real_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);
  this->electric_field_imag_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);
  this->magnetic_field_real_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);
  this->magnetic_field_imag_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // Allocate physical properties
  this->electric_conductivity      = std::vector<double>(this->n_q_points);
  this->magnetic_permeability_imag = std::vector<double>(this->n_q_points);
  this->electric_permittivity_imag = std::vector<double>(this->n_q_points);

  // Define extractors for the time-harmonic Maxwell fields
  this->extractor_E_real = FEValuesExtractors::Vector(0);
  this->extractor_E_imag = FEValuesExtractors::Vector(dim);
  this->extractor_H_real = FEValuesExtractors::Vector(2 * dim);
  this->extractor_H_imag = FEValuesExtractors::Vector(3 * dim);
}

template <int dim>
void
HeatTransferScratchData<dim>::calculate_physical_properties()
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

              if (gather_time_harmonic_maxwell)
                {
                  const auto electric_conductivity_model =
                    properties_manager.get_electric_conductivity();
                  const auto magnetic_permeability_model =
                    properties_manager.get_magnetic_permeability_imag();
                  const auto electric_permittivity_model =
                    properties_manager.get_electric_permittivity_imag();

                  electric_conductivity_model->vector_value(
                    fields, electric_conductivity);
                  magnetic_permeability_model->vector_value(
                    fields, magnetic_permeability_imag);
                  electric_permittivity_model->vector_value(
                    fields, electric_permittivity_imag);
                }

              break;
            }
          case 2:
            {
              const auto density_models =
                properties_manager.get_density_vector();
              const auto specific_heat_models =
                properties_manager.get_specific_heat_vector();
              const auto thermal_conductivity_models =
                properties_manager.get_thermal_conductivity_vector();
              const auto rheology_models =
                properties_manager.get_rheology_vector();

              density_models[0]->vector_value(fields, density_0);
              specific_heat_models[0]->vector_value(fields, specific_heat_0);
              thermal_conductivity_models[0]->vector_value(
                fields, thermal_conductivity_0);
              rheology_models[0]->get_dynamic_viscosity_vector(
                density_models[0]->get_density_ref(),
                fields,
                dynamic_viscosity_0);
              specific_heat_models[0]->vector_jacobian(
                fields, field::temperature, grad_specific_heat_temperature_0);

              density_models[1]->vector_value(fields, density_1);
              specific_heat_models[1]->vector_value(fields, specific_heat_1);
              thermal_conductivity_models[1]->vector_value(
                fields, thermal_conductivity_1);
              rheology_models[1]->get_dynamic_viscosity_vector(
                density_models[1]->get_density_ref(),
                fields,
                dynamic_viscosity_1);
              specific_heat_models[1]->vector_jacobian(
                fields, field::temperature, grad_specific_heat_temperature_1);

              // Blend the physical properties using the CLS field
              for (unsigned int q = 0; q < this->n_q_points; ++q)
                {
                  auto filtered_phase_value = this->filtered_phase_values[q];

                  density[q] = calculate_point_property(filtered_phase_value,
                                                        this->density_0[q],
                                                        this->density_1[q]);

                  specific_heat[q] =
                    calculate_point_property(filtered_phase_value,
                                             this->specific_heat_0[q],
                                             this->specific_heat_1[q]);

                  thermal_conductivity[q] =
                    calculate_point_property(filtered_phase_value,
                                             this->thermal_conductivity_0[q],
                                             this->thermal_conductivity_1[q]);

                  dynamic_viscosity[q] =
                    calculate_point_property(filtered_phase_value,
                                             this->dynamic_viscosity_0[q],
                                             this->dynamic_viscosity_1[q]);

                  grad_specific_heat_temperature[q] = calculate_point_property(
                    filtered_phase_value,
                    this->grad_specific_heat_temperature_0[q],
                    this->grad_specific_heat_temperature_1[q]);
                }
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

      if (gather_time_harmonic_maxwell)
        {
          const auto electric_conductivity_model =
            properties_manager.get_electric_conductivity(0, material_id);
          const auto magnetic_permeability_model =
            properties_manager.get_magnetic_permeability_imag(0, material_id);
          const auto electric_permittivity_model =
            properties_manager.get_electric_permittivity_imag(0, material_id);

          electric_conductivity_model->vector_value(fields,
                                                    electric_conductivity);
          magnetic_permeability_model->vector_value(fields,
                                                    magnetic_permeability_imag);
          electric_permittivity_model->vector_value(fields,
                                                    electric_permittivity_imag);
        }
    }
}

template class HeatTransferScratchData<2>;
template class HeatTransferScratchData<3>;
