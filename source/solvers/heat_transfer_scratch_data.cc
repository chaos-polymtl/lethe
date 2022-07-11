#include <core/bdf.h>
#include <core/sdirk.h>
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
  this->JxW = std::vector<double>(n_q_points);

  // Forcing term array
  this->source = std::vector<double>(n_q_points);

  // Initialize arrays related to velocity and pressure
  this->velocities.first_vector_component = 0;
  // Velocity
  this->velocity_values          = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_gradient_values = std::vector<Tensor<2, dim>>(n_q_points);
  this->shear_rate_values        = std::vector<double>(n_q_points);


  // Initialize arrays related to temperature
  this->present_temperature_values = std::vector<double>(n_q_points);
  this->temperature_gradients      = std::vector<Tensor<1, dim>>(n_q_points);
  this->present_temperature_laplacians = std::vector<double>(n_q_points);

  // Initial arrays for physical properties
  this->specific_heat                  = std::vector<double>(n_q_points);
  this->grad_specific_heat_temperature = std::vector<double>(n_q_points);
  this->thermal_conductivity           = std::vector<double>(n_q_points);
  this->density                        = std::vector<double>(n_q_points);
  this->viscosity                      = std::vector<double>(n_q_points);

  // Velocity for BDF schemes
  this->previous_temperature_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));
  this->previous_temperature_gradients =
    std::vector<std::vector<Tensor<1, dim>>>(
      maximum_number_of_previous_solutions(),
      std::vector<Tensor<1, dim>>(n_q_points));

  // Velocity for SDIRK schemes
  this->stages_temperature_values =
    std::vector<std::vector<double>>(max_number_of_intermediary_stages(),
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
    std::pair<field, std::vector<double>>(field::temperature, n_q_points));
  fields.insert(
    std::pair<field, std::vector<double>>(field::previous_temperature,
                                          n_q_points));
  fields.insert(
    std::pair<field, std::vector<double>>(field::shear_rate, n_q_points));
  specific_heat        = std::vector<double>(n_q_points);
  density              = std::vector<double>(n_q_points);
  thermal_conductivity = std::vector<double>(n_q_points);
  viscosity            = std::vector<double>(n_q_points);
}

template <int dim>
void
HeatTransferScratchData<dim>::enable_vof(const FiniteElement<dim> &fe,
                                         const Quadrature<dim> &   quadrature,
                                         const Mapping<dim> &      mapping)
{
  gather_vof    = true;
  fe_values_vof = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Allocate VOF values
  phase_values = std::vector<double>(this->n_q_points);

  // Allocate physical properties
  specific_heat_0                  = std::vector<double>(n_q_points);
  density_0                        = std::vector<double>(n_q_points);
  thermal_conductivity_0           = std::vector<double>(n_q_points);
  viscosity_0                      = std::vector<double>(n_q_points);
  grad_specific_heat_temperature_0 = std::vector<double>(n_q_points);

  specific_heat_1                  = std::vector<double>(n_q_points);
  density_1                        = std::vector<double>(n_q_points);
  thermal_conductivity_1           = std::vector<double>(n_q_points);
  viscosity_1                      = std::vector<double>(n_q_points);
  grad_specific_heat_temperature_1 = std::vector<double>(n_q_points);

  viscous_dissipation_coefficient = std::vector<double>(n_q_points);
}


template <int dim>
void
HeatTransferScratchData<dim>::calculate_physical_properties()
{
  set_field_vector(field::temperature,
                   this->present_temperature_values,
                   this->fields);

  if (properties_manager.field_is_required(field::previous_temperature))
    set_field_vector(field::previous_temperature,
                     this->previous_temperature_values[0],
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
          specific_heat_model->vector_jacobian(fields,
                                               field::temperature,
                                               grad_specific_heat_temperature);
          thermal_conductivity_model->vector_value(fields,
                                                   thermal_conductivity);
          rheology_model->vector_value(fields, viscosity);

          break;
        }
      case 2:
        {
          const auto density_models = properties_manager.get_density_vector();
          const auto specific_heat_models =
            properties_manager.get_specific_heat_vector();
          const auto thermal_conductivity_models =
            properties_manager.get_thermal_conductivity_vector();
          const auto rheology_models = properties_manager.get_rheology_vector();

          density_models[0]->vector_value(fields, density_0);
          specific_heat_models[0]->vector_value(fields, specific_heat_0);
          thermal_conductivity_models[0]->vector_value(fields,
                                                       thermal_conductivity_0);
          rheology_models[0]->vector_value(fields, viscosity_0);
          specific_heat_models[0]->vector_jacobian(
            fields, field::temperature, grad_specific_heat_temperature_0);


          density_models[1]->vector_value(fields, density_1);
          specific_heat_models[1]->vector_value(fields, specific_heat_1);
          thermal_conductivity_models[1]->vector_value(fields,
                                                       thermal_conductivity_1);
          rheology_models[1]->vector_value(fields, viscosity_1);
          specific_heat_models[1]->vector_jacobian(
            fields, field::temperature, grad_specific_heat_temperature_1);

          // Blend the physical properties using the VOF field
          for (unsigned int q = 0; q < this->n_q_points; ++q)
            {
              density[q] = calculate_point_property(this->phase_values[q],
                                                    this->density_0[q],
                                                    this->density_1[q]);

              specific_heat[q] =
                calculate_point_property(this->phase_values[q],
                                         this->specific_heat_0[q],
                                         this->specific_heat_1[q]);

              thermal_conductivity[q] =
                calculate_point_property(this->phase_values[q],
                                         this->thermal_conductivity_0[q],
                                         this->thermal_conductivity_1[q]);

              viscosity[q] = calculate_point_property(this->phase_values[q],
                                                      this->viscosity_0[q],
                                                      this->viscosity_1[q]);

              grad_specific_heat_temperature[q] = calculate_point_property(
                this->phase_values[q],
                this->grad_specific_heat_temperature_0[q],
                this->grad_specific_heat_temperature_1[q]);

              // Coefficient used to neglect viscous dissipation in the fluid
              // which is more than 10 times less dense than the other fluid
              if (this->viscosity_1[q] > 10. * this->viscosity_0[q])
                {
                  // if phase = 0, no viscous dissipation
                  // if phase = 1, maximum viscous dissipation
                  this->viscous_dissipation_coefficient[q] =
                    this->phase_values[q];
                }
              else
                {
                  // if phase = 1, no viscous dissipation
                  // if phase = 0, maximum viscous dissipation
                  this->viscous_dissipation_coefficient[q] =
                    1 - this->phase_values[q];
                }
            }
          break;
        }
      default:
        throw std::runtime_error("Unsupported number of fluids (>2)");
    }
}


template class HeatTransferScratchData<2>;
template class HeatTransferScratchData<3>;
