/*---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

#include <core/bdf.h>

#include <solvers/reactive_species_scratch_data.h>
template <int dim>
void
ReactiveSpeciesScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_reactive_species.get_quadrature().size();
  this->n_dofs     = fe_values_reactive_species.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW = std::vector<double>(n_q_points);

  // Forcing term array
  this->source = std::vector<std::vector<double>>(4);
  // Variables for the Reactive species equations
  this->reactive_species_values = std::vector<std::vector<double>>(4);
  this->reactive_species_gradients =
    std::vector<std::vector<Tensor<1, dim>>>(4);
  this->reactive_species_laplacians = std::vector<std::vector<double>>(4);
  // Reactive species values for BDF schemes
  this->previous_reactive_species_values =
    std::vector<std::vector<std::vector<double>>>(4);
  // Arrays related to shape functions
  this->phi      = std::vector<std::vector<std::vector<double>>>(4);
  this->grad_phi = std::vector<std::vector<std::vector<Tensor<1, dim>>>>(4);
  this->hess_phi = std::vector<std::vector<std::vector<Tensor<2, dim>>>>(4);
  this->laplacian_phi = std::vector<std::vector<std::vector<double>>>(4);
  // TODO Make flexible number of species
  for (unsigned int i = 0; i < 4; i++)
    {
      this->source[i] = std::vector<double>(n_q_points);
      // Initialize arrays related to phase order and chemical potential
      this->fe_values_extractors[i].component = i;

      this->reactive_species_values[i] = std::vector<double>(n_q_points);
      this->reactive_species_gradients[i] =
        std::vector<Tensor<1, dim>>(n_q_points);
      this->reactive_species_laplacians[i] = std::vector<double>(n_q_points);

      this->previous_reactive_species_values[i] =
        std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                         std::vector<double>(n_q_points));

      this->phi[i] =
        std::vector<std::vector<double>>(n_q_points,
                                         std::vector<double>(n_dofs));
      this->grad_phi[i] = std::vector<std::vector<Tensor<1, dim>>>(
        n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
      this->hess_phi[i] = std::vector<std::vector<Tensor<2, dim>>>(
        n_q_points, std::vector<Tensor<2, dim>>(n_dofs));
      this->laplacian_phi[i] =
        std::vector<std::vector<double>>(n_q_points,
                                         std::vector<double>(n_dofs));
    }

  // Velocity values
  this->velocity_values = std::vector<Tensor<1, dim>>(this->n_q_points);
}

template <int dim>
void
ReactiveSpeciesScratchData<dim>::calculate_physical_properties()
{
  switch (properties_manager.get_number_of_fluids())
    {
      case 1:
        {
          // TODO We need to compute physical properties
          /*throw std::runtime_error("Unsupported number of fluids (<2)");*/
        }
      case 2:
        {
          /*
            set_field_vector(field::phase_order_reactive_species,
                           this->phase_order_values,
                           this->fields);

          // Gather properties from material interactions
          const auto material_interaction_id =
            properties_manager.get_material_interaction_id(
              material_interactions_type::fluid_fluid, 0, 1);

          // Gather surface tension
          const auto surface_tension_model =
            properties_manager.get_surface_tension(material_interaction_id);
          surface_tension_model->vector_value(fields, surface_tension);

          // Gather mobility
          const auto mobility_reactive_species_model =
            properties_manager.get_mobility_reactive_species(
              material_interaction_id);
          mobility_reactive_species_model->vector_value(
            fields, mobility_reactive_species);
          mobility_reactive_species_model->vector_jacobian(
            fields,
            phase_order_reactive_species,
            mobility_reactive_species_gradient);
*/
          break;
        }
      default:
        throw std::runtime_error("Unsupported number of fluids (>2)");
    }
}

template class ReactiveSpeciesScratchData<2>;
template class ReactiveSpeciesScratchData<3>;
