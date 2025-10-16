// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>

#include <solvers/cahn_hilliard_scratch_data.h>
template <int dim>
void
CahnHilliardScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_cahn_hilliard.get_quadrature().size();
  this->n_dofs     = fe_values_cahn_hilliard.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW.reinit(n_q_points);

  // Forcing term array
  this->source_phase_order        = std::vector<double>(n_q_points);
  this->source_chemical_potential = std::vector<double>(n_q_points);

  // Initialize arrays related to phase order and chemical potential
  this->phase_order.component        = 0;
  this->chemical_potential.component = 1;

  // Variables for the Cahn-Hilliard equations
  this->phase_order_values            = std::vector<double>(n_q_points);
  this->phase_order_gradients         = std::vector<Tensor<1, dim>>(n_q_points);
  this->phase_order_laplacians        = std::vector<double>(n_q_points);
  this->chemical_potential_values     = std::vector<double>(n_q_points);
  this->chemical_potential_gradients  = std::vector<Tensor<1, dim>>(n_q_points);
  this->chemical_potential_laplacians = std::vector<double>(n_q_points);

  // Phase order for BDF schemes
  this->previous_phase_order_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));

  this->previous_chemical_potential_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));

  // Initialize arrays related to shape functions
  // Phase-order shape functions
  this->phi_phase.reinit(n_q_points, n_dofs);
  this->grad_phi_phase.reinit(n_q_points, n_dofs);
  this->hess_phi_phase.reinit(n_q_points, n_dofs);
  this->laplacian_phi_phase.reinit(n_q_points, n_dofs);

  // Chemical potential shape functions
  this->phi_potential.reinit(n_q_points, n_dofs);
  this->grad_phi_potential.reinit(n_q_points, n_dofs);
  this->hess_phi_potential.reinit(n_q_points, n_dofs);
  this->laplacian_phi_potential.reinit(n_q_points, n_dofs);

  // Velocity values
  this->velocity_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  this->previous_velocity_values = std::vector<std::vector<Tensor<1, dim>>>(
    maximum_number_of_previous_solutions(),
    std::vector<Tensor<1, dim>>(this->n_q_points));
  this->velocity_gradient_values =
    std::vector<Tensor<2, dim>>(this->n_q_points);

  // Allocate physical properties
  this->surface_tension                 = std::vector<double>(n_q_points);
  this->mobility_cahn_hilliard          = std::vector<double>(n_q_points);
  this->mobility_cahn_hilliard_gradient = std::vector<double>(n_q_points);

  // Physical properties
  fields.insert(
    std::pair<field, std::vector<double>>(field::phase_order_cahn_hilliard,
                                          n_q_points));
  fields.insert(std::pair<field, std::vector<double>>(
    field::phase_order_cahn_hilliard_filtered, n_q_points));
}

template <int dim>
void
CahnHilliardScratchData<dim>::calculate_physical_properties()
{
  switch (properties_manager.get_number_of_fluids())
    {
      case 1:
        {
          throw std::runtime_error("Unsupported number of fluids (<2)");
        }
      case 2:
        {
          set_field_vector(field::phase_order_cahn_hilliard,
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
          const auto mobility_cahn_hilliard_model =
            properties_manager.get_mobility_cahn_hilliard(
              material_interaction_id);
          mobility_cahn_hilliard_model->vector_value(fields,
                                                     mobility_cahn_hilliard);
          mobility_cahn_hilliard_model->vector_jacobian(
            fields, phase_order_cahn_hilliard, mobility_cahn_hilliard_gradient);

          break;
        }
      default:
        throw std::runtime_error("Unsupported number of fluids (>2)");
    }
}

template class CahnHilliardScratchData<2>;
template class CahnHilliardScratchData<3>;
