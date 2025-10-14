// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>

#include <solvers/tracer_scratch_data.h>

template <int dim>
void
TracerScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_tracer.get_quadrature().size();
  this->n_dofs     = fe_values_tracer.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW = std::vector<double>(n_q_points);

  // Forcing term array
  this->source = std::vector<double>(n_q_points);

  // Initialize arrays related to velocity
  this->velocities.first_vector_component = 0;
  // Velocity
  this->velocity_values      = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_divergences = std::vector<double>(n_q_points);

  // Tracer
  this->tracer_values    = std::vector<double>(n_q_points);
  this->tracer_gradients = std::vector<Tensor<1, dim>>(n_q_points);
  this->previous_tracer_gradients =
    std::vector<Tensor<1, dim>>(this->n_q_points);
  this->tracer_laplacians              = std::vector<double>(n_q_points);
  this->tracer_diffusivity             = std::vector<double>(n_q_points);
  this->tracer_diffusivity_face        = std::vector<double>(n_q_points);
  this->tracer_diffusivity_0           = std::vector<double>(n_q_points);
  this->tracer_diffusivity_1           = std::vector<double>(n_q_points);
  this->tracer_reaction_prefactor      = std::vector<double>(n_q_points);
  this->grad_tracer_reaction_prefactor = std::vector<double>(n_q_points);

  // Solid signed distance function
  if (properties_manager.field_is_required(field::levelset))
    this->sdf_values = std::vector<double>(n_q_points);

  // Velocity for BDF schemes
  this->previous_tracer_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));

  // Initialize arrays related to shape functions
  // Velocity shape functions
  this->phi =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));
  this->grad_phi = std::vector<std::vector<Tensor<1, dim>>>(
    n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
  this->hess_phi = std::vector<std::vector<Tensor<2, dim>>>(
    n_q_points, std::vector<Tensor<2, dim>>(n_dofs));
  this->laplacian_phi =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));

  if (properties_manager.field_is_required(field::levelset))
    fields.insert(
      std::pair<field, std::vector<double>>(field::levelset, n_q_points));
  if (properties_manager.field_is_required(field::tracer_concentration))
    fields.insert(
      std::pair<field, std::vector<double>>(field::tracer_concentration,
                                            n_q_points));
}


template <int dim>
void
TracerScratchData<dim>::calculate_physical_properties()
{
  if (properties_manager.field_is_required(field::levelset))
    set_field_vector(field::levelset, this->sdf_values, this->fields);

  if (properties_manager.field_is_required(field::tracer_concentration))
    set_field_vector(field::tracer_concentration,
                     this->tracer_values,
                     this->fields);

  switch (properties_manager.get_number_of_fluids())
    {
      case 1:
        {
          const auto diffusivity_model =
            properties_manager.get_tracer_diffusivity();
          diffusivity_model->vector_value(fields, tracer_diffusivity);
          const auto reaction_prefactor_model =
            properties_manager.get_tracer_reaction_prefactor();
          reaction_prefactor_model->vector_value(fields,
                                                 tracer_reaction_prefactor);
          reaction_prefactor_model->vector_jacobian(
            fields,
            field::tracer_concentration,
            grad_tracer_reaction_prefactor);

          break;
        }
      case 2:
        {
          throw std::runtime_error(
            "Unsupported number of fluids (2) with the tracer physics.");
          // In this case, we need both density and viscosity
          const auto diffusivity_models =
            properties_manager.get_tracer_diffusivity_vector();

          diffusivity_models[0]->vector_value(fields, tracer_diffusivity_0);
          diffusivity_models[1]->vector_value(fields, tracer_diffusivity_1);

          // TODO Incomplete at the present time because the tracer VOF
          // complete is not finished Blend the physical properties
          // using the VOF field
          // for (unsigned int q = 0; q < this->n_q_points; ++q)
          // {
          //          tracer_diffusivity[q] =
          //            calculate_point_property(this->phase_values[q],
          //                                     this->density_0[q],
          //                                     this->density_1[q]);
          // }
          break;
        }
      default:
        throw std::runtime_error("Unsupported number of fluids (>2)");
    }
}


template class TracerScratchData<2>;
template class TracerScratchData<3>;
