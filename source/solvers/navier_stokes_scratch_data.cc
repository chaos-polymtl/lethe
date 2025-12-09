// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/utilities.h>

#include <solvers/navier_stokes_scratch_data.h>

template <int dim>
void
NavierStokesScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values.get_quadrature().size();
  this->n_dofs     = fe_values.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW.reinit(n_q_points);

  // Initialize component array
  this->components = std::vector<unsigned int>(n_dofs);

  // Forcing term array
  this->rhs_force =
    std::vector<Vector<double>>(n_q_points, Vector<double>(dim + 1));
  this->force       = std::vector<Tensor<1, dim>>(n_q_points);
  this->mass_source = std::vector<double>(n_q_points);

  // Initialize arrays related to velocity and pressure
  this->velocities.first_vector_component = 0;
  this->pressure.component                = dim;

  // Velocity
  this->velocity_values            = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_divergences       = std::vector<double>(n_q_points);
  this->velocity_gradients         = std::vector<Tensor<2, dim>>(n_q_points);
  this->velocity_laplacians        = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_hessians          = std::vector<Tensor<3, dim>>(n_q_points);
  this->shear_rate                 = std::vector<double>(n_q_points);
  this->velocity_for_stabilization = std::vector<Tensor<1, dim>>(n_q_points);

  // For SDIRK method: sum(a_ij * k_j)
  if (this->simulation_control->is_sdirk())
    {
      this->sdirk_stage_sum = std::vector<Tensor<1, dim>>(n_q_points);
    }

  // Velocity for BDF schemes
  this->previous_velocity_values = std::vector<std::vector<Tensor<1, dim>>>(
    maximum_number_of_previous_solutions(),
    std::vector<Tensor<1, dim>>(n_q_points));

  // Pressure
  this->pressure_values         = std::vector<double>(n_q_points);
  this->pressure_gradients      = std::vector<Tensor<1, dim>>(n_q_points);
  this->pressure_scaling_factor = 1;

  // Pressure for BDF schemes (compressible Navier-Stokes)
  this->previous_pressure_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));
  // Initialize arrays related to shape functions
  // Velocity shape functions
  this->phi_u.reinit(n_q_points, n_dofs);
  this->grad_phi_u.reinit(n_q_points, n_dofs);
  this->div_phi_u.reinit(n_q_points, n_dofs);
  this->hess_phi_u.reinit(n_q_points, n_dofs);
  this->laplacian_phi_u.reinit(n_q_points, n_dofs);

  // Pressure shape functions
  this->phi_p.reinit(n_q_points, n_dofs);
  this->grad_phi_p.reinit(n_q_points, n_dofs);

  // Physical properties
  fields.insert(
    std::pair<field, std::vector<double>>(field::pressure, n_q_points));
  fields.insert(
    std::pair<field, std::vector<double>>(field::shear_rate, n_q_points));

  density                               = std::vector<double>(n_q_points);
  dynamic_viscosity                     = std::vector<double>(n_q_points);
  kinematic_viscosity                   = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization   = std::vector<double>(n_q_points);
  kinematic_viscosity_for_stabilization = std::vector<double>(n_q_points);
  thermal_expansion                     = std::vector<double>(n_q_points);
  grad_kinematic_viscosity_shear_rate   = std::vector<double>(n_q_points);

  previous_density =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_vof(
  const FiniteElement<dim>          &fe,
  const Quadrature<dim>             &quadrature,
  const Mapping<dim>                &mapping,
  const Parameters::VOF_PhaseFilter &phase_filter_parameters)
{
  gather_vof    = true;
  fe_values_vof = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Allocate VOF values
  phase_values          = std::vector<double>(this->n_q_points);
  filtered_phase_values = std::vector<double>(this->n_q_points);
  previous_phase_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
  // For STF calculation
  phase_gradient_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  filtered_phase_gradient_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // Allocate physical properties
  density_0                             = std::vector<double>(n_q_points);
  density_1                             = std::vector<double>(n_q_points);
  dynamic_viscosity_0                   = std::vector<double>(n_q_points);
  dynamic_viscosity_1                   = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_0 = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_1 = std::vector<double>(n_q_points);
  thermal_expansion_0                   = std::vector<double>(n_q_points);
  thermal_expansion_1                   = std::vector<double>(n_q_points);
  surface_tension                       = std::vector<double>(n_q_points);
  surface_tension_gradient              = std::vector<double>(n_q_points);
  compressibility_multiplier            = std::vector<double>(n_q_points);

  // Create filter
  filter = VolumeOfFluidFilterBase::model_cast(phase_filter_parameters);
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_vof(
  const FiniteElement<dim>                       &fe,
  const Quadrature<dim>                          &quadrature,
  const Mapping<dim>                             &mapping,
  const std::shared_ptr<VolumeOfFluidFilterBase> &filter)
{
  gather_vof    = true;
  fe_values_vof = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Allocate VOF values
  phase_values          = std::vector<double>(this->n_q_points);
  filtered_phase_values = std::vector<double>(this->n_q_points);
  previous_phase_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
  // For STF calculation
  phase_gradient_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  filtered_phase_gradient_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // Allocate physical properties
  density_0                             = std::vector<double>(n_q_points);
  density_1                             = std::vector<double>(n_q_points);
  dynamic_viscosity_0                   = std::vector<double>(n_q_points);
  dynamic_viscosity_1                   = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_0 = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_1 = std::vector<double>(n_q_points);
  thermal_expansion_0                   = std::vector<double>(n_q_points);
  thermal_expansion_1                   = std::vector<double>(n_q_points);
  surface_tension                       = std::vector<double>(n_q_points);
  surface_tension_gradient              = std::vector<double>(n_q_points);
  compressibility_multiplier            = std::vector<double>(n_q_points);

  // Create filter
  this->filter = filter;
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_cahn_hilliard(
  const FiniteElement<dim>       &fe,
  const Quadrature<dim>          &quadrature,
  const Mapping<dim>             &mapping,
  const Parameters::CahnHilliard &cahn_hilliard_parameters)
{
  gather_cahn_hilliard    = true;
  fe_values_cahn_hilliard = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Allocate CahnHilliard values
  phase_order_cahn_hilliard_values = std::vector<double>(this->n_q_points);
  chemical_potential_cahn_hilliard_values =
    std::vector<double>(this->n_q_points);

  // Allocate CahnHilliard gradients
  phase_order_cahn_hilliard_gradients =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // For STF calculation
  filtered_phase_order_cahn_hilliard_values =
    std::vector<double>(this->n_q_points);

  fields.insert(
    std::pair<field, std::vector<double>>(field::phase_order_cahn_hilliard,
                                          n_q_points));
  fields.insert(std::pair<field, std::vector<double>>(
    field::phase_order_cahn_hilliard_filtered, n_q_points));

  // Allocate physical properties
  density_0                             = std::vector<double>(n_q_points);
  density_1                             = std::vector<double>(n_q_points);
  dynamic_viscosity_0                   = std::vector<double>(n_q_points);
  dynamic_viscosity_1                   = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_0 = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_1 = std::vector<double>(n_q_points);
  thermal_expansion_0                   = std::vector<double>(n_q_points);
  thermal_expansion_1                   = std::vector<double>(n_q_points);
  surface_tension                       = std::vector<double>(n_q_points);
  surface_tension_gradient              = std::vector<double>(n_q_points);

  // Create filter
  cahn_hilliard_filter =
    CahnHilliardFilterBase::model_cast(cahn_hilliard_parameters);
}

template <int dim>
void
NavierStokesScratchData<dim>::enable_cahn_hilliard(
  const FiniteElement<dim>                      &fe,
  const Quadrature<dim>                         &quadrature,
  const Mapping<dim>                            &mapping,
  const std::shared_ptr<CahnHilliardFilterBase> &cahn_hilliard_filter)
{
  gather_cahn_hilliard    = true;
  fe_values_cahn_hilliard = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Allocate CahnHilliard values
  phase_order_cahn_hilliard_values = std::vector<double>(this->n_q_points);
  chemical_potential_cahn_hilliard_values =
    std::vector<double>(this->n_q_points);

  // Allocate CahnHilliard gradients
  phase_order_cahn_hilliard_gradients =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // For STF calculation
  filtered_phase_order_cahn_hilliard_values =
    std::vector<double>(this->n_q_points);

  fields.insert(
    std::pair<field, std::vector<double>>(field::phase_order_cahn_hilliard,
                                          n_q_points));

  fields.insert(std::pair<field, std::vector<double>>(
    field::phase_order_cahn_hilliard_filtered, n_q_points));

  // Allocate physical properties
  density_0                             = std::vector<double>(n_q_points);
  density_1                             = std::vector<double>(n_q_points);
  dynamic_viscosity_0                   = std::vector<double>(n_q_points);
  dynamic_viscosity_1                   = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_0 = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_1 = std::vector<double>(n_q_points);
  thermal_expansion_0                   = std::vector<double>(n_q_points);
  thermal_expansion_1                   = std::vector<double>(n_q_points);
  surface_tension                       = std::vector<double>(n_q_points);
  surface_tension_gradient              = std::vector<double>(n_q_points);

  // Create filter
  this->cahn_hilliard_filter = cahn_hilliard_filter;
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_projected_phase_fraction_gradient(
  const FiniteElement<dim> &fe_projected_phase_fraction_gradient,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_projected_phase_fraction_gradient = true;
  fe_values_projected_phase_fraction_gradient =
    std::make_shared<FEValues<dim>>(mapping,
                                    fe_projected_phase_fraction_gradient,
                                    quadrature,
                                    update_values | update_gradients);

  // phase fraction gradient (PFG)
  projected_phase_fraction_gradient_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_curvature(
  const FiniteElement<dim> &fe_curvature,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_curvature    = true;
  fe_values_curvature = std::make_shared<FEValues<dim>>(
    mapping, fe_curvature, quadrature, update_values | update_gradients);

  // curvature
  curvature_values = std::vector<double>(this->n_q_points);
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_void_fraction(
  const FiniteElement<dim> &fe,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_void_fraction = true;

  // Contrary to the other physics, we enable the calculation of the JxW values
  // on the void fraction
  fe_values_void_fraction =
    std::make_shared<FEValues<dim>>(mapping,
                                    fe,
                                    quadrature,
                                    update_values | update_gradients |
                                      update_JxW_values);

  // Void Fraction
  void_fraction_values = std::vector<double>(this->n_q_points);
  previous_void_fraction_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
  void_fraction_gradient_values = std::vector<Tensor<1, dim>>(this->n_q_points);
}

template <int dim>
void
NavierStokesScratchData<dim>::enable_particle_field_projection(
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping,
  const FiniteElement<dim> &fe_particle_drag_proj,
  const FiniteElement<dim> &fe_particle_two_way_coupling_force_proj,
  const FiniteElement<dim> &fe_particle_velocity_proj,
  const FiniteElement<dim> &fe_particle_momentum_transfer_coefficient)
{
  gather_particle_field_project = true;

  fe_values_particle_drag = std::make_shared<FEValues<dim>>(
    mapping, fe_particle_drag_proj, quadrature, update_values);

  fe_values_particle_two_way_coupling_force =
    std::make_shared<FEValues<dim>>(mapping,
                                    fe_particle_two_way_coupling_force_proj,
                                    quadrature,
                                    update_values);
  fe_values_particle_velocity = std::make_shared<FEValues<dim>>(
    mapping, fe_particle_velocity_proj, quadrature, update_values);
  fe_values_particle_momentum_transfer_coefficient =
    std::make_shared<FEValues<dim>>(mapping,
                                    fe_particle_momentum_transfer_coefficient,
                                    quadrature,
                                    update_values);

  particle_drag_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  particle_two_way_coupling_force_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);
  particle_velocity_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  particle_momentum_transfer_coefficient_values =
    std::vector<double>(this->n_q_points);
}

template <int dim>
void
NavierStokesScratchData<dim>::enable_particle_fluid_interactions(
  const unsigned int n_global_max_particles_per_cell,
  const bool         enable_void_fraction_interpolation)
{
  gather_particles_information     = true;
  max_number_of_particles_per_cell = n_global_max_particles_per_cell;
  interpolated_void_fraction       = enable_void_fraction_interpolation;

  // Reinitialize vectors used to store flow information at the particle
  // location that do not rely on a quadrature
  particle_velocity =
    std::vector<Tensor<1, dim>>(n_global_max_particles_per_cell);
  fluid_particle_relative_velocity_at_particle_location =
    std::vector<Tensor<1, dim>>(n_global_max_particles_per_cell);
  Re_particle = std::vector<double>(n_global_max_particles_per_cell);

  // This table is not used within an FeValues so we can avoid resizing it
  // dynamically
  density_at_particle_location.resize(n_global_max_particles_per_cell);
  kinematic_viscosity_at_particle_location.resize(
    n_global_max_particles_per_cell);
}

template <int dim>
void
NavierStokesScratchData<dim>::enable_mortar()
{
  gather_mortar = true;
}

template <int dim>
void
NavierStokesScratchData<dim>::reinit_mortar(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const Parameters::Mortar<dim>                        &mortar_parameters,
  const double                                         &radius)
{
  const auto cell_center = cell->center();

  // Get updated rotor angular velocity
  mortar_parameters.rotor_angular_velocity->set_time(
    this->simulation_control->get_current_time());
  const double rotor_angular_velocity =
    mortar_parameters.rotor_angular_velocity->value(Point<dim>());

  // Compute radius between center of rotation and current cell center
  const double radius_current =
    cell_center.distance(mortar_parameters.center_of_rotation);

  // Use prescribed rotor angular velocity only if cell is part of the rotor
  double cell_rotor_angular_velocity;
  if (radius_current > radius)
    cell_rotor_angular_velocity = 0.0;
  else
    cell_rotor_angular_velocity = rotor_angular_velocity;

  // Compute rotor linear velocity at quadrature points
  rotor_linear_velocity_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Assumption in 3D case: rotation axis is in z
      // TODO generalize rotation axis
      const auto x                       = fe_values.quadrature_point(q)[0];
      const auto y                       = fe_values.quadrature_point(q)[1];
      rotor_linear_velocity_values[q][0] = -cell_rotor_angular_velocity * y;
      rotor_linear_velocity_values[q][1] = cell_rotor_angular_velocity * x;

      // Update velocity for stabilization
      this->velocity_for_stabilization[q] -= rotor_linear_velocity_values[q];
    }
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_heat_transfer(
  const FiniteElement<dim> &fe,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_temperature    = true;
  fe_values_temperature = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  temperature_values = std::vector<double>(this->n_q_points);
  previous_temperature_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
  temperature_gradients = std::vector<Tensor<1, dim>>(this->n_q_points);
  previous_temperature_gradients = std::vector<std::vector<Tensor<1, dim>>>(
    maximum_number_of_previous_solutions(),
    std::vector<Tensor<1, dim>>(this->n_q_points));
  fields.insert(
    std::pair<field, std::vector<double>>(field::temperature, n_q_points));
}


template <int dim>
void
NavierStokesScratchData<dim>::calculate_physical_properties()
{
  if (properties_manager.field_is_required(field::temperature) &&
      gather_temperature)
    {
      set_field_vector(field::temperature,
                       this->temperature_values,
                       this->fields);
    }

  if (properties_manager.field_is_required(field::pressure))
    {
      set_field_vector(field::pressure, this->pressure_values, this->fields);
    }

  if (properties_manager.field_is_required(field::shear_rate))
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // Calculate shear rate (at each q)
          const Tensor<2, dim> shear_rate_tensor =
            velocity_gradients[q] + transpose(velocity_gradients[q]);

          // Calculate the shear rate magnitude
          shear_rate[q] = calculate_shear_rate_magnitude(shear_rate_tensor);
        }

      set_field_vector(field::shear_rate, shear_rate, this->fields);
    }

  switch (properties_manager.get_number_of_fluids())
    {
      case 1:
        {
          // In the case of an incompressible flow, viscosity is the only
          // required property
          const auto rheology_model = properties_manager.get_rheology();
          rheology_model->vector_value(fields, kinematic_viscosity);

          rheology_model->get_kinematic_viscosity_for_stabilization_vector(
            fields, kinematic_viscosity_for_stabilization);

          kinematic_viscosity_scale =
            rheology_model->get_kinematic_viscosity_scale();

          if (!properties_manager.density_is_constant())
            {
              // For a weakly compressible flow, density variations will play a
              // role
              const auto density_model = properties_manager.get_density();
              density_model->vector_value(fields, density);
              density_psi = density_model->get_psi();
              density_ref = density_model->get_density_ref();
              // Dynamic viscosity is also necessary for compressible flows
              rheology_model->get_dynamic_viscosity_vector(density_ref,
                                                           fields,
                                                           dynamic_viscosity);
              rheology_model->get_dynamic_viscosity_for_stabilization_vector(
                density_ref, fields, dynamic_viscosity_for_stabilization);
            }
          else
            {
              density_scale = properties_manager.get_density_scale();
            }

          if (properties_manager.is_non_newtonian())
            {
              // Calculate derivative of viscosity with respect to shear rate
              rheology_model->vector_jacobian(
                fields, field::shear_rate, grad_kinematic_viscosity_shear_rate);
            }

          if (gather_temperature)
            {
              const auto thermal_expansion_model =
                properties_manager.get_thermal_expansion();
              thermal_expansion_model->vector_value(fields, thermal_expansion);
            }
          break;
        }
      case 2:
        {
          // We need both density and viscosity
          const auto density_model_0  = properties_manager.get_density(0);
          const auto rheology_model_0 = properties_manager.get_rheology(0);
          const auto density_model_1  = properties_manager.get_density(1);
          const auto rheology_model_1 = properties_manager.get_rheology(1);

          density_ref_0 = density_model_0->get_density_ref();
          density_ref_1 = density_model_1->get_density_ref();

          kinematic_viscosity_scale_0 =
            rheology_model_0->get_kinematic_viscosity_scale();
          kinematic_viscosity_scale_1 =
            rheology_model_1->get_kinematic_viscosity_scale();

          kinematic_viscosity_scale =
            std::max(kinematic_viscosity_scale_0, kinematic_viscosity_scale_1);

          // Gather properties from material interactions if necessary
          if (properties_manager.get_number_of_material_interactions() > 0)
            {
              const auto material_interaction_id =
                properties_manager.get_material_interaction_id(
                  material_interactions_type::fluid_fluid, 0, 1);
              // Gather surface tension
              const auto surface_tension_model =
                properties_manager.get_surface_tension(material_interaction_id);
              surface_tension_model->vector_value(fields, surface_tension);
              // Gather surface tension gradient only if necessary
              if (!properties_manager.surface_tension_is_constant())
                surface_tension_model->vector_jacobian(
                  fields, field::temperature, surface_tension_gradient);
            }

          density_model_0->vector_value(fields, density_0);
          rheology_model_0->get_dynamic_viscosity_vector(density_ref_0,
                                                         fields,
                                                         dynamic_viscosity_0);
          rheology_model_0->get_dynamic_viscosity_for_stabilization_vector(
            density_ref_0, fields, dynamic_viscosity_for_stabilization_0);

          density_model_1->vector_value(fields, density_1);
          rheology_model_1->get_dynamic_viscosity_vector(density_ref_1,
                                                         fields,
                                                         dynamic_viscosity_1);
          rheology_model_1->get_dynamic_viscosity_for_stabilization_vector(
            density_ref_1, fields, dynamic_viscosity_for_stabilization_1);
          if (gather_temperature)
            {
              const auto thermal_expansion_model_0 =
                properties_manager.get_thermal_expansion(0);
              const auto thermal_expansion_model_1 =
                properties_manager.get_thermal_expansion(1);
              thermal_expansion_model_0->vector_value(fields,
                                                      thermal_expansion_0);
              thermal_expansion_model_1->vector_value(fields,
                                                      thermal_expansion_1);
            }

          if (gather_vof && !gather_cahn_hilliard)
            {
              for (unsigned int q = 0; q < this->n_q_points; ++q)
                {
                  double filtered_phase_value = this->filtered_phase_values[q];

                  density[q] = calculate_point_property(filtered_phase_value,
                                                        this->density_0[q],
                                                        this->density_1[q]);

                  dynamic_viscosity[q] =
                    calculate_point_property(filtered_phase_value,
                                             this->dynamic_viscosity_0[q],
                                             this->dynamic_viscosity_1[q]);

                  kinematic_viscosity[q] = dynamic_viscosity[q] / density[q];

                  dynamic_viscosity_for_stabilization[q] =
                    calculate_point_property(
                      filtered_phase_value,
                      this->dynamic_viscosity_for_stabilization_0[q],
                      this->dynamic_viscosity_for_stabilization_1[q]);

                  thermal_expansion[q] =
                    calculate_point_property(filtered_phase_value,
                                             this->thermal_expansion_0[q],
                                             this->thermal_expansion_1[q]);
                }

              // Gather density_psi for isothermal compressible NS equations
              if (!properties_manager.density_is_constant())
                {
                  density_psi_0 = density_model_0->get_psi();
                  density_psi_1 = density_model_1->get_psi();

                  for (unsigned int q = 0; q < this->n_q_points; ++q)
                    {
                      double filtered_phase_value =
                        this->filtered_phase_values[q];
                      // To avoid calculating twice in the assemblers
                      compressibility_multiplier[q] =
                        filtered_phase_value *
                          (density_0[q] / density_1[q] * density_psi_1 -
                           density_psi_0) +
                        density_psi_0;
                    }
                }

              for (unsigned p = 0; p < previous_phase_values.size(); ++p)
                {
                  // Blend the physical properties using the VOF field
                  for (unsigned int q = 0; q < this->n_q_points; ++q)
                    {
                      // Calculate previous density (right now assumes constant
                      // density model per phase)
                      previous_density[p][q] = calculate_point_property(
                        filter->filter_phase(this->previous_phase_values[p][q]),
                        this->density_0[q],
                        this->density_1[q]);
                    }
                }
              break;
            }
          else if (gather_cahn_hilliard && !gather_vof)
            {
              // Blend the physical properties using the CahnHilliard field
              for (unsigned int q = 0; q < this->n_q_points; ++q)
                {
                  double phase_order_cahn_hilliard_value =
                    this->filtered_phase_order_cahn_hilliard_values[q];

                  density[q] = calculate_point_property_cahn_hilliard(
                    phase_order_cahn_hilliard_value,
                    this->density_0[q],
                    this->density_1[q]);

                  dynamic_viscosity[q] = calculate_point_property_cahn_hilliard(
                    phase_order_cahn_hilliard_value,
                    this->dynamic_viscosity_0[q],
                    this->dynamic_viscosity_1[q]);

                  dynamic_viscosity_for_stabilization[q] =
                    calculate_point_property_cahn_hilliard(
                      phase_order_cahn_hilliard_value,
                      this->dynamic_viscosity_for_stabilization_0[q],
                      this->dynamic_viscosity_for_stabilization_1[q]);
                }
              break;
            }
          break;
        }
      default:
        throw std::runtime_error("Unsupported number of fluids (>2)");
    }
}

template <int dim>
void
NavierStokesScratchData<dim>::reinit_particle_fluid_forces()
{
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();
      // Set the particle_fluid_interactions properties and vectors to 0
      for (int d = 0; d < dim; ++d)
        {
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                fem_force_two_way_coupling_x +
                              d] = 0.;
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                                fem_force_one_way_coupling_x +
                              d] = 0.;
          particle_properties
            [DEM::CFDDEMProperties::PropertiesIndex::fem_drag_x + d] = 0.;
          particle_properties
            [DEM::CFDDEMProperties::PropertiesIndex::fem_torque_x + d] = 0.;
          explicit_particle_volumetric_acceleration_on_fluid[d]        = 0.;
        }
    }
}

template <int dim>
void
NavierStokesScratchData<dim>::extract_particle_properties()
{
  this->average_particle_velocity = 0;
  // Loop over particles in cell
  this->total_particle_volume = 0;
  unsigned int i_particle     = 0;

  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();
      // Stores the values of particle velocity in a tensor
      particle_velocity[i_particle][0] =
        particle_properties[DEM::CFDDEMProperties::PropertiesIndex::v_x];
      particle_velocity[i_particle][1] =
        particle_properties[DEM::CFDDEMProperties::PropertiesIndex::v_y];
      if constexpr (dim == 3)
        particle_velocity[i_particle][2] =
          particle_properties[DEM::CFDDEMProperties::PropertiesIndex::v_z];

      if (!interpolated_void_fraction)
        total_particle_volume +=
          M_PI *
          pow(particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp],
              dim) /
          (2 * dim);

      average_particle_velocity += particle_velocity[i_particle];
      i_particle++;
    }
  number_of_particles = i_particle;
  if (number_of_particles != 0)
    { // Calculate the average particle velocity within the cell
      average_particle_velocity =
        average_particle_velocity / number_of_particles;
    }
}

template <int dim>
void
NavierStokesScratchData<dim>::calculate_cell_void_fraction(
  const double &total_particle_volume)
{
  cell_volume = compute_cell_measure_with_JxW(
    this->fe_values_void_fraction->get_JxW_values());

  // The cell_void_fraction vector needs to be resized to the correct size
  // since the quadrature used in the FEValues for the void fraction will
  // be tailored exactly to the number of particles.
  // We also initialize it to zero and it will be filled if a cell-averaged
  // void fraction is used.
  cell_void_fraction.resize(number_of_particles, 0);

  if (!this->interpolated_void_fraction)
    {
      double cell_void_fraction_bulk = 0;
      cell_void_fraction_bulk =
        (cell_volume - total_particle_volume) / cell_volume;

      for (unsigned int j = 0; j < number_of_particles; ++j)
        cell_void_fraction[j] = cell_void_fraction_bulk;
    }
}

template <int dim>
Quadrature<dim>
NavierStokesScratchData<dim>::gather_particles_reference_location()
{
  // Create local vector that will be used to spawn an in-situ quadrature to
  // interpolate at the locations of the particles
  std::vector<Point<dim>> particle_reference_location(number_of_particles);
  std::vector<double>     particle_weights(number_of_particles, 1);
  unsigned int            i_particle = 0;

  // Loop over particles in cell and cache their reference location
  for (auto &particle : pic)
    {
      // Store particle positions and weights
      // Reference location of the particle
      particle_reference_location[i_particle] =
        particle.get_reference_location();
      i_particle++;
    }

  // Return a quadrature for the Navier-Stokes equations that is based on the
  // particle reference location
  return Quadrature<dim>(particle_reference_location, particle_weights);
}
template class NavierStokesScratchData<2>;
template class NavierStokesScratchData<3>;
