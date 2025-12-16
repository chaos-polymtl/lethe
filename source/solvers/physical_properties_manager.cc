// SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/physical_properties_manager.h>

void
PhysicalPropertiesManager::establish_fields_required_by_model(
  PhysicalPropertyModel &model)
{
  // Loop through the map. The use of or (||) is there to ensure
  // that if a field is already required, it won't be erased.
  for (auto &f : required_fields)
    {
      f.second = f.second || model.depends_on(f.first);
    }
}

void
PhysicalPropertiesManager::establish_fields_required_by_model(
  InterfacePropertyModel &model)
{
  // Loop through the map. The use of or (||) is there to ensure
  // that if a field is already required, it won't be erased.
  for (auto &f : required_fields)
    {
      f.second = f.second || model.depends_on(f.first);
    }
}

void
PhysicalPropertiesManager::initialize(
  Parameters::PhysicalProperties physical_properties)
{
  // Keep a copy of the physical properties used to build the physical property
  // manager
  physical_properties_parameters = physical_properties;

  is_initialized = true;

  number_of_fluids = physical_properties.number_of_fluids;
  number_of_solids = physical_properties.number_of_solids;
  number_of_material_interactions =
    physical_properties.number_of_material_interactions;
  fluid_fluid_interactions_with_material_interaction_ids =
    physical_properties.fluid_fluid_interactions_with_material_interaction_ids;
  fluid_solid_interactions_with_material_interaction_ids =
    physical_properties.fluid_solid_interactions_with_material_interaction_ids;
  reference_temperature = physical_properties.reference_temperature;

  non_newtonian_flow       = false;
  phase_change             = false;
  constant_density         = true;
  constant_surface_tension = true;

  required_fields[field::temperature]               = false;
  required_fields[field::temperature_p1]            = false;
  required_fields[field::temperature_p2]            = false;
  required_fields[field::shear_rate]                = false;
  required_fields[field::pressure]                  = false;
  required_fields[field::phase_order_cahn_hilliard] = false;
  required_fields[field::levelset]                  = false;
  required_fields[field::tracer_concentration]      = false;

  // For each fluid, declare the physical properties
  for (unsigned int f = 0; f < number_of_fluids; ++f)
    {
      density.push_back(
        DensityModel::model_cast(physical_properties.fluids[f]));
      establish_fields_required_by_model(*density[f]);

      // Store an indicator for the density to indicate if it is not constant
      // This indicator is used elsewhere in the code to throw assertions
      // if non-constant density is not implemented in a post-processing utility
      if (!density.back()->is_constant_density_model())
        constant_density = false;

      specific_heat.push_back(
        SpecificHeatModel::model_cast(physical_properties.fluids[f]));
      establish_fields_required_by_model(*specific_heat[f]);

      // Store an indicator that a phase change model is present
      phase_change =
        phase_change || physical_properties.fluids[f].specific_heat_model ==
                          Parameters::Material::SpecificHeatModel::phase_change;

      thermal_conductivity.push_back(
        ThermalConductivityModel::model_cast(physical_properties.fluids[f]));
      establish_fields_required_by_model(*thermal_conductivity[f]);

      rheology.push_back(
        RheologicalModel::model_cast(physical_properties.fluids[f]));
      this->establish_fields_required_by_model(*rheology[f]);
      rheology[f]->set_dynamic_viscosity(density[f]->get_density_ref());

      tracer_diffusivity.push_back(
        TracerDiffusivityModel::model_cast(physical_properties.fluids[f]));
      establish_fields_required_by_model(*tracer_diffusivity[f]);

      tracer_reaction_prefactor.push_back(
        TracerReactionPrefactorModel::model_cast(
          physical_properties.fluids[f]));
      establish_fields_required_by_model(*tracer_reaction_prefactor[f]);

      tracer_reaction_order.push_back(
        physical_properties.fluids[f].tracer_reaction_order);

      thermal_expansion.push_back(
        ThermalExpansionModel::model_cast(physical_properties.fluids[f]));
      establish_fields_required_by_model(*thermal_expansion[f]);

      phase_change_parameters.push_back(
        physical_properties.fluids[f].phase_change_parameters);

      if (rheology.back()->is_non_newtonian_rheological_model())
        non_newtonian_flow = true;

      electric_conductivity.push_back(
        ElectricConductivityModel::model_cast(physical_properties.fluids[f]));

      electric_permittivity_real.push_back(
        ElectricPermittivityModel::model_cast_real(
          physical_properties.fluids[f]));
      electric_permittivity_imag.push_back(
        ElectricPermittivityModel::model_cast_imag(
          physical_properties.fluids[f]));

      magnetic_permeability_real.push_back(
        MagneticPermeabilityModel::model_cast_real(
          physical_properties.fluids[f]));
      magnetic_permeability_imag.push_back(
        MagneticPermeabilityModel::model_cast_imag(
          physical_properties.fluids[f]));
    }

  // For each solid, append the physical properties
  for (unsigned int s = 0; s < number_of_solids; ++s)
    {
      density.push_back(
        DensityModel::model_cast(physical_properties.solids[s]));
      establish_fields_required_by_model(*density[s + number_of_fluids]);

      // Store an indicator for the density to indicate if it is not constant
      // This indicator is used elsewhere in the code to throw assertions
      // if non-constant density is not implemented in a post-processing utility
      if (!density.back()->is_constant_density_model())
        constant_density = false;

      specific_heat.push_back(
        SpecificHeatModel::model_cast(physical_properties.solids[s]));
      establish_fields_required_by_model(*specific_heat[s + number_of_fluids]);

      thermal_conductivity.push_back(
        ThermalConductivityModel::model_cast(physical_properties.solids[s]));
      establish_fields_required_by_model(
        *thermal_conductivity[s + number_of_fluids]);

      rheology.push_back(
        RheologicalModel::model_cast(physical_properties.solids[s]));
      this->establish_fields_required_by_model(*rheology[s + number_of_fluids]);
      rheology.back()->set_dynamic_viscosity(density.back()->get_density_ref());

      tracer_diffusivity.push_back(
        TracerDiffusivityModel::model_cast(physical_properties.solids[s]));
      establish_fields_required_by_model(
        *tracer_diffusivity[s + number_of_fluids]);

      tracer_reaction_prefactor.push_back(
        TracerReactionPrefactorModel::model_cast(
          physical_properties.solids[s]));
      establish_fields_required_by_model(
        *tracer_reaction_prefactor[s + number_of_fluids]);

      tracer_reaction_order.push_back(
        physical_properties.fluids[s].tracer_reaction_order);

      thermal_expansion.push_back(
        ThermalExpansionModel::model_cast(physical_properties.solids[s]));
      establish_fields_required_by_model(
        *thermal_expansion[s + number_of_fluids]);

      electric_conductivity.push_back(
        ElectricConductivityModel::model_cast(physical_properties.solids[s]));

      electric_permittivity_real.push_back(
        ElectricPermittivityModel::model_cast_real(
          physical_properties.solids[s]));
      electric_permittivity_imag.push_back(
        ElectricPermittivityModel::model_cast_imag(
          physical_properties.solids[s]));

      magnetic_permeability_real.push_back(
        MagneticPermeabilityModel::model_cast_real(
          physical_properties.solids[s]));
      magnetic_permeability_imag.push_back(
        MagneticPermeabilityModel::model_cast_imag(
          physical_properties.solids[s]));
    }

  // For each pair of interacting materials (fluid-fluid or fluid-solid), append
  // physical properties
  for (unsigned int i = 0; i < number_of_material_interactions; ++i)
    {
      surface_tension.push_back(SurfaceTensionModel::model_cast(
        physical_properties.material_interactions[i]));
      establish_fields_required_by_model(*surface_tension[i]);
      if (!surface_tension.back()->is_constant_surface_tension_model())
        constant_surface_tension = false;

      mobility_cahn_hilliard.push_back(MobilityCahnHilliardModel::model_cast(
        physical_properties.material_interactions[i]));
      establish_fields_required_by_model(*mobility_cahn_hilliard[i]);
    }
}
