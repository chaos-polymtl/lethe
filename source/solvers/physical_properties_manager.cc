/*---------------------------------------------------------------------
 *
 * Copyright (C) 2021 - by the Lethe authors
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
  is_initialized = true;

  number_of_fluids = physical_properties.number_of_fluids;
  number_of_solids = physical_properties.number_of_solids;
  number_of_material_interactions =
    physical_properties.number_of_material_interactions;
  fluid_fluid_interactions_with_material_interaction_ids =
    physical_properties.fluid_fluid_interactions_with_material_interaction_ids;
  fluid_solid_interactions_with_material_interaction_ids =
    physical_properties.fluid_solid_interactions_with_material_interaction_ids;

  viscosity_scale = physical_properties.fluids[0].viscosity;
  density_scale   = physical_properties.fluids[0].density;

  non_newtonian_flow = false;
  constant_density   = true;

  required_fields[field::temperature]          = false;
  required_fields[field::previous_temperature] = false;
  required_fields[field::shear_rate]           = false;
  required_fields[field::pressure]             = false;

  // For each fluid, declare the physical properties
  for (unsigned int f = 0; f < number_of_fluids; ++f)
    {
      density.push_back(
        DensityModel::model_cast(physical_properties.fluids[f]));
      establish_fields_required_by_model(*density[f]);

      // Store an indicator for the density to indicate if it is not constant
      // This indicator is used elsewhere in the code to throw assertions
      // if non-constant density is not implemented in a post-processing utility
      if (physical_properties.fluids[f].density_model !=
          Parameters::Material::DensityModel::constant)
        constant_density = false;

      specific_heat.push_back(
        SpecificHeatModel::model_cast(physical_properties.fluids[f]));
      establish_fields_required_by_model(*specific_heat[f]);

      thermal_conductivity.push_back(
        ThermalConductivityModel::model_cast(physical_properties.fluids[f]));
      establish_fields_required_by_model(*thermal_conductivity[f]);

      rheology.push_back(
        RheologicalModel::model_cast(physical_properties.fluids[f]));
      this->establish_fields_required_by_model(*rheology[f]);

      tracer_diffusivity.push_back(
        TracerDiffusivityModel::model_cast(physical_properties.fluids[f]));
      establish_fields_required_by_model(*tracer_diffusivity[f]);

      thermal_expansion.push_back(
        ThermalExpansionModel::model_cast(physical_properties.fluids[f]));
      establish_fields_required_by_model(*thermal_expansion[f]);

      if (physical_properties.fluids[f].rheological_model !=
            Parameters::Material::RheologicalModel::newtonian &&
          physical_properties.fluids[f].rheological_model !=
            Parameters::Material::RheologicalModel::phase_change)
        non_newtonian_flow = true;
    }

  // For each solid, append the physical properties
  for (unsigned int s = 0; s < number_of_solids; ++s)
    {
      density.push_back(
        DensityModel::model_cast(physical_properties.solids[s]));
      establish_fields_required_by_model(*density[s]);

      // Store an indicator for the density to indicate if it is not constant
      // This indicator is used elsewhere in the code to throw assertions
      // if non-constant density is not implemented in a post-processing utility
      if (physical_properties.fluids[s].density_model !=
          Parameters::Material::DensityModel::constant)
        constant_density = false;

      specific_heat.push_back(
        SpecificHeatModel::model_cast(physical_properties.solids[s]));
      establish_fields_required_by_model(*specific_heat[s]);

      thermal_conductivity.push_back(
        ThermalConductivityModel::model_cast(physical_properties.solids[s]));
      establish_fields_required_by_model(*thermal_conductivity[s]);

      rheology.push_back(
        RheologicalModel::model_cast(physical_properties.solids[s]));
      this->establish_fields_required_by_model(*rheology[s]);

      tracer_diffusivity.push_back(
        TracerDiffusivityModel::model_cast(physical_properties.solids[s]));
      establish_fields_required_by_model(*tracer_diffusivity[s]);

      thermal_expansion.push_back(
        ThermalExpansionModel::model_cast(physical_properties.solids[s]));
      establish_fields_required_by_model(*thermal_expansion[s]);
    }

  // For each pair of interacting materials (fluid-fluid or fluid-solid), append
  // physical properties
  for (unsigned int i = 0; i < number_of_material_interactions; ++i)
    {
      surface_tension.push_back(SurfaceTensionModel::model_cast(
        physical_properties.material_interactions[i]
          .surface_tension_parameters));
      establish_fields_required_by_model(*surface_tension[i]);
    }
}
