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


PhysicalPropertiesManager::PhysicalPropertiesManager(
  Parameters::PhysicalProperties physical_properties)
  : is_initialized(true)
{
  initialize(physical_properties);
}

void
PhysicalPropertiesManager::initialize(
  Parameters::PhysicalProperties physical_properties)
{
  number_of_fluids = physical_properties.number_of_fluids;

  non_newtonian_flow = false;
  // For each fluid, declare the physical properties
  for (unsigned int f = 0; f < number_of_fluids; ++f)
    {
      density.push_back(
        DensityModel::model_cast(physical_properties.fluids[f]));

      specific_heat.push_back(
        SpecificHeatModel::model_cast(physical_properties.fluids[f]));

      thermal_conductivity.push_back(
        ThermalConductivityModel::model_cast(physical_properties.fluids[f]));

      rheology.push_back(
        RheologicalModel::model_cast(physical_properties.fluids[f]));
      if (physical_properties.fluids[f].rheology_model !=
          Parameters::Fluid::RheologyModel::newtonian)
        non_newtonian_flow = true;
    }
}
