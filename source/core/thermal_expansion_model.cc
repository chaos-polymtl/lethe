/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
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

#include <core/thermal_expansion_model.h>

std::shared_ptr<ThermalExpansionModel>
ThermalExpansionModel::model_cast(const Parameters::Fluid &fluid_properties)
{
  if (fluid_properties.thermal_expansion_model ==
      Parameters::Fluid::ThermalExpansionModel::phase_change)
    return std::make_shared<ThermalExpansionPhaseChange>(
      fluid_properties.phase_change_parameters);
  else
    return std::make_shared<ConstantThermalExpansion>(
      fluid_properties.thermal_expansion);
}
