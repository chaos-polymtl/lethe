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

#include <core/thermal_conductivity_model.h>

std::shared_ptr<ThermalConductivityModel>
ThermalConductivityModel::model_cast(
  const Parameters::Material &material_properties)
{
  if (material_properties.thermal_conductivity_model ==
      Parameters::Material::ThermalConductivityModel::linear)
    return std::make_shared<ThermalConductivityLinear>(
      material_properties.k_A0, material_properties.k_A1);
  else if (material_properties.thermal_conductivity_model ==
           Parameters::Material::ThermalConductivityModel::phase_change)
    return std::make_shared<ThermalConductivityPhaseChange>(
      material_properties.phase_change_parameters);
  else
    return std::make_shared<ConstantThermalConductivity>(
      material_properties.thermal_conductivity);
}
