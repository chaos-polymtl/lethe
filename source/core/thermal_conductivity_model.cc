// SPDX-FileCopyrightText: Copyright (c) 2022-2023 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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
