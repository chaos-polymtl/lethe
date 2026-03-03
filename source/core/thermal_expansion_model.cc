// SPDX-FileCopyrightText: Copyright (c) 2022-2023 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/thermal_expansion_model.h>

std::shared_ptr<ThermalExpansionModel>
ThermalExpansionModel::model_cast(
  const Parameters::Material &material_properties)
{
  if (material_properties.thermal_expansion_model ==
      Parameters::Material::ThermalExpansionModel::phase_change)
    return std::make_shared<ThermalExpansionPhaseChange>(
      material_properties.phase_change_parameters);
  else
    return std::make_shared<ConstantThermalExpansion>(
      material_properties.thermal_expansion);
}
