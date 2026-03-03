// SPDX-FileCopyrightText: Copyright (c) 2022-2023 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/specific_heat_model.h>

std::shared_ptr<SpecificHeatModel>
SpecificHeatModel::model_cast(const Parameters::Material &material_properties)
{
  if (material_properties.specific_heat_model ==
      Parameters::Material::SpecificHeatModel::phase_change)
    return std::make_shared<PhaseChangeSpecificHeat>(
      material_properties.phase_change_parameters);

  else
    return std::make_shared<ConstantSpecificHeat>(
      material_properties.specific_heat);
}
