// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/magnetic_permeability_model.h>

std::shared_ptr<MagneticPermeabilityModel>
MagneticPermeabilityModel::model_cast_real(
  const Parameters::Material &material_properties)
{
  if (material_properties.magnetic_permeability_model ==
      Parameters::Material::MagneticPermeabilityModel::constant)
    return std::make_shared<ConstantMagneticPermeability>(
      material_properties.magnetic_permeability_real);
  else
    AssertThrow(
      false,
      ExcMessage(
        "Invalid magnetic permeability model. The choices are <constant>"));
}

std::shared_ptr<MagneticPermeabilityModel>
MagneticPermeabilityModel::model_cast_imag(
  const Parameters::Material &material_properties)
{
  if (material_properties.magnetic_permeability_model ==
      Parameters::Material::MagneticPermeabilityModel::constant)
    return std::make_shared<ConstantMagneticPermeability>(
      material_properties.magnetic_permeability_imag);
  else
    AssertThrow(
      false,
      ExcMessage(
        "Invalid magnetic permeability model. The choices are <constant>"));
}
