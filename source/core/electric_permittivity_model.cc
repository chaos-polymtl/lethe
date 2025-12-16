// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/electric_permittivity_model.h>

std::shared_ptr<ElectricPermittivityModel>
ElectricPermittivityModel::model_cast_real(
  const Parameters::Material &material_properties)
{
  if (material_properties.electric_permittivity_model ==
      Parameters::Material::ElectricPermittivityModel::constant)
    return std::make_shared<ConstantElectricPermittivity>(
      material_properties.electric_permittivity_real);
  else
    AssertThrow(
      false,
      ExcMessage(
        "Invalid electric permittivity model. The choices are <constant>"));
}

std::shared_ptr<ElectricPermittivityModel>
ElectricPermittivityModel::model_cast_imag(
  const Parameters::Material &material_properties)
{
  if (material_properties.electric_permittivity_model ==
      Parameters::Material::ElectricPermittivityModel::constant)
    return std::make_shared<ConstantElectricPermittivity>(
      material_properties.electric_permittivity_imag);
  else
    AssertThrow(
      false,
      ExcMessage(
        "Invalid electric permittivity model. The choices are <constant>"));
}
