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
  else if (material_properties.magnetic_permeability_model ==
           Parameters::Material::MagneticPermeabilityModel::polynomial)
    return std::make_shared<PolynomialMagneticPermeability>(
      material_properties.magnetic_permeability_real_polynomial_coefficients);
  else
    AssertThrow(
      false,
      ExcMessage(
        "Invalid magnetic permeability model. The choices are <constant|polynomial>"));
}

std::shared_ptr<MagneticPermeabilityModel>
MagneticPermeabilityModel::model_cast_imag(
  const Parameters::Material &material_properties)
{
  if (material_properties.magnetic_permeability_model ==
      Parameters::Material::MagneticPermeabilityModel::constant)
    return std::make_shared<ConstantMagneticPermeability>(
      material_properties.magnetic_permeability_imag);
  else if (material_properties.magnetic_permeability_model ==
           Parameters::Material::MagneticPermeabilityModel::polynomial)
    return std::make_shared<PolynomialMagneticPermeability>(
      material_properties.magnetic_permeability_imag_polynomial_coefficients);
  else
    AssertThrow(
      false,
      ExcMessage(
        "Invalid magnetic permeability model. The choices are <constant|polynomial>"));
}
