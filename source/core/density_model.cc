// SPDX-FileCopyrightText: Copyright (c) 2022-2023 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/density_model.h>

std::shared_ptr<DensityModel>
DensityModel::model_cast(const Parameters::Material &material_properties)
{
  if (material_properties.density_model ==
      Parameters::Material::DensityModel::isothermal_ideal_gas)
    {
      return std::make_shared<DensityIsothermalIdealGas>(
        material_properties.isothermal_ideal_gas_density_parameters.density_ref,
        material_properties.isothermal_ideal_gas_density_parameters.R,
        material_properties.isothermal_ideal_gas_density_parameters.T);
    }

  return std::make_shared<DensityConstant>(material_properties.density);
}
