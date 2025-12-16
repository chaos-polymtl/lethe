// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/electromagnetic_permittivity_model.h>

std::shared_ptr<ElectroMagneticPermittivityModel>
ElectroMagneticPermittivityModel::model_cast(
  const Parameters::Material &material_properties)
{
  return std::make_shared<ConstantElectroMagneticConductivity>(
    material_properties.electromagnetic_permittivity);
}
