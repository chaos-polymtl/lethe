// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/electromagnetic_conductivity_model.h>

std::shared_ptr<ElectroMagneticConductivityModel>
ElectroMagneticConductivityModel::model_cast(
  const Parameters::Material &material_properties)
{
  return std::make_shared<ConstantElectroMagneticConductivity>(
    material_properties.electromagnetic_conductivity);
}
