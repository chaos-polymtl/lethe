// SPDX-FileCopyrightText: Copyright (c) 2023 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/mobility_cahn_hilliard_model.h>

std::shared_ptr<MobilityCahnHilliardModel>
MobilityCahnHilliardModel::model_cast(
  const Parameters::MaterialInteractions &material_interaction_parameters)
{
  if (material_interaction_parameters.mobility_cahn_hilliard_model ==
      Parameters::MaterialInteractions::MobilityCahnHilliardModel::quartic)
    {
      return std::make_shared<MobilityCahnHilliardModelQuartic>(
        material_interaction_parameters.mobility_cahn_hilliard_parameters
          .mobility_cahn_hilliard_constant);
    }
  else
    {
      return std::make_shared<MobilityCahnHilliardModelConstant>(
        material_interaction_parameters.mobility_cahn_hilliard_parameters
          .mobility_cahn_hilliard_constant);
    }
}
