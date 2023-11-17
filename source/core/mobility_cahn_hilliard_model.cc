/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

#include <core/mobility_cahn_hilliard_model.h>

/*
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
*/
