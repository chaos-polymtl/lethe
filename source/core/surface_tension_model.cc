// SPDX-FileCopyrightText: Copyright (c) 2023 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/surface_tension_model.h>

std::shared_ptr<SurfaceTensionModel>
SurfaceTensionModel::model_cast(
  const Parameters::MaterialInteractions &material_interaction_parameters)
{
  if (material_interaction_parameters.surface_tension_model ==
      Parameters::MaterialInteractions::SurfaceTensionModel::linear)
    return std::make_shared<SurfaceTensionLinear>(
      material_interaction_parameters.surface_tension_parameters);
  else if (material_interaction_parameters.surface_tension_model ==
           Parameters::MaterialInteractions::SurfaceTensionModel::phase_change)
    return std::make_shared<SurfaceTensionPhaseChange>(
      material_interaction_parameters.surface_tension_parameters);
  else
    return std::make_shared<SurfaceTensionConstant>(
      material_interaction_parameters.surface_tension_parameters
        .surface_tension_coefficient);
}
