// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/tracer_reaction_model.h>

std::shared_ptr<TracerReactionPrefactorModel>
TracerReactionPrefactorModel::model_cast(
  const Parameters::Material &material_properties)
{
  if (material_properties.tracer_reaction_constant_model ==
      Parameters::Material::TracerReactionConstantModel::immersed_boundary_tanh)
    return std::make_shared<TanhLevelsetTracerReactionPrefactor>(
      material_properties.immersed_solid_tanh_parameters
        .tracer_reaction_constant_outside,
      material_properties.immersed_solid_tanh_parameters
        .tracer_reaction_constant_inside,
      material_properties.immersed_solid_tanh_parameters.thickness,
      material_properties.tracer_reaction_order);
  else if (material_properties.tracer_reaction_constant_model ==
           Parameters::Material::TracerReactionConstantModel::constant)
    return std::make_shared<ConstantTracerReactionPrefactor>(
      material_properties.tracer_reaction_constant,
      material_properties.tracer_reaction_order);
  else
    return std::make_shared<NoneTracerReactionPrefactor>();
}
