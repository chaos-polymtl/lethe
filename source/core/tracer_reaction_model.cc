// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/tracer_reaction_model.h>

std::shared_ptr<TracerReactionPrefactorModel>
TracerReactionPrefactorModel::model_cast(
  const Parameters::Material &material_properties)
{
  if (material_properties.tracer_reaction_prefactor_model ==
      Parameters::Material::TracerReactionPrefactorModel::
        immersed_boundary_tanh)
    return std::make_shared<TanhLevelsetTracerReactionPrefactor>(
      material_properties.immersed_solid_tanh_parameters
        .tracer_reaction_constant_outside,
      material_properties.immersed_solid_tanh_parameters
        .tracer_reaction_constant_inside,
      material_properties.immersed_solid_tanh_parameters.thickness,
      material_properties.tracer_reaction_order,
      material_properties.tracer_reaction_epsilon);
  else if (material_properties.tracer_reaction_prefactor_model ==
           Parameters::Material::TracerReactionPrefactorModel::
             immersed_boundary_gaussian)
    return std::make_shared<GaussianLevelsetTracerReactionPrefactor>(
      material_properties.immersed_solid_gaussian_parameters
        .tracer_reaction_constant_interface,
      material_properties.immersed_solid_gaussian_parameters
        .tracer_reaction_constant_bulk,
      material_properties.immersed_solid_gaussian_parameters.thickness,
      material_properties.tracer_reaction_order,
      material_properties.tracer_reaction_epsilon);
  else if (material_properties.tracer_reaction_prefactor_model ==
           Parameters::Material::TracerReactionPrefactorModel::constant)
    return std::make_shared<ConstantTracerReactionPrefactor>(
      material_properties.tracer_reaction_constant,
      material_properties.tracer_reaction_order,
      material_properties.tracer_reaction_epsilon);
  else
    return std::make_shared<NoneTracerReactionPrefactor>();
}
