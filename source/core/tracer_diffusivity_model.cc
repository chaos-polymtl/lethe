// SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/tracer_diffusivity_model.h>

std::shared_ptr<TracerDiffusivityModel>
TracerDiffusivityModel::model_cast(
  const Parameters::Material &material_properties)
{
  if (material_properties.tracer_diffusivity_model ==
      Parameters::Material::TracerDiffusivityModel::immersed_boundary_tanh)
    return std::make_shared<TanhLevelsetTracerDiffusivity>(
      material_properties.immersed_solid_tanh_parameters
        .tracer_diffusivity_outside,
      material_properties.immersed_solid_tanh_parameters
        .tracer_diffusivity_inside,
      material_properties.immersed_solid_tanh_parameters.thickness);
  else if (material_properties.tracer_diffusivity_model ==
           Parameters::Material::TracerDiffusivityModel::
             immersed_boundary_gaussian)
    return std::make_shared<GaussianLevelsetTracerDiffusivity>(
      material_properties.immersed_solid_gaussian_parameters
        .tracer_diffusivity_interface,
      material_properties.immersed_solid_gaussian_parameters
        .tracer_diffusivity_bulk,
      material_properties.immersed_solid_gaussian_parameters.thickness);
  else
    return std::make_shared<ConstantTracerDiffusivity>(
      material_properties.tracer_diffusivity);
}
