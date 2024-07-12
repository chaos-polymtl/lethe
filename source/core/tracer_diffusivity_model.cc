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
  else
    return std::make_shared<ConstantTracerDiffusivity>(
      material_properties.tracer_diffusivity);
}
