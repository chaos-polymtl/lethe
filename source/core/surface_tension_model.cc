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

#include <core/surface_tension_model.h>

std::shared_ptr<SurfaceTensionModel>
SurfaceTensionModel::model_cast(
  const Parameters::SurfaceTensionParameters &surface_tension_parameters)
{
  return std::make_shared<SurfaceTensionConstant>(
    surface_tension_parameters.surface_tension_coefficient);
}
