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

#include <core/evaporation_model.h>

std::shared_ptr<EvaporationModel>
EvaporationModel::model_cast(
  const Parameters::Evaporation &evaporation_parameters)
{
  if (evaporation_parameters.evaporative_mass_flux_model_type ==
      Parameters::Evaporation::EvaporativeMassFluxModelType::
        temperature_dependent)
    return std::make_shared<EvaporationModelTemperature>(
      evaporation_parameters);
  else
    return std::make_shared<EvaporationModelConstant>(evaporation_parameters);
}
