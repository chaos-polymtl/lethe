// SPDX-FileCopyrightText: Copyright (c) 2023 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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
