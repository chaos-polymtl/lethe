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

#include <solvers/vof_filter.h>

std::shared_ptr<VolumeOfFluidFilterBase>
VolumeOfFluidFilterBase::model_cast(
  const Parameters::VOF_PhaseFilter &phase_filter_parameters)
{
  if (phase_filter_parameters.type == Parameters::FilterType::tanh)
    return std::make_shared<VolumeOfFluidTanhFilter>(
      phase_filter_parameters.beta);
  else
    return std::make_shared<VolumeOfFluidNoFilter>();
}
