// SPDX-FileCopyrightText: Copyright (c) 2023 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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
