// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/cahn_hilliard_filter.h>

std::shared_ptr<CahnHilliardFilterBase>
CahnHilliardFilterBase::model_cast(
  const Parameters::CahnHilliard &cahn_hilliard_parameters)
{
  if (cahn_hilliard_parameters.cahn_hilliard_phase_filter.type ==
      Parameters::FilterType::tanh)
    return std::make_shared<CahnHilliardTanhFilter>(
      cahn_hilliard_parameters.cahn_hilliard_phase_filter.beta);
  else if (cahn_hilliard_parameters.cahn_hilliard_phase_filter.type ==
           Parameters::FilterType::clip)
    return std::make_shared<CahnHilliardClipFilter>();
  else
    return std::make_shared<CahnHilliardNoFilter>();
}
