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

#include <solvers/cahn_hilliard_filter.h>

std::shared_ptr<CahnHilliardFilterBase>
CahnHilliardFilterBase::model_cast(
  const Parameters::CahnHilliard &cahn_hilliard_parameters)
{
  if (cahn_hilliard_parameters.cahn_hilliard_phase_filter.type ==
      Parameters::FilterType::tanh)
    return std::make_shared<CahnHilliardTanhFilter>(
      cahn_hilliard_parameters.cahn_hilliard_phase_filter.beta,
      cahn_hilliard_parameters.well_height,
      cahn_hilliard_parameters.epsilon);
  else
    return std::make_shared<CahnHilliardNoFilter>();
}
