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

#ifndef phase_change_h
#define phase_change_h

#include <core/physical_property_model.h>

using namespace dealii;

/**
 * @brief calculate_liquid_fraction Calculates the liquid fraction of a phase change material at
 * a temperature T
 *
 * @param T temperature at which to calculate the solid fraction
 * @return value of the liquid_fraction
 *
 */
inline double
calculate_liquid_fraction(
  const double                  &T,
  const Parameters::PhaseChange &phase_change_parameters)
{
  return std::min(std::max((T - phase_change_parameters.T_solidus) /
                             (phase_change_parameters.T_liquidus -
                              phase_change_parameters.T_solidus),
                           0.),
                  1.);
}

#endif
