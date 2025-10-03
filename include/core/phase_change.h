// SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_phase_change_h
#define lethe_phase_change_h

#include <core/parameters.h>

using namespace dealii;

/**
 * @brief calculate_liquid_fraction Calculates the liquid fraction of a phase change material at
 * a temperature T
 *
 * @param T temperature at which to calculate the liquid fraction
 * @param phase_change_parameters
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

/**
 * @brief calculate_liquid_fraction Calculates the liquid fraction of a phase change material at
 * a temperature T
 *
 * @param T temperature at which to calculate the liquid fraction
 * @param T_solidus temperature
 * @param T_liquidus temperature
 * @return value of the liquid_fraction
 *
 */
inline double
calculate_liquid_fraction(const double &T,
                          const double  T_solidus,
                          const double  T_liquidus)
{
  return std::min(std::max((T - T_solidus) / (T_liquidus - T_solidus), 0.), 1.);
}

#endif
