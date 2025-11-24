// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_time_integration_utilities_h
#define lethe_time_integration_utilities_h

#include <core/parameters.h>

/**
 * @brief Determines if the time integration method is a steady-state method
 *
 * @param method A time integration method
 */
inline bool
is_steady(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method == Parameters::SimulationControl::TimeSteppingMethod::steady);
}

/**
 * @brief Determines if the time integration method is within the bdf-1 family
 *
 * @param method A time integration method
 */
inline bool
time_stepping_is_bdf1(
  const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method ==
            Parameters::SimulationControl::TimeSteppingMethod::steady_bdf ||
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf1);
}

/**
 * @brief Determines if the time integration method is within the bdf family
 *
 * @param method A time integration method
 *
 * @return true if the method is BDF1, BDF2, BDF3, or steady BDF, false otherwise.
 */
inline bool
time_stepping_is_bdf(
  const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method ==
            Parameters::SimulationControl::TimeSteppingMethod::steady_bdf ||
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf2 ||
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf3);
}

/**
 * @brief Determines if the time integration method is within the sdirk family
 *
 * @param method A time integration method
 *
 * @return true if the method is SDIRK, false otherwise.
 */
inline bool
time_stepping_is_sdirk(
  const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk22 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk33 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk43);
}

/**
 * @brief Determines if the time integration method is within the high-order (>=2) bdf family
 *
 * @param method A time integration method
 */
inline bool
time_stepping_is_bdf_high_order(
  const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method == Parameters::SimulationControl::TimeSteppingMethod::bdf2 ||
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf3);
}


/**
 * @brief Determines if the time integration method requires an additional array
 *
 * @param method A time integration method
 */
inline bool
time_stepping_method_uses_two_previous_solutions(
  const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method == Parameters::SimulationControl::TimeSteppingMethod::bdf2 ||
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf3);
}

/**
 * @brief Determines if the time integration method requires two additional arrays
 *
 * @param method A time integration method
 */
inline bool
time_stepping_method_uses_three_previous_solutions(
  const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method == Parameters::SimulationControl::TimeSteppingMethod::bdf3);
}

#endif
