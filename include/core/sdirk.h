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

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019 -
 */


#ifndef lethe_sdirk_h
#define lethe_sdirk_h

#include <core/parameters.h>

#include <deal.II/lac/full_matrix.h>

#include <vector>

using namespace dealii;

/**
 * @brief Calculates the SDIRK coefficient
 *
 * @param order Order of the SDIRK method to be used
 *
 * @param time_step Value of the time step used for time integration.
 *
 * SDIRK 22 -  Coefficients from Kennedy and Carpenter 2016
 * Intermediary step is a t+ (1-sqrt(2)/2) * dt
 *
 * #SDIRK 33 - Coefficients from Kennedy and Carpenter 2016
 * gamma = 0.435866521508458999416019
 * b = 1.20849664917601007033648
 * c =0.717933260754229499708010
 *
 * Butcher's Tableau for SDIRK33
 * --------------------------
 * | gamma      0        0  |
 * |c-gamma   gamma      0  |
 * |   b    1-b-gamma  gamma|
 * --------------------------
 *
 * Intermediary 1 step is at
 * Intermediary 2 step is at
 *
 */
FullMatrix<double>
sdirk_coefficients(const unsigned int order, const double time_step);


/**
 * @brief Returns the number of previous time steps supposed by the BDF schemes implemented in Lethe.
 *  At the moment this is hardcoded to 3, but eventually this could be made
 * larger or smaller depending on the methods used.
 *
 * @param Time-stepping scheme used
 *
 */
inline unsigned int
number_of_intermediary_stages(
  Parameters::SimulationControl::TimeSteppingMethod method)
{
  unsigned int n_stages = 0;

  if (method == Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1 ||
      method == Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2 ||
      method == Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3 ||
      method == Parameters::SimulationControl::TimeSteppingMethod::sdirk33||
      method == Parameters::SimulationControl::TimeSteppingMethod::bdf3
      )
    n_stages = 2;

  if (method == Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1 ||
      method == Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2 ||
      method == Parameters::SimulationControl::TimeSteppingMethod::sdirk22||
      method == Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    n_stages = 1;

  return n_stages;
}

inline unsigned int
intermediary_stages(Parameters::SimulationControl::TimeSteppingMethod method)
{
  if (method == Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1 ||
      method == Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1)
    return 0;

  else if (method ==
             Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2 ||
           method ==
             Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2)
    return 1;

  else // if (method ==
       // Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3)
    return 2;
}

inline unsigned int
max_number_of_intermediary_stages()
{
  return 3;
}


/**
 * @brief Determines if this is the first step of an sdirk method
 *
 * @param method A time integration method
 */
inline bool
is_sdirk_step1(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1);
}

/**
 * @brief Determines if this is the second step of an sdirk method
 *
 * @param method A time integration method
 */
inline bool
is_sdirk_step2(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2 ||
    method == Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2);
}

/**
 * @brief Determines if this is the third step of an sdirk method
 *
 * @param method A time integration method
 */
inline bool
is_sdirk_step3(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  return (method ==
          Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3);
}


inline unsigned int
sdirk_step(Parameters::SimulationControl::TimeSteppingMethod method)
{
  if (is_sdirk_step1(method))
    return 1;
  else if (is_sdirk_step2(method))
    return 2;
  else if (is_sdirk_step3(method))
    return 3;
  else
    return 0;
}



#endif
