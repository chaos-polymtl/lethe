/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
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

#ifndef lethe_cahn_hilliard_filter_h
#define lethe_cahn_hilliard_filter_h


#include <core/parameters_multiphysics.h>

/**
 * @brief CahnHilliardFilterBase class allows phase fraction filtering
 */
class CahnHilliardFilterBase
{
public:
  CahnHilliardFilterBase()
  {}

  /**
   * @brief Instantiates and returns a pointer to a CahnHilliardFilterBase
   * object by casting it to the proper child class
   *
   * @param phase_filter_parameters CahnHilliard model parameters
   */
  static std::shared_ptr<CahnHilliardFilterBase>
  model_cast(const Parameters::CahnHilliard &cahn_hilliard_parameters);

  /**
   * @brief filter_phase calculates the value of the filtered phase fraction
   * @param unfiltered_phase Value of the phase fraction before applying the filter
   * @return Value of the phase fraction after applying the filter
   */
  virtual double
  filter_phase(const double &unfiltered_phase) = 0;
};

/**
 * @brief CahnHilliardNoFilter class is used as a default when no
 * filter is applied to the phase fraction
 */
class CahnHilliardNoFilter : public CahnHilliardFilterBase
{
public:
  CahnHilliardNoFilter()
  {}

  /**
   * @brief CahnHilliardNoFilter class is used as a default when no
   * filter is applied to the phase fraction. In this case, a simple clamping is
   * performed on the phase fraction parameter for it to remain in the [-1,1]
   * interval.
   */
  virtual double
  filter_phase(const double &unfiltered_phase) override
  {
    return (std::abs(unfiltered_phase) < 1) ? unfiltered_phase :
                                              sgn(unfiltered_phase);
  }
};

/**
 * @brief CahnHilliardTanhFilter class is used to calculate
 * a filtered phase for Cahn-Hilliard simulations. The filtered phase is defined
 * as
 * $$\phi_f = \tanh(\beta(\phi-0.5)) $$.
 */
class CahnHilliardTanhFilter : public CahnHilliardFilterBase
{
public:
  CahnHilliardTanhFilter(const double beta)
    : beta(beta)
  {}

  /**
   * @brief filter_phase calculates the value of the filtered phase fraction
   * @param unfiltered_phase Value of the phase fraction before applying the filter
   * @return Value of the phase fraction after applying the tanh filter
   */
  virtual double
  filter_phase(const double &unfiltered_phase) override
  {
    return tanh(beta * unfiltered_phase);
  }

private:
  // User-defined parameter that influences how sharp the filtering is.
  const double beta;
};

#endif