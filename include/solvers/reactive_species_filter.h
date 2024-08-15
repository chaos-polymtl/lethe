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

#ifndef lethe_reactive_species_filter_h
#define lethe_reactive_species_filter_h

#include <core/parameters_multiphysics.h>

/**
 * @brief Filter phase fraction according to the selected filtering method.
 */
class ReactiveSpeciesFilterBase
{
public:
  ReactiveSpeciesFilterBase()
  {}

  /**
   * @brief Instantiate and return a pointer to a ReactiveSpeciesFilterBase
   * object by casting it to the proper child class.
   *
   * @param[in] reactive_species_parameters ReactiveSpecies model parameters.
   */
  static std::shared_ptr<ReactiveSpeciesFilterBase>
  model_cast(const Parameters::ReactiveSpecies &reactive_species_parameters);

  /**
   * @brief Calculate the value of the filtered phase fraction.
   * @param[in] unfiltered_phase Value of the phase fraction before applying
   * the filter.
   * @return Value of the phase fraction after applying the filter.
   */
  virtual double
  filter_phase(const double &unfiltered_phase) = 0;
};

/**
 * @brief Used as a default when no filter is applied to the phase fraction.
 *
 */
class ReactiveSpeciesNoFilter : public ReactiveSpeciesFilterBase
{
public:
  ReactiveSpeciesNoFilter()
  {}

  /**
   * @brief Return the phase fraction with no modification.
   */
  virtual double
  filter_phase(const double &unfiltered_phase) override
  {
    return unfiltered_phase;
  }
};

/**
 * @brief Calculate a filtered phase fraction for Reactive species simulations.
 * In this case, a simple clamping is performed on the phase fraction parameter
 * for it to remain in the [-1,1] interval.
 */
class ReactiveSpeciesClipFilter : public ReactiveSpeciesFilterBase
{
public:
  ReactiveSpeciesClipFilter()
  {}

  /**
   * @brief Calculate the value of the filtered phase fraction.
   * @param[in] unfiltered_phase Value of the phase fraction before applying the
   * filter.
   * @return Value of the phase fraction after applying the filter.
   */
  virtual double
  filter_phase(const double &unfiltered_phase) override
  {
    return (std::abs(unfiltered_phase) < 1) ? unfiltered_phase :
                                              sgn(unfiltered_phase);
  }
};

/**
 * @brief Calculate a filtered phase fraction for Reactive species simulations.
 *
 * The filtered phase is defined as \f$\phi_f = \tanh(\beta \phi)\f$ with
 * \f$ \phi \f$ the unfiltered phase field.
 */
class ReactiveSpeciesTanhFilter : public ReactiveSpeciesFilterBase
{
public:
  ReactiveSpeciesTanhFilter(const double beta)
    : beta(beta)
  {}

  /**
   * @brief Calculate the value of the filtered phase fraction.
   * @param[in] unfiltered_phase Value of the phase fraction before applying the
   * filter.
   * @return Value of the phase fraction after applying the filter.
   */
  virtual double
  filter_phase(const double &unfiltered_phase) override
  {
    return tanh(beta * unfiltered_phase);
  }

private:
  /**
   * User-defined parameter that influences how sharp the filtering is.
   */
  const double beta;
};

#endif
