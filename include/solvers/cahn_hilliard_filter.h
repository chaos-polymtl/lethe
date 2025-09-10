// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_cahn_hilliard_filter_h
#define lethe_cahn_hilliard_filter_h

#include <core/parameters_multiphysics.h>

/**
 * @brief Filter phase fraction according to the selected filtering method.
 */
class CahnHilliardFilterBase
{
public:
  CahnHilliardFilterBase()
  {}

  /**
   * @brief Default destructor.
   */
  virtual ~CahnHilliardFilterBase() = default;

  /**
   * @brief Instantiate and return a pointer to a CahnHilliardFilterBase
   * object by casting it to the proper child class.
   *
   * @param[in] cahn_hilliard_parameters CahnHilliard model parameters.
   */
  static std::shared_ptr<CahnHilliardFilterBase>
  model_cast(const Parameters::CahnHilliard &cahn_hilliard_parameters);

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
class CahnHilliardNoFilter : public CahnHilliardFilterBase
{
public:
  CahnHilliardNoFilter()
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
 * @brief Calculate a filtered phase fraction for Cahn-Hilliard simulations.
 * In this case, a simple clamping is performed on the phase fraction parameter
 * for it to remain in the [-1,1] interval.
 */
class CahnHilliardClipFilter : public CahnHilliardFilterBase
{
public:
  CahnHilliardClipFilter()
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
 * @brief Calculate a filtered phase fraction for Cahn-Hilliard simulations.
 *
 * The filtered phase is defined as \f$\phi_f = \tanh(\beta \phi)\f$ with
 * \f$ \phi \f$ the unfiltered phase field.
 */
class CahnHilliardTanhFilter : public CahnHilliardFilterBase
{
public:
  CahnHilliardTanhFilter(const double beta)
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
