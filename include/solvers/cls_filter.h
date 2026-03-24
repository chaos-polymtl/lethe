// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_cls_filter_h
#define lethe_cls_filter_h


#include <core/parameters_multiphysics.h>

/**
 * @brief Abstract class for filtering phase fraction in Conservative Level Set (CLS)
 * simulations.
 */
class ConservativeLevelSetFilterBase
{
public:
  /**
   * @brief Default destructor.
   */
  virtual ~ConservativeLevelSetFilterBase() = default;

  /**
   * @brief Instantiate and return a pointer to a ConservativeLevelSetFilterBase
   * object by casting it to the proper child class.
   *
   * @param[in] phase_filter_parameters CLS filtration parameters.
   *
   * @return Casted ConservativeLevelSetFilterBase child class object.
   */
  static std::shared_ptr<ConservativeLevelSetFilterBase>
  model_cast(const Parameters::CLS_PhaseFilter &phase_filter_parameters);

  /**
   * @brief Compute the value of the filtered phase fraction.
   *
   * @param[in] unfiltered_phase Value of the phase fraction before applying the
   * filter.
   *
   * @return Value of the computed phase fraction after applying the filter.
   */
  virtual double
  filter_phase(const double &unfiltered_phase) = 0;
};

/**
 * @brief Default filter applied to phase fraction. The ConservativeLevelSetNoFilter
 * filter does not modify the phase fraction value and returns it as it is.
 */
class ConservativeLevelSetNoFilter : public ConservativeLevelSetFilterBase
{
public:
  /**
   * @brief Constructor of the transparent filter.
   *
   * No filter is applied to the computed filtered phase fraction and the
   * unfiltered phase fraction is returned.
   */
  ConservativeLevelSetNoFilter()
  {}

  /**
   * @brief Computes the value of the filtered phase fraction.
   *
   * Here, the filtered phase fraction value corresponds to the unfiltered one
   * (@p unfiltered_phase).
   *
   * @param[in] unfiltered_phase Value of the phase fraction before applying the
   * filter.
   *
   * @return Value of @p unfiltered_phase.
   */
  virtual double
  filter_phase(const double &unfiltered_phase) override
  {
    return unfiltered_phase;
  }
};

/**
 * @brief Filter phase fraction with a hyperbolic tangent function
 * for CLS simulations.
 *
 * The filtered phase is defined as:
 *
 * \f$\phi' = 0.5 \tanh(\beta(\phi-0.5)) + 0.5\f$
 *
 * where \f$\phi'\f$ is the filtered phase fraction value and \f$\beta\f$
 * is a model parameter that enables sharper definition when increased.
 *
 * The filter allows a sharper definition of the interface and clamps the value
 * of the phase fraction between 0 and 1.
 */
class ConservativeLevelSetTanhFilter : public ConservativeLevelSetFilterBase
{
public:
  /**
   * @brief Constructor of the hyperbolic tangent phase fraction filter.
   *
   * @param[in] beta Value of the \f$\beta\f$ parameter modulating the
   * interface sharpness.
   */
  ConservativeLevelSetTanhFilter(const double beta)
    : beta(beta)
  {}

  /**
   * @brief Computes the value of the filtered phase fraction.
   *
   * @param[in] unfiltered_phase Value of the phase fraction before applying the
   * filter.
   *
   * @return Value of the computed phase fraction after applying the filter.
   */
  virtual double
  filter_phase(const double &unfiltered_phase) override
  {
    return 0.5 * tanh(beta * (unfiltered_phase - 0.5)) + 0.5;
  }

private:
  /**
   * User-defined parameter that influences the sharpness of the interface
   * between fluids. A greater value leads to a sharper interface. The
   * recommended value is \f$\beta=20\f$.
   */
  const double beta;
};

#endif
