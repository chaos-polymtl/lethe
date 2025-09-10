// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_filter_h
#define lethe_vof_filter_h


#include <core/parameters_multiphysics.h>

/**
 * @brief Abstract class for filtering phase fraction in volume of fluid (VOF)
 * simulations.
 */
class VolumeOfFluidFilterBase
{
public:
  /**
   * @brief Default destructor.
   */
  virtual ~VolumeOfFluidFilterBase() = default;

  /**
   * @brief Instantiate and return a pointer to a VolumeOfFluidFilterBase
   * object by casting it to the proper child class.
   *
   * @param[in] phase_filter_parameters VOF filtration parameters.
   *
   * @return Casted VolumeOfFluidFilterBase child class object.
   */
  static std::shared_ptr<VolumeOfFluidFilterBase>
  model_cast(const Parameters::VOF_PhaseFilter &phase_filter_parameters);

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
 * @brief Default filter applied to phase fraction. The VolumeOfFluidNoFilter
 * filter does not modify the phase fraction value and returns it as it is.
 */
class VolumeOfFluidNoFilter : public VolumeOfFluidFilterBase
{
public:
  /**
   * @brief Constructor of the transparent filter.
   *
   * No filter is applied to the computed filtered phase fraction and the
   * unfiltered phase fraction is returned.
   */
  VolumeOfFluidNoFilter()
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
 * for VOF simulations.
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
class VolumeOfFluidTanhFilter : public VolumeOfFluidFilterBase
{
public:
  /**
   * @brief Constructor of the hyperbolic tangent phase fraction filter.
   *
   * @param[in] beta Value of the \f$\beta\f$ parameter modulating the
   * interface sharpness.
   */
  VolumeOfFluidTanhFilter(const double beta)
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
