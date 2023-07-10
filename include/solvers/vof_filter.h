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

#ifndef lethe_vof_filter_h
#define lethe_vof_filter_h


#include <core/parameters_multiphysics.h>

/**
 * @brief VolumeOfFluidFilterBase class allows phase fraction filtering
 */
class VolumeOfFluidFilterBase
{
public:
  VolumeOfFluidFilterBase()
  {}

  /**
   * @brief Instantiates and returns a pointer to a VolumeOfFluidFilterBase
   * object by casting it to the proper child class
   *
   * @param phase_filter_parameters VOF model parameters
   */
  static std::shared_ptr<VolumeOfFluidFilterBase>
  model_cast(const Parameters::VOF_PhaseFilter &phase_filter_parameters);

  /**
   * @brief filter_phase calculates the value of the filtered phase fraction
   * @param unfiltered_phase Value of the phase fraction before applying the filter
   * @return Value of the phase fraction after applying the filter
   */
  virtual double
  filter_phase(const double &unfiltered_phase) = 0;
};

/**
 * @brief VolumeOfFluidNoFilter class is used as a default when no
 * filter is applied to the phase fraction
 */
class VolumeOfFluidNoFilter : public VolumeOfFluidFilterBase
{
public:
  VolumeOfFluidNoFilter()
  {}

  /**
   * @brief filter_phase calculates the value of the filtered phase fraction
   * @param unfiltered_phase Value of the phase fraction before applying the filter
   * @return unfiltered_phase
   */
  virtual double
  filter_phase(const double &unfiltered_phase) override
  {
    return unfiltered_phase;
  }
};

/**
 * @brief VolumeOfFluidTanhFilter class is used to calculate
 * a filtered phase for VOF simulations. The filtered phase is defined as
 * $$\phi_f = 0.5 \tanh(\beta(\phi-0.5)) + 0.5$$.
 */
class VolumeOfFluidTanhFilter : public VolumeOfFluidFilterBase
{
public:
  VolumeOfFluidTanhFilter(const double beta)
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
    return 0.5 * tanh(beta * (unfiltered_phase - 0.5)) + 0.5;
  }

private:
  // User-defined parameter that influences the definition of the interface
  const double beta;
};

#endif
