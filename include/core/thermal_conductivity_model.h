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

#ifndef lethe_thermal_conductivity_model_h
#define lethe_thermal_conductivity_model_h

#include <core/parameters.h>

/**
 * @brief ThermalConductivityModel. Abstract class that allows to calculate the
 * thermal conductivity on each quadrature point using the temperature of the
 * fluid.
 */
class ThermalConductivityModel
{
public:
  /**
   * @brief Returns the the thermal conductivity
   */
  inline virtual double
  get_thermal_conductivity(const double temperature) = 0;
};


/**
 * @brief Constant specific heat. Returns a constant specific
 * heat for a fluid
 */
class ThermalConductivityConstant : public ThermalConductivityModel
{
public:
  /**
   * @brief Default constructor
   */
  ThermalConductivityConstant(const double p_thermal_conductivity)
    : thermal_conductivity(p_thermal_conductivity)
  {}

  /**
   * @brief Returns the specific heat
   *
   * @param temperature Temperature at time t+dt
   *
   * @param previous_temperature Temperature at time t
   */
  inline virtual double
  get_thermal_conductivity(const double /*temperature*/) override
  {
    return thermal_conductivity;
  }

private:
  const double thermal_conductivity;
};

#endif
