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

#ifndef lethe_density_model_h
#define lethe_density_model_h

/**
 * @brief DensityModel. Abstract class that allows to calculate the
 * density on each quadrature point using the temperature of the fluid.
 * DensityModel::get_density() is a pure virtual method,
 * since it can only be calculated knowing the model for the density has been
 * specified
 */
class DensityModel
{
public:
  /**
   * @brief Returns the density
   *
   * @param temperature Temperature
   */
  virtual double inline get_density(const double temperature) = 0;
};


/**
 * @brief Constant specific heat. Returns a constant specific
 * heat for a fluid
 */
class DensityConstant : public DensityModel
{
public:
  /**
   * @brief Default constructor
   */
  DensityConstant(const double p_density)
    : density(p_density)
  {}

  /**
   * @brief Returns the specific heat
   *
   * @param temperature Temperature at time t+dt
   *
   * @param previous_temperature Temperature at time t
   */
  virtual double inline get_density(const double /*temperature*/) override
  {
    return density;
  }

private:
  const double density;
};

#endif
