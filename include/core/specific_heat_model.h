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

#ifndef lethe_specific_heat_model_h
#define lethe_specific_heat_model_h

#include <core/parameters.h>

using namespace dealii;

/**
 * @brief SpecificHeatModel. Abstract class that allows to calculate the
 * specific heat on each quadrature point using the temperature of the fluid.
 * magnitude. SpecficiHeat::get_specific_heat() is a pure virtual method,
 * since it can only be calculated knowing the model for the specific
 * heat that has been specifid
 */
class SpecificHeatModel
{
public:
  /**
   * @brief Returns the specific heat
   *
   * @param temperature Temperature at time t+dt
   *
   * @param previous_temperature Temperature at time t
   */
  virtual double inline get_specific_heat(
    const double temperature,
    const double previous_temperature) = 0;
};


/**
 * @brief Constant specific heat. Returns a constant specific
 * heat for a fluid
 */
class SpecificHeatConstant : public SpecificHeatModel
{
public:
  /**
   * @brief Default constructor
   */
  SpecificHeatConstant(const double p_specific_heat)
    : specific_heat(p_specific_heat)
  {}

  /**
   * @brief Returns the specific heat
   *
   * @param temperature Temperature at time t+dt
   *
   * @param previous_temperature Temperature at time t
   */
  virtual double inline get_specific_heat(
    const double /*temperature*/,
    const double /*previous_temperature*/) override
  {
    return specific_heat;
  }

  const double specific_heat;
};


/**
 * @brief Constant specific heat. Returns a constant specific
 * heat for a fluid
 */
class PhaseChangeSpecificHeat : public SpecificHeatModel
{
public:
  /**
   * @brief Default constructor
   */
  PhaseChangeSpecificHeat(
    const Parameters::PhaseChange p_phase_change_parameters)
    : param(p_phase_change_parameters)
  {}

  /**
   * @brief Returns the specific heat
   *
   * @param temperature Temperature at time t+dt
   *
   * @param previous_temperature Temperature at time t
   */
  virtual double inline get_specific_heat(
    const double temperature,
    const double previous_temperature) override
  {
    if (temperature > previous_temperature)
      {
        const double dT = std::max(temperature - previous_temperature, 1e-6);
        return (enthalpy(temperature) - enthalpy(previous_temperature)) / dT;
      }
    else
      {
        const double dT = std::max(previous_temperature - temperature, 1e-6);
        return (enthalpy(previous_temperature) - enthalpy(temperature)) / dT;
      }
  }

  /**
   * @brief solid_fraction Calculates the solid fraction of a phase change material at
   * a temperature T
   *
   * @param T temperature at which to calculate the solid fraction
   *
   */
  inline double
  solid_fraction(const double T)
  {
    return std::min(std::max((T - param.T_solidus) /
                               (param.T_liquidus - param.T_solidus),
                             0.),
                    1.);
  }


  /**
   * @brief enthalpy Calculates the enthalpy of a phase change material for a temperature T
   * The enthalpy is defined as :
   *
   * !! Pure liquid !!
   * if (T>T_liquidus) : H = cp_solid * T_solidus + 0.5*(cp_solid+cp_liquid) *
   * (T_liquidus-T_solidus) + latent_enthalpy + cp_liquid * (T-T_liquidus)
   *
   * !! Liquid-solid mix !!
   * else if (T>T_solidus) : cp_solid * T_solidus + 0.5*(cp_solid+cp_liquid) *
   * (T-T_solidus) + liquid_fraction * latent_enthalpy
   *
   * !! Pure solid !!
   * else : cp_solid * T
   *
   * @param T temperature at which to calculate the enthalpy
   * @return Value of the enthalpy
   */

  inline double
  enthalpy(const double T)
  {
    if (T > param.T_liquidus)
      return (param.cp_s * param.T_solidus +
              0.5 * (param.cp_l + param.cp_s) *
                (param.T_liquidus - param.T_solidus) +
              param.latent_enthalpy + (T - param.T_liquidus) * param.cp_l);
    else if (T > param.T_solidus)
      {
        const double s_frac = solid_fraction(T);
        return (param.cp_s * param.T_solidus +
                0.5 * (param.cp_l + param.cp_s) * (T - param.T_solidus) +
                param.latent_enthalpy * s_frac);
      }
    else
      return param.cp_s * T;
  }

private:
  const Parameters::PhaseChange param;
};


#endif
