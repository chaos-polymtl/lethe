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

#ifndef lethe_dimensionality_h
#define lethe_dimensionality_h

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

/**
 * @brief Dimensionality - This is used to rescale physical properties and other elements for cases where the main dimensions
 * of the problems (Length, Mass, Time) are not expressed in SI, but are
 * expressed in another set of unit. For example, this can be used to carry
 * out a simulation where the mesh is expressed in millimeter instead of
 * meters
 */
namespace Parameters
{
  class Dimensionality
  {
  public:
    Dimensionality()
      : length(1)
      , mass(1)
      , time(1)
      , temperature(1)
    {
      define_all_scales();
    }

    void
    define_all_scales()
    {
      // Define fundamental quantities using Wikipedia's notation
      const double L     = length;
      const double M     = mass;
      const double theta = temperature;
      const double T     = time;

      density_scaling              = 1 * L * L * L / M;
      viscosity_scaling            = 1 * L * L / T;
      specific_heat_scaling        = 1. / L / L * T * T * theta;
      thermal_conductivity_scaling = 1. / M / L * T * T * T * theta;
      enthalpy_scaling             = 1. / M / L * T * T * T;
      diffusivity_scaling          = 1. / L / L * T;
    }

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);

    double density_scaling;
    double viscosity_scaling;
    double specific_heat_scaling;
    double thermal_conductivity_scaling;
    double enthalpy_scaling;
    double diffusivity_scaling;


    double length;
    double mass;
    double time;
    double temperature;
  };
} // namespace Parameters

#endif
