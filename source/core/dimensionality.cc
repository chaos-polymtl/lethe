// SPDX-FileCopyrightText: Copyright (c) 2022, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dimensionality.h>


using namespace dealii;
namespace Parameters
{
  void
  Dimensionality::define_all_scales()
  {
    // Define fundamental quantities using Wikipedia's notation
    const double L     = length;
    const double M     = mass;
    const double theta = temperature;
    const double T     = time;

    density_scaling                  = 1. * L * L * L / M;
    specific_gas_constant_scaling    = 1. / L / L * T * T * theta;
    viscosity_scaling                = 1. / L / L * T;
    specific_heat_scaling            = 1. / L / L * T * T * theta;
    thermal_conductivity_scaling     = 1. / M / L * T * T * T * theta;
    enthalpy_scaling                 = 1. / M / L / L * T * T;
    diffusivity_scaling              = 1. / L / L * T;
    thermal_expansion_scaling        = 1. * theta;
    surface_tension_scaling          = 1. * T * T / M;
    surface_tension_gradient_scaling = 1. * theta * T * T / M;
    cahn_hilliard_mobility_scaling   = 1. * M / L / L / L / T;
    cahn_hilliard_epsilon_scaling    = 1. / L;
  }

  void
  Dimensionality::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("dimensionality");
    {
      prm.declare_entry(
        "length",
        "1",
        Patterns::Double(),
        "Length scale. The default value assumed is meters. If the simulation is carried out in centimeters,"
        " a length of 0.01 should be specified. If the simulation is carried out in kilometers, a length of 1000 should be specified.");

      prm.declare_entry(
        "time",
        "1",
        Patterns::Double(),
        "Time scale. The default value assumed is seconds. If the simulation is carried out in milliseconds,"
        " a time of 0.001 should be specified.");

      prm.declare_entry(
        "mass",
        "1",
        Patterns::Double(),
        "Mass scale. The default value assumed is kilograms. If the simulation is carried out in grams,"
        " a mass of 0.001 should be specified.");

      prm.declare_entry(
        "temperature",
        "1",
        Patterns::Double(),
        "Temperature scale. The default value assumed is Kelvin. If the simulation is carried out in kiloKelvin,"
        " a temperature of 1000 should be specified.");
    }
    prm.leave_subsection();
  }


  void
  Dimensionality::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("dimensionality");
    {
      length      = prm.get_double("length");
      time        = prm.get_double("time");
      mass        = prm.get_double("mass");
      temperature = prm.get_double("temperature");
    }
    prm.leave_subsection();

    // Define all scales using these dimensions
    define_all_scales();
  }
} // namespace Parameters
