/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

/*
 * This file defines the parameters in the parameter namespace
 * that pertain to multiphysics simulations
 */


#ifndef lethe_parameters_multiphysics_h
#define lethe_parameters_multiphysics_h

#include <core/parameters.h>

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

namespace Parameters
{
  /**
   * @brief Defines the subparameters for free surface monitoring.
   * Has to be declared before member creation in VOF structure.
   */
  struct VOF_MassConservation
  {
    bool skip_mass_conservation_fluid_0;
    bool skip_mass_conservation_fluid_1;
    bool monitoring;
    int  id_fluid_monitored;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief VOF_InterfaceSharpening - Defines the parameters for
   * interface sharpening in the VOF solver.
   */
  struct VOF_InterfaceSharpening
  {
    // Interface sharpening parameters. The sharpening method and parameters are
    // explained in the dam break VOF example:
    // https://github.com/lethe-cfd/lethe/wiki/Dam-break-VOF
    // sharpening_threshold is the phase fraction threshold for sharpening. It
    // should be chosen in the range of (0,1), but generally it is equal to 0.5
    // interface_sharpness is a parameter which defines the sharpness of the
    // interface. It should be chosen in the range of (1,2] sharpening_frequency
    // (integer) is the frequency at which the interface sharpneing is called.
    // Users may set this variable to 1 to call interface sharpening at every
    // step, but it could be chosen in the range of [1-20]

    bool interface_sharpening;

    double sharpening_threshold;
    double interface_sharpness;
    int    sharpening_frequency;
    // Type of verbosity for the interface sharpening calculation
    Parameters::Verbosity verbosity;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief VOF - Defines the parameters for free surface simulations
   * using the VOF method.
   * Has to be declared before member creation in Multiphysics structure.
   */
  struct VOF
  {
    bool   continuum_surface_force;
    bool   peeling_wetting;
    double diffusion;

    Parameters::VOF_MassConservation    conservation;
    Parameters::VOF_InterfaceSharpening sharpening;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Multiphysics - the parameters for multiphysics simulations
   * and handles sub-physics parameters.
   */
  struct Multiphysics
  {
    bool fluid_dynamics;
    bool heat_transfer;
    bool tracer;
    bool VOF;

    // subparameters for heat_transfer
    bool viscous_dissipation;
    bool buoyancy_force;

    Parameters::VOF vof_parameters;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };
} // namespace Parameters
#endif
