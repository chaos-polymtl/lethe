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
  /** @brief Class to account for different sharpening types:
   *  - constant: the sharpening threshold is the same throughout the
   * simulation,
   *  - adaptive: the sharpening threshold is determined by binary search, to
   * ensure mass conservation of the monitored phase
   */
  enum class SharpeningType
  {
    constant,
    adaptive
  };

  /** @brief Class to account for different phase fraction filtering types:
   * - none: no filter wil be applied on the calculated phase fraction
   * - tanh: the tanh filter function will be applied to the phase fraction,
   * a $$\beta$$ parameter influencing the interface definition must be defined
   */
  enum class FilterType
  {
    none,
    tanh
  };

  enum class EpsilonSetStrategy
  {
    automatic,
    manual
  };

  /**
   * @brief Defines the sub-parameters for free surface mass conservation.
   * Has to be declared before member creation in VOF structure.
   */
  struct VOF_MassConservation
  {
    bool monitoring;

    // Conservation tolerance on the fluid monitored,
    // used with adaptive Sharpening
    double tolerance;

    Parameters::FluidIndicator monitored_fluid;

    // Type of verbosity for the mass conservation algorithm
    Parameters::Verbosity verbosity;

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
    // https://lethe-cfd.github.io/lethe/examples/multiphysics/dam-break-VOF/dam-break-VOF.html

    bool enable;

    Parameters::SharpeningType type;

    // Parameters for constant sharpening
    double threshold;

    // Parameters for adaptive sharpening
    double threshold_max_deviation;
    int    max_iterations;

    // Other sharpening parameters
    double interface_sharpness;
    int    frequency;

    // Type of verbosity for the interface sharpening calculation
    Parameters::Verbosity verbosity;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief SurfaceTensionForce - Defines the parameters for
   * the calculation of surface tension force in the VOF solver.
   */
  struct VOF_SurfaceTensionForce
  {
    bool enable;

    double phase_fraction_gradient_diffusion_factor;
    double curvature_diffusion_factor;

    bool output_vof_auxiliary_fields;

    // Type of verbosity for the surface tension force calculation
    Parameters::Verbosity verbosity;

    // Enable marangoni effect
    bool enable_marangoni_effect;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief VOF_PhaseFilter - Defines the parameters for the phase filtration
   */
  struct VOF_PhaseFilter
  {
    // Type of filter
    Parameters::FilterType type;

    // $$\beta$$ value for the tanh filter
    double beta;

    // Type of verbosity for the phase filter
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
    Parameters::VOF_MassConservation    conservation;
    Parameters::VOF_InterfaceSharpening sharpening;
    Parameters::VOF_SurfaceTensionForce surface_tension_force;
    Parameters::VOF_PhaseFilter         phase_filter;

    Parameters::FluidIndicator viscous_dissipative_fluid;

    // artificial diffusivity (diffusion coefficient) (in L^2/s) added to the
    // VOF transport equation. This parameter is zero by default, and can be
    // increased to improve the wetting of the phases in the vicinity of
    // boundaries
    double diffusivity;

    bool compressible;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct CahnHilliard
  {
    // Well height value (W) in the Cahn-Hilliard equations
    double well_height;

    // Smoothing parameter Xi in the Cahn-Hilliard equations (potential
    // equation)
    double potential_smoothing_coefficient;

    // Epsilon set strategy (automatic|manual)
    Parameters::EpsilonSetStrategy epsilon_set_method;

    // Epsilon value in the Cahn-Hilliard equations
    double epsilon;

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
    bool cahn_hilliard;

    bool use_time_average_velocity_field;

    // subparameters for heat_transfer
    bool viscous_dissipation;
    bool buoyancy_force;

    Parameters::VOF          vof_parameters;
    Parameters::CahnHilliard cahn_hilliard_parameters;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };
} // namespace Parameters
#endif
