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

namespace Parameters
{
  /**
   * @brief Class used to rescale physical properties for cases where the main dimensions
   * of the problems (Length, Mass, Time and Temperature) are not expressed in
   * SI, but are expressed in another set of units. For example, this can be
   * used to carry out a simulation where the mesh is expressed in millimeters
   * instead of meters. This class is currently only used to rescale the
   * physical properties in the FEM-based physics (e.g. all of those that are
   * solvers). Support for DEM, CFD-DEM and rescaling of boundary conditions is
   * planned in the mid-term.
   */
  class Dimensionality
  {
  public:
    /**
     * @brief
     * Default constructor that assumes that the units used in the simulation
     * are consistent.
     */
    Dimensionality()
      : length(1)
      , mass(1)
      , time(1)
      , temperature(1)
    {
      define_all_scales();
    }

    /**
     * @brief
     * Defines all the sub scales that are relevant to simulations. This is used
     * to define scaling laws for kinematic viscosity or other quantities which
     * appear multiple times in the simulation. This prevents the recalculation
     * of these scales (and thus, human errors).
     */
    void
    define_all_scales();

    /**
     * @brief
     * Declares the fundamental dimensions of the parameter file
     *
     * @param prm The parameter handler being used
     */
    void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief
     * Parses the fundamental dimensions from the parameter file
     *
     * @param prm The parameter handler being used
     */
    void
    parse_parameters(ParameterHandler &prm);

    /*
     * Fundamental dimensions. These are the dimensions that are read from the
     * parameter file
     */
    double length;
    double mass;
    double time;
    double temperature;

    /*
     * Scaling laws. These are constants that are used to rescale quantities
     * arising frequently. Some of these quantities are not trivial to calculate
     * (e.g. thermal conductivity). Thus, it is better to pre-calculate these
     * scalings and to reuse them instead of calculating them on the fly
     * everywhere.
     */
    double density_scaling;
    double specific_gas_constant_scaling;
    double viscosity_scaling;
    double specific_heat_scaling;
    double thermal_conductivity_scaling;
    double enthalpy_scaling;
    double diffusivity_scaling;
    double thermal_expansion_scaling;
  };
} // namespace Parameters

#endif
