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

#ifndef lethe_analytical_solutions_h
#define lethe_analytical_solutions_h

#include <core/parameters.h>

#include <deal.II/base/function.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;


/**
 * @brief Contains all the classes related to the calculation of the analytical solutions
 * for the various physics
 **/

namespace AnalyticalSolutions
{

  /**
  * @brief Implements analytical solutions for the Navier-Stokes equations
  * and the other physics supported by lethe
  *
  * @tparam dim An integer that denotes the dimension of the space in which
  * the flow is solved.

  * The analytical solution class provides an all common
  * element required for the calculation of analytical solution for all the
  * physics supported by Lethe.
  **/
  template <int dim>
  class AnalyticalSolution
  {
  private:
    /**
     * Establishes if the calculation of the analytical solution is enabled or
     * not
     */
    bool enable;

    /**
     * Filename used to store the L2 norm of the error
     */
    std::string filename;


    /**
     * @brief Construct a new AnalyticalSolution object.
     *
     * The constructor automatically sets the number of components
     * of each equation manually since all equations currently supported have
     * a set number of equations.
     */
  public:
    AnalyticalSolution()
      : enable(false)
      , uvwp(dim + 1)
      , temperature(1)
      , tracer(1)
      , cahn_hilliard(2)
    {}

    /**
     * @brief Declares the parameters required by the analytical solution
     * within a parameter file
     *
     * @param prm ParameterHandler used to declare the parameters.
     *
     */
    virtual void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Declares the parameters required by the analytical solution
     * within a parameter file
     *
     * @param prm ParameterHandler used to parse the parameters.
     *
     */
    virtual void
    parse_parameters(ParameterHandler &prm);

    /**
     * @brief Checks if the calculation of the error using the analytical
     * solution is enabled. This option is there to prevent unnecessary
     * calculations when no analytical solutions are required.
     */
    bool
    calculate_error()
    {
      return enable;
    }


    /**
     * @brief Enables the usage of the analytical solution
     */
    void
    set_enable(bool is_enable)
    {
      enable = is_enable;
    }

    /**
     * @brief Get the file name associated with the output of the
     * L2 norm of the error calculated using the analytical solution.
     *
     * @return filename A string that contains the filename used to output
     * the L2 norm of the error.
     */
    std::string
    get_filename()
    {
      return filename;
    }


    /**
     * Controls if the L2 norm of the error is printed to the terminal
     * during the simulation
     */
    Parameters::Verbosity verbosity;

    /**
     * ParsedFunction that contains  the analytical solution for the velocity
     * components and the pressure
     */
    Functions::ParsedFunction<dim> uvwp;

    /**
     * ParsedFunction that contains the analytical solution for the temperature
     */
    Functions::ParsedFunction<dim> temperature;

    /**
     * ParsedFunction that contains the analytical solution for the tracer
     */
    Functions::ParsedFunction<dim> tracer;

    /**
     * ParsedFunction that contains the analytical solution for the phase
     */
    Functions::ParsedFunction<dim> phase;

    /**
     * ParsedFunction that contains the analytical solution for the phase order
     * and the chemical potential
     */
    Functions::ParsedFunction<dim> cahn_hilliard;
  };
} // namespace AnalyticalSolutions

#endif
