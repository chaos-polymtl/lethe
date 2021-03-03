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

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019 -
 */

// TODO : Refactor so the class itself is not a pointer, but contains a pointer
//        to a function. This would be a lot more coherent...

#ifndef LETHE_ANALYTICALSOLUTIONS_H
#define LETHE_ANALYTICALSOLUTIONS_H

#include <deal.II/base/function.h>
#include <deal.II/base/parsed_function.h>

#include <core/parameters.h>

using namespace dealii;
/**
 * The analytical solution class provides an interface for all common
 * element required for the calculation of analytical solution
 * All equation-specific analytical solution should derive
 * from the base class but also call it's declare_parameters and
 *parse_parameters routine. This allows specialize class to focus on their
 *specificity and forget about other non-specific elements that are generic to
 *the calculation of analytical solutions
 **/

namespace AnalyticalSolutions
{
  template <int dim>
  class AnalyticalSolution
  {
  protected:
    bool        enable;
    std::string filename;

  public:
    AnalyticalSolution()
      : enable(false)
      , velocity(dim + 1)
      , temperature(1)
      , tracer(1)
    {}

    virtual void
    declare_parameters(ParameterHandler &prm);
    virtual void
    parse_parameters(ParameterHandler &prm);

    bool
    calculate_error()
    {
      return enable;
    }

    std::string
    get_filename()
    {
      return filename;
    }

    Parameters::Verbosity verbosity;
    // Velocity components
    Functions::ParsedFunction<dim> velocity;

    Functions::ParsedFunction<dim> temperature;

    Functions::ParsedFunction<dim> tracer;
  };
} // namespace AnalyticalSolutions

#endif
