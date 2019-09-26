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
 * Element required for the calculation of analytical solution
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
  public:
    AnalyticalSolution()
      : enable(false)
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

  protected:
    bool        enable;
    std::string filename;
  };

  template <int dim>
  void
  AnalyticalSolution<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("analytical solution");
    prm.declare_entry(
      "enable",
      "false",
      Patterns::Bool(),
      "Enable the calculation of the analytical solution and L2 error");
    prm.declare_entry(
      "verbosity",
      "quiet",
      Patterns::Selection("quiet|verbose"),
      "State whether from the post-processing values should be printed "
      "Choices are <quiet|verbose>.");

    prm.declare_entry(
      "filename",
      "L2Error",
      Patterns::FileName(),
      "File name for the output for the L2Error table with respect to time or mesh ");
    prm.leave_subsection();
  }

  template <int dim>
  void
  AnalyticalSolution<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("analytical solution");
    enable               = prm.get_bool("enable");
    filename             = prm.get("filename");
    const std::string op = prm.get("verbosity");
    if (op == "verbose")
      verbosity = Parameters::verbose;
    if (op == "quiet")
      verbosity = Parameters::quiet;
    prm.leave_subsection();
  }

  template <int dim>
  class NSAnalyticalSolution : public AnalyticalSolution<dim>
  {
  public:
    NSAnalyticalSolution()
      : velocity(dim + 1)
    {}

    // Velocity components
    Functions::ParsedFunction<dim> velocity;

    virtual void
    declare_parameters(ParameterHandler &prm);
    virtual void
    parse_parameters(ParameterHandler &prm);
  };

  template <int dim>
  void
  NSAnalyticalSolution<dim>::declare_parameters(ParameterHandler &prm)
  {
    this->AnalyticalSolution<dim>::declare_parameters(prm);
    prm.enter_subsection("analytical solution");
    prm.enter_subsection("uvw");
    velocity.declare_parameters(prm, dim);
    if (dim == 2)
      prm.set("Function expression", "0; 0; 0;");
    if (dim == 3)
      prm.set("Function expression", "0; 0; 0; 0;");
    prm.leave_subsection();
    prm.leave_subsection();
  }

  template <int dim>
  void
  NSAnalyticalSolution<dim>::parse_parameters(ParameterHandler &prm)
  {
    this->AnalyticalSolution<dim>::parse_parameters(prm);
    prm.enter_subsection("analytical solution");
    prm.enter_subsection("uvw");
    velocity.parse_parameters(prm);
    prm.leave_subsection();
    prm.leave_subsection();
  }
} // namespace AnalyticalSolutions

#endif
