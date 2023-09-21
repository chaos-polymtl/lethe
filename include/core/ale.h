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
#ifndef lethe_ale_h
#define lethe_ale_h

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <deal.II/lac/vector.h>

using namespace dealii;
/**
 * The source term class provides an interface for all common
 * Element required for the introduction of a source term in equations
 * All equation-specific source term should derive
 * from the base class but also call it's declare_parameters and
 *parse_parameters routine. This allows specialize class to focus on their
 *specificity and forget about other non-specific elements that are generic to
 *the calculation of analytical solutions
 **/

namespace Parameters
{
  template <int dim>
  class ALE
  {
  public:
    ALE()
    {
      velocity = std::make_shared<Functions::ParsedFunction<dim>>(dim);
    }

    virtual void
    declare_parameters(ParameterHandler &prm);
    virtual void
    parse_parameters(ParameterHandler &prm);


    // ALE-Velocity function
    std::shared_ptr<Functions::ParsedFunction<dim>> velocity;

    bool
    enabled()
    {
      return enable;
    }

  private:
    bool enable;
  };

  template <int dim>
  void
  ALE<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("ALE");
    prm.declare_entry(
      "enable",
      "false",
      Patterns::Bool(),
      "Enable simulations in a Galilean moving frame of reference");

    prm.enter_subsection("velocity");
    velocity->declare_parameters(prm, dim);
    prm.leave_subsection();

    prm.leave_subsection();
  }

  template <int dim>
  void
  ALE<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("ALE");
    enable = prm.get_bool("enable");

    prm.enter_subsection("velocity");
    velocity->parse_parameters(prm);
    prm.leave_subsection();

    prm.leave_subsection();
  }
} // namespace Parameters



#endif
