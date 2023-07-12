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

#ifndef lethe_source_terms_h
#define lethe_source_terms_h

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

namespace SourceTerms
{
  template <int dim>
  class SourceTerm
  {
  public:
    SourceTerm()
      : navier_stokes_source(dim + 1)
      , heat_transfer_source(1)
      , tracer_source(1)
      , cahn_hilliard_source(2)
    {}

    virtual void
    declare_parameters(ParameterHandler &prm);
    virtual void
    parse_parameters(ParameterHandler &prm);


    bool
    source_term()
    {
      return enable;
    }

    // Velocity-pressure components
    Functions::ParsedFunction<dim> navier_stokes_source;

    // Heat transfer source
    Functions::ParsedFunction<dim> heat_transfer_source;

    // Tracer source
    Functions::ParsedFunction<dim> tracer_source;

    // Cahn-Hilliard source
    Functions::ParsedFunction<dim> cahn_hilliard_source;

  protected:
    bool enable;
  };

  template <int dim>
  void
  SourceTerm<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("source term");
    prm.declare_entry("enable",
                      "false",
                      Patterns::Bool(),
                      "Enable the calculation of a source term");

    prm.enter_subsection("xyz");
    navier_stokes_source.declare_parameters(prm, dim);
    if (dim == 2)
      prm.set("Function expression", "0; 0; 0");
    if (dim == 3)
      prm.set("Function expression", "0; 0; 0; 0;");
    prm.leave_subsection();


    prm.enter_subsection("heat transfer");
    heat_transfer_source.declare_parameters(prm);
    prm.set("Function expression", "0");
    prm.leave_subsection();

    prm.enter_subsection("tracer");
    tracer_source.declare_parameters(prm);
    prm.set("Function expression", "0");
    prm.leave_subsection();

    prm.enter_subsection("cahn hilliard");
    cahn_hilliard_source.declare_parameters(prm);
    prm.set("Function expression", "0; 0;");
    prm.leave_subsection();

    prm.leave_subsection();
  }

  template <int dim>
  void
  SourceTerm<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("source term");
    enable = prm.get_bool("enable");

    prm.enter_subsection("xyz");
    navier_stokes_source.parse_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("heat transfer");
    heat_transfer_source.parse_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("tracer");
    tracer_source.parse_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("cahn hilliard");
    cahn_hilliard_source.parse_parameters(prm);
    prm.leave_subsection();

    prm.leave_subsection();
  }
} // namespace SourceTerms

template <int dim>
class NoForce : public Function<dim>
{
public:
  NoForce()
    : Function<dim>(3){};
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override;
};
template <int dim>
void
NoForce<dim>::vector_value(const Point<dim> & /*p*/,
                           Vector<double> &values) const
{
  values(0) = 0.;
  values(1) = 0.;
  if (dim == 3)
    values(2) = 0.;
}


#endif
