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

#ifndef LETHE_SOURCETERMS_H
#define LETHE_SOURCETERMS_H

#include <deal.II/base/function.h>
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

  protected:
    bool enable;
  };

  template <int dim>
  void
  SourceTerm<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("source term");
    prm.declare_entry(
      "enable",
      "true",
      Patterns::Bool(),
      "Enable the calculation of the analytical solution and L2 error");
    prm.leave_subsection();
  }

  template <int dim>
  void
  SourceTerm<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("source term");
    enable = prm.get_bool("enable");
    prm.leave_subsection();
  }



  template <int dim>
  class NSSourceTerm : public SourceTerm<dim>
  {
  public:
    NSSourceTerm()
      : source(dim + 1)
    {}

    // Velocity components
    Functions::ParsedFunction<dim> source;

    virtual void
    declare_parameters(ParameterHandler &prm);
    virtual void
    parse_parameters(ParameterHandler &prm);
  };

  template <int dim>
  void
  NSSourceTerm<dim>::declare_parameters(ParameterHandler &prm)
  {
    this->SourceTerm<dim>::declare_parameters(prm);
    prm.enter_subsection("source term");
    prm.enter_subsection("xyz");
    source.declare_parameters(prm, dim);
    if (dim == 2)
      prm.set("Function expression", "0; 0; 0");
    if (dim == 3)
      prm.set("Function expression", "0; 0; 0; 0;");
    prm.leave_subsection();
    prm.leave_subsection();
    //    if (this->enable=false) source=NULL;
  }

  template <int dim>
  void
  NSSourceTerm<dim>::parse_parameters(ParameterHandler &prm)
  {
    this->SourceTerm<dim>::parse_parameters(prm);
    prm.enter_subsection("source term");
    prm.enter_subsection("xyz");
    source.parse_parameters(prm);
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
  vector_value(const Point<dim> &p, Vector<double> &values) const;
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
