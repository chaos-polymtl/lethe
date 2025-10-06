// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// TODO : Refactor so the class itself is not a pointer, but contains a pointer
//        to a function. This would be a lot more coherent...

#ifndef lethe_source_terms_h
#define lethe_source_terms_h

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <memory>

using namespace dealii;
/**
 * The source term class provides an interface for all common
 * Element required for the introduction of a source term in equations
 * All equation-specific source term should derive
 * from the base class but also call it's declare_parameters and
 * parse_parameters routine. This allows specialized classes to focus on their
 * specificity and forget about other non-specific elements that are generic to
 * the calculation of analytical solutions
 **/

namespace SourceTerms
{
  template <int dim>
  class SourceTerm
  {
  public:
    SourceTerm()
    {
      navier_stokes_source =
        std::make_shared<Functions::ParsedFunction<dim>>(dim + 1);
      heat_transfer_source =
        std::make_shared<Functions::ParsedFunction<dim>>(1);
      tracer_source = std::make_shared<Functions::ParsedFunction<dim>>(1);
      cahn_hilliard_source =
        std::make_shared<Functions::ParsedFunction<dim>>(2);
    }

    virtual void
    declare_parameters(ParameterHandler &prm);
    virtual void
    parse_parameters(ParameterHandler &prm);

    /// Enable the Navier-Stokes source term
    bool enable;

    /// Velocity-pressure components
    std::shared_ptr<Functions::ParsedFunction<dim>> navier_stokes_source;

    /// Heat transfer source
    std::shared_ptr<Functions::ParsedFunction<dim>> heat_transfer_source;

    /// Tracer source
    std::shared_ptr<Functions::ParsedFunction<dim>> tracer_source;

    /// Cahn-Hilliard source
    std::shared_ptr<Functions::ParsedFunction<dim>> cahn_hilliard_source;
  };

  template <int dim>
  void
  SourceTerm<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("source term");

    prm.enter_subsection("fluid dynamics");
    prm.declare_entry("enable",
                      "true",
                      Patterns::Bool(),
                      "Enable the usage of a source term for the fluid solver");
    navier_stokes_source->declare_parameters(prm, dim + 1);
    prm.leave_subsection();

    prm.enter_subsection("heat transfer");
    heat_transfer_source->declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("tracer");
    tracer_source->declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("cahn hilliard");
    cahn_hilliard_source->declare_parameters(prm, 2);
    prm.leave_subsection();

    prm.leave_subsection();
  }

  template <int dim>
  void
  SourceTerm<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("source term");

    prm.enter_subsection("fluid dynamics");
    enable = prm.get_bool("enable");
    navier_stokes_source->parse_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("heat transfer");
    heat_transfer_source->parse_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("tracer");
    tracer_source->parse_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("cahn hilliard");
    cahn_hilliard_source->parse_parameters(prm);
    prm.leave_subsection();

    prm.leave_subsection();
  }
} // namespace SourceTerms


#endif
