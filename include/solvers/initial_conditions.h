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

#ifndef lethe_initial_conditions_h
#define lethe_initial_conditions_h

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;

namespace Parameters
{
  // Type of initial conditions
  enum class InitialConditionType
  {
    none,
    L2projection,
    viscous,
    nodal,
    ramp
  };

  struct Ramp_n
  {
    double n_init;
    int    n_iter;
    double alpha;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct Ramp_viscosity
  {
    double viscosity_init;
    int    n_iter;
    double alpha;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct Ramp
  {
    Ramp_n         ramp_n;
    Ramp_viscosity ramp_viscosity;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  template <int dim>
  class InitialConditions
  {
  public:
    InitialConditions()
      : uvwp(dim + 1)
    {}

    InitialConditionType type;

    // Velocity components
    Functions::ParsedFunction<dim> uvwp;

    // Artificial viscosity
    double viscosity;

    // Temperature
    Functions::ParsedFunction<dim> temperature;

    // Tracer
    Functions::ParsedFunction<dim> tracer;

    // VOF
    Functions::ParsedFunction<dim> VOF;
    // Bool to apply a Galerkin projection (with a diffusion term) to the VOF
    // initial condition
    bool   enable_projection_step;
    double projection_step_diffusion_factor;

    // Non-Newtonian
    Ramp ramp;

    // Cahn-Hilliard
    Functions::ParsedFunction<dim> cahn_hilliard =
      Functions::ParsedFunction<dim>(2);

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  template <int dim>
  void
  InitialConditions<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("initial conditions");
    {
      prm.declare_entry("type",
                        "nodal",
                        Patterns::Selection("L2projection|viscous|nodal|ramp"),
                        "Type of initial condition"
                        "Choices are <L2projection|viscous|nodal|ramp>.");
      prm.enter_subsection("uvwp");
      uvwp.declare_parameters(prm, dim);
      if (dim == 2)
        prm.set("Function expression", "0; 0; 0");
      if (dim == 3)
        prm.set("Function expression", "0; 0; 0; 0");
      prm.leave_subsection();

      prm.declare_entry("viscosity",
                        "1",
                        Patterns::Double(),
                        "viscosity for viscous initial conditions");


      prm.enter_subsection("temperature");
      temperature.declare_parameters(prm);
      prm.set("Function expression", "0");
      prm.leave_subsection();

      prm.enter_subsection("tracer");
      tracer.declare_parameters(prm);
      prm.set("Function expression", "0");
      prm.leave_subsection();

      prm.enter_subsection("VOF");
      VOF.declare_parameters(prm);
      prm.set("Function expression", "0");
      prm.enter_subsection("projection step");
      prm.declare_entry(
        "enable",
        "false",
        Patterns::Bool(),
        "Apply a projection step with diffusion to smooth the VOF initial condition");
      prm.declare_entry(
        "diffusion factor",
        "1",
        Patterns::Double(),
        "Factor applied to the diffusion term in the projection step");
      prm.leave_subsection();
      prm.leave_subsection();

      prm.enter_subsection("cahn hilliard");
      cahn_hilliard.declare_parameters(prm, 2);
      prm.set("Function expression", "0; 0;");
      prm.leave_subsection();

      ramp.declare_parameters(prm);
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  InitialConditions<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("initial conditions");
    {
      const std::string op = prm.get("type");
      if (op == "L2projection")
        type = InitialConditionType::L2projection;
      else if (op == "viscous")
        type = InitialConditionType::viscous;
      else if (op == "nodal")
        type = InitialConditionType::nodal;
      else if (op == "ramp")
        type = InitialConditionType::ramp;

      viscosity = prm.get_double("viscosity");
      prm.enter_subsection("uvwp");
      uvwp.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("temperature");
      temperature.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("tracer");
      tracer.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("VOF");
      {
        VOF.parse_parameters(prm);
        prm.enter_subsection("projection step");
        {
          enable_projection_step           = prm.get_bool("enable");
          projection_step_diffusion_factor = prm.get_double("diffusion factor");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      ramp.parse_parameters(prm);

      prm.enter_subsection("cahn hilliard");
      cahn_hilliard.parse_parameters(prm);
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
} // namespace Parameters

#endif
