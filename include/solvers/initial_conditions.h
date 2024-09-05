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
    ramp,
    average_velocity_profile
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
    double kinematic_viscosity_init;
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

    // Artificial kinematic viscosity
    double kinematic_viscosity;

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

    // Reactive species
    Functions::ParsedFunction<dim> reactive_species =
      Functions::ParsedFunction<dim>(
        4); // TODO Change to a flexible number of species

    // Path to the checkpointed average velocity profile
    std::string average_velocity_folder;
    std::string average_velocity_file_name;

    void
    declare_parameters(ParameterHandler  &prm,
                       const unsigned int reactive_species_count = 0);
    void
    parse_parameters(ParameterHandler &prm);
  };

  template <int dim>
  void
  InitialConditions<dim>::declare_parameters(
    ParameterHandler  &prm,
    const unsigned int reactive_species_count)
  {
    prm.enter_subsection("initial conditions");
    {
      prm.declare_entry(
        "type",
        "nodal",
        Patterns::Selection(
          "L2projection|viscous|nodal|ramp|average_velocity_profile"),
        "Type of initial condition"
        "Choices are <L2projection|viscous|nodal|ramp|average_velocity_profile>.");
      prm.enter_subsection("uvwp");
      uvwp.declare_parameters(prm, dim + 1);
      prm.leave_subsection();

      prm.declare_entry("kinematic viscosity",
                        "1",
                        Patterns::Double(),
                        "Kinematic viscosity for viscous initial conditions");


      prm.enter_subsection("temperature");
      temperature.declare_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("tracer");
      tracer.declare_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("VOF");
      VOF.declare_parameters(prm);
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
      prm.leave_subsection();

      prm.enter_subsection("reactive species");
      reactive_species.declare_parameters(
        prm, 4); // TODO Change to more flexible number of species
      prm.leave_subsection();

      prm.enter_subsection("average velocity profile");
      prm.declare_entry(
        "checkpoint folder",
        "./",
        Patterns::FileName(),
        "the path leading to the checkpointed average velocity profile");
      prm.declare_entry("checkpoint file name",
                        "restart",
                        Patterns::FileName(),
                        "checkpoint file name");
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
      else if (op == "average_velocity_profile")
        type = InitialConditionType::average_velocity_profile;

      kinematic_viscosity = prm.get_double("kinematic viscosity");
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

      prm.enter_subsection("reactive species");
      reactive_species.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("average velocity profile");
      average_velocity_folder    = prm.get("checkpoint folder");
      average_velocity_file_name = prm.get("checkpoint file name");

      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
} // namespace Parameters

#endif
