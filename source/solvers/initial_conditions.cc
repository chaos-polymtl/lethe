// SPDX-FileCopyrightText: Copyright (c) 2019, 2021-2023, 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/utilities.h>

#include <solvers/initial_conditions.h>

namespace Parameters
{
  void
  Ramp_n::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("n");
    {
      prm.declare_entry(
        "initial n",
        Patterns::Tools::Convert<double>::to_string(n_init),
        Patterns::Double(),
        "First n value with which to start the initial condition");

      prm.declare_entry(
        "iterations",
        Patterns::Tools::Convert<int>::to_string(n_iter),
        Patterns::Integer(0),
        "Number of iterations used in the ramp before reaching the final n value");

      prm.declare_entry("alpha",
                        Patterns::Tools::Convert<double>::to_string(alpha),
                        Patterns::Double(),
                        "Coefficient used for n-spacing.");
    }
    prm.leave_subsection();
  }

  void
  Ramp_n::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("n");
    {
      n_init = prm.get_double("initial n");
      n_iter = prm.get_integer("iterations");
      alpha  = prm.get_double("alpha");
    }
    prm.leave_subsection();
  }

  void
  Ramp_viscosity::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("kinematic viscosity");
    {
      prm.declare_entry(
        "initial kinematic viscosity",
        Patterns::Tools::Convert<double>::to_string(kinematic_viscosity_init),
        Patterns::Double(),
        "First kinematic viscosity value with which to start the initial condition");

      prm.declare_entry(
        "iterations",
        Patterns::Tools::Convert<int>::to_string(n_iter),
        Patterns::Integer(),
        "Number of iterations used in the ramp before reaching the final kinematic viscosity value");

      prm.declare_entry("alpha",
                        Patterns::Tools::Convert<double>::to_string(alpha),
                        Patterns::Double(),
                        "Coefficient used for kinematic viscosity-spacing.");
    }
    prm.leave_subsection();
  }

  void
  Ramp_viscosity::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("kinematic viscosity");
    {
      kinematic_viscosity_init = prm.get_double("initial kinematic viscosity");
      n_iter                   = prm.get_integer("iterations");
      alpha                    = prm.get_double("alpha");
    }
    prm.leave_subsection();
  }

  void
  Ramp::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("ramp");
    {
      ramp_n.declare_parameters(prm);
      ramp_viscosity.declare_parameters(prm);
    }
    prm.leave_subsection();
  }

  void
  Ramp::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("ramp");
    {
      ramp_n.parse_parameters(prm);
      ramp_viscosity.parse_parameters(prm);
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  InitialConditions<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("initial conditions");
    {
      declare_enum_entry(prm, "type", type, "Type of initial condition.");
      prm.enter_subsection("uvwp");
      uvwp.declare_parameters(prm, dim + 1);
      prm.leave_subsection();

      prm.declare_entry("kinematic viscosity",
                        Patterns::Tools::Convert<double>::to_string(
                          kinematic_viscosity),
                        Patterns::Double(),
                        "Kinematic viscosity for viscous initial conditions");


      prm.enter_subsection("temperature");
      temperature.declare_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("tracer");
      tracer.declare_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("CLS");
      CLS.declare_parameters(prm);
      declare_enum_entry(
        prm,
        "smoothing type",
        cls_initial_condition_smoothing,
        "Apply a projection step with diffusion to smooth the CLS initial condition.");

      prm.declare_entry(
        "diffusion factor",
        Patterns::Tools::Convert<double>::to_string(
          projection_step_diffusion_factor),
        Patterns::Double(),
        "Factor applied to the diffusion term in the projection step");
      prm.leave_subsection();

      prm.enter_subsection("cahn hilliard");
      cahn_hilliard.declare_parameters(prm, 2);
      prm.leave_subsection();

      prm.enter_subsection("electromagnetics");
      electromagnetics.declare_parameters(prm, 4 * dim);
      prm.leave_subsection();

      prm.enter_subsection("average velocity profile");
      prm.declare_entry(
        "checkpoint folder",
        average_velocity_folder,
        Patterns::FileName(),
        "the path leading to the checkpointed average velocity profile");
      prm.declare_entry("checkpoint file name",
                        average_velocity_file_name,
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
      type = string_to_enum<FluidDynamicsInitialConditionType>(prm.get("type"));

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

      prm.enter_subsection("CLS");
      {
        CLS.parse_parameters(prm);
        cls_initial_condition_smoothing =
          string_to_enum<CLSInitialConditionType>(prm.get("smoothing type"));

        projection_step_diffusion_factor = prm.get_double("diffusion factor");
      }
      prm.leave_subsection();

      ramp.parse_parameters(prm);

      prm.enter_subsection("cahn hilliard");
      cahn_hilliard.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("electromagnetics");
      electromagnetics.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("average velocity profile");
      average_velocity_folder    = prm.get("checkpoint folder");
      average_velocity_file_name = prm.get("checkpoint file name");

      prm.leave_subsection();
    }
    prm.leave_subsection();
  }


  template class InitialConditions<2>;
  template class InitialConditions<3>;
} // namespace Parameters
