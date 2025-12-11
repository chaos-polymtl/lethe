// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_initial_conditions_h
#define lethe_initial_conditions_h

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;

namespace Parameters
{
  // Type of initial conditions for fluid dynamics
  enum class FluidDynamicsInitialConditionType
  {
    none,
    L2projection,
    viscous,
    nodal,
    ramp,
    average_velocity_profile
  };

  // Type of initial conditions for VOF
  enum class VOFInitialConditionType
  {
    none,      // Uses the function as is
    diffusive, // Uses a L2 projection with a smoothing coefficient
    geometric  // Uses geometric redistanciation
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

    FluidDynamicsInitialConditionType type;

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


    VOFInitialConditionType vof_initial_condition_smoothing;
    double                  projection_step_diffusion_factor;

    // Non-Newtonian
    Ramp ramp;

    // Cahn-Hilliard
    Functions::ParsedFunction<dim> cahn_hilliard =
      Functions::ParsedFunction<dim>(2);

    // Time-Harmonic Maxwells
    Functions::ParsedFunction<dim> electromagnetics =
      Functions::ParsedFunction<dim>(4 * dim);

    // Path to the checkpointed average velocity profile
    std::string average_velocity_folder;
    std::string average_velocity_file_name;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };
} // namespace Parameters

#endif
