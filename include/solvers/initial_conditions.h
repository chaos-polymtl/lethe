// SPDX-FileCopyrightText: Copyright (c) 2019-2026 The Lethe Authors
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
    L2projection,
    viscous,
    nodal,
    ramp,
    average_velocity_profile
  };

  // Type of initial conditions for CLS
  enum class CLSInitialConditionType
  {
    none,      // Uses the function as is
    diffusive, // Uses a L2 projection with a smoothing coefficient
    geometric  // Uses geometric redistanciation
  };

  struct Ramp_n
  {
    double n_init = 1.0;
    int    n_iter = 0;
    double alpha  = 0.5;

    /**
     * @brief Declares parameters for the ramped power-law exponent.
     *
     * @param[in,out] prm The parameter handler object.
     */
    void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parses parameters for the ramped power-law exponent.
     *
     * @param[in] prm The parameter handler object.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct Ramp_viscosity
  {
    double kinematic_viscosity_init = 1.0;
    int    n_iter                   = 0;
    double alpha                    = 0.5;

    /**
     * @brief Declares parameters for the ramped kinematic viscosity.
     *
     * @param[in,out] prm The parameter handler object.
     */
    void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parses parameters for the ramped kinematic viscosity.
     *
     * @param[in] prm The parameter handler object.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct Ramp
  {
    Ramp_n         ramp_n;
    Ramp_viscosity ramp_viscosity;

    /**
     * @brief Declares parameters for all ramp functions.
     *
     * @param[in,out] prm The parameter handler object.
     */
    void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parses parameters for all ramp functions.
     *
     * @param[in] prm The parameter handler object.
     */
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

    FluidDynamicsInitialConditionType type =
      FluidDynamicsInitialConditionType::nodal;

    // Velocity components
    Functions::ParsedFunction<dim> uvwp;

    // Artificial kinematic viscosity
    double kinematic_viscosity = 1;

    // Temperature
    Functions::ParsedFunction<dim> temperature;

    // Tracer
    Functions::ParsedFunction<dim> tracer;

    // CLS
    Functions::ParsedFunction<dim> CLS;


    CLSInitialConditionType cls_initial_condition_smoothing =
      CLSInitialConditionType::none;
    double projection_step_diffusion_factor = 1;

    // Non-Newtonian
    Ramp ramp;

    // Cahn-Hilliard
    Functions::ParsedFunction<dim> cahn_hilliard =
      Functions::ParsedFunction<dim>(2);

    // Time-Harmonic Maxwell
    Functions::ParsedFunction<dim> electromagnetics =
      Functions::ParsedFunction<dim>(4 * dim);

    // Path to the checkpointed average velocity profile
    std::string average_velocity_folder    = "./";
    std::string average_velocity_file_name = "restart";

    /**
     * @brief Declares initial-condition parameters.
     *
     * @param[in,out] prm The parameter handler object.
     */
    void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parses initial-condition parameters.
     *
     * @param[in] prm The parameter handler object.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };
} // namespace Parameters

#endif
