// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_sub_simulation_control_h
#define lethe_sub_simulation_control_h

#include <deal.II/base/exceptions.h>

using namespace dealii;

/**
 * @file sub_simulation_control.h
 * @brief Controls smaller simulation blocks within simulations.
 *
 * Lethe simulations are controlled using the SimulationControl class and its
 * derivative. However, some models contain a sub model with it's own
 * time-stepping. This is the case for the CFD-DEM modules which needs to do a
 * number of DEM iterations per CFD time-step. Whenever this is required, a sub
 * simulation control is used. This file contains classes which are used for sub
 * simulation control. The sub simulation control do not all depend from a
 * base class because they are all very different and don't share a common
 * architecture.
 *
 */

/**
 * @brief Controls DEM sub iteration logic.
 *
 * The SubSimulationControlDEM class is responsible for managing the progression
 * of the DEM component of a CFD-DEM simulation.
 *
 * This class supports two logic for the DEM iteration:
 * - Fixed number of iterations: In this case, a fixed number of DEM iterations
 * is carried out for each CFD time step. If the CFD uses adaptive time
 * stepping, then the DEM simulation time step may change dynamically during the
 * simulation.
 * - Fixed fraction of Rayleigh characteristic time: In this case, the DEM time
 * step is fixed to a fraction of the Rayleigh characteristic time step. This
 * means that the number of DEM iterations per CFD time step is adjusted
 * dynamically.
 *
 */

class SubSimulationControlDEM
{
public:
  /**
   * @brief Type of iteration methods for DEM sub simulation control
   *
   * Declares the types of sub simulation control that are available for CFD-DEM
   * simulations.
   */

  enum class DEMSubIterationLogic
  {
    fixed_number_of_iterations           = 0,
    fixed_fraction_of_rayleigh_time_step = 1
  };

  SubSimulationControlDEM(
    const DEMSubIterationLogic iteration_logic,
    const double               time_interval,
    const unsigned int         coupling_frequency,
    const double               rayleigh_characteristic_time,
    const double               fraction_of_rayleigh_characteristic_time);


  bool
  iterate();

  double
  get_time_step_value() const
  {
    return time_step;
  }

private:
  /// Logic used for the DEM iterations
  const DEMSubIterationLogic iteration_logic;

  /// Time interval over which the time integration is to be carried out.
  const double time_interval;

  /// Coupling frequency. This is only used when
  /// iteration_logic=fixed_number_of_iterations
  const unsigned int coupling_frequency;

  /// Rayleigh characteristic time for the particles
  const double rayleigh_characteristic_time;

  /// Targetted frction of the Rayleigh characteristic time
  const double fraction_of_rayleigh_characteristic_time;

  /// Current value of the sub iteration number. This is initialized to zero.
  unsigned int iteration_number;

  /// Total number of iterations. This is the number of iterations that is
  /// required until the completion of the sub simulation control
  unsigned int total_number_of_iterations;

  /// Time step value of the sub time stepping
  double time_step;
};

#endif
