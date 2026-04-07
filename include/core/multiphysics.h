// SPDX-FileCopyrightText: Copyright (c) 2019, 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file multiphysics.h
 * @brief Identifiers for multiphysics modules and their solution fields.
 *
 * This file defines the enumerations used to identify the different physics
 * modules available in Lethe (fluid dynamics, heat transfer, CLS, etc.) and
 * the solution fields they produce. It also provides a utility function to
 * convert a physics name string into its corresponding PhysicsID.
 */

#ifndef lethe_multiphysics_h
#define lethe_multiphysics_h

/**
 * @brief Unique identifiers for each physics module in Lethe.
 *
 * These identifiers are used to index and refer to the different physics
 * solvers throughout the multiphysics framework.
 */
enum PhysicsID : unsigned int
{
  /// Fluid Dynamics (Velocity,Pressure)
  fluid_dynamics = 0,
  /// Heat transfer (Temperature)
  heat_transfer = 1,
  /// Passive scalar tracer (Tracer)
  tracer = 2,
  /// Conservative level-set / Volume-of-Fluid (Phase indicator / phase
  /// fraction)
  CLS = 3,
  /// Cahn-Hilliard (Phase indicator, Potential)
  cahn_hilliard = 4,
  /// Void fraction
  void_fraction = 5,
  /// Time-harmonic electromagnetics (Electric field, Magnetic field)
  electromagnetics = 6
};

/**
 * @brief Solution fields of the different physics that are used as an indicator
 * for multiple purposes (e.g. adaptive mesh refinement, solid domain
 * constraints).
 */
enum class Variable : unsigned int
{ /// Velocity vector field from fluid dynamics
  velocity = 0,
  /// Pressure scalar field from fluid dynamics
  pressure = 1,
  /// Phase indicator scalar field from CLS
  phase = 2,
  /// Temperature scalar field from heat transfer
  temperature = 3,
  /// Phase order scalar field from Cahn Hilliard
  phase_cahn_hilliard = 4,
  /// Chemical potential scalar field from Cahn Hilliard
  chemical_potential_cahn_hilliard = 5,
  /// Tracer scalar field
  tracer = 6,
  /// Electric field vector field from electromagnetics
  electric_field = 7,
  /// Magnetic field vector field from electromagnetics
  magnetic_field = 8,
  /// Combination of both electric and magnetic fields from electromagnetics
  electromagnetic_fields = 9
};

/**
 * @brief Utility function used for parsing physics-based
 * parameters
 *
 */
inline PhysicsID
get_physics_id(std::string physics_name)
{
  if (physics_name == "fluid dynamics")
    return PhysicsID::fluid_dynamics;
  else if (physics_name == "heat transfer")
    return PhysicsID::heat_transfer;
  else if (physics_name == "tracer")
    return PhysicsID::tracer;
  else if (physics_name == "CLS")
    return PhysicsID::CLS;
  else if (physics_name == "cahn hilliard")
    return PhysicsID::cahn_hilliard;
  else if (physics_name == "void fraction")
    return PhysicsID::void_fraction;
  else if (physics_name == "electromagnetics")
    return PhysicsID::electromagnetics;
  else
    AssertThrow(false,
                dealii::StandardExceptions::ExcMessage(
                  "An unknown Physics name was requested"));
}

#endif
