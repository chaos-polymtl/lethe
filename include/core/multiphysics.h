// SPDX-FileCopyrightText: Copyright (c) 2019, 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_multiphysics_h
#define lethe_multiphysics_h

enum PhysicsID : unsigned int
{
  fluid_dynamics = 0,
  heat_transfer  = 1,
  tracer         = 2,
  VOF            = 3,
  cahn_hilliard  = 4,
  void_fraction  = 5
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
  /// Phase fraction scalar field from VOF
  phase = 2,
  /// Temperature scalar field from heat transfer
  temperature = 3,
  /// Phase order scalar field from Cahn Hilliard
  phase_cahn_hilliard = 4,
  /// Chemical potential scalar field from Cahn Hilliard
  chemical_potential_cahn_hilliard = 5,
  /// Tracer scalar field
  tracer = 6
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
  else if (physics_name == "VOF")
    return PhysicsID::VOF;
  else if (physics_name == "cahn hilliard")
    return PhysicsID::cahn_hilliard;
  else if (physics_name == "void fraction")
    return PhysicsID::void_fraction;
  else
    AssertThrow(false,
                dealii::StandardExceptions::ExcMessage(
                  "An unknown Physics name was requested"));
}

#endif
