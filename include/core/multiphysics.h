/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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
 * This file defines a small enum which is used to identify
 * the various physics solved within Lethe.
 */


#ifndef lethe_multiphysics_h
#define lethe_multiphysics_h

enum PhysicsID : unsigned int
{
  fluid_dynamics   = 0,
  heat_transfer    = 1,
  tracer           = 2,
  VOF              = 3,
  cahn_hilliard    = 4,
  reactive_species = 5
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
  /// Reactive species concentration field from Reactive species
  concentration_reactive_species = 6
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
  else
    return PhysicsID::reactive_species;
}

#endif
