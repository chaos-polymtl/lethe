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
 * Author: Bruno Blais, Polytechnique Montreal, 2021-
 */

/*
 * This file defines a small enum which is used to identify
 * the various physics solved within Lethe.
 */


#ifndef lethe_multiphysics_h
#define lethe_multiphysics_h

enum PhysicsID : unsigned int
{
  fluid_dynamics = 0,
  heat_transfer  = 1,
  tracer         = 2,
  VOF            = 3,
  cahn_hilliard  = 4
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
  else
    return PhysicsID::cahn_hilliard;
}

#endif
