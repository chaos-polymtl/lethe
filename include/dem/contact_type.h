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
 */

#ifndef lethe_contact_type_h
#define lethe_contact_type_h

/**
 * @brief Label for particle-object contact type.
 *
 */

enum ContactType
{
  local_particle_particle,
  ghost_particle_particle,
  particle_wall,
  particle_floating_wall,
  particle_floating_mesh,
  particle_point,
  particle_line
};

#endif // lethe_contact_type_h
