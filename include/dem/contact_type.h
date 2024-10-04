// SPDX-FileCopyrightText: Copyright (c) 2022, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_contact_type_h
#define lethe_contact_type_h

/**
 * @brief Label for particle-object contact type.
 */
enum ContactType
{
  local_particle_particle,
  ghost_particle_particle,
  local_periodic_particle_particle,
  ghost_periodic_particle_particle,
  ghost_local_periodic_particle_particle,
  particle_wall,
  particle_floating_wall,
  particle_floating_mesh,
  particle_point,
  particle_line
};

#endif
