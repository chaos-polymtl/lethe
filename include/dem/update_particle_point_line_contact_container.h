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
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */
#include <dem/particle_point_line_contact_info_struct.h>

using namespace dealii;

#ifndef UPDATEPARTICLEPOINTLINECONTACTCONTAINER_H_
#  define UPDATEPARTICLEPOINTLINECONTACTCONTAINER_H_

/**
 * Updates the iterators to particles in particle_points_in_contact and
 * particle_lines_in_contact (output of particle point line fine search)
 *
 * @param particle_points_in_contact Output of particle-point fine search
 * @param particle_lines_in_contact Output of particle-line fine search
 * @param particle_container Output of update_particle_container function
 */

template <int dim>
void
update_particle_point_line_contact_container_iterators(
  std::map<int, particle_point_line_contact_info_struct<dim>>
    &particle_points_in_contact,
  std::map<int, particle_point_line_contact_info_struct<dim>>
    &particle_lines_in_contact,
  const std::unordered_map<int, Particles::ParticleIterator<dim>>
    &particle_container);

#endif /* UPDATEPARTICLEPOINTLINECONTACTCONTAINER_H_ */
