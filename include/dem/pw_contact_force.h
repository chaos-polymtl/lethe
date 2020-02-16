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

#include "dem/physical_info_struct.h"
#include "dem/pw_contact_info_struct.h"

using namespace dealii;

#ifndef PWCONTACTFORCE_H_
#define PWCONTACTFORCE_H_

/**
 * Base interface for classes that carry out the calculation of particle-wall
 * contact force
 */

template <int dim, int spacedim = dim> class PWContactForce {
public:
  PWContactForce() {}

  virtual ~PWContactForce() {}

  /**
   * Carries out the calculation of the particle-wall contact force using the
   * contact pair information obtained from the particle-wall fine search and
   * physical properties of particles and walls
   *
   * @param pw_pairs_in_contact Required information for calculation of the
   * particle-wall contact force
   * @param physical_properties Physical properties of particles and walls
   */
  virtual void calculate_pw_contact_force(
      std::vector<std::map<int, pw_contact_info_struct<dim, spacedim>>>
          &pw_pairs_in_contact,
      const physical_info_struct<dim> &physical_info_struct) = 0;
};

#endif /* PWCONTACTFORCE_H_ */
