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

#include <dem/dem_solver_parameters.h>
#include <dem/pw_contact_info_struct.h>

using namespace dealii;

#ifndef PWCONTACTFORCE_H_
#define PWCONTACTFORCE_H_

/**
 * Base interface for classes that carry out the calculation of particle-wall
 * contact force
 */

template <int dim> class PWContactForce {
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
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  virtual void calculate_pw_contact_force(
      const std::vector<std::map<int, pw_contact_info_struct<dim>>>
          *pw_pairs_in_contact,
      const DEMSolverParameters<dim> &dem_parameters) = 0;
};

#endif /* PWCONTACTFORCE_H_ */
