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

#include <deal.II/particles/particle_handler.h>

#include <dem/dem_solver_parameters.h>
#include <dem/pp_contact_info_struct.h>

using namespace dealii;

#ifndef PPCONTACTFORCE_H_
#  define PPCONTACTFORCE_H_

/**
 * Base interface for classes that carry out the calculation of particle-paricle
 * contact force
 */

template <int dim>
class PPContactForce
{
public:
  PPContactForce()
  {}

  virtual ~PPContactForce()
  {}

  /**
   * Carries out the calculation of the contact force using the contact pair
   * information
   * obtained in the fine search and physical properties of particles
   *
   * @param pairs_in_contact_info Required information for calculation of the
   * particle-particle contact force, these information were obtained in the
   * fine search
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  virtual void
  calculate_pp_contact_force(
    const std::vector<std::map<int, pp_contact_info_struct<dim>>>
      &                             pairs_in_contact_info,
    const DEMSolverParameters<dim> &dem_parameters) = 0;
};

#endif /* PPCONTACTFORCE_H_ */
