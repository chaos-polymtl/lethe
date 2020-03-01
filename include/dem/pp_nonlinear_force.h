/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>
#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/pp_contact_force.h>
#include <dem/pp_contact_info_struct.h>
#include <iostream>
#include <math.h>
#include <vector>

using namespace dealii;

#ifndef PPNONLINEARFORCE_H_
#define PPNONLINEARFORCE_H_

/**
 * Calculation of the non-linear particle-particle contact force using the
 * information obtained from the fine search and physical properties of
 * particles
 *
 * @note
 *
 * @author Shahab Golshan, Bruno Blais, Polytechnique Montreal 2019-
 */

template <int dim> class PPNonLinearForce : public PPContactForce<dim> {
public:
  PPNonLinearForce() {}

  /**
   * Carries out the calculation of the particle-particle contact force using
   * non-linear (Hertzian) model
   *
   * @param pairs_in_contact_info Required information for calculation of the
   * particle-particle contact force, these information were obtained in the
   * fine search
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  virtual void calculate_pp_contact_force(
      const std::vector<std::map<int, pp_contact_info_struct<dim>>>
          &pairs_in_contact_info,
      const DEMSolverParameters<dim> &dem_parameters) override;
};

#endif
