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
#include <math.h>

#include <iostream>
#include <vector>

#include "dem/dem_solver_parameters.h"
#include "dem/physical_info_struct.h"
#include "dem/pp_contact_force.h"
#include "dem/pp_contact_info_struct.h"

using namespace dealii;

#ifndef PPLINEARFORCE_H_
#  define PPLINEARFORCE_H_

/**
 * Calculation of the linear particle-particle contact force using the
 * information obtained from the fine search and physical properties of
 * particles
 *
 * @note
 *
 * @author Shahab Golshan, Bruno Blais, Polytechnique Montreal 2019-
 */

template <int dim, int spacedim = dim>
class PPLinearForce : public PPContactForce<dim, spacedim>
{
public:
  PPLinearForce()
  {}

  /**
   * Carries out the calculation of the particle-particle contact force using
   * linear (Hookean) model
   *
   * @param pairs_in_contact_info Required information for calculation of the
   * particle-particle contact force, these information were obtained in the
   * fine search
   * @param physical_properties Physical properties of particles
   */
  virtual void
  calculate_pp_contact_force(
    const std::vector<std::map<int, pp_contact_info_struct<dim, spacedim>>>
      &                              pairs_in_contact_info,
    const physical_info_struct<dim> &physical_properties) override;
};

#endif
