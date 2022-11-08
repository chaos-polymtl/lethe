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
#include <dem/particle_particle_contact_force.h>

using namespace std;

#ifndef set_particle_particle_contact_force_model_h
#  define set_particle_particle_contact_force_model_h

/**
 * Sets the selected particle-particle contact force model in the parameter
 * handler file.
 *
 * @param dem_parameters DEM parameters
 * @return A pointer to the particle-particle contact force object
 */
template <int dim>
std::shared_ptr<ParticleParticleContactForceBase<dim>>
set_particle_particle_contact_force_model(
  const DEMSolverParameters<dim> &dem_parameters);

#endif /* set_particle_particle_contact_force_model_h */
