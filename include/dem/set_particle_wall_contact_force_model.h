/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
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
 */

#ifndef lethe_set_particle_wall_contact_force_model_h
#define lethe_set_particle_wall_contact_force_model_h

#include <dem/dem_solver_parameters.h>
#include <dem/particle_wall_dmt_force.h>
#include <dem/particle_wall_jkr_force.h>
#include <dem/particle_wall_linear_force.h>
#include <dem/particle_wall_nonlinear_force.h>

#include <deal.II/distributed/tria.h>

/**
 * @brief Set the selected particle-wall contact force model in the parameter
 * handler file.
 *
 * @param dem_parameters DEM parameters
 * @param triangulation Triangulation
 * @return A pointer to the particle-wall contact force object
 */
template <int dim>
std::shared_ptr<ParticleWallContactForce<dim>>
set_particle_wall_contact_force_model(
  const DEMSolverParameters<dim>                  &dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation);

#endif
