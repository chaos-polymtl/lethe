// SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam solve_type Type of solver used for the DEM.
 *
 * @param dem_parameters DEM parameters
 * @param triangulation Triangulation
 * @return A pointer to the particle-wall contact force object
 */
template <int dim, typename PropertiesIndex>
std::shared_ptr<ParticleWallContactForce<dim, PropertiesIndex>>
set_particle_wall_contact_force_model(
  const DEMSolverParameters<dim>                  &dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation);

#endif
