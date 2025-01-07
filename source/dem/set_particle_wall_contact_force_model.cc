// SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/set_particle_wall_contact_force_model.h>

using namespace dealii;

template <int dim, typename PropertiesIndex>
std::shared_ptr<ParticleWallContactForce<dim, PropertiesIndex>>
set_particle_wall_contact_force_model(
  const DEMSolverParameters<dim>                  &dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation)
{
  std::shared_ptr<ParticleWallContactForce<dim, PropertiesIndex>>
    particle_wall_contact_force_object;


  std::vector<types::boundary_id> boundary_index =
    triangulation.get_boundary_ids();
  if (dem_parameters.model_parameters.particle_wall_contact_force_method ==
      Parameters::Lagrangian::ModelParameters::ParticleWallContactForceModel::
        linear)
    {
      particle_wall_contact_force_object =
        std::make_shared<ParticleWallLinearForce<dim, PropertiesIndex>>(
          dem_parameters, boundary_index);
    }
  else if (dem_parameters.model_parameters.particle_wall_contact_force_method ==
           Parameters::Lagrangian::ModelParameters::
             ParticleWallContactForceModel::nonlinear)
    {
      particle_wall_contact_force_object =
        std::make_shared<ParticleWallNonLinearForce<dim, PropertiesIndex>>(
          dem_parameters, boundary_index);
    }
  else if (dem_parameters.model_parameters.particle_wall_contact_force_method ==
           Parameters::Lagrangian::ModelParameters::
             ParticleWallContactForceModel::JKR)
    {
      particle_wall_contact_force_object =
        std::make_shared<ParticleWallJKRForce<dim, PropertiesIndex>>(
          dem_parameters, boundary_index);
    }
  else if (dem_parameters.model_parameters.particle_wall_contact_force_method ==
           Parameters::Lagrangian::ModelParameters::
             ParticleWallContactForceModel::DMT)
    {
      particle_wall_contact_force_object =
        std::make_shared<ParticleWallDMTForce<dim, PropertiesIndex>>(
          dem_parameters, boundary_index);
    }
  else
    {
      throw "The chosen particle-wall contact force model is invalid";
    }
  return particle_wall_contact_force_object;
}

template std::shared_ptr<ParticleWallContactForce<2, DEM::DEMProperties::PropertiesIndex>>
set_particle_wall_contact_force_model(
  const DEMSolverParameters<2>                  &dem_parameters,
  const parallel::distributed::Triangulation<2> &triangulation);

template std::shared_ptr<ParticleWallContactForce<3, DEM::DEMProperties::PropertiesIndex>>
set_particle_wall_contact_force_model(
  const DEMSolverParameters<3>                  &dem_parameters,
  const parallel::distributed::Triangulation<3> &triangulation);

template std::shared_ptr<ParticleWallContactForce<2, DEM::CFDDEMProperties::PropertiesIndex>>
set_particle_wall_contact_force_model(
  const DEMSolverParameters<2>                  &dem_parameters,
  const parallel::distributed::Triangulation<2> &triangulation);

template std::shared_ptr<ParticleWallContactForce<3, DEM::CFDDEMProperties::PropertiesIndex>>
set_particle_wall_contact_force_model(
  const DEMSolverParameters<3>                  &dem_parameters,
  const parallel::distributed::Triangulation<3> &triangulation);
