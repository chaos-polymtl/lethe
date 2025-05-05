// SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


#include <dem/set_particle_wall_contact_force_model.h>


using namespace dealii;

using namespace Parameters::Lagrangian;

template <int dim, typename PropertiesIndex>
std::shared_ptr<ParticleWallContactForceBase<dim, PropertiesIndex>>
set_particle_wall_contact_force_model(
  const DEMSolverParameters<dim>                  &dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation)
{
  ParticleWallContactForceModel particle_wall_contact_force_model =
    dem_parameters.model_parameters.particle_wall_contact_force_model;

  std::shared_ptr<ParticleWallContactForceBase<dim, PropertiesIndex>>
    particle_wall_contact_force_object;

  const std::vector<types::boundary_id> boundary_index =
    triangulation.get_boundary_ids();

  switch (particle_wall_contact_force_model)
    {
      case ParticleWallContactForceModel::linear:
        {
          set_rolling_resistance_model<dim,
                                       PropertiesIndex,
                                       ParticleWallContactForceModel::linear>(
            dem_parameters, boundary_index, particle_wall_contact_force_object);
          break;
        }
      case ParticleWallContactForceModel::nonlinear:
        {
          set_rolling_resistance_model<
            dim,
            PropertiesIndex,
            ParticleWallContactForceModel::nonlinear>(
            dem_parameters, boundary_index, particle_wall_contact_force_object);
          break;
        }
      case ParticleWallContactForceModel::JKR:
        {
          set_rolling_resistance_model<dim,
                                       PropertiesIndex,
                                       ParticleWallContactForceModel::JKR>(
            dem_parameters, boundary_index, particle_wall_contact_force_object);
          break;
        }
      case ParticleWallContactForceModel::DMT:
        {
          set_rolling_resistance_model<dim,
                                       PropertiesIndex,
                                       ParticleWallContactForceModel::DMT>(
            dem_parameters, boundary_index, particle_wall_contact_force_object);
          break;
        }
      default:
        {
          throw std::runtime_error(
            "The chosen particle-wall contact force model is invalid");
          break;
        }
    }
  return particle_wall_contact_force_object;
}


template <int dim,
          typename PropertiesIndex,
          Parameters::Lagrangian::ParticleWallContactForceModel
            particle_wall_contact_force_model>
void
set_rolling_resistance_model(
  const DEMSolverParameters<dim>       &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<ParticleWallContactForceBase<dim, PropertiesIndex>>
    &particle_wall_contact_force_object)
{
  RollingResistanceMethod rolling_resistance_method =
    dem_parameters.model_parameters.rolling_resistance_method;

  switch (rolling_resistance_method)
    {
      case RollingResistanceMethod::no_resistance:
        {
          particle_wall_contact_force_object = std::make_shared<
            ParticleWallContactForce<dim,
                                     PropertiesIndex,
                                     particle_wall_contact_force_model,
                                     RollingResistanceMethod::no_resistance>>(
            dem_parameters, boundary_index);
          break;
        }
      case RollingResistanceMethod::constant_resistance:
        {
          particle_wall_contact_force_object =
            std::make_shared<ParticleWallContactForce<
              dim,
              PropertiesIndex,
              particle_wall_contact_force_model,
              RollingResistanceMethod::constant_resistance>>(dem_parameters,
                                                             boundary_index);
          break;
        }
      case RollingResistanceMethod::viscous_resistance:
        {
          particle_wall_contact_force_object =
            std::make_shared<ParticleWallContactForce<
              dim,
              PropertiesIndex,
              particle_wall_contact_force_model,
              RollingResistanceMethod::viscous_resistance>>(dem_parameters,
                                                            boundary_index);
          break;
        }
      case RollingResistanceMethod::epsd_resistance:
        {
          particle_wall_contact_force_object = std::make_shared<
            ParticleWallContactForce<dim,
                                     PropertiesIndex,
                                     particle_wall_contact_force_model,
                                     RollingResistanceMethod::epsd_resistance>>(
            dem_parameters, boundary_index);
          break;
        }
      default:
        throw std::runtime_error("Invalid rolling resistance method");
    }
}


template std::shared_ptr<
  ParticleWallContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
set_particle_wall_contact_force_model<2, DEM::DEMProperties::PropertiesIndex>(
  const DEMSolverParameters<2>                    &dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation);
template std::shared_ptr<
  ParticleWallContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
set_particle_wall_contact_force_model<3, DEM::DEMProperties::PropertiesIndex>(
  const DEMSolverParameters<3>                    &dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation);

template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);

//////////////////////////

template std::shared_ptr<
  ParticleWallContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
set_particle_wall_contact_force_model<2,
                                      DEM::CFDDEMProperties::PropertiesIndex>(
  const DEMSolverParameters<2>                    &dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation);
template std::shared_ptr<
  ParticleWallContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
set_particle_wall_contact_force_model<3,
                                      DEM::CFDDEMProperties::PropertiesIndex>(
  const DEMSolverParameters<3>                    &dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation);

template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);


//////////////////////////
template std::shared_ptr<
  ParticleWallContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
set_particle_wall_contact_force_model<2, DEM::DEMMPProperties::PropertiesIndex>(
  const DEMSolverParameters<2>                    &dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation);
template std::shared_ptr<
  ParticleWallContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
set_particle_wall_contact_force_model<3, DEM::DEMMPProperties::PropertiesIndex>(
  const DEMSolverParameters<3>                    &dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation);

template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<2>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<3>         &dem_parameters,
  const std::vector<types::boundary_id> boundary_index,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
