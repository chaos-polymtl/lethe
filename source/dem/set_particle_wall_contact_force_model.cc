// SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


#include <dem/set_particle_wall_contact_force_model.h>


using namespace dealii;

using namespace Parameters::Lagrangian;

template <int dim, typename PropertiesIndex>
std::shared_ptr<ParticleWallContactForceBase<dim, PropertiesIndex>>
set_particle_wall_contact_force_model(
  const DEMSolverParameters<dim> &dem_parameters)
{
  ParticleWallContactForceModel particle_wall_contact_force_model =
    dem_parameters.model_parameters.particle_wall_contact_force_method;

  std::shared_ptr<ParticleWallContactForceBase<dim, PropertiesIndex>>
    particle_wall_contact_force_object;

  switch (particle_wall_contact_force_model)
    {
      case ParticleWallContactForceModel::linear:
        {
          set_rolling_resistance_model<dim,
                                       PropertiesIndex,
                                       ParticleWallContactForceModel::linear>(
            dem_parameters, particle_wall_contact_force_object);
          break;
        }
      case ParticleWallContactForceModel::nonlinear:
        {
          set_rolling_resistance_model<
            dim,
            PropertiesIndex,
            ParticleWallContactForceModel::nonlinear>(
            dem_parameters, particle_wall_contact_force_object);
          break;
        }
      case ParticleWallContactForceModel::JKR:
        {
          set_rolling_resistance_model<dim,
                                       PropertiesIndex,
                                       ParticleWallContactForceModel::JKR>(
            dem_parameters, particle_wall_contact_force_object);
          break;
        }
      case ParticleWallContactForceModel::DMT:
        {
          set_rolling_resistance_model<dim,
                                       PropertiesIndex,
                                       ParticleWallContactForceModel::DMT>(
            dem_parameters, particle_wall_contact_force_object);
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
  const DEMSolverParameters<dim> &dem_parameters,
  std::shared_ptr<ParticleWallContactForceBase<dim, PropertiesIndex>>
    &particle_wall_contact_force_object)
{
  RollingResistanceMethod rolling_resistance_method =
    dem_parameters.model_parameters.rolling_resistance_method;

  switch (rolling_resistance_method)
    {
      case RollingResistanceMethod::none:
        {
          particle_wall_contact_force_object = std::make_shared<
            ParticleWallContactForce<dim,
                                     PropertiesIndex,
                                     particle_wall_contact_force_model,
                                     RollingResistanceMethod::none>>(
            dem_parameters);
          break;
        }
      case RollingResistanceMethod::constant:
        {
          particle_wall_contact_force_object = std::make_shared<
            ParticleWallContactForce<dim,
                                     PropertiesIndex,
                                     particle_wall_contact_force_model,
                                     RollingResistanceMethod::constant>>(
            dem_parameters);
          break;
        }
      case RollingResistanceMethod::viscous:
        {
          particle_wall_contact_force_object = std::make_shared<
            ParticleWallContactForce<dim,
                                     PropertiesIndex,
                                     particle_wall_contact_force_model,
                                     RollingResistanceMethod::viscous>>(
            dem_parameters);
          break;
        }
      case RollingResistanceMethod::epsd:
        {
          particle_wall_contact_force_object = std::make_shared<
            ParticleWallContactForce<dim,
                                     PropertiesIndex,
                                     particle_wall_contact_force_model,
                                     RollingResistanceMethod::epsd>>(
            dem_parameters);
          break;
        }
      default:
        throw std::runtime_error("Invalid rolling resistance method");
    }
}

// dem
template std::shared_ptr<
  ParticleWallContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
set_particle_wall_contact_force_model<2, DEM::DEMProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<
  ParticleWallContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
set_particle_wall_contact_force_model<3, DEM::DEMProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &dem_parameters);

template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);


// cfd_dem
template std::shared_ptr<
  ParticleWallContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
set_particle_wall_contact_force_model<2,
                                      DEM::CFDDEMProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<
  ParticleWallContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
set_particle_wall_contact_force_model<3,
                                      DEM::CFDDEMProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &dem_parameters);

template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);


// dem_mp
template std::shared_ptr<
  ParticleWallContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
set_particle_wall_contact_force_model<2, DEM::DEMMPProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<
  ParticleWallContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
set_particle_wall_contact_force_model<3, DEM::DEMMPProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &dem_parameters);

template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::DMT>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::JKR>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::linear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleWallContactForceModel::nonlinear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleWallContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_wall_contact_force_object);
