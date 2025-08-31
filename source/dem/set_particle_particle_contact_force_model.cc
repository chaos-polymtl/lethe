// SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/force_chains_visualization.h>
#include <dem/set_particle_particle_contact_force_model.h>

using namespace dealii;

using namespace Parameters::Lagrangian;


template <int dim, typename PropertiesIndex>
std::shared_ptr<ParticleParticleContactForceBase<dim, PropertiesIndex>>
set_particle_particle_contact_force_model(
  const DEMSolverParameters<dim> &dem_parameters)
{
  ParticleParticleContactForceModel particle_particle_contact_force_model =
    dem_parameters.model_parameters.particle_particle_contact_force_model;

  std::shared_ptr<ParticleParticleContactForceBase<dim, PropertiesIndex>>
    particle_particle_contact_force_object;

  switch (particle_particle_contact_force_model)
    {
      case ParticleParticleContactForceModel::linear:
        {
          set_rolling_resistance_model<
            dim,
            PropertiesIndex,
            ParticleParticleContactForceModel::linear>(
            dem_parameters, particle_particle_contact_force_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz_mindlin_limit_overlap:
        {
          set_rolling_resistance_model<
            dim,
            PropertiesIndex,
            ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
            dem_parameters, particle_particle_contact_force_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz_mindlin_limit_force:
        {
          set_rolling_resistance_model<
            dim,
            PropertiesIndex,
            ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
            dem_parameters, particle_particle_contact_force_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz:
        {
          set_rolling_resistance_model<
            dim,
            PropertiesIndex,
            ParticleParticleContactForceModel::hertz>(
            dem_parameters, particle_particle_contact_force_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz_JKR:
        {
          set_rolling_resistance_model<
            dim,
            PropertiesIndex,
            ParticleParticleContactForceModel::hertz_JKR>(
            dem_parameters, particle_particle_contact_force_object);
          break;
        }
      case ParticleParticleContactForceModel::DMT:
        {
          set_rolling_resistance_model<dim,
                                       PropertiesIndex,
                                       ParticleParticleContactForceModel::DMT>(
            dem_parameters, particle_particle_contact_force_object);
          break;
        }
      default:
        {
          throw std::runtime_error(
            "The chosen particle-particle contact force model is invalid");
          break;
        }
    }
  return particle_particle_contact_force_object;
}


template <int dim,
          typename PropertiesIndex,
          Parameters::Lagrangian::ParticleParticleContactForceModel
            particle_particle_contact_force_model>
void
set_rolling_resistance_model(
  const DEMSolverParameters<dim> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<dim, PropertiesIndex>>
    &particle_particle_contact_force_object)
{
  RollingResistanceMethod rolling_resistance_method =
    dem_parameters.model_parameters.rolling_resistance_method;

  switch (rolling_resistance_method)
    {
      case RollingResistanceMethod::none:
        {
          particle_particle_contact_force_object = std::make_shared<
            ParticleParticleContactForce<dim,
                                         PropertiesIndex,
                                         particle_particle_contact_force_model,
                                         RollingResistanceMethod::none>>(
            dem_parameters);
          break;
        }
      case RollingResistanceMethod::constant:
        {
          particle_particle_contact_force_object = std::make_shared<
            ParticleParticleContactForce<dim,
                                         PropertiesIndex,
                                         particle_particle_contact_force_model,
                                         RollingResistanceMethod::constant>>(
            dem_parameters);
          break;
        }
      case RollingResistanceMethod::viscous:
        {
          particle_particle_contact_force_object = std::make_shared<
            ParticleParticleContactForce<dim,
                                         PropertiesIndex,
                                         particle_particle_contact_force_model,
                                         RollingResistanceMethod::viscous>>(
            dem_parameters);
          break;
        }
      case RollingResistanceMethod::epsd:
        {
          particle_particle_contact_force_object = std::make_shared<
            ParticleParticleContactForce<dim,
                                         PropertiesIndex,
                                         particle_particle_contact_force_model,
                                         RollingResistanceMethod::epsd>>(
            dem_parameters);
          break;
        }
      default:
        throw std::runtime_error("Invalid rolling resistance method");
    }
}


template <int dim, typename PropertiesIndex>
std::shared_ptr<ParticlesForceChainsBase<dim, PropertiesIndex>>
set_force_chains_contact_force_model(
  const DEMSolverParameters<dim> &dem_parameters)
{
  ParticleParticleContactForceModel particle_particle_contact_force_model =
    dem_parameters.model_parameters.particle_particle_contact_force_model;

  std::shared_ptr<ParticlesForceChainsBase<dim, PropertiesIndex>>
    particles_force_chains_object;

  switch (particle_particle_contact_force_model)
    {
      case ParticleParticleContactForceModel::linear:
        {
          set_rolling_resistance_model<
            dim,
            PropertiesIndex,
            ParticleParticleContactForceModel::linear>(
            dem_parameters, particles_force_chains_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz_mindlin_limit_overlap:
        {
          set_rolling_resistance_model<
            dim,
            PropertiesIndex,
            ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
            dem_parameters, particles_force_chains_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz_mindlin_limit_force:
        {
          set_rolling_resistance_model<
            dim,
            PropertiesIndex,
            ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
            dem_parameters, particles_force_chains_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz:
        {
          set_rolling_resistance_model<
            dim,
            PropertiesIndex,
            ParticleParticleContactForceModel::hertz>(
            dem_parameters, particles_force_chains_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz_JKR:
        {
          set_rolling_resistance_model<
            dim,
            PropertiesIndex,
            ParticleParticleContactForceModel::hertz_JKR>(
            dem_parameters, particles_force_chains_object);
          break;
        }
      case ParticleParticleContactForceModel::DMT:
        {
          set_rolling_resistance_model<dim,
                                       PropertiesIndex,
                                       ParticleParticleContactForceModel::DMT>(
            dem_parameters, particles_force_chains_object);
          break;
        }
      default:
        {
          throw std::runtime_error(
            "The chosen particle-particle contact force model is invalid");
          break;
        }
    }
  return particles_force_chains_object;
}


template <int dim,
          typename PropertiesIndex,
          Parameters::Lagrangian::ParticleParticleContactForceModel
            particle_particle_contact_force_model>
void
set_rolling_resistance_model(
  const DEMSolverParameters<dim> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<dim, PropertiesIndex>>
    &particles_force_chains_object)
{
  RollingResistanceMethod rolling_resistance_method =
    dem_parameters.model_parameters.rolling_resistance_method;

  switch (rolling_resistance_method)
    {
      case RollingResistanceMethod::none:
        {
          particles_force_chains_object = std::make_shared<
            ParticlesForceChains<dim,
                                 PropertiesIndex,
                                 particle_particle_contact_force_model,
                                 RollingResistanceMethod::none>>(
            dem_parameters);
          break;
        }
      case RollingResistanceMethod::constant:
        {
          particles_force_chains_object = std::make_shared<
            ParticlesForceChains<dim,
                                 PropertiesIndex,
                                 particle_particle_contact_force_model,
                                 RollingResistanceMethod::constant>>(
            dem_parameters);
          break;
        }
      case RollingResistanceMethod::viscous:
        {
          particles_force_chains_object = std::make_shared<
            ParticlesForceChains<dim,
                                 PropertiesIndex,
                                 particle_particle_contact_force_model,
                                 RollingResistanceMethod::viscous>>(
            dem_parameters);
          break;
        }
      case RollingResistanceMethod::epsd:
        {
          particles_force_chains_object = std::make_shared<
            ParticlesForceChains<dim,
                                 PropertiesIndex,
                                 particle_particle_contact_force_model,
                                 RollingResistanceMethod::epsd>>(
            dem_parameters);
          break;
        }
      default:
        throw std::runtime_error("Invalid rolling resistance method");
    }
}


template std::shared_ptr<
  ParticleParticleContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
set_particle_particle_contact_force_model<2,
                                          DEM::DEMProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<
  ParticleParticleContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
set_particle_particle_contact_force_model<3,
                                          DEM::DEMProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &dem_parameters);

template std::shared_ptr<
  ParticlesForceChainsBase<2, DEM::DEMProperties::PropertiesIndex>>
set_force_chains_contact_force_model<2, DEM::DEMProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<
  ParticlesForceChainsBase<3, DEM::DEMProperties::PropertiesIndex>>
set_force_chains_contact_force_model<3, DEM::DEMProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &dem_parameters);

template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMProperties::PropertiesIndex>>
    &particles_force_chains_object);

//////////////////////////
template std::shared_ptr<
  ParticleParticleContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
set_particle_particle_contact_force_model<
  2,
  DEM::CFDDEMProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<
  ParticleParticleContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
set_particle_particle_contact_force_model<
  3,
  DEM::CFDDEMProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &dem_parameters);

template std::shared_ptr<
  ParticlesForceChainsBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
set_force_chains_contact_force_model<2, DEM::CFDDEMProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<
  ParticlesForceChainsBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
set_force_chains_contact_force_model<3, DEM::CFDDEMProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &dem_parameters);

template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::CFDDEMProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::CFDDEMProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::CFDDEMProperties::PropertiesIndex>>
    &particles_force_chains_object);

//////////////////////////
template std::shared_ptr<
  ParticleParticleContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
set_particle_particle_contact_force_model<
  2,
  DEM::DEMMPProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<
  ParticleParticleContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
set_particle_particle_contact_force_model<
  3,
  DEM::DEMMPProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &dem_parameters);

template std::shared_ptr<
  ParticlesForceChainsBase<2, DEM::DEMMPProperties::PropertiesIndex>>
set_force_chains_contact_force_model<2, DEM::DEMMPProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<
  ParticlesForceChainsBase<3, DEM::DEMMPProperties::PropertiesIndex>>
set_force_chains_contact_force_model<3, DEM::DEMMPProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &dem_parameters);

template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticleParticleContactForceBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<2, DEM::DEMMPProperties::PropertiesIndex>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::DEMMPProperties::PropertiesIndex,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<
    ParticlesForceChainsBase<3, DEM::DEMMPProperties::PropertiesIndex>>
    &particles_force_chains_object);
