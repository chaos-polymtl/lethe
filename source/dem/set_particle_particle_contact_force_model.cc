// SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/force_chains_visualization.h>
#include <dem/set_particle_particle_contact_force_model.h>

using namespace dealii;

using namespace Parameters::Lagrangian;


template <int dim, DEM::SolverType solver_type>
std::shared_ptr<ParticleParticleContactForceBase<dim, solver_type>>
set_particle_particle_contact_force_model(
  const DEMSolverParameters<dim> &dem_parameters)
{
  ParticleParticleContactForceModel particle_particle_contact_force_model =
    dem_parameters.model_parameters.particle_particle_contact_force_model;

  std::shared_ptr<ParticleParticleContactForceBase<dim, solver_type>>
    particle_particle_contact_force_object;

  switch (particle_particle_contact_force_model)
    {
      case ParticleParticleContactForceModel::linear:
        {
          set_rolling_resistance_model<
            dim,
            solver_type,
            ParticleParticleContactForceModel::linear>(
            dem_parameters, particle_particle_contact_force_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz_mindlin_limit_overlap:
        {
          set_rolling_resistance_model<
            dim,
            solver_type,
            ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
            dem_parameters, particle_particle_contact_force_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz_mindlin_limit_force:
        {
          set_rolling_resistance_model<
            dim,
            solver_type,
            ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
            dem_parameters, particle_particle_contact_force_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz:
        {
          set_rolling_resistance_model<
            dim,
            solver_type,
            ParticleParticleContactForceModel::hertz>(
            dem_parameters, particle_particle_contact_force_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz_JKR:
        {
          set_rolling_resistance_model<
            dim,
            solver_type,
            ParticleParticleContactForceModel::hertz_JKR>(
            dem_parameters, particle_particle_contact_force_object);
          break;
        }
      case ParticleParticleContactForceModel::DMT:
        {
          set_rolling_resistance_model<dim,
                                       solver_type,
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


template <int             dim,
          DEM::SolverType solver_type,
          Parameters::Lagrangian::ParticleParticleContactForceModel
            particle_particle_contact_force_model>
void
set_rolling_resistance_model(
  const DEMSolverParameters<dim> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<dim, solver_type>>
    &particle_particle_contact_force_object)
{
  RollingResistanceMethod rolling_resistance_method =
    dem_parameters.model_parameters.rolling_resistance_method;

  switch (rolling_resistance_method)
    {
      case RollingResistanceMethod::no_resistance:
        {
          particle_particle_contact_force_object =
            std::make_shared<ParticleParticleContactForce<
              dim,
              solver_type,
              particle_particle_contact_force_model,
              RollingResistanceMethod::no_resistance>>(dem_parameters);
          break;
        }
      case RollingResistanceMethod::constant_resistance:
        {
          particle_particle_contact_force_object =
            std::make_shared<ParticleParticleContactForce<
              dim,
              solver_type,
              particle_particle_contact_force_model,
              RollingResistanceMethod::constant_resistance>>(dem_parameters);
          break;
        }
      case RollingResistanceMethod::viscous_resistance:
        {
          particle_particle_contact_force_object =
            std::make_shared<ParticleParticleContactForce<
              dim,
              solver_type,
              particle_particle_contact_force_model,
              RollingResistanceMethod::viscous_resistance>>(dem_parameters);
          break;
        }
      default:
        throw std::runtime_error("Invalid rolling resistance method");
    }
}


template <int dim, DEM::SolverType solver_type>
std::shared_ptr<ParticlesForceChainsBase<dim, solver_type>>
set_force_chains_contact_force_model(
  const DEMSolverParameters<dim> &dem_parameters)
{
  ParticleParticleContactForceModel particle_particle_contact_force_model =
    dem_parameters.model_parameters.particle_particle_contact_force_model;

  std::shared_ptr<ParticlesForceChainsBase<dim, solver_type>>
    particles_force_chains_object;

  switch (particle_particle_contact_force_model)
    {
      case ParticleParticleContactForceModel::linear:
        {
          set_rolling_resistance_model<
            dim,
            solver_type,
            ParticleParticleContactForceModel::linear>(
            dem_parameters, particles_force_chains_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz_mindlin_limit_overlap:
        {
          set_rolling_resistance_model<
            dim,
            solver_type,
            ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
            dem_parameters, particles_force_chains_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz_mindlin_limit_force:
        {
          set_rolling_resistance_model<
            dim,
            solver_type,
            ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
            dem_parameters, particles_force_chains_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz:
        {
          set_rolling_resistance_model<
            dim,
            solver_type,
            ParticleParticleContactForceModel::hertz>(
            dem_parameters, particles_force_chains_object);
          break;
        }
      case ParticleParticleContactForceModel::hertz_JKR:
        {
          set_rolling_resistance_model<
            dim,
            solver_type,
            ParticleParticleContactForceModel::hertz_JKR>(
            dem_parameters, particles_force_chains_object);
          break;
        }
      case ParticleParticleContactForceModel::DMT:
        {
          set_rolling_resistance_model<dim,
                                       solver_type,
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


template <int             dim,
          DEM::SolverType solver_type,
          Parameters::Lagrangian::ParticleParticleContactForceModel
            particle_particle_contact_force_model>
void
set_rolling_resistance_model(
  const DEMSolverParameters<dim> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<dim, solver_type>>
    &particles_force_chains_object)
{
  RollingResistanceMethod rolling_resistance_method =
    dem_parameters.model_parameters.rolling_resistance_method;

  switch (rolling_resistance_method)
    {
      case RollingResistanceMethod::no_resistance:
        {
          particles_force_chains_object = std::make_shared<
            ParticlesForceChains<dim,
                                 solver_type,
                                 particle_particle_contact_force_model,
                                 RollingResistanceMethod::no_resistance>>(
            dem_parameters);
          break;
        }
      case RollingResistanceMethod::constant_resistance:
        {
          particles_force_chains_object = std::make_shared<
            ParticlesForceChains<dim,
                                 solver_type,
                                 particle_particle_contact_force_model,
                                 RollingResistanceMethod::constant_resistance>>(
            dem_parameters);
          break;
        }
      case RollingResistanceMethod::viscous_resistance:
        {
          particles_force_chains_object = std::make_shared<
            ParticlesForceChains<dim,
                                 solver_type,
                                 particle_particle_contact_force_model,
                                 RollingResistanceMethod::viscous_resistance>>(
            dem_parameters);
          break;
        }
      default:
        throw std::runtime_error("Invalid rolling resistance method");
    }
}


template std::shared_ptr<
  ParticleParticleContactForceBase<2, DEM::SolverType::dem>>
set_particle_particle_contact_force_model<2, DEM::SolverType::dem>(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<
  ParticleParticleContactForceBase<3, DEM::SolverType::dem>>
set_particle_particle_contact_force_model<3, DEM::SolverType::dem>(
  const DEMSolverParameters<3> &dem_parameters);

template std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::dem>>
set_force_chains_contact_force_model<2, DEM::SolverType::dem>(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::dem>>
set_force_chains_contact_force_model<3, DEM::SolverType::dem>(
  const DEMSolverParameters<3> &dem_parameters);

template void
set_rolling_resistance_model<2,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<2, DEM::SolverType::dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<3, DEM::SolverType::dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  2,
  DEM::SolverType::dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<2, DEM::SolverType::dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  3,
  DEM::SolverType::dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<3, DEM::SolverType::dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  2,
  DEM::SolverType::dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<2, DEM::SolverType::dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  3,
  DEM::SolverType::dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<3, DEM::SolverType::dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<2, DEM::SolverType::dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<3, DEM::SolverType::dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<2, DEM::SolverType::dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<3, DEM::SolverType::dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<2, DEM::SolverType::dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<3, DEM::SolverType::dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  2,
  DEM::SolverType::dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  3,
  DEM::SolverType::dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  2,
  DEM::SolverType::dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  3,
  DEM::SolverType::dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::dem,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::dem>>
    &particles_force_chains_object);

//////////////////////////
template std::shared_ptr<
  ParticleParticleContactForceBase<2, DEM::SolverType::cfd_dem>>
set_particle_particle_contact_force_model<2, DEM::SolverType::cfd_dem>(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<
  ParticleParticleContactForceBase<3, DEM::SolverType::cfd_dem>>
set_particle_particle_contact_force_model<3, DEM::SolverType::cfd_dem>(
  const DEMSolverParameters<3> &dem_parameters);

template std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::cfd_dem>>
set_force_chains_contact_force_model<2, DEM::SolverType::cfd_dem>(
  const DEMSolverParameters<2> &dem_parameters);
template std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::cfd_dem>>
set_force_chains_contact_force_model<3, DEM::SolverType::cfd_dem>(
  const DEMSolverParameters<3> &dem_parameters);

template void
set_rolling_resistance_model<2,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<2, DEM::SolverType::cfd_dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<3, DEM::SolverType::cfd_dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  2,
  DEM::SolverType::cfd_dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<2, DEM::SolverType::cfd_dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  3,
  DEM::SolverType::cfd_dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<3, DEM::SolverType::cfd_dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  2,
  DEM::SolverType::cfd_dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<2, DEM::SolverType::cfd_dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<
  3,
  DEM::SolverType::cfd_dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<3, DEM::SolverType::cfd_dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<2, DEM::SolverType::cfd_dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<3, DEM::SolverType::cfd_dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<2, DEM::SolverType::cfd_dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<3, DEM::SolverType::cfd_dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<2, DEM::SolverType::cfd_dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<3, DEM::SolverType::cfd_dem>>
    &particle_particle_contact_force_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::cfd_dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::linear>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::cfd_dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  2,
  DEM::SolverType::cfd_dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::cfd_dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  3,
  DEM::SolverType::cfd_dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::cfd_dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  2,
  DEM::SolverType::cfd_dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::cfd_dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<
  3,
  DEM::SolverType::cfd_dem,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::cfd_dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::cfd_dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::hertz>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::cfd_dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::cfd_dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::hertz_JKR>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::cfd_dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<2,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<2> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<2, DEM::SolverType::cfd_dem>>
    &particles_force_chains_object);
template void
set_rolling_resistance_model<3,
                             DEM::SolverType::cfd_dem,
                             ParticleParticleContactForceModel::DMT>(
  const DEMSolverParameters<3> &dem_parameters,
  std::shared_ptr<ParticlesForceChainsBase<3, DEM::SolverType::cfd_dem>>
    &particles_force_chains_object);
