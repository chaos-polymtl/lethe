// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters_lagrangian.h>

#include <dem/particle_wall_contact_force_new.h>

using namespace DEM;
using namespace Parameters::Lagrangian;

template <int dim,
          typename PropertiesIndex,
          ParticleWallContactForceModel contact_model,
          RollingResistanceMethod       rolling_friction_model>
ParticleWallContactForce<dim,
                         PropertiesIndex,
                         contact_model,
                         rolling_friction_model>::
  ParticleWallContactForce(const DEMSolverParameters<dim> &dem_parameters)
  : dmt_cut_off_threshold(dem_parameters.model_parameters.dmt_cut_off_threshold)
  , f_coefficient_epsd(dem_parameters.model_parameters.f_coefficient_epsd)
{
  set_effective_wall_properties(dem_parameters);
}

template <int dim,
          typename PropertiesIndex,
          ParticleWallContactForceModel contact_model,
          RollingResistanceMethod       rolling_friction_model>
void
ParticleWallContactForce<dim,
                         PropertiesIndex,
                         contact_model,
                         rolling_friction_model>::
  calculate_particle_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                &particle_wall_pairs_in_contact,
    const double dt,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome)
{}

template <int dim,
          typename PropertiesIndex,
          ParticleWallContactForceModel contact_model,
          RollingResistanceMethod       rolling_friction_model>
void
ParticleWallContactForce<dim,
                         PropertiesIndex,
                         contact_model,
                         rolling_friction_model>::
  calculate_particle_floating_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
                &particle_floating_mesh_in_contact,
    const double dt,
    const std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solids,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome)
{}

// dem
// No resistance
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;

// EPSD resistance
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;

// cfd_dem
// No resistance
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;

// EPSD resistance
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;

// dem_mp
//  No resistance
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;

// EPSD resistance
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
