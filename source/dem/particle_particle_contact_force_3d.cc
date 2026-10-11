// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dem_properties.h>
#include <core/parameters_lagrangian.h>

#include <dem/particle_particle_contact_force.templates.h>

using namespace Parameters::Lagrangian;

// Explicit instantiations of ParticleParticleContactForce in 3D. The 2D ones
// are in particle_particle_contact_force_2d.cc: the instantiations are
// distributed over two source files so that they are compiled in parallel.

// dem
// No resistance
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::shift,
  RollingResistanceMethod::none>;

// Constant resistance
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant>;

// Viscous resistance
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous>;

// EPSD resistance
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::epsd>;

// cfd_dem
// No resistance
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::none>;

// Constant resistance
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant>;

// Viscous resistance
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous>;

// EPSD resistance
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::epsd>;

// dem_mp
// No resistance
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::none>;

// Constant resistance
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant>;

// Viscous resistance
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous>;

// EPSD resistance
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::epsd>;

// cfd_dem_mp
// No resistance
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::none>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::none>;

// Constant resistance
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant>;

// Viscous resistance
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous>;

// EPSD resistance
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::epsd>;
